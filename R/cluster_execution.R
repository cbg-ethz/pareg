#' @noRd
#' @importFrom logger log_trace log_debug log_error
#' @importFrom tibble tibble add_row
#' @importFrom dplyr bind_rows
#' @importFrom glue glue
#' @importFrom stringr str_match
#' @importFrom progress progress_bar
cluster_apply <- function(
  df_iter,
  func,
  .bsub_params = c("-n", "2", "-W", "24:00", "-R", "rusage[mem=10000]"),
  .tempdir = ".",
  ...
) {
  dir.create(.tempdir, showWarnings = FALSE, recursive = TRUE)

  # submit jobs
  df_jobs <- tibble(
    index = numeric(),
    job_id = character(),
    result_path = character()
  )
  for (index in seq_len(nrow(df_iter))) {
    # save environment needed to execute function
    image_path <- file.path(.tempdir, glue("image_{index}.RData"))
    log_trace("[index={index}] Saving image to {image_path}")

    current_df_row <- df_iter[index, ]

    with(
      c(
        as.list(environment()), # current env for function and argument
        as.list(parent.frame()) # env of function to deal with globals
      ), {
        save(
          list = ls(),
          file = image_path
        )
      }
    )

    # create script
    script_path <- file.path(.tempdir, glue("script_{index}.R"))
    result_path <- file.path(.tempdir, glue("result_{index}.rds"))
    code <- glue("
      load('{image_path}')
      result <- do.call(func, current_df_row)
      saveRDS(result, '{result_path}')
    ")
    log_trace("[index={index}] Writing script to {script_path}")
    writeLines(code, script_path)

    # submit script
    bsub_param_str <- paste(
      c(.bsub_params, "Rscript", script_path),
      sep = "",
      collapse = " "
    )
    log_trace("[index={index}] Executing 'bsub {bsub_param_str}'")
    stdout <- system2(
      "bsub",
      c(.bsub_params, "Rscript", script_path),
      stdout = TRUE,
      stderr = FALSE
    )

    job_id <- str_match(stdout, "Job <(.*?)>")[1, 2]
    log_debug("[index={index}] Submitted job {job_id}")

    # finalize
    df_jobs <- add_row(
      df_jobs,
      index = index,
      job_id = job_id,
      result_path = result_path
    )
  }

  # check job status and retrieve results
  successful_job_list <- c()
  result_list <- list()

  pb <- progress_bar$new(total = index)
  pb$tick(0)

  while (length(successful_job_list) < nrow(df_jobs)) {
    for (i in seq_len(nrow(df_jobs))) {
      row <- df_jobs[i, ]

      # skip processed jobs
      if (row$job_id %in% successful_job_list) {
        next
      }

      # check job
      status <- system2(
        "bjobs",
        c("-o", "stat", "-noheader", row$job_id),
        stdout = TRUE
      )
      log_trace("Job {row$job_id} has status {status}")

      if (status == "DONE") {
        log_debug("Job {row$job_id} at index {row$index} is done")
        successful_job_list <- c(successful_job_list, row$job_id)

        result <- readRDS(row$result_path)
        result_list[[row$index]] <- result

        pb$tick()
      } else if (status == "EXIT") {
        log_error("Job {row$job_id} crashed, killing all other jobs")
        for (i in seq_len(nrow(df_jobs))) {
          system2("bkill", df_jobs[i, ]$job_id)
        }
        stop("Cluster job crashed")
      }

      Sys.sleep(0.1)
    }

    Sys.sleep(10)
  }

  # finalize
  return(bind_rows(result_list))
}
