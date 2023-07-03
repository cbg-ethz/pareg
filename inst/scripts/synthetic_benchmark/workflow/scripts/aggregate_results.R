library(tidyverse)


# parameters
fname_list_enr <- snakemake@input$fname_list_enr
fname_list_study <- snakemake@input$fname_list_study
fname_list_benchmark <- snakemake@input$fname_list_benchmark

fname_out <- snakemake@output$fname
fname_benchmark_out <- snakemake@output$fname_benchmark

# aggregate
purrr::transpose(list(
  fname_enr = fname_list_enr,
  fname_study = fname_list_study
)) %>%
  purrr::map_dfr(function(x) {
    parts <- (strsplit(x$fname_enr, "/")[[1]])
    replicate <- parts[[length(parts) - 2]]

    study <- readRDS(x$fname_study)
    on_terms <- study$on_terms

    tmp <- read_csv(x$fname_enr)

    if (dim(tmp)[[1]] == 0) {
      # method execution was skipped
      return(tmp)
    }

    tmp %>%
      mutate(
        replicate = replicate,
        is_on_term = term %in% on_terms
      )
  }) %>%
  write_csv(fname_out)

df_time <- fname_list_benchmark %>%
  purrr::map_dfr(function (x) {
    # parse file name
    filename <- basename(x)
    filename_noext <- substr(filename, 1, nchar(filename) - 14)
    parts <- strsplit(filename_noext, "__")[[1]]

    method <- parts[[1]]
    replicate <- parts[[3]]

    # parse benchmark file
    tmp <- read_tsv(x)

    duration_time <- tmp[["h:m:s"]]
    duration_seconds <- as.numeric(lubridate::seconds(duration_time))

    data.frame(method = method, replicate = replicate, duration_time = duration_time, duration_seconds = duration_seconds)
  }) %>%
  write_csv(fname_benchmark_out)
