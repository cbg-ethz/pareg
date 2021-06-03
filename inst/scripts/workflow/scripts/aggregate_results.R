library(tidyverse)


# parameters
fname_list_enr <- snakemake@input$fname_list_enr
fname_list_study <- snakemake@input$fname_list_study

fname_out <- snakemake@output$fname

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

    read_csv(x$fname_enr) %>%
      mutate(
        replicate = replicate,
        is_on_term = term %in% on_terms
      )
  }) %>%
  write_csv(fname_out)
