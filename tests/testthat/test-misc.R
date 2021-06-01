library(purrr)
library(ggplot2)

test_that("idea works", {
  # setup terms
  df_terms <- rbind(
    data.frame(
      term = "foo",
      gene = paste("g", 100:200, sep = "")
    ),
    data.frame(
      term = "control",
      gene = c(paste("g", 10:20, sep = ""), "g90")
    )
  )

  # compute enrichments
  set.seed(13)

  window_size <- 10
  genes <- paste("g", 90:109, sep = "")

  df_enr <- map_df(0:10, function(x) {
    df_genes <- data.frame(
      gene = genes,
      pvalue = c(
        rbeta(x, 1, 1), # uniform
        rbeta(window_size, .1, 1), # "significant" p-values
        rbeta(20 - x - window_size, 1, 1) # uniform
      )
    )

    df_res <- pareg(df_genes, df_terms)

    cbind(
      df_res,
      data.frame(x = x)
    )
  })

  # validation plot
  ggplot(df_enr, aes(x = x, y = enrichment)) +
    geom_line() +
    facet_grid(~ term) +
    theme_minimal()

  # assertions
  expect_true(
    df_enr %>%
      filter(x == 0 & term == "foo") %>%
      pull("enrichment") > 0
  )
  expect_true(
    df_enr %>%
    filter(x == 10 & term == "foo") %>%
    pull("enrichment") < 0
  )
})
