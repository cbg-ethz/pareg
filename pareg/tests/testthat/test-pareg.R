test_that("package doesn't crash for trivial case", {
  # create synthetic data
  set.seed(42)

  df_genes <- data.frame(
    gene = paste("g", 1:20, sep = ""),
    pvalue = c(
      rbeta(10, .1, 1),
      rbeta(10, 1.1, 1)
    )
  )

  df_terms <- rbind(
    data.frame(
      name = "foo",
      gene = paste("g", 1:10, sep = "")
    ),
    data.frame(
      name = "bar",
      gene = paste("g", 11:20, sep = "")
    )
  )

  # run model
  res <- pareg(df_genes, df_terms)

  # check results
  expect_lt(
    res %>% filter(name == "foo") %>% pull(enrichment),
    res %>% filter(name == "bar") %>% pull(enrichment)
  )
})
