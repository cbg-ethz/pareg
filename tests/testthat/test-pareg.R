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
      term = "foo",
      gene = paste("g", 1:10, sep = "")
    ),
    data.frame(
      term = "bar",
      gene = paste("g", 11:20, sep = "")
    )
  )

  # run model
  res <- pareg(df_genes, df_terms)

  # check results
  expect_lt(
    res %>% as.data.frame() %>% filter(term == "foo") %>% pull(enrichment),
    res %>% as.data.frame() %>% filter(term == "bar") %>% pull(enrichment)
  )
})

test_that("Bernoulli family works", {
  # create synthetic data
  set.seed(42)

  df_genes <- data.frame(
    gene = paste("g", 1:200, sep = ""),
    pvalue = c(
      rbeta(100, .1, 1),
      rbeta(100, 1, 1)
    )
  )

  df_terms <- rbind(
    data.frame(
      term = "foo",
      gene = paste("g", 1:100, sep = "")
    ),
    data.frame(
      term = "bar",
      gene = paste("g", 101:200, sep = "")
    )
  )

  # run model
  res <- pareg(
    df_genes, df_terms,
    family = bernoulli,
    response_column_name = "pvalue_notsig"
  )

  # check results
  expect_lt(
    res %>% as.data.frame() %>% filter(term == "foo") %>% pull(enrichment),
    res %>% as.data.frame() %>% filter(term == "bar") %>% pull(enrichment)
  )
})

test_that("term input network mismatch leads to crash", {
  df_genes <- data.frame(gene = c("g1", "g2"), pvalue = c(0.1, 0.01))
  df_terms <- data.frame(term = c("A", "B"), gene = c("g1", "g2"))

  network <- matrix(c(1, 0.3, 0.3, 1), 2, 2)

  rownames(network) <- colnames(network) <- c("A", "B")
  res <- pareg(df_genes, df_terms, term_network = network)

  rownames(network) <- colnames(network) <- c("A", "C")
  expect_error(
    pareg(df_genes, df_terms, term_network = network),
    "The following covariates do not appear in term network: B"
  )
})
