test_that("model creation works with only common genes", {
  df.genes <- data.frame(
    gene=c("g1", "g2"),
    pvalue=c(0.1, 0.2)
  )
  df.terms <- data.frame(
    name=c("A", "A", "B", "B", "C"),
    gene=c("g1", "g2", "g1", "g2", "g2")
  )

  df.res <- create_model_df(df.genes, df.terms)
  df.expected <- data.frame(
    gene=c("g1", "g2"),
    pvalue=c(0.1, 0.2),
    A.member=c(TRUE, TRUE),
    B.member=c(TRUE, TRUE),
    C.member=c(FALSE, TRUE)
  )

  expect_equal(df.res, df.expected)
})

test_that("model creation works with gene which is in no term", {
  df.genes <- data.frame(
    gene=c("g1", "g2", "g3"),
    pvalue=c(0.1, 0.2, 0.3)
  )
  df.terms <- data.frame(
    name=c("A", "A", "B", "B", "C"),
    gene=c("g1", "g2", "g1", "g2", "g2")
  )

  df.res <- create_model_df(df.genes, df.terms)
  df.expected <- data.frame(
    gene=c("g1", "g2", "g3"),
    pvalue=c(0.1, 0.2, 0.3),
    A.member=c(TRUE, TRUE, FALSE),
    B.member=c(TRUE, TRUE, FALSE),
    C.member=c(FALSE, TRUE, FALSE)
  )

  expect_equal(df.res, df.expected)
})
