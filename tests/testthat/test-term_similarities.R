test_that("term similarities can be computed", {
  df_terms <- data.frame(
    term = c("A", "A", "B", "B", "B", "C", "C", "C"),
    gene = c("a", "b", "a", "b", "c", "a", "c", "d")
  )

  mat <- compute_term_similarities(df_terms)

  mat_expected <- matrix(c(1, 0.667, 0.25, 0.667, 1, 0.5, 0.25, 0.5, 1), 3, 3)
  rownames(mat_expected) <- colnames(mat_expected) <- c("A", "B", "C")

  expect_equal(mat, mat_expected, tolerance = 0.001)
})


test_that("term similarities can be computed on tibble", {
  df_terms <- tibble(
    term = c("A", "A", "B", "B", "B", "C", "C", "C"),
    gene = c("a", "b", "a", "b", "c", "a", "c", "d")
  )

  mat <- compute_term_similarities(df_terms)

  mat_expected <- matrix(c(1, 0.667, 0.25, 0.667, 1, 0.5, 0.25, 0.5, 1), 3, 3)
  rownames(mat_expected) <- colnames(mat_expected) <- c("A", "B", "C")

  expect_equal(mat, mat_expected, tolerance = 0.001)
})
