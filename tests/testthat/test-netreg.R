test_that("Gaussian case works", {
  skip_on_bioc()

  set.seed(42)

  alpha <- 2.3
  beta <- c(1.7, -0.4, 6)

  X <- matrix(rnorm(1000 * 3), nrow = 1000)
  Y <- rnorm(1000, X %*% beta + alpha)

  fit <- pareg::edgenet(X, Y, family = pareg::gaussian)

  expect_equal(as.vector(coef(fit)), c(alpha, beta), tolerance = 0.02)
})


test_that("Beta case works", {
  skip_on_bioc()

  set.seed(42)

  N <- 1000
  alpha <- 2.3
  beta <- c(1.7, 4)

  X <- matrix(rnorm(N * 2), nrow = N)
  eta <- X %*% beta + alpha
  mu <- pareg::beta_phi_var()$linkinv(eta)$numpy()
  phi <- 1
  Y <- rbeta(N, mu * phi, (1 - mu) * phi)

  Y <- transform_y(Y)

  fit_br <- betareg::betareg(Y ~ X)
  fit <- pareg::edgenet(X, Y, family = pareg::beta_phi_var())

  expect_equal(
    c(as.vector(coef(fit)), fit$dispersion),
    as.vector(coef(fit_br)),
    tolerance = 0.05
  )
  # expect_equal(as.vector(coef(fit)), c(alpha, beta), tolerance = 0.05)
})


test_that("Cross-validation works", {
  skip_on_bioc()

  set.seed(42)

  alpha <- 2.3
  beta <- c(1.7, -0.4, 6)

  X <- matrix(rnorm(1000 * 3), nrow = 1000)
  Y <- rnorm(1000, X %*% beta + alpha)

  # small iteration count to reduce runtime
  fit <- pareg::cv_edgenet(
    X,
    Y,
    family = pareg::gaussian,
    maxit = 10,
    lambda_range = c(0, 1)
  )

  # expect_equal(as.vector(coef(fit)), c(alpha, beta), tolerance = 0.02)
  expect_equal(dim(coef(fit)), c(4, 1))
})


test_that("Cross-validation with network regularization works", {
  skip_on_bioc()

  set.seed(42)

  alpha <- 2.3
  beta <- c(1.7, -0.4, 6)

  X <- matrix(rnorm(1000 * 3), nrow = 1000)
  Y <- rnorm(1000, X %*% beta + alpha)

  amat <- matrix(runif(3 * 3), 3, 3)

  # small iteration count to reduce runtime
  fit <- pareg::cv_edgenet(
    X,
    Y,
    G.X = amat,
    family = pareg::gaussian,
    maxit = 10,
    lambda_range = c(0, 1),
    psigx_range = c(0, 5, 10)
  )

  # expect_equal(as.vector(coef(fit)), c(alpha, beta), tolerance = 0.02)
  expect_equal(dim(coef(fit)), c(4, 1))
})


test_that("laplacian computation works", {
  mat_adj <- matrix(c(1, 0.5, 0, 0.5, 1, 0.8, 0, 0.8, 1), 3, 3)
  mat_lap <- compute_norm_laplacian(mat_adj)

  expect_equal(
    mat_lap,
    matrix(c(
      0.3333333, -0.2691910, 0,
      -0.2691910, 0.5652174, -0.3931785,
      0, -0.3931785, 0.4444444
    ), 3, 3),
    tolerance = 1e-7
  )
})
