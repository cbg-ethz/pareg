test_that("Gaussian case works", {
  set.seed(42)

  alpha <- 2.3
  beta <- c(1.7, -0.4, 6)

  X <- matrix(rnorm(1000 * 3), nrow = 1000)
  Y <- rnorm(1000, X %*% beta + alpha)

  fit <- pareg::edgenet(X, Y, family = pareg::gaussian)

  expect_equal(as.vector(coef(fit)), c(alpha, beta), tolerance = 0.02)
})


test_that("Beta case works", {
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
  set.seed(42)

  alpha <- 2.3
  beta <- c(1.7, -0.4, 6)

  X <- matrix(rnorm(1000 * 3), nrow = 1000)
  Y <- rnorm(1000, X %*% beta + alpha)

  # small iteration count to reduce runtime
  fit <- pareg::cv.edgenet(X, Y, family = pareg::gaussian, maxit = 10)

  # expect_equal(as.vector(coef(fit)), c(alpha, beta), tolerance = 0.02)
})


test_that("Cross-validation with network regularization works", {
  set.seed(42)

  alpha <- 2.3
  beta <- c(1.7, -0.4, 6)

  X <- matrix(rnorm(1000 * 3), nrow = 1000)
  Y <- rnorm(1000, X %*% beta + alpha)

  amat <- matrix(runif(3 * 3), 3, 3)

  # small iteration count to reduce runtime
  fit <- pareg::cv.edgenet(
    X,
    Y,
    G.X = amat,
    family = pareg::gaussian,
    maxit = 10
  )

  # expect_equal(as.vector(coef(fit)), c(alpha, beta), tolerance = 0.02)
})
