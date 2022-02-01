test_that("Gaussian case works", {
  set.seed(42)

  alpha <- 2.3
  beta <- c(1.7, -0.4, 6)

  X <- matrix(rnorm(1000 * 3), nrow = 1000)
  Y <- rnorm(1000, X %*% beta + alpha)

  fit <- netReg::edgenet(X, Y, family = netReg::gaussian())

  expect_equal(as.vector(coef(fit)), c(alpha, beta), tolerance = 0.02)
})


test_that("Beta case works", {
  set.seed(42)

  N <- 1000
  alpha <- 2.3
  beta <- c(1.7, 4)

  X <- matrix(rnorm(N * 2), nrow = N)
  eta <- X %*% beta + alpha
  mu <- netReg::beta_phi_var()$linkinv(eta)$numpy()
  phi <- 1
  Y <- rbeta(N, mu * phi, (1 - mu) * phi)

  Y <- pareg::transform_y(Y)

  # betareg::betareg(Y ~ X)
  fit <- netReg::edgenet(X, Y, family = netReg::beta_phi_var())

  # expect_equal(as.vector(coef(fit)), c(alpha, beta), tolerance = 0.05)
})
