######
# Code adapted from https://github.com/dirmeier/netReg
######


### model definition


#' @noRd
#' @import tensorflow
linear.predictor <- function(alpha, beta, x) {
  eta <- tf$matmul(x, beta) + alpha
  eta
}


#' @noRd
model <- function(p, q, family) {
  keras::keras_model_custom(function(self) {
    self$alpha <- init_vector(q)
    self$beta <- init_matrix(p, q)
    if (family$family %in% c("beta_phi_lm")) {
      self$gamma <- init_matrix(p, q)
    }
    if (family$family %in% c("beta_phi_var")) {
      initializer <- tf$keras.initializers$Ones()
      self$dispersion <- tf$Variable(initializer(tf$shape(1), tensorflow::tf$float32))
    }
    self$family <- family
    self$linkinv <- family$linkinv
    self$init_weights <- self$get_weights()

    self$reinit <- function() {
      self$set_weights(self$init_weights)
    }


    function(x, mask = NULL, training = FALSE) {
      if (self$family$family %in% c("beta_phi_lm")) {
        eta <- linear.predictor(self$alpha, self$beta, x)
        eta_phi <- linear.predictor(self$alpha, self$gamma, x)

        return(list(mu_hat = self$linkinv(eta), phi_hat = self$family$secondary_linkinv(eta_phi)))
      }

      eta <- linear.predictor(self$alpha, self$beta, x)
      if (self$family$family %in% c("gamma", "inverse.gaussian")) {
        eta <- tf$exp(eta)
      }
      self$linkinv(eta)
    }
  })
}


### link functions ###


#' @noRd
identity <- function(x) x


#' @noRd
#' @importFrom tensorflow tf
exp <- function(x) tf$maximum(tf$exp(x), tf$float32$min)


#' @noRd
#' @importFrom tensorflow tf
inverse <- function(x) 1 / x


#' @noRd
#' @importFrom tensorflow tf
inverse.exp <- function(x) tf$exp(-x)


#' @noRd
#' @importFrom tensorflow tf
logistic <- function(x) 1 / (1 + cast_float(tf$exp(-x)))


#' @noRd
#' @import tensorflow
gcdf <- function(x) {
  std <- tfp$distributions$Normal(0, 1)
  std$cdf(x)
}


#' @noRd
#' @importFrom tensorflow tf
inverse.sqrt <- function(x) 1 / tf$sqrt(x)


### utils ###


#' @noRd
intercept <- function(Y, X, B, n) {
  (t(Y - X %*% B) %*% rep(1, n)) / n
}


#' @noRd
intercept.matrix <- function(n, alpha) {
  rep(1, n) %*% t(alpha)
}


#' @noRd
rss <- function(Y, Y.hat) {
  sum((Y - Y.hat)**2)
}


#' @noRd
check.matrices <- function(X, Y) {
  stopifnot(is.matrix(X), is.matrix(Y))
}


#' @noRd
check.graphs <- function(X, Y, G.X, G.Y, psigx, psigy) {
  if (psigx != 0 & any(dim(G.X) != dim(X)[2])) {
    stop("ncol(X) and dim(G.X) do not fit!")
  }
  if (psigy != 0 & any(dim(G.Y) != dim(Y)[2])) {
    stop("ncol(Y) and dim(G.Y) do not fit!")
  }
  if (is.matrix(G.X) & any(G.X < 0)) {
    stop("Some elements G.X<0; please use non-negative matrix!")
  }
  if (is.matrix(G.Y) & any(G.Y < 0)) {
    stop("Some elements G.Y<0; please use non-negative matrix!")
  }
}


#' @noRd
check.dimensions <- function(X, Y, n, p) {
  if (dim(X)[1] != n) {
    stop("X and Y have not same number of observations!")
  }
  if (dim(X)[1] != dim(Y)[1]) {
    stop("X and Y have not same number of observations!")
  }
  if (n != dim(Y)[1]) {
    stop("X and Y have not same number of observations!")
  }
  if (p < 2) {
    stop("Pls use a X matrix with at least 2 covariables!")
  }
}


#' @noRd
is.positive.numeric <- function(d) {
  is.numeric(d) && d > 0
}


#' @noRd
check.param <- function(param, comp, op, replace.with) {
  if (!is.na(param) & op(param, comp)) {
    warning(sprintf("%s < 0, setting to 0!", deparse(substitute(param))))
    param <- replace.with
  }

  param
}


# shamelessly copied from stats::glm
#' @noRd
get.family <- function(family) {
  if (is.character(family)) {
    family <- get(family, mode = "function")
  }
  if (is.function(family)) {
    family <- family()
  }
  if (is.null(family$family)) {
    stop("'family' not recognized", call. = FALSE)
  }

  family
}


#' @noRd
not.supported.yet <- function(family) {
  err <- sprintf(
    "family '%s' not supported yet. choose 'gaussian'/'binomial' please.",
    family
  )
  stop(err, call. = FALSE)
}


warn.experimental <- function(family) {
  warning(paste("family", family, "is still experimental. enjoy with care."),
    call. = FALSE
  )
}


#' @noRd
#' @import tensorflow
cast_float <- function(x) {
  tensorflow::tf$cast(x, tensorflow::tf$float32)
}


#' @noRd
#' @import tensorflow
constant_float <- function(x) {
  tensorflow::tf$constant(x, tensorflow::tf$float32)
}


#' @noRd
#' @import tensorflow
init_matrix <- function(m, n, trainable = TRUE) {
  initializer <- tf$keras.initializers$glorot_normal(23L)
  tensorflow::tf$Variable(initializer(shape(m, n), tensorflow::tf$float32),
    trainable = trainable
  )
}


#' @noRd
#' @import tensorflow
init_vector <- function(m, trainable = TRUE) {
  initializer <- tf$keras.initializers$glorot_normal(23L)
  tensorflow::tf$Variable(initializer(shape(m), tensorflow::tf$float32),
    trainable = trainable
  )
}

#' @noRd
#' @import tensorflow
init_zero_scalar <- function(trainable = TRUE) {
  tensorflow::tf$Variable(tf$zeros(shape(), tensorflow::tf$float32),
    trainable = trainable
  )
}


### edgenet methods ###


#' @export
#' @method coef edgenet
coef.edgenet <- function(object, ...) {
  alpha <- object$alpha
  beta <- object$beta
  coefs <- rbind(alpha, beta)
  rownames(coefs) <- c("(Intercept)", sprintf("x[%s]", seq(nrow(beta))))
  colnames(coefs) <- sprintf("y[%s]", seq(ncol(beta)))
  coefs
}


### response distribution families ###


#' @title Family objects for models
#'
#' @export
#' @docType methods
#' @rdname family-methods
#'
#' @description Family objects provide a convenient way to specify the details
#'  of the models used by \code{pareg}.
#'  See also \code{\link[stats:family]{stats::family}} for more details.
#'
#' @param link  name of a link function
#' @param object  a object for which the family shoulr be retured
#'  (e.g. \code{edgenet})
#' @param ... further arguments passed to methods
#'
#' @return An object of class \code{pareg.family}
#'  \item{family }{ name of the family}
#'  \item{link }{ name of the link function}
#'  \item{linkinv }{ inverse link function}
#'  \item{loss }{ loss function}
#' @examples
#' gaussian()
#' binomial("probit")$link
#' poisson()$linkinv
#' gamma()$linkinv
#' beta()$loss
#' inverse.gaussian()$loss
family <- function(object, ...) UseMethod("family")


#' @export
#' @rdname family-methods
gaussian <- function(link = c("identity")) {
  link <- match.arg(link)
  linkinv <- switch(
    link,
    "identity" = identity,
    stop("did not recognize link function", call. = FALSE)
  )

  .as.family(
    "gaussian",
    link,
    NULL,
    linkinv,
    gaussian.loss
  )
}


#' @export
#' @rdname family-methods
bernoulli <- function(link = c("logit", "probit", "log")) {
  link <- match.arg(link)
  linkinv <- switch(
    link,
    "logit" = logistic,
    "log" = exp,
    "probit" = gcdf,
    stop("did not recognize link function", call. = FALSE)
  )

  .as.family(
    "bernoulli",
    link,
    NULL,
    linkinv,
    bernoulli.loss
  )
}


#' @export
#' @rdname family-methods
beta <- function(link = c("logit", "probit", "log")) {
  warn.experimental("beta")
  link <- match.arg(link)
  linkfun <- switch(
    link,
    "logit" = stats::binomial(link = "logit")$linkfun,
    "log" = stats::binomial(link = "log")$linkfun,
    "probit" = stats::binomial(link = "probit")$linkfun,
  )
  linkinv <- switch(
    link,
    "logit" = logistic,
    "log" = exp,
    "probit" = gcdf,
    stop("did not recognize link function", call. = FALSE)
  )

  loss_constant_phi <- function(y, mu_hat, phi_hat, ...) {
    beta.loss(y, mu_hat, tf$ones(tf$shape(mu_hat)), ...)
  }

  .as.family(
    "beta",
    link,
    linkfun,
    linkinv,
    loss_constant_phi
  )
}


#' @export
#' @rdname family-methods
beta_phi_lm <- function(link = c("logit", "probit", "log")) {
  warn.experimental("beta_phi_lm")
  link <- match.arg(link)
  linkfun <- switch(
    link,
    "logit" = stats::binomial(link = "logit")$linkfun,
    "log" = stats::binomial(link = "log")$linkfun,
    "probit" = stats::binomial(link = "probit")$linkfun,
  )
  linkinv <- switch(
    link,
    "logit" = logistic,
    "log" = exp,
    "probit" = gcdf,
    stop("did not recognize link function", call. = FALSE)
  )

  .as.family(
    "beta_phi_lm",
    link,
    linkfun,
    linkinv,
    beta.loss,
    exp
  )
}


#' @export
#' @rdname family-methods
beta_phi_var <- function(link = c("logit", "probit", "log")) {
  warn.experimental("beta_phi_var")
  link <- match.arg(link)
  linkfun <- switch(
    link,
    "logit" = stats::binomial(link = "logit")$linkfun,
    "log" = stats::binomial(link = "log")$linkfun,
    "probit" = stats::binomial(link = "probit")$linkfun,
  )
  linkinv <- switch(
    link,
    "logit" = logistic,
    "log" = exp,
    "probit" = gcdf,
    stop("did not recognize link function", call. = FALSE)
  )

  .as.family(
    "beta_phi_var",
    link,
    linkfun,
    linkinv,
    beta.loss
  )
}


#' @noRd
.as.family <- function(family, link, linkfun, linkinv, loss, secondary_linkinv = NULL) {
  structure(
    list(
      family = family,
      link = link,
      linkfun = linkfun,
      linkinv = linkinv,
      loss = loss,
      secondary_linkinv = secondary_linkinv
    ),
    class = "pareg.family"
  )
}



### regularization terms ###


#' @noRd
#' @importFrom tensorflow tf
.lasso.penalty <- function(lambda, beta) {
  lambda * tf$reduce_sum(tf$abs(beta))
}


#' @noRd
#' @import tensorflow
.edgenet.x.penalty <- function(gx, beta) {
  tf$linalg$trace(tf$matmul(tf$transpose(beta), tf$matmul(gx, beta)))
}


#' @noRd
#' @import tensorflow
.edgenet.y.penalty <- function(gy, beta) {
  tf$linalg$trace(tf$matmul(beta, tf$matmul(gy, tf$transpose(beta))))
}


### loss terms ###


#' @noRd
#' @importFrom tensorflow tf
gaussian.loss <- function(y, mu.hat, ...) {
  obj <- tf$reduce_mean(tf$square(y - mu.hat))
  obj
}


#' @noRd
#' @importFrom tensorflow tf
#' @importFrom tfprobability tfd_bernoulli
bernoulli.loss <- function(y, mu.hat, ...) {
  obj <- 0
  for (j in seq(ncol(y))) {
    prob <- tfprobability::tfd_bernoulli(probs = mu.hat[, j])
    obj <- obj + tf$reduce_sum(prob$log_prob(y[, j]))
  }

  -obj
}


#' @noRd
#' @importFrom tensorflow tf
#' @importFrom tfprobability tfd_beta
beta.loss <- function(y, mu.hat, phi_hat, ...) {
  obj <- 0
  eps <- .Machine$double.eps * 1e9
  for (j in seq(ncol(y))) {
    mu <- mu.hat[, j]
    phi <- phi_hat[, j]

    # reparametrize
    # concentration1 := alpha = mu * phi
    p <- mu * phi
    # concentration0 := beta = (1. - mu) * phi
    q <- (1 - mu) * phi

    # deal with numerical instabilities
    p.trans <- tf$math$maximum(p, eps)
    q.trans <- tf$math$maximum(q, eps)

    prob <- tfprobability::tfd_beta(
      concentration1 = p.trans, concentration0 = q.trans
    )
    obj <- obj + tf$reduce_sum(prob$log_prob(y[, j]))
  }

  -obj
}


#' @noRd
#' @import tensorflow
edgenet.loss <- function(lambda, psigx, psigy, gx, gy, family) {
  invlink <- family$linkinv
  loss.function <- family$loss

  loss <- function(mod, x, y) {
    ret <- list()

    if (family$family %in% c("beta_phi_lm")) {
      res <- mod(x)
      obj <- loss.function(y, res$mu_hat, res$phi_hat)
    } else if (family$family %in% c("beta_phi_var")) {
      mu_hat <- mod(x)
      obj <- loss.function(
        y,
        mu_hat,
        tf$ones(tf$shape(mu_hat)) * mod$dispersion
      )
    } else {
      mu.hat <- mod(x)
      obj <- loss.function(y, mu.hat)
    }
    ret$likelihood <- obj$numpy()

    obj <- obj + .lasso.penalty(lambda, mod$beta)
    ret$lasso <- .lasso.penalty(lambda, mod$beta)$numpy()

    if (family$family %in% c("beta_phi_lm")) {
      obj <- obj + .lasso.penalty(lambda, mod$gamma)
      ret$lasso_phi <- .lasso.penalty(lambda, mod$gamma)$numpy()
    }

    if (!is.null(gx)) {
      obj <- obj + psigx * .edgenet.x.penalty(gx, mod$beta)
      ret$nf_x <- psigx * .edgenet.x.penalty(gx, mod$beta)$numpy()

      if (family$family %in% c("beta_phi_lm")) {
        obj <- obj + psigx * .edgenet.x.penalty(gx, mod$gamma)
        ret$nf_x_phi <- psigx * .edgenet.x.penalty(gx, mod$gamma)$numpy()
      }
    }
    if (!is.null(gy)) {
      obj <- obj + psigy * .edgenet.y.penalty(gy, mod$beta)
      ret$nf_y <- psigy * .edgenet.y.penalty(gy, mod$beta)$numpy()

      if (family$family %in% c("beta_phi_lm")) {
        obj <- obj + psigy * .edgenet.y.penalty(gy, mod$gamma)
        ret$nf_y_phi <- psigy * .edgenet.y.penalty(gy, mod$gamma)$numpy()
      }
    }

    ret$obj <- obj
    ret$total_loss <- obj$numpy()
    ret
  }

  loss
}


### fit function ###


#' @noRd
#' @import tensorflow
#' @importFrom purrr transpose
fit <- function(mod, loss, x, y, maxit = 1000, learning.rate = 0.03, thresh = 1e-4) {
  optimizer <- keras::optimizer_adam(learning.rate)
  lo.old <- Inf
  loss_hist <- vector("list", length = maxit)
  stopping_reason <- "max_iterations"

  for (step in seq_len(maxit)) {
    with(tf$GradientTape() %as% t, {
      loss_obj <- loss(mod, x, y)
      lo <- loss_obj$obj
    })

    loss_hist[[step]] <- loss_obj
    loss_hist[[step]][["obj"]] <- NULL

    gradients <- t$gradient(lo, mod$trainable_variables)

    if (any(is.na(gradients[[1]]$numpy()))) {
      stopping_reason <- "nan_gradient"
      break
    }
    if (any(is.na(gradients[[2]]$numpy()))) {
      stopping_reason <- "nan_gradient"
      break
    }

    optimizer$apply_gradients(purrr::transpose(list(
      gradients, mod$trainable_variables
    )))

    if (step %% 25 == 0) {
      if (sum(abs(lo$numpy() - lo.old)) < thresh) {
        stopping_reason <- "loss_threshold"
        break
      }
      lo.old <- lo$numpy()
    }
  }

  ret <- list(
    beta = mod$beta$numpy(),
    alpha = mod$alpha$numpy(),
    loss_hist = loss_hist[!sapply(loss_hist, is.null)],
    stopping_reason = stopping_reason
  )

  if (mod$family$family %in% c("beta_phi_lm")) {
    ret$gamma <- mod$gamma$numpy()
  }
  if (mod$family$family %in% c("beta_phi_var")) {
    ret$dispersion <- mod$dispersion$numpy()
  }

  ret
}


### main entrypoint ###


#' @title Fit a graph-regularized linear regression model using
#'  edge-based regularization. Adapted from https://github.com/dirmeier/netReg.
#'
#' @export
#' @docType methods
#' @rdname edgenet-methods
#'
#' @importFrom stats gaussian binomial
#'
#' @description  Fit a graph-regularized linear regression model using
#'  edge-penalization. The coefficients are computed using graph-prior
#'  knowledge in the form of one/two affinity matrices. Graph-regularization is
#'  an extension to previously introduced regularization techniques,
#'  such as the LASSO. See the vignette for details on the objective function of
#'  the model: \href{../doc/edgenet.html}{\code{vignette("edgenet", package="netReg")}}
#'
#' @param X  input matrix, of dimension (\code{n} x \code{p})
#' where \code{n} is the number of observations and \code{p} is the number
#' of covariables. Each row is an observation vector.
#' @param Y  output matrix, of dimension (\code{n} x \code{q})
#' where \code{n} is the number of observations and \code{q} is the number
#' of response variables. Each row is an observation vector.
#' @param G.X  non-negativ affinity matrix for \code{X}, of dimensions
#' (\code{p} x \code{p}) where \code{p} is the number of covariables
#' @param G.Y  non-negativ affinity matrix for \code{Y}, of dimensions
#' (\code{q} x \code{q}) where \code{q} is the number of responses
#' @param lambda  \code{numerical} shrinkage parameter for LASSO.
#' @param psigx  \code{numerical} shrinkage parameter for graph-regularization
#'  of \code{G.X}
#' @param psigy  \code{numerical} shrinkage parameter for graph-regularization
#'  of \code{G.Y}
#' @param thresh  \code{numerical} threshold for optimizer
#' @param maxit  maximum number of iterations for optimizer
#'  (\code{integer})
#' @param learning.rate   step size for Adam optimizer (\code{numerical})
#' @param family  family of response, e.g. \emph{gaussian} or \emph{binomial}
#'
#' @return An object of class \code{edgenet}
#' \item{beta }{ the estimated (\code{p} x \code{q})-dimensional
#'  coefficient matrix B.hat}
#' \item{alpha }{ the estimated (\code{q} x \code{1})-dimensional
#'  vector of intercepts}
#' \item{parameters }{ regularization parameters}
#' \item{lambda }{ regularization parameter lambda)}
#' \item{psigx }{ regularization parameter psigx}
#' \item{psigy }{ regularization parameter psigy}
#' \item{family }{ a description of the error distribution and link function
#'    to be used. Can be a \code{\link[pareg:family]{pareg::family}} function or a character string
#'    naming a family function, e.g. \code{gaussian} or "gaussian".}
#' \item{call }{ the call that produced the object}
#'
#' @examples
#' X <- matrix(rnorm(100 * 10), 100, 10)
#' b <- matrix(rnorm(100), 10)
#' G.X <- abs(rWishart(1, 10, diag(10))[, , 1])
#' G.Y <- abs(rWishart(1, 10, diag(10))[, , 1])
#' diag(G.X) <- diag(G.Y) <- 0
#'
#' # estimate the parameters of a Gaussian model
#' Y <- X %*% b + matrix(rnorm(100 * 10), 100)
#' ## dont use affinity matrices
#' fit <- edgenet(X = X, Y = Y, family = gaussian, maxit = 10)
#' ## only provide one matrix
#' fit <- edgenet(X = X, Y = Y, G.X = G.X, psigx = 1, family = gaussian, maxit = 10)
#' ## use two matrices
#' fit <- edgenet(X = X, Y = Y, G.X = G.X, G.Y, family = gaussian, maxit = 10)
#' ## if Y is vectorial, we cannot use an affinity matrix for Y
#' fit <- edgenet(X = X, Y = Y[, 1], G.X = G.X, family = gaussian, maxit = 10)
#' @references
#'  Cheng, Wei and Zhang, Xiang and Guo, Zhishan and Shi, Yu and Wang, Wei (2014),
#'  Graph-regularized dual Lasso for robust eQTL mapping. \cr
#'  \emph{Bioinformatics}
#'
setGeneric(
  "edgenet",
  function(X, Y, G.X = NULL, G.Y = NULL,
           lambda = 0, psigx = 0, psigy = 0,
           thresh = 1e-5, maxit = 1e5, learning.rate = 0.01,
           family = gaussian) {
    standardGeneric("edgenet")
  },
  package = "pareg"
)


#' @rdname edgenet-methods
setMethod(
  "edgenet",
  signature = signature(X = "matrix", Y = "numeric"),
  function(X, Y, G.X = NULL, G.Y = NULL,
           lambda = 0, psigx = 0, psigy = 0,
           thresh = 1e-5, maxit = 1e5, learning.rate = 0.01,
           family = gaussian) {
    edgenet(
      X, as.matrix(Y), G.X, G.Y,
      lambda, psigx, psigy,
      thresh, maxit, learning.rate,
      family
    )
  }
)


#' @rdname edgenet-methods
setMethod(
  "edgenet",
  signature = signature(X = "matrix", Y = "matrix"),
  function(X, Y, G.X = NULL, G.Y = NULL,
           lambda = 0, psigx = 0, psigy = 0,
           thresh = 1e-5, maxit = 1e5, learning.rate = 0.01,
           family = gaussian) {
    stopifnot(
      is.numeric(maxit), is.numeric(thresh),
      is.numeric(learning.rate)
    )

    if (is.null(G.X)) psigx <- 0
    if (is.null(G.Y)) psigy <- 0

    check.matrices(X, Y)
    check.graphs(X, Y, G.X, G.Y, psigx, psigy)
    check.dimensions(X, Y, nrow(X), ncol(X))
    lambda <- check.param(lambda, 0, `<`, 0)
    psigx <- check.param(psigx, 0, `<`, 0)
    psigy <- check.param(psigy, 0, `<`, 0)
    maxit <- check.param(maxit, 0, `<`, 1e5)
    thresh <- check.param(thresh, 0, `<`, 1e-5)
    family <- get.family(family)

    if (ncol(Y) == 1) {
      psigy <- 0
      G.Y <- NULL
    }

    # estimate coefficients
    ret <- .edgenet(
      x = X, y = Y, gx = G.X, gy = G.Y,
      lambda = lambda, psigx = psigx, psigy = psigy,
      thresh = thresh, maxit = maxit,
      learning.rate = learning.rate, family = family
    )

    ret$call <- match.call()
    class(ret) <- c(class(ret), "edgenet")

    ret
  }
)


#' @noRd
.edgenet <- function(x, y, gx, gy,
                     lambda, psigx, psigy,
                     thresh, maxit, learning.rate, family) {
  cols.x <- colnames(x)
  cols.y <- colnames(y)
  x <- cast_float(x)
  y <- cast_float(y)

  if (!is.null(gx)) {
    gx <- cast_float(laplacian_(gx))
  }
  if (!is.null(gy)) {
    gy <- cast_float(laplacian_(gy))
  }

  mod <- model(ncol(x), ncol(y), family)
  loss <- edgenet.loss(lambda, psigx, psigy, gx, gy, family)
  res <- fit(mod, loss, x, y, maxit, learning.rate, thresh)

  # finalize output
  beta <- res$beta
  alpha <- res$alpha
  rownames(beta) <- cols.x
  colnames(beta) <- cols.y

  gamma <- NULL
  if (family$family %in% c("beta_phi_lm")) {
    gamma <- res$gamma
    rownames(gamma) <- cols.x
    colnames(gamma) <- cols.y
  }

  mse <- NULL
  if (family$family %in% c("beta", "beta_phi_lm", "beta_phi_var")) {
    residuals <- y$numpy() - family$linkinv(x$numpy() %*% beta)$numpy()
    residual_degrees_of_freedom <- nrow(x) - ncol(x) - 1
    mse <- sum(residuals^2) / residual_degrees_of_freedom
  }

  ret <- list(
    beta = beta,
    alpha = alpha,
    gamma = gamma,
    dispersion = res$dispersion,
    parameters = c("lambda" = lambda, "psigx" = psigx, "psigy" = psigy),
    lambda = lambda,
    psigx = psigx,
    psigy = psigy,
    loss_hist = res$loss_hist,
    stopping_reason = res$stopping_reason,
    pseudo_r_squared = if (is.null(family$linkfun)) NA else cor(x$numpy() %*% beta, family$linkfun(y$numpy()))^2,
    mse = mse
  )

  ret$family <- family
  class(ret) <- paste0(family$family, ".edgenet")

  ret
}
