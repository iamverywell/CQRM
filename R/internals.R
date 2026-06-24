# internals.R — Score-test primitive helpers and core fitting wrappers.
# These functions are NOT exported; they are called by the four user-facing
# CQRM functions defined in CQRM.R.
#
# Required packages (declared in DESCRIPTION Imports):
#   lqmm     — lqmControl(), lqm.fit.gs(), addnoise()
#   quantreg — (loaded via CQRM.approx; not needed here directly)


# ── Quantile check function ────────────────────────────────────────────────────

phi.tau <- function(tau, u) {
  tau - as.numeric(u < 0)
}


# ── Rank-score statistic for predictor j under null theta_j = 0 ───────────────

s.tau <- function(tau, logy, n, theta_null, x_null, xj) {
  sum(phi.tau(tau, logy - x_null %*% theta_null) *
        (xj - x_null %*% solve(t(x_null) %*% x_null) %*% t(x_null) %*% xj)) /
    sqrt(n)
}


# ── Inverse-jitter-variance weighted rank-score statistic ─────────────────────

s.tau.wi <- function(tau, logy, n, theta_null, x_null, xj, wi) {
  sum(wi * phi.tau(tau, logy - x_null %*% theta_null) *
        (xj - x_null %*% solve(t(x_null) %*% x_null) %*% t(x_null) %*% xj)) /
    sqrt(n)
}


# ── Variance of rank-score statistic at a single quantile ─────────────────────

v.tau <- function(tau, n, x_null, xj) {
  xjstar <- xj - x_null %*% solve(t(x_null) %*% x_null) %*% t(x_null) %*% xj
  tau * (1 - tau) * (t(xjstar) %*% xjstar) / n
}


# ── Covariance of rank-score statistics at two quantile levels ────────────────

cov.tau <- function(tau1, tau2, n, x_null, xj) {
  xjstar <- xj - x_null %*% solve(t(x_null) %*% x_null) %*% t(x_null) %*% xj
  (min(tau1, tau2) - tau1 * tau2) * (t(xjstar) %*% xjstar) / n
}


# ── Single-quantile score test for count data (jitter-based) ──────────────────
# Tests each predictor separately at each tau in tau_vec via M jittered
# replicates using lqm.fit.gs (Machado–Santos-Silva transformation).
# Returns a list with elements: theta, score, vscore, chi2, p.

qrcount.score.tau <- function(formula, data,
                               tau_vec  = seq(0.1, 0.9, 0.1),
                               M        = 10,
                               weights  = NULL,
                               offset   = NULL,
                               constrasts = NULL,
                               zeta     = 1e-5,
                               B        = 0.999,
                               cn       = NULL,
                               control  = list()) {
  mf <- match.call(expand.dots = FALSE)
  m  <- match(c("formula", "data", "weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf   <- eval(mf, parent.frame())
  mt   <- attr(mf, "terms")
  y    <- model.response(mf, "numeric")
  w    <- as.vector(model.weights(mf))
  if (!is.null(w) && !is.numeric(w)) stop("'weights' must be a numeric vector")
  if (is.null(w)) w <- rep(1, length(y))
  x  <- model.matrix(mt, mf, contrasts)
  p  <- ncol(x); n <- nrow(x)
  if (is.null(offset)) offset <- rep(0, n)

  control <- .make_control(control, y)

  theta_0  <- glm.fit(x = x, y = y, offset = offset,
                       family = poisson())$coefficients
  Z        <- replicate(M, lqmm::addnoise(y, centered = FALSE, B = B))
  K        <- length(tau_vec)
  theta_mat <- matrix(NA, K, p)
  v_mat <- s_mat <- p_mat <- chi2_mat <- matrix(NA, K, p)
  rownames(v_mat) <- rownames(s_mat) <-
    rownames(p_mat) <- rownames(chi2_mat) <- tau_vec

  for (k in seq_len(K)) {
    tau <- tau_vec[k]
    TZ  <- apply(Z, 2, function(z, off, tau, zeta)
      log(ifelse((z - tau) > zeta, z - tau, zeta)) - off,
      off = offset, tau = tau, zeta = zeta)
    fit <- apply(TZ, 2, function(yy, xx, wts, tau, ctrl, theta)
      lqmm::lqm.fit.gs(theta = theta, x = xx, y = yy, weights = wts,
                        tau = tau, control = ctrl),
      xx = x, wts = w, tau = tau, ctrl = control, theta = theta_0)
    theta_mat[k, ] <- rowMeans(sapply(fit, function(obj) obj$theta))

    for (j in seq_len(p)) {
      s_mat[k, j]    <- s.tau(tau, rowMeans(TZ), n,
                               as.matrix(theta_mat[k, -j]),
                               as.matrix(x[, -j]), x[, j])
      v_mat[k, j]    <- v.tau(tau, n, as.matrix(x[, -j]), x[, j])
      chi2_mat[k, j] <- s_mat[k, j]^2 / v_mat[k, j]
      p_mat[k, j]    <- 1 - pchisq(chi2_mat[k, j], df = 1)
    }
  }
  list(theta = theta_mat, score = s_mat, vscore = v_mat,
       chi2 = chi2_mat, p = p_mat)
}


# ── Composite (joint) score test for count data (jitter-based) ────────────────
# Aggregates rank-score statistics across all tau_vec quantiles via a joint
# chi-squared test with the analytic covariance matrix Sigma.
# Returns a list with elements: theta, score, sigma, chi2, p.

qrcount.score.composite <- function(formula, data,
                                     tau_vec  = seq(0.1, 0.9, 0.1),
                                     M        = 10,
                                     weights  = NULL,
                                     offset   = NULL,
                                     constrasts = NULL,
                                     zeta     = 1e-5,
                                     B        = 0.999,
                                     cn       = NULL,
                                     control  = list()) {
  mf <- match.call(expand.dots = FALSE)
  m  <- match(c("formula", "data", "weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf   <- eval(mf, parent.frame())
  mt   <- attr(mf, "terms")
  y    <- model.response(mf, "numeric")
  w    <- as.vector(model.weights(mf))
  if (!is.null(w) && !is.numeric(w)) stop("'weights' must be a numeric vector")
  if (is.null(w)) w <- rep(1, length(y))
  x  <- model.matrix(mt, mf, contrasts)
  p  <- ncol(x); n <- nrow(x)
  if (is.null(offset)) offset <- rep(0, n)

  control <- .make_control(control, y)

  theta_0  <- glm.fit(x = x, y = y, offset = offset,
                       family = poisson())$coefficients
  Z        <- replicate(M, lqmm::addnoise(y, centered = FALSE, B = B))
  K        <- length(tau_vec)
  theta_mat <- matrix(NA, K, p)
  TZ        <- NULL

  for (k in seq_len(K)) {
    tau <- tau_vec[k]
    TZ  <- apply(Z, 2, function(z, off, tau, zeta)
      log(ifelse((z - tau) > zeta, z - tau, zeta)) - off,
      off = offset, tau = tau, zeta = zeta)
    fit <- apply(TZ, 2, function(yy, xx, wts, tau, ctrl, theta)
      lqmm::lqm.fit.gs(theta = theta, x = xx, y = yy, weights = wts,
                        tau = tau, control = ctrl),
      xx = x, wts = w, tau = tau, ctrl = control, theta = theta_0)
    theta_mat[k, ] <- rowMeans(sapply(fit, function(obj) obj$theta))
  }

  s_mat      <- matrix(NA, K, p)
  sigma_list <- vector("list", p)
  pvec <- chi2vec <- rep(NA, p)

  for (j in seq_len(p)) {
    Sigma <- matrix(NA, K, K)
    for (k in seq_len(K)) {
      tau <- tau_vec[k]
      s_mat[k, j] <- s.tau(tau, rowMeans(TZ), n,
                            as.matrix(theta_mat[k, -j]),
                            as.matrix(x[, -j]), x[, j])
    }
    for (k1 in seq_len(K)) {
      for (k2 in seq_len(K)) {
        Sigma[k1, k2] <- if (k1 == k2)
          v.tau(tau_vec[k1], n, as.matrix(x[, -j]), x[, j])
        else
          cov.tau(tau_vec[k1], tau_vec[k2], n, as.matrix(x[, -j]), x[, j])
      }
    }
    if (!isSymmetric(Sigma)) Sigma[lower.tri(Sigma)] <- t(Sigma)[lower.tri(Sigma)]
    sigma_list[[j]] <- Sigma
    chi2vec[j]      <- t(s_mat[, j]) %*% solve(Sigma) %*% s_mat[, j]
    pvec[j]         <- 1 - pchisq(chi2vec[j], df = K)
  }
  list(theta = theta_mat, score = s_mat, sigma = sigma_list,
       chi2 = chi2vec, p = pvec)
}


# ── Internal: build lqmm control list with sensible defaults ──────────────────

.make_control <- function(control, y) {
  ctrl <- lqmm::lqmControl()
  if (length(names(control)) > 0) {
    nm <- intersect(names(control), names(ctrl))
    ctrl[nm] <- control[nm]
  }
  if (is.null(ctrl$loop_step)) ctrl$loop_step <- sd(as.numeric(y))
  if (ctrl$beta  > 1 || ctrl$beta  < 0) stop("Beta must be in (0,1)")
  if (ctrl$gamma < 1)                    stop("Gamma must be >= 1")
  ctrl
}
