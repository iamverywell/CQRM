# CQRM.R — Four user-facing differential expression test functions.
#
# All four functions return a numeric vector of length G (one calibrated
# p-value per gene).  Apply p.adjust(..., method = "BH") for FDR control.
#
# In all functions, solve(subsigma) for the null distribution is precomputed
# once per weight pattern outside the B-loop.


# ── CQRM ──────────────────────────────────────────────────────────────────────

#' Composite quantile regression MinP test (all 2^K-1 subsets)
#'
#' Permutation-calibrated MinP test for differential expression of count data.
#' Tests all \eqn{2^K - 1} non-empty subsets of \code{tau_vec} quantiles and
#' takes the minimum p-value (MinP). The MinP is calibrated against draws from
#' the asymptotic multivariate normal null distribution of the joint score
#' vector. Fitting uses \code{M} jittered replicates of the Machado–Santos-
#' Silva transformed response via \code{lqm.fit.gs}.
#'
#' @param count1 Gene \eqn{\times} sample integer count matrix, group 1.
#'   Genes in rows, samples in columns.
#' @param count2 Gene \eqn{\times} sample integer count matrix, group 2.
#'   Genes in rows, samples in columns.
#' @param tau_vec Numeric vector of quantile levels to aggregate.
#'   Default \code{c(0.2, 0.4, 0.6, 0.8)}.
#' @param B Integer. Number of draws from the null distribution for
#'   permutation calibration. Default \code{10000}.
#' @param seed Integer. Random seed for reproducibility. Default \code{12345}.
#' @param M Integer. Number of jitter replicates per gene per quantile.
#'   Default \code{50}.
#'
#' @return Numeric vector of length \code{nrow(count1)} containing one
#'   permutation-calibrated p-value per gene.
#'
#' @details
#' The score statistics follow a multivariate normal distribution under the
#' null. For each gene, the function computes a chi-squared statistic for
#' every non-empty subset of quantiles and takes the minimum p-value. The null
#' distribution of the minimum p-value is approximated by drawing \code{B}
#' realisations from \eqn{N(0, \Sigma)} where \eqn{\Sigma} is the analytic
#' covariance matrix. \code{solve(subsigma)} is precomputed once per weight
#' pattern outside the B-loop for efficiency.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' G  <- 100; n <- 10
#' count1 <- matrix(rnbinom(G * n, mu = 5,  size = 0.5), G, n)
#' count2 <- matrix(rnbinom(G * n, mu = 10, size = 0.5), G, n)
#' pv  <- CQRM(count1, count2, B = 500, M = 10)
#' fdr <- p.adjust(pv, method = "BH")
#' }
#'
#' @importFrom MASS mvrnorm
#' @importFrom lqmm lqmControl lqm.fit.gs addnoise
#' @export
CQRM <- function(count1, count2,
                 tau_vec = c(0.2, 0.4, 0.6, 0.8),
                 B = 10000, seed = 12345, M = 50) {

  n1 <- ncol(count1); n2 <- ncol(count2)
  cond_idx <- factor(c(rep(1, n1), rep(2, n2)))
  counts   <- as.matrix(cbind(count1, count2))
  rownames(counts) <- seq_len(nrow(counts))
  offset   <- log(colSums(counts))
  G        <- nrow(counts)
  K        <- length(tau_vec)

  wt_all <- expand.grid(rep(list(0:1), K))
  wt_all <- wt_all[-which(rowSums(wt_all) == 0), , drop = FALSE]
  L      <- nrow(wt_all)

  score_mat <- matrix(NA, G, K)
  ch2_mat   <- matrix(NA, G, L)
  p_mat     <- matrix(NA, G, L)
  sigma     <- NULL

  for (j in seq_len(G)) {
    y   <- counts[j, ]
    dat <- data.frame(y = y, x = cond_idx)
    out <- qrcount.score.composite(y ~ x, data = dat, offset = offset,
                                   tau_vec = tau_vec, M = M)
    score          <- out$score[, 2]
    sigma          <- out$sigma[[2]]
    score_mat[j, ] <- score

    for (l in seq_len(L)) {
      wt       <- as.logical(wt_all[l, ])
      subscore <- score[wt]
      subsigma <- sigma[wt, wt, drop = FALSE]
      chi2_wt  <- drop(t(subscore) %*% solve(subsigma) %*% subscore)
      ch2_mat[j, l] <- chi2_wt
      p_mat[j, l]   <- 1 - pchisq(chi2_wt, df = sum(wt))
    }
    if (j %% 100 == 0) cat(sprintf("gene %d / %d\n", j, G))
  }

  wqcountp_min    <- apply(p_mat, 1, min)
  wqcountp_min_df <- apply(p_mat, 1, function(x) sum(wt_all[which.min(x), ]))

  # Precompute solve(subsigma) once per weight pattern outside the B-loop
  set.seed(seed)
  score_null   <- MASS::mvrnorm(B, mu = rep(0, K), Sigma = sigma)
  wt_logical   <- lapply(seq_len(L), function(l) as.logical(wt_all[l, ]))
  subsigma_inv <- lapply(seq_len(L), function(l) {
    wt <- wt_logical[[l]]
    solve(sigma[wt, wt, drop = FALSE])
  })
  wt_df        <- apply(wt_all, 1, sum)
  wq_minp_null <- matrix(NA, B, L)
  subp         <- numeric(L)

  for (b in seq_len(B)) {
    bscore <- score_null[b, ]
    for (l in seq_len(L)) {
      wt      <- wt_logical[[l]]
      subscore <- bscore[wt]
      subchi2  <- drop(t(subscore) %*% subsigma_inv[[l]] %*% subscore)
      subp[l]  <- 1 - pchisq(subchi2, df = wt_df[l])
    }
    wq_minp_null[b, ] <- subp
  }

  vapply(seq_len(G), function(j) {
    null_indx <- wt_df == wqcountp_min_df[j]
    p_null_df <- if (sum(null_indx) > 1)
      apply(wq_minp_null[, null_indx, drop = FALSE], 1, min)
    else
      wq_minp_null[, null_indx]
    (sum(p_null_df < wqcountp_min[j]) + 1) / (B + 1)
  }, numeric(1))
}


# ── CQRM.best ─────────────────────────────────────────────────────────────────

#' Single-quantile MinP test (best-tau selection)
#'
#' Tests each quantile in \code{tau_vec} independently (one chi-squared test
#' per quantile, scalar variance) and takes the minimum p-value across
#' quantiles. The null distribution is generated analytically from independent
#' normals, making this faster than \code{\link{CQRM}}.
#'
#' @inheritParams CQRM
#'
#' @return Numeric vector of length \code{nrow(count1)} containing one
#'   calibrated p-value per gene.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' G  <- 100; n <- 10
#' count1 <- matrix(rnbinom(G * n, mu = 5,  size = 0.5), G, n)
#' count2 <- matrix(rnbinom(G * n, mu = 10, size = 0.5), G, n)
#' pv <- CQRM.best(count1, count2, B = 500, M = 10)
#' }
#'
#' @importFrom MASS mvrnorm
#' @importFrom lqmm lqmControl lqm.fit.gs addnoise
#' @export
CQRM.best <- function(count1, count2,
                      tau_vec = c(0.2, 0.4, 0.6, 0.8),
                      B = 10000, seed = 12345, M = 50) {

  n1 <- ncol(count1); n2 <- ncol(count2)
  cond_idx <- factor(c(rep(1, n1), rep(2, n2)))
  counts   <- as.matrix(cbind(count1, count2))
  rownames(counts) <- seq_len(nrow(counts))
  offset   <- log(colSums(counts))
  G        <- nrow(counts)
  K        <- length(tau_vec)

  qcountp <- matrix(NA, G, K)
  sigma2  <- NULL

  for (j in seq_len(G)) {
    y   <- counts[j, ]
    dat <- data.frame(y = y, x = cond_idx)
    out <- qrcount.score.tau(y ~ x, data = dat, tau_vec = tau_vec,
                             offset = offset, M = M)
    qcountp[j, ] <- out$p[, 2]
    sigma2       <- out$vscore[, 2]
    if (j %% 100 == 0) cat(sprintf("gene %d / %d\n", j, G))
  }

  qcount_minp <- apply(qcountp, 1, min)

  # Scalar null — no matrix solve needed
  set.seed(seed)
  p_null <- matrix(NA, B, K)
  for (k in seq_len(K)) {
    s_null     <- rnorm(B, mean = 0, sd = sqrt(sigma2[k]))
    chi2_null  <- s_null^2 / sigma2[k]
    p_null[, k] <- 1 - pchisq(chi2_null, df = 1)
  }
  minp_null <- apply(p_null, 1, min)

  vapply(seq_len(G), function(j) {
    (sum(minp_null < qcount_minp[j]) + 1) / (B + 1)
  }, numeric(1))
}


# ── CQRM.region ───────────────────────────────────────────────────────────────

#' Contiguous-window quantile MinP test
#'
#' Tests all contiguous windows of \code{k} consecutive quantiles from
#' \code{tau_vec} and takes the MinP across windows. Intermediate between
#' \code{\link{CQRM.best}} (\code{k = 1}) and \code{\link{CQRM}} (all
#' subsets). \code{solve(subsigma)} is precomputed once per window outside
#' the B-loop.
#'
#' @inheritParams CQRM
#' @param k Integer. Window width (number of consecutive quantiles per test).
#'   Default \code{2}.
#'
#' @return Numeric vector of length \code{nrow(count1)} containing one
#'   calibrated p-value per gene.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' G  <- 100; n <- 10
#' count1 <- matrix(rnbinom(G * n, mu = 5,  size = 0.5), G, n)
#' count2 <- matrix(rnbinom(G * n, mu = 10, size = 0.5), G, n)
#' pv <- CQRM.region(count1, count2, k = 2, B = 500, M = 10)
#' }
#'
#' @importFrom MASS mvrnorm
#' @importFrom lqmm lqmControl lqm.fit.gs addnoise
#' @export
CQRM.region <- function(count1, count2,
                         tau_vec = c(0.2, 0.4, 0.6, 0.8), k = 2,
                         B = 10000, seed = 12345, M = 50) {

  n1 <- ncol(count1); n2 <- ncol(count2)
  cond_idx <- factor(c(rep(1, n1), rep(2, n2)))
  counts   <- as.matrix(cbind(count1, count2))
  rownames(counts) <- seq_len(nrow(counts))
  offset   <- log(colSums(counts))
  G        <- nrow(counts)
  K        <- length(tau_vec)

  n_windows <- K - k + 1
  wt_all    <- matrix(0, nrow = n_windows, ncol = K)
  for (i in seq_len(n_windows)) wt_all[i, i:(i + k - 1)] <- 1
  L <- nrow(wt_all)

  score_mat <- matrix(NA, G, K)
  ch2_mat   <- matrix(NA, G, L)
  p_mat     <- matrix(NA, G, L)
  sigma     <- NULL

  for (j in seq_len(G)) {
    y   <- counts[j, ]
    dat <- data.frame(y = y, x = cond_idx)
    out <- qrcount.score.composite(y ~ x, data = dat, offset = offset,
                                   tau_vec = tau_vec, M = M)
    score          <- out$score[, 2]
    sigma          <- out$sigma[[2]]
    score_mat[j, ] <- score

    for (l in seq_len(L)) {
      wt       <- as.logical(wt_all[l, ])
      subscore <- score[wt]
      subsigma <- sigma[wt, wt, drop = FALSE]
      chi2_wt  <- drop(t(subscore) %*% solve(subsigma) %*% subscore)
      ch2_mat[j, l] <- chi2_wt
      p_mat[j, l]   <- 1 - pchisq(chi2_wt, df = sum(wt))
    }
    if (j %% 100 == 0) cat(sprintf("gene %d / %d\n", j, G))
  }

  wqcountp_min <- apply(p_mat, 1, min)

  # Precompute solve(subsigma) once per window outside the B-loop
  set.seed(seed)
  score_null   <- MASS::mvrnorm(B, mu = rep(0, K), Sigma = sigma)
  wt_logical   <- lapply(seq_len(L), function(l) as.logical(wt_all[l, ]))
  subsigma_inv <- lapply(seq_len(L), function(l) {
    wt <- wt_logical[[l]]
    solve(sigma[wt, wt, drop = FALSE])
  })
  wt_df        <- apply(wt_all, 1, sum)
  wq_minp_null <- numeric(B)
  subp         <- numeric(L)

  for (b in seq_len(B)) {
    bscore <- score_null[b, ]
    for (l in seq_len(L)) {
      wt      <- wt_logical[[l]]
      subscore <- bscore[wt]
      subchi2  <- drop(t(subscore) %*% subsigma_inv[[l]] %*% subscore)
      subp[l]  <- 1 - pchisq(subchi2, df = wt_df[l])
    }
    wq_minp_null[b] <- min(subp)
  }

  vapply(seq_len(G), function(j) {
    (sum(wq_minp_null < wqcountp_min[j]) + 1) / (B + 1)
  }, numeric(1))
}


# ── CQRM.approx ───────────────────────────────────────────────────────────────

#' Fast approximation of CQRM (no jittering)
#'
#' Approximation of \code{\link{CQRM}} without jittering. Replaces the
#' M-replicate \code{lqm.fit.gs} loop with a single \code{rq.wfit} call
#' (C-compiled LP solver from \pkg{quantreg}) on the un-jittered
#' Machado–Santos-Silva transformed response. Each quantile uses its own
#' tau-specific transformation in the score computation. The permutation null
#' and MinP calibration proceed identically to \code{CQRM}.
#' \code{solve(subsigma)} is precomputed once per weight pattern outside the
#' B-loop.
#'
#' @inheritParams CQRM
#' @param zeta Numeric. Minimum value for the log transformation (numerical
#'   stability). Default \code{1e-5}.
#'
#' @return Numeric vector of length \code{nrow(count1)} containing one
#'   calibrated p-value per gene.
#'
#' @note This function requires \pkg{quantreg} instead of \pkg{lqmm} and does
#'   not have an \code{M} argument (no jittering is performed).
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' G  <- 100; n <- 10
#' count1 <- matrix(rnbinom(G * n, mu = 5,  size = 0.5), G, n)
#' count2 <- matrix(rnbinom(G * n, mu = 10, size = 0.5), G, n)
#' pv  <- CQRM.approx(count1, count2, B = 500)
#' fdr <- p.adjust(pv, method = "BH")
#' }
#'
#' @importFrom MASS mvrnorm
#' @importFrom quantreg rq.wfit
#' @export
CQRM.approx <- function(count1, count2,
                          tau_vec = c(0.2, 0.4, 0.6, 0.8),
                          B = 10000, seed = 12345, zeta = 1e-5) {

  n1 <- ncol(count1); n2 <- ncol(count2)
  cond_idx <- factor(c(rep(1, n1), rep(2, n2)))
  counts   <- as.matrix(cbind(count1, count2))
  rownames(counts) <- seq_len(nrow(counts))
  offset   <- log(colSums(counts))
  G        <- nrow(counts)
  K        <- length(tau_vec)

  wt_all <- expand.grid(rep(list(0:1), K))
  wt_all <- wt_all[-which(rowSums(wt_all) == 0), , drop = FALSE]
  L      <- nrow(wt_all)

  score_mat <- matrix(NA, G, K)
  ch2_mat   <- matrix(NA, G, L)
  p_mat     <- matrix(NA, G, L)
  sigma     <- NULL

  for (j in seq_len(G)) {
    y   <- counts[j, ]
    dat <- data.frame(y = y, x = cond_idx)
    x   <- model.matrix(~ x, data = dat)
    n   <- nrow(x)

    theta_mat <- matrix(NA, K, ncol(x))
    tz_list   <- vector("list", K)

    for (k in seq_len(K)) {
      tau            <- tau_vec[k]
      tz             <- log(ifelse((y - tau) > zeta, y - tau, zeta)) - offset
      theta_mat[k, ] <- quantreg::rq.wfit(x, tz, weights = rep(1, n),
                                           tau = tau)$coefficients
      tz_list[[k]]  <- tz
    }

    s_vec <- numeric(K)
    Sigma <- matrix(NA, K, K)

    for (k in seq_len(K)) {
      tau      <- tau_vec[k]
      s_vec[k] <- s.tau(tau       = tau,
                        logy      = tz_list[[k]],
                        n         = n,
                        theta_null = as.matrix(theta_mat[k, -2]),
                        x_null    = as.matrix(x[, -2]),
                        xj        = x[, 2])
    }
    for (k1 in seq_len(K)) {
      for (k2 in seq_len(K)) {
        Sigma[k1, k2] <- if (k1 == k2)
          v.tau(tau_vec[k1], n, as.matrix(x[, -2]), x[, 2])
        else
          cov.tau(tau_vec[k1], tau_vec[k2], n, as.matrix(x[, -2]), x[, 2])
      }
    }
    if (!isSymmetric(Sigma)) Sigma[lower.tri(Sigma)] <- t(Sigma)[lower.tri(Sigma)]

    score_mat[j, ] <- s_vec
    sigma          <- Sigma

    for (l in seq_len(L)) {
      wt       <- as.logical(wt_all[l, ])
      subscore <- s_vec[wt]
      subsigma <- Sigma[wt, wt, drop = FALSE]
      chi2_wt  <- drop(t(subscore) %*% solve(subsigma) %*% subscore)
      ch2_mat[j, l] <- chi2_wt
      p_mat[j, l]   <- 1 - pchisq(chi2_wt, df = sum(wt))
    }
    if (j %% 100 == 0) cat(sprintf("gene %d / %d\n", j, G))
  }

  wqcountp_min    <- apply(p_mat, 1, min)
  wqcountp_min_df <- apply(p_mat, 1, function(x) sum(wt_all[which.min(x), ]))

  # Precompute solve(subsigma) once per weight pattern outside the B-loop
  set.seed(seed)
  score_null   <- MASS::mvrnorm(B, mu = rep(0, K), Sigma = sigma)
  wt_logical   <- lapply(seq_len(L), function(l) as.logical(wt_all[l, ]))
  subsigma_inv <- lapply(seq_len(L), function(l) {
    wt <- wt_logical[[l]]
    solve(sigma[wt, wt, drop = FALSE])
  })
  wt_df        <- apply(wt_all, 1, sum)
  wq_minp_null <- matrix(NA, B, L)
  subp         <- numeric(L)

  for (b in seq_len(B)) {
    bscore <- score_null[b, ]
    for (l in seq_len(L)) {
      wt      <- wt_logical[[l]]
      subscore <- bscore[wt]
      subchi2  <- drop(t(subscore) %*% subsigma_inv[[l]] %*% subscore)
      subp[l]  <- 1 - pchisq(subchi2, df = wt_df[l])
    }
    wq_minp_null[b, ] <- subp
  }

  vapply(seq_len(G), function(j) {
    null_indx <- wt_df == wqcountp_min_df[j]
    p_null_df <- if (sum(null_indx) > 1)
      apply(wq_minp_null[, null_indx, drop = FALSE], 1, min)
    else
      wq_minp_null[, null_indx]
    (sum(p_null_df < wqcountp_min[j]) + 1) / (B + 1)
  }, numeric(1))
}
