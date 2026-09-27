#' Panel Unit Root Test Based on Recursive Detrending
#'
#' Implements the t-REC and t-RREC panel unit root tests of Westerlund (2015).
#' The t-REC assumes iid errors; the robust t-RREC accounts for serial
#' correlation, cross-sectional dependence, and heteroskedasticity.
#'
#' @param data A data frame in long format containing panel data.
#' @param var Character. Name of the variable to test.
#' @param panel_id Character. Name of the panel identifier variable.
#' @param time_id Character. Name of the time variable.
#' @param trend Integer. Polynomial trend degree for recursive detrending.
#'   \code{0} (default) = constant only; \code{1} = constant + linear trend;
#'   \code{2} = constant + linear + quadratic trend.
#' @param robust Logical. If \code{TRUE}, compute the robust t-RREC test that
#'   removes cross-sectional dependence via cross-section averaging and
#'   augments with BIC-selected lags. Default is \code{FALSE}.
#' @param maxlag Integer. Maximum lag for BIC selection in the robust test.
#'   If \code{-1} (default), automatically set as
#'   \code{floor(4 * (T/100)^(2/9))}.
#'
#' @return An object of class \code{"xtrec"} containing:
#' \describe{
#'   \item{statistic}{Numeric. The t-REC or t-RREC test statistic.}
#'   \item{pvalue}{Numeric. Left-tail p-value from the standard normal distribution.}
#'   \item{N}{Integer. Number of panel units.}
#'   \item{TT}{Integer. Number of time periods.}
#'   \item{trend}{Integer. Polynomial trend degree used.}
#'   \item{robust}{Logical. Whether robust t-RREC was computed.}
#'   \item{sigma2}{Numeric. Estimated error variance (t-REC only; \code{NA} for t-RREC).}
#'   \item{a_p}{Numeric. First-order bias coefficient from Westerlund (2015) Table 1.}
#'   \item{b_p}{Numeric. Second-order bias coefficient from Westerlund (2015) Table 1.}
#'   \item{kappa}{Numeric. Power neighbourhood exponent (0.5 for trend <= 0; 0.25 otherwise).}
#'   \item{maxlag_used}{Integer. Lag order used (t-RREC only).}
#'   \item{panel_stats}{Named numeric vector of summary statistics of per-panel
#'     test statistics: mean, median, sd, min, max.}
#'   \item{test_name}{Character. "t-REC" or "t-RREC".}
#'   \item{var}{Character. Name of the tested variable.}
#'   \item{panel_id}{Character. Name of the panel identifier.}
#'   \item{time_id}{Character. Name of the time variable.}
#' }
#'
#' @references
#' Westerlund, J. (2015). The effect of recursive detrending on panel unit root
#' tests. \emph{Journal of Econometrics}, 185(2), 453--467.
#' \doi{10.1016/j.jeconom.2014.06.015}
#'
#' @examples
#' dat <- grunfeld_data()
#' # t-REC (iid, intercept only)
#' res <- xtrec(dat, var = "invest", panel_id = "firm",
#'              time_id = "year", trend = 0L, robust = FALSE)
#' print(res)
#'
#' # t-RREC (robust, linear trend)
#' \donttest{
#' res2 <- xtrec(dat, var = "invest", panel_id = "firm",
#'               time_id = "year", trend = 1L, robust = TRUE)
#' summary(res2)
#' }
#'
#' @export
xtrec <- function(data, var, panel_id, time_id,
                  trend = 0L, robust = FALSE, maxlag = -1L) {

  # --- Validate inputs ---
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  for (v in c(var, panel_id, time_id)) {
    if (!v %in% names(data)) stop(sprintf("Variable '%s' not found in data.", v))
  }
  trend <- as.integer(trend)
  if (trend < 0L) stop("'trend' must be >= 0.")
  maxlag <- as.integer(maxlag)

  # --- Sort and check balanced panel ---
  data <- data[order(data[[panel_id]], data[[time_id]]), ]
  panels <- sort(unique(data[[panel_id]]))
  N <- length(panels)
  times  <- sort(unique(data[[time_id]]))
  TT <- length(times)

  if (N < 3L) stop("At least 3 panel units are required.")
  if (TT < (trend + 5L)) {
    stop(sprintf("Insufficient time periods. Need at least %d, have %d.",
                 trend + 5L, TT))
  }

  # Check balance
  obs_per_panel <- tapply(data[[time_id]], data[[panel_id]], length)
  if (any(obs_per_panel != TT)) stop("Panel must be strongly balanced (no gaps).")
  if (nrow(data) != N * TT) stop("Panel must be strongly balanced.")

  # --- Extract Y matrix (T x N) ---
  Y <- matrix(data[[var]], nrow = TT, ncol = N)

  if (maxlag == -1L) {
    maxlag <- max(1L, floor(4 * (TT / 100)^(2 / 9)))
  }

  # --- Build bias coefficients ---
  coefs <- .xtrec_get_coefficients(trend)
  a_p <- coefs[1]
  b_p <- coefs[2]
  kappa <- if (trend < 1L) 0.5 else 0.25

  # --- Compute ---
  if (!robust) {
    res <- .xtrec_trec(Y, trend)
    test_name <- "t-REC"
    maxlag_used <- 0L
    sigma2 <- res$sigma2
  } else {
    res <- .xtrec_trrec(Y, trend, maxlag)
    test_name <- "t-RREC"
    maxlag_used <- res$maxlag_used
    sigma2 <- NA_real_
  }

  # Panel-level summary statistics
  panel_stats_vec <- .xtrec_panel_stats(Y, trend)

  structure(
    list(
      statistic   = res$statistic,
      pvalue      = stats::pnorm(res$statistic),
      N           = N,
      TT          = TT,
      trend       = trend,
      robust      = robust,
      sigma2      = sigma2,
      a_p         = a_p,
      b_p         = b_p,
      kappa       = kappa,
      maxlag_used = maxlag_used,
      panel_stats = panel_stats_vec,
      test_name   = test_name,
      var         = var,
      panel_id    = panel_id,
      time_id     = time_id
    ),
    class = "xtrec"
  )
}


# ============================================================
# Internal helpers
# ============================================================

#' @keywords internal
.xtrec_get_coefficients <- function(p) {
  if (p <= 0L) {
    a_p <- 0.5
    b_p <- 1 / 3
  } else {
    a_p <- .xtrec_compute_ap(p)
    b_p <- .xtrec_compute_bp(p)
  }
  c(a_p, b_p)
}

#' @keywords internal
.xtrec_hmn <- function(m, n, p) {
  ((-1)^(m + n)) * (m + n - 1) *
    choose(p - 1 + m, p - n) *
    choose(p - 1 + n, p - m) *
    (choose(m + n - 2, m - 1))^2
}

#' @keywords internal
.xtrec_hmp <- function(m, p) {
  if (p < 1L) return(0)
  sum(vapply(seq_len(p), function(n) .xtrec_hmn(m, n, p), numeric(1)))
}

#' @keywords internal
.xtrec_compute_ap <- function(p) {
  if (p < 1L) return(0.5)
  ap <- 0.5
  for (k in seq_len(p)) {
    hk <- .xtrec_hmp(k, p)
    ap <- ap - hk / k
  }
  for (k in seq_len(p)) {
    for (m in seq_len(p)) {
      hk <- .xtrec_hmp(k, p)
      hm <- .xtrec_hmp(m, p)
      ap <- ap + hk * hm *
        ((k + 1) * (m + 2) + m * (m + 1)) /
        (2 * m * (m + k) * (m + 1) * (k + 1))
    }
  }
  ap
}

#' @keywords internal
.xtrec_compute_bp <- function(p) {
  bp <- 1 / 3
  for (k in seq_len(p)) {
    hk <- .xtrec_hmp(k, p)
    bp <- bp - hk *
      ((2 * k + 1) * (k + 2) + k * (2 * k + 5) + 4) /
      (6 * k * (k + 1) * (k + 2))
  }
  for (k in seq_len(p)) {
    for (m in seq_len(p)) {
      hk <- .xtrec_hmp(k, p)
      hm <- .xtrec_hmp(m, p)
      bp <- bp + hk * hm *
        (k * (k + 2) * (k + m + 1) + m) /
        (3 * m * k * (k + 1) * (k + 2) * (k + m + 1))
    }
  }
  bp
}

#' @keywords internal
.xtrec_build_dt <- function(t_val, p) {
  vapply(seq_len(p), function(k) t_val^k - (t_val - 1)^k, numeric(1))
}

#' @keywords internal
.xtrec_recursive_detrend <- function(y_vec, p) {
  TT <- length(y_vec)
  if (p <= 0L) return(y_vec)
  y_det <- rep(NA_real_, TT)
  DtDt_sum <- matrix(0, p, p)
  for (t_val in seq(2, TT)) {
    d_t <- .xtrec_build_dt(t_val, p)
    DtDt_sum <- DtDt_sum + outer(d_t, d_t)
    if (t_val <= p) next
    DtDt_inv <- solve(DtDt_sum)
    s_val <- 0
    for (k in seq(2, t_val)) {
      d_k <- .xtrec_build_dt(k, p)
      s_val <- s_val + y_vec[k] * as.numeric(d_k %*% DtDt_inv %*% d_t)
    }
    y_det[t_val] <- y_vec[t_val] - s_val
  }
  y_det
}

#' @keywords internal
.xtrec_trec <- function(Y, p) {
  TT_orig <- nrow(Y)
  N <- ncol(Y)
  # First difference
  dY <- Y[seq(2, TT_orig), ] - Y[seq(1, TT_orig - 1), ]
  TT_fd <- nrow(dY)
  t_start <- if (p <= 0L) 1L else p + 1L

  # Recursively detrend each panel
  dY_rd <- matrix(NA_real_, TT_fd, N)
  for (i in seq_len(N)) {
    dY_rd[, i] <- .xtrec_recursive_detrend(dY[, i], p)
  }

  # sigma^2
  valid_vals <- dY_rd[seq(t_start, TT_fd), ]
  valid_vals <- valid_vals[!is.na(valid_vals)]
  sigma2 <- sum(valid_vals^2) / length(valid_vals)

  # A_NT and B_NT
  A_NT <- 0
  B_NT <- 0
  for (i in seq_len(N)) {
    R_i <- numeric(TT_fd)
    cumsum_val <- 0
    for (t in seq(t_start, TT_fd)) {
      if (!is.na(dY_rd[t, i])) cumsum_val <- cumsum_val + dY_rd[t, i]
      R_i[t] <- cumsum_val
    }
    for (t in seq(t_start + 1L, TT_fd)) {
      if (!is.na(dY_rd[t, i])) {
        A_NT <- A_NT + R_i[t - 1] * dY_rd[t, i]
        B_NT <- B_NT + R_i[t - 1]^2
      }
    }
  }

  T_eff <- TT_orig
  A_NT <- A_NT / (sqrt(N) * T_eff)
  B_NT <- B_NT / (N * T_eff^2)
  statistic <- A_NT / (sqrt(sigma2) * sqrt(B_NT))

  list(statistic = statistic, sigma2 = sigma2)
}

#' @keywords internal
.xtrec_trrec <- function(Y, p, maxlag) {
  TT_orig <- nrow(Y)
  N <- ncol(Y)
  dY <- Y[seq(2, TT_orig), ] - Y[seq(1, TT_orig - 1), ]
  TT_fd <- nrow(dY)
  t_start <- if (p <= 0L) 1L else p + 1L

  dY_rd <- matrix(NA_real_, TT_fd, N)
  for (i in seq_len(N)) {
    dY_rd[, i] <- .xtrec_recursive_detrend(dY[, i], p)
  }

  # BIC lag selection per panel
  q_selected <- integer(N)
  for (i in seq_len(N)) {
    best_bic <- Inf
    best_q <- 0L
    for (qq in seq(0L, maxlag)) {
      t_reg_start <- max(t_start, qq + 2L)
      obs_idx <- seq(t_reg_start, TT_fd)
      n_obs <- length(obs_idx)
      if (n_obs < qq + 2L) next
      y_dep <- dY_rd[obs_idx, i]
      if (any(is.na(y_dep))) next
      if (qq == 0L) {
        rss <- sum(y_dep^2)
      } else {
        X_lags <- matrix(0, n_obs, qq)
        for (ll in seq_len(qq)) {
          lag_idx <- obs_idx - ll
          valid <- lag_idx >= 1L
          X_lags[valid, ll] <- dY_rd[lag_idx[valid], i]
          if (any(is.na(X_lags[, ll]))) {
            rss <- Inf
            break
          }
        }
        if (any(is.na(X_lags))) next
        beta <- tryCatch(
          solve(crossprod(X_lags), crossprod(X_lags, y_dep)),
          error = function(e) NULL
        )
        if (is.null(beta)) next
        rss <- sum((y_dep - X_lags %*% beta)^2)
      }
      bic_val <- n_obs * log(rss / n_obs) + qq * log(n_obs)
      if (bic_val < best_bic) {
        best_bic <- bic_val
        best_q <- as.integer(qq)
      }
    }
    q_selected[i] <- best_q
  }

  q_max_used <- max(q_selected)
  h_val <- max(p, q_max_used + 2L)
  if (h_val < 1L) h_val <- 1L
  if ((TT_fd - h_val + 1L) < 5L) stop("Insufficient time periods after lag augmentation.")

  # Cross-section average (factor proxy)
  f_hat <- numeric(TT_fd)
  for (t in seq_len(TT_fd)) {
    vals <- dY_rd[t, ]
    vals <- vals[!is.na(vals)]
    if (length(vals) > 0) f_hat[t] <- mean(vals)
  }

  # Defactor per panel
  r_defact <- matrix(NA_real_, TT_fd, N)
  for (i in seq_len(N)) {
    qi <- q_selected[i]
    t_reg_s <- h_val
    obs_idx <- seq(t_reg_s, TT_fd)
    n_obs <- length(obs_idx)
    if (n_obs < 3L) next
    y_vec <- dY_rd[obs_idx, i]
    if (any(is.na(y_vec))) next
    # Regressors: f_hat + qi lags
    X_aug <- matrix(f_hat[obs_idx], nrow = n_obs, ncol = 1)
    if (qi > 0L) {
      for (ll in seq_len(qi)) {
        lag_idx <- obs_idx - ll
        x_lag <- ifelse(lag_idx >= 1L, dY_rd[pmax(lag_idx, 1L), i], 0)
        X_aug <- cbind(X_aug, x_lag)
      }
    }
    beta_def <- tryCatch(
      solve(crossprod(X_aug), crossprod(X_aug, y_vec)),
      error = function(e) NULL
    )
    if (!is.null(beta_def)) {
      r_defact[obs_idx, i] <- y_vec - X_aug %*% beta_def
    }
  }

  # t-RREC statistic
  num_sum <- 0
  den_sum <- 0
  for (i in seq_len(N)) {
    vals_i <- r_defact[seq(h_val, TT_fd), i]
    vals_i_clean <- vals_i[!is.na(vals_i)]
    n_valid <- length(vals_i_clean)
    if (n_valid == 0L) next
    sig2_ei <- max(sum(vals_i_clean^2) / n_valid, 1e-15)

    R_rob <- numeric(TT_fd)
    cs_val <- 0
    for (t in seq(h_val, TT_fd)) {
      if (!is.na(r_defact[t, i])) cs_val <- cs_val + r_defact[t, i]
      R_rob[t] <- cs_val
    }
    for (t in seq(h_val + 1L, TT_fd)) {
      if (!is.na(r_defact[t, i])) {
        num_sum <- num_sum + (1 / sig2_ei) * R_rob[t - 1] * r_defact[t, i]
        den_sum <- den_sum + (1 / sig2_ei) * R_rob[t - 1]^2
      }
    }
  }

  statistic <- if (den_sum > 0) num_sum / sqrt(den_sum) else 0
  list(statistic = statistic, maxlag_used = q_max_used)
}

#' @keywords internal
.xtrec_panel_stats <- function(Y, p) {
  TT_orig <- nrow(Y)
  N <- ncol(Y)
  dY <- Y[seq(2, TT_orig), ] - Y[seq(1, TT_orig - 1), ]
  TT_fd <- nrow(dY)
  t_start <- if (p <= 0L) 1L else p + 1L

  panel_ts <- numeric(N)
  for (i in seq_len(N)) {
    dYrd_i <- .xtrec_recursive_detrend(dY[, i], p)
    sig2_i <- max(mean(dYrd_i[seq(t_start, TT_fd)]^2, na.rm = TRUE), 1e-15)
    R_i <- numeric(TT_fd)
    cs <- 0
    for (t in seq(t_start, TT_fd)) {
      if (!is.na(dYrd_i[t])) cs <- cs + dYrd_i[t]
      R_i[t] <- cs
    }
    A_i <- 0; B_i <- 0
    for (t in seq(t_start + 1L, TT_fd)) {
      if (!is.na(dYrd_i[t])) {
        A_i <- A_i + R_i[t - 1] * dYrd_i[t]
        B_i <- B_i + R_i[t - 1]^2
      }
    }
    panel_ts[i] <- if (B_i > 0 && sig2_i > 0) A_i / (sqrt(sig2_i) * sqrt(B_i)) else 0
  }

  c(
    mean   = mean(panel_ts),
    median = stats::median(panel_ts),
    sd     = stats::sd(panel_ts),
    min    = min(panel_ts),
    max    = max(panel_ts)
  )
}

#' Print Method for xtrec Objects
#'
#' @param x An object of class \code{"xtrec"}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{x}.
#' @export
print.xtrec <- function(x, ...) {
  trend_lab <- switch(as.character(x$trend),
    "0" = "Constant only",
    "1" = "Constant + linear trend",
    "2" = "Constant + linear + quadratic trend",
    sprintf("Polynomial degree %d", x$trend)
  )
  stars <- if (x$pvalue < 0.01) "***" else if (x$pvalue < 0.05) "**" else
            if (x$pvalue < 0.10) "*" else ""
  decision <- if (x$pvalue < (1 - 0.95)) "Reject H0" else "Fail to reject H0"

  message(strrep("-", 65))
  message("Westerlund (2015) Panel Unit Root Test")
  message(sprintf("Test: %s | Trend: %s", x$test_name, trend_lab))
  message(sprintf("N = %d panels, T = %d periods", x$N, x$TT))
  message(strrep("-", 65))
  message(sprintf("%-20s %12.4f", x$test_name, x$statistic))
  message(sprintf("%-20s %12.4f%s", "p-value", x$pvalue, stars))
  message(sprintf("%-20s %12s", "Decision (5%)", decision))
  message(strrep("-", 65))
  message(sprintf("a_p = %.5f, b_p = %.5f, kappa = %.2f", x$a_p, x$b_p, x$kappa))
  message("H0: All panels contain a unit root")
  message("H1: At least some panels are stationary")
  message("*** p<0.01, ** p<0.05, * p<0.10")
  invisible(x)
}

#' Summary Method for xtrec Objects
#'
#' @param object An object of class \code{"xtrec"}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{object}.
#' @export
summary.xtrec <- function(object, ...) {
  print(object)
  message("")
  message("Critical values (standard normal, left tail):")
  message(sprintf("  1%%: %.4f  5%%: %.4f  10%%: %.4f",
    stats::qnorm(0.01), stats::qnorm(0.05), stats::qnorm(0.10)))
  message("")
  message("Individual panel statistics:")
  message(sprintf("  Mean=%.4f  Median=%.4f  SD=%.4f  Min=%.4f  Max=%.4f",
    object$panel_stats["mean"], object$panel_stats["median"],
    object$panel_stats["sd"], object$panel_stats["min"], object$panel_stats["max"]))
  invisible(object)
}
