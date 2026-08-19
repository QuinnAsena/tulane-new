## ocfs_helpers.r
##
## Shared utilities for the OCFS TVARSS analysis. source() this AFTER
## kalman-smoother/TVARSS_12Aug26.r.
##
## Contents
##   1. transform          pw() / ipw(), the lambda = 0.25 power transform
##   2. measurement error  build_ocfs_me(), per-bin observation variance
##   3. data assembly      build_ocfs_data()
##   4. fitting            fit_ocfs(), fit_ms(), ascend()
##   5. inference          aicc(), drop_one(), longrun()
##   6. caching            cached()
##
## Design notes that matter for correctness:
##
##  * Predictors enter through c (the process equation), so their effect
##    accumulates through the autoregression. The interpretable quantity is the
##    long-run effect c / (1 - sum(b)), not c itself. longrun() computes it and
##    every reporting function carries both.
##
##  * TVARSS's optimiser is Nelder-Mead and it stops early on this problem. A
##    single fit is never trusted: fit_ms() multi-starts, and ascend() lifts each
##    reduced solution back into the full parameter space and re-optimises. Both
##    are required; multi-start alone was not enough on the log-scale version of
##    this problem (it reached -383.30 where the true optimum was -381.16).
##
##  * Every likelihood-ratio deviance is asserted non-negative. A negative
##    deviance between nested models means the optimiser failed, not that the
##    variable is unimportant, and it is the failure mode that silently produced
##    a table of P = 1 in the original script.

OCFS_LAMBDA <- 0.25   # power transform exponent; see reports/ for the derivation
BIN_WIDTH   <- 200    # years, must match tulane.R

# ---------------------------------------------------------------------------
# The degenerate branch, and why fits are filtered on se
# ---------------------------------------------------------------------------
# This likelihood has a second, higher-likelihood family of solutions in which
# the process-error SD collapses to ~0 while the autoregression sits close to the
# unit circle. With se = 0 the latent series is a DETERMINISTIC AR(2) recursion,
# and near the unit circle such a recursion is flexible enough to trace almost
# any smooth curve - so it chases the data and wins on likelihood by 15-25 log
# units. On this dataset it appears once the observation variance is inflated
# (tau >~ 1.25) and it is reached only from some starting values, which makes the
# tau profile non-reproducible and produces negative drop-one deviances.
#
# It is rejected on scientific grounds, not numerical ones: se = 0 asserts that a
# 60-kyr record of megaherbivore activity has NO process noise at all - that every
# fluctuation is measurement error and the underlying signal is exactly
# predictable from its own two previous values. That is not an admissible model of
# an ecological time series.
#
# Fits are therefore required to have se >= se.min, a fraction of sd(response).
#
# This filter is a GUARD, not a selection device, and the distinction matters. At
# tau = 0 - the primary observation model, where su is fixed entirely from the
# measured counting variance - se comes out around 0.31 of sd(response) and the
# guard never binds. It only starts to bind once the observation variance is
# inflated, which is exactly the regime the analysis avoids.
#
# It is deliberately set low (2% of sd) so that it excludes only genuine
# degeneracy. An earlier draft used 5%, which on this dataset sat right between
# two competing optima (se/sd of 0.0395 and 0.0522) and so silently decided which
# fit was reported - a threshold doing real inferential work is a fudge, not a
# guard. If you see the "all fits fell below se.min" warning, do not simply lower
# the floor: that means the model is unidentified at that tau and the fix is to
# reduce tau, not to admit the degenerate branch.
#
# Set se.min = 0 to disable the filter and see the degenerate solutions yourself.
SE_MIN_FRAC <- 0.02

# TRUE if a fit is usable: it exists, has a finite likelihood, and its process
# noise is not in the degenerate corner.
se_admissible <- function(m, se.min = 0) {
  !is.null(m) && !inherits(m, "try-error") && is.finite(m$logLik) &&
    (se.min <= 0 || (is.finite(m$se) && m$se >= se.min))
}

# =============================================================================
# 1. Transform
# =============================================================================

# lambda = 0.25 was chosen by asking which scale makes the *measurement error*
# homoscedastic and symmetric, using the per-sample quantiles in
# ocfs_uncertainty.csv. The optimum sits at 0.24-0.30 depending on the criterion
# (CV of the per-sample SD, max/min variance ratio, correlation between SD and
# value, interval asymmetry); 0.25 is a round number inside that window. Raw
# scale and sqrt both leave the error growing with the value; log overshoots and
# makes it shrink. Zero maps to zero, so the seven zero bins need no offset.
pw  <- function(x, lambda = OCFS_LAMBDA) x^lambda
ipw <- function(y, lambda = OCFS_LAMBDA) pmax(y, 0)^(1 / lambda)

# =============================================================================
# 2. Measurement error
# =============================================================================

# The 200-yr bin grid, reconstructed exactly as tulane.R builds it (lines
# 199-203 and 253-257). all_composite$age is pollen_wide_binned$age, so the grid
# can be rebuilt from the saved csv without re-running the pollen pipeline.
ocfs_bin_breaks <- function(all_composite, bin_width = BIN_WIDTH) {
  seq(from = min(all_composite$age, na.rm = TRUE),
      to   = max(all_composite$age + bin_width, na.rm = TRUE),
      by   = bin_width)
}

# Per-bin observation variance for OCFS on the lambda scale.
#
# Three sources, recorded per bin in the `source` column:
#   quantiles     - delta method on the stored Q2.5/Q97.5 (the normal case)
#   poisson_zero  - OCFS_total == 0, so no quantiles were computed. The exact
#                   Poisson 97.5% upper bound for an observed count of zero is
#                   qgamma(.975, 1) = 3.689 grains; convert to a concentration
#                   with that sample's own Lycopodium and Tracer_Added, then to
#                   the lambda scale, and take sd = upper / 1.96. The observation
#                   itself stays at zero. This is deliberately conservative: it
#                   down-weights the zero-count samples rather than inventing a
#                   value for them.
#   imputed       - neither available; median of the rest. Should be empty.
#
# CAVEAT worth knowing: for a zero count the lambda-scale interval is [0, upper],
# and describing that with a Gaussian centred on 0 puts half its mass below zero.
# That is the fundamental limitation of a Gaussian observation model on
# transformed counts, and the reason a binomial/Poisson model is the principled
# alternative. Inflating the SD is the safe direction to be wrong in.
build_ocfs_me <- function(all_composite,
                          unc_file   = "data/TULA20_age-depth_files/ocfs_uncertainty.csv",
                          chron_file = "data/TULA20_age-depth_files/TULA20_compsiteCore_w-ages.csv",
                          lambda = OCFS_LAMBDA, bin_width = BIN_WIDTH, verbose = TRUE) {

  u  <- read.csv(unc_file)
  ch <- read.csv(chron_file)

  # tulane.R bins ocfs on chron's median_age (renamed cov_age), joined by depth,
  # so do the same here rather than using the uncertainty file's own median_age.
  ch2 <- data.frame(depth_core = ch$depth_composite, cov_age = ch$median_age)
  if (!all(u$depth_core %in% ch2$depth_core))
    stop("build_ocfs_me: some spore depths are absent from the chronology")
  m <- merge(u, ch2, by = "depth_core", all.x = TRUE)

  n <- nrow(m)
  sd_s <- rep(NA_real_, n)
  src  <- rep(NA_character_, n)
  z975 <- qnorm(0.975)

  haveq <- !is.na(m$Q2.5) & !is.na(m$Q97.5)
  sd_s[haveq] <- (pw(m$Q97.5[haveq], lambda) - pw(m$Q2.5[haveq], lambda)) / (2 * z975)
  src[haveq]  <- "quantiles"

  zero <- !haveq & !is.na(m$OCFS_total) & m$OCFS_total == 0 &
          !is.na(m$Lycopodium) & m$Lycopodium > 0 & !is.na(m$Tracer_Added)
  if (any(zero)) {
    up_conc <- qgamma(0.975, shape = 1) / m$Lycopodium[zero] * m$Tracer_Added[zero]
    sd_s[zero] <- pw(up_conc, lambda) / z975
    src[zero]  <- "poisson_zero"
  }

  if (any(is.na(sd_s))) {
    src[is.na(sd_s)] <- "imputed"
    sd_s[is.na(sd_s)] <- median(sd_s, na.rm = TRUE)
  }
  if (any(!is.finite(sd_s)) || any(sd_s <= 0))
    stop("build_ocfs_me: non-finite or non-positive per-sample SD")

  # bin on the same grid as tulane.R
  br <- ocfs_bin_breaks(all_composite, bin_width)
  m$bins <- cut(m$cov_age, breaks = br, include.lowest = TRUE, labels = FALSE)
  m$var_s <- sd_s^2
  m$src <- src
  mm <- m[!is.na(m$bins), ]

  # Variance of the binned value. all_composite takes the MEAN of the samples in
  # a bin, so the variance of that mean is mean(var)/k under independence. Only
  # one bin has k > 1, so the approximation is nearly vacuous here, but it is
  # recorded honestly rather than ignored.
  agg <- aggregate(cbind(var_s = mm$var_s) ~ bins, data = mm, FUN = mean)
  k   <- aggregate(cbind(k = mm$var_s) ~ bins, data = mm, FUN = length)
  s   <- aggregate(cbind(src = mm$src) ~ bins, data = mm,
                   FUN = function(x) paste(sort(unique(x)), collapse = "+"))
  out <- merge(merge(agg, k, by = "bins"), s, by = "bins")
  out$var_bin <- out$var_s / out$k
  out$sd_bin  <- sqrt(out$var_bin)

  if (verbose) {
    cat("build_ocfs_me: lambda =", lambda, "\n")
    cat("  samples:", n, " binned:", nrow(mm), " bins:", nrow(out), "\n")
    cat("  per-sample SD source:", paste(names(table(src)), table(src),
                                         sep = " = ", collapse = ", "), "\n")
    cat("  per-bin SD  min/med/max:", paste(round(range(out$sd_bin), 3), collapse = " / "),
        "median", round(median(out$sd_bin), 3), "\n")
    cat("  variance ratio max/min:", round(max(out$var_bin) / min(out$var_bin), 1), "\n")
  }
  out[, c("bins", "k", "src", "var_bin", "sd_bin")]
}

# =============================================================================
# 3. Data assembly
# =============================================================================

# Builds X (response), U (predictors), ME (relative observation variance) and a
# lookup table, all trimmed to start at the first real OCFS observation.
#
# Trimming rather than inventing X[1] matters for two reasons: TVARSS's
# initial-updating block runs with no NA check so X[1] must be non-missing, and
# with a small su an invented first value is treated as near-exact data.
build_ocfs_data <- function(all_composite, me_tab,
                            predictors = c("char_acc", "d18O", "heinrich", "mean_co2", "PrDens"),
                            lambda = OCFS_LAMBDA, verbose = TRUE) {

  # all_composite as saved is already sorted bins 309 -> 1, i.e. oldest first,
  # so time runs forward down the rows. arrange(desc(bins)) is a no-op; assert it.
  if (!all(diff(all_composite$bins) < 0))
    stop("build_ocfs_data: expected all_composite sorted by descending bins")

  d <- all_composite
  d$ocfs_t <- pw(d$ocfs, lambda)

  first <- min(which(!is.na(d$ocfs_t)))
  d <- d[first:nrow(d), ]
  if (verbose)
    cat("build_ocfs_data: trimmed", first - 1, "leading bins with no OCFS;",
        nrow(d), "rows,", sum(!is.na(d$ocfs_t)), "observations\n")

  # predictors: interpolate the two with internal gaps, then standardise
  # everything except the 0/1 Heinrich indicator, exactly as before
  interp <- intersect(c("char_acc", "d18O", "mean_co2"), predictors)
  filled <- list()
  for (v in interp) {
    filled[[v]] <- d$bins[is.na(d[[v]])]
    d[[v]] <- as.numeric(forecast::na.interp(d[[v]]))
  }
  for (v in setdiff(predictors, "heinrich")) d[[v]] <- as.numeric(scale(d[[v]]))

  X <- matrix(d$ocfs_t, ncol = 1, dimnames = list(NULL, "ocfs"))
  U <- as.matrix(d[, predictors, drop = FALSE])
  if (anyNA(U)) stop("build_ocfs_data: U contains NA")

  # ME is a multiplier on su^2. Normalising to mean 1 lets su carry the overall
  # scale, so fitting su.fixed = NA and comparing the result with
  # sqrt(mean(var_bin)) is a calibration check on the stated uncertainties.
  v <- me_tab$var_bin[match(d$bins, me_tab$bins)]
  obs <- !is.na(X[, 1])
  if (any(is.na(v[obs]))) {
    miss <- d$bins[obs & is.na(v)]
    warning("build_ocfs_data: no measurement variance for bins ",
            paste(miss, collapse = ", "), "; using the median")
    v[obs & is.na(v)] <- median(v[obs], na.rm = TRUE)
  }
  v_mean <- mean(v[obs])
  v[!obs] <- v_mean          # never used by the filter, but must not be NA
  var_obs <- v               # counting variance per bin, in lambda-scale units
  ME <- v / v_mean           # normalised: use with su.fixed = sqrt(v_mean)
  if (anyNA(ME) || any(ME <= 0)) stop("build_ocfs_data: bad ME")
  if (length(ME) != length(X)) stop("build_ocfs_data: length(ME) != length(X)")

  if (verbose) {
    cat("  response range", paste(round(range(X, na.rm = TRUE), 3), collapse = " to "),
        " sd", round(sd(X, na.rm = TRUE), 3), "\n")
    cat("  ME range", paste(round(range(ME[obs]), 3), collapse = " to "),
        "; mean per-bin observation SD", round(sqrt(v_mean), 3), "\n")
  }

  list(X = X, U = U, ME = ME, var_obs = var_obs, bins = d$bins, age = d$age,
       data = d, n_obs = sum(obs), obs = obs, var_mean = v_mean, filled = filled,
       lambda = lambda)
}

# =============================================================================
# 4. Fitting
# =============================================================================

# One TVARSS fit with the predictors in c (the process equation) and d pinned at
# zero. This is TVARSS's own default parameterisation; the original script
# reversed it.
fit_ocfs <- function(X, U, ME, p = 2, su.fixed = 1, c.start = NULL,
                     b0.start = NA, b.start = array(NA, dim = p),
                     sb0.fixed = 0, sb.fixed = matrix(0, 1, p)) {
  args <- list(X = X, p = p, ME = ME, Tsamplefract = .9, show.fig = FALSE,
               annealing = FALSE, initial.points = "stationary",
               sb0.fixed = sb0.fixed, sb.fixed = sb.fixed, su.fixed = su.fixed,
               b0.start = b0.start, b.start = b.start)
  if (!is.null(U) && ncol(U) > 0) {
    U <- as.matrix(U)
    args <- c(args, list(U = U,
                         c.fixed = rep(NA, ncol(U)),
                         d.fixed = rep(0,  ncol(U)),
                         c.start = if (is.null(c.start)) rep(0.05, ncol(U)) else c.start))
  }
  do.call(TVARSS, args)
}

# Multi-start over c.start; returns the best fit by logLik.
#
# The grid includes a start derived from arima(): an arima xreg coefficient beta
# is a long-run effect, and the long-run effect of c is c/(1 - sum(b)), so
# c ~ beta * (1 - sum(b)). Note TVARSS.r has a commented-out c.init that DIVIDES
# by (1 - sum(b)); that goes the wrong way.
fit_ms <- function(X, U, ME, p = 2, su.fixed = 1,
                   starts = c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5),
                   b0.start = NA, b.start = array(NA, dim = p),
                   extra.starts = list(), arima.start = TRUE, quiet = TRUE,
                   se.min = 0, ...) {

  if (is.null(U) || ncol(U) == 0) {
    # no predictors: still multi-start over b, because the degenerate branch
    # exists here too
    cand <- list(try(fit_ocfs(X, NULL, ME, p = p, su.fixed = su.fixed,
                              b0.start = b0.start, b.start = b.start, ...),
                     silent = TRUE))
    if (!se_admissible(cand[[1]], se.min))
      cand <- c(cand, list(try(fit_ocfs(X, NULL, ME, p = p, su.fixed = su.fixed,
                                        b0.start = mean(X, na.rm = TRUE),
                                        b.start = rep(0.3, p), ...), silent = TRUE)))
    ok <- Filter(function(m) se_admissible(m, se.min), cand)
    if (length(ok)) return(ok[[which.max(vapply(ok, function(m) m$logLik, 0))]])
    ok <- Filter(function(m) se_admissible(m, 0), cand)
    if (!length(ok)) stop("fit_ms: null-model fit failed")
    warning("fit_ms: no fit met se.min for the null model; returning the best available")
    return(ok[[which.max(vapply(ok, function(m) m$logLik, 0))]])
  }

  U <- as.matrix(U)
  grid <- lapply(starts, function(s) rep(s, ncol(U)))

  if (arima.start) {
    a <- try(arima(X, xreg = U, order = c(p, 0, 0)), silent = TRUE)
    if (!inherits(a, "try-error")) {
      shrink <- 1 - sum(a$coef[1:p])
      beta <- a$coef[(p + 2):(p + 1 + ncol(U))]
      if (all(is.finite(beta)) && is.finite(shrink))
        grid <- c(grid, list(as.numeric(beta * shrink)))
    }
  }
  grid <- c(grid, extra.starts)

  fits <- lapply(grid, function(g)
    try(fit_ocfs(X, U, ME, p = p, su.fixed = su.fixed, c.start = g,
                 b0.start = b0.start, b.start = b.start, ...), silent = TRUE))

  # Choose the best ADMISSIBLE fit. Without this filter the optimiser sometimes
  # returns the se ~ 0 deterministic-trend solution, which beats every honest fit
  # on likelihood and then produces negative drop-one deviances downstream.
  ok <- Filter(function(m) se_admissible(m, se.min), fits)
  n_degen <- sum(vapply(fits, function(m) se_admissible(m, 0), TRUE)) - length(ok)

  if (!length(ok)) {
    any_ok <- Filter(function(m) se_admissible(m, 0), fits)
    if (!length(any_ok)) stop("fit_ms: every start failed")
    warning(sprintf(paste("fit_ms: all %d fits fell below se.min = %.4g (degenerate",
                          "zero-process-noise branch); returning the best available.",
                          "Treat its estimates as unreliable."), length(any_ok), se.min))
    return(any_ok[[which.max(vapply(any_ok, function(m) m$logLik, 0))]])
  }

  best <- ok[[which.max(vapply(ok, function(m) m$logLik, 0))]]
  if (!quiet)
    cat(sprintf("    fit_ms: %d starts, best admissible logLik = %.4f (se = %.4f)%s\n",
                length(grid), best$logLik, best$se,
                if (n_degen > 0) sprintf(", %d degenerate fit(s) rejected", n_degen) else ""))
  best
}

# Lift-back ascent. Given a full model and the list of one-term-dropped reduced
# models, insert a zero for each dropped coefficient, re-optimise the full model
# from there, and keep any improvement. Iterate until no reduced model beats the
# full one. This is what makes the drop-one deviances non-negative.
ascend <- function(full, reduced, X, U, ME, p = 2, su.fixed = 1,
                   max_pass = 4, quiet = FALSE, se.min = 0) {
  U <- as.matrix(U)
  for (pass in seq_len(max_pass)) {
    improved <- FALSE
    for (i in seq_along(reduced)) {
      if (is.null(reduced[[i]])) next
      cs <- numeric(ncol(U))
      cs[-i] <- as.vector(reduced[[i]]$c)
      cs[i]  <- 0
      m <- try(fit_ocfs(X, U, ME, p = p, su.fixed = su.fixed, c.start = cs,
                        b0.start = reduced[[i]]$b0, b.start = reduced[[i]]$b),
               silent = TRUE)
      if (se_admissible(m, se.min) && m$logLik > full$logLik) {
        if (!quiet) cat(sprintf("    ascend pass %d: full model %.4f -> %.4f (via %s)\n",
                                pass, full$logLik, m$logLik, colnames(U)[i]))
        full <- m
        improved <- TRUE
      }
    }
    if (!improved) break
  }
  full
}

# =============================================================================
# 5. Inference
# =============================================================================

# AICc on the number of OBSERVED time points, not the number of rows. TVARSS's
# own $AIC field uses a log-likelihood constant based on Tmax, so it is only
# comparable within one response definition - never across transforms or
# different trimmings.
aicc <- function(mod, n_obs) {
  k <- mod$npar
  if (n_obs - k - 1 <= 0) return(NA_real_)
  -2 * mod$logLik + 2 * k + 2 * k * (k + 1) / (n_obs - k - 1)
}

# Long-run effect of a c coefficient: a sustained one-SD change in the predictor
# eventually moves the latent series by c / (1 - sum(b)).
#
# TREAT WITH CARE on this dataset. sum(b) sits near 1 (0.85 with an SE of 0.06),
# so 1 - sum(b) is a small number with large relative uncertainty and the
# multiplier 1/(1 - sum(b)) ranges from about 4 to 23 across its own confidence
# interval. The c coefficients are far better determined than their long-run
# counterparts, so quote c as the primary estimate and use longrun_ci() - never a
# bare longrun() point estimate - when the accumulated effect is discussed.
longrun <- function(mod) as.vector(mod$c) / (1 - sum(mod$b))

# Delta-method SE and CI for the long-run effects, propagating the uncertainty in
# both c and b. With g = c_j / s and s = 1 - sum(b):
#   dg/dc_j = 1/s        dg/db_i = c_j / s^2
longrun_ci <- function(mod, V, level = 0.95) {
  if (is.null(V) || is.null(mod$c)) return(NULL)
  p <- mod$p
  s <- 1 - sum(mod$b)
  bn <- paste0("b", seq_len(p))
  cn <- rownames(V)[!(rownames(V) %in% c("b0", bn, "se"))]
  z <- qnorm(1 - (1 - level) / 2)
  cc <- as.vector(mod$c)
  out <- do.call(rbind, lapply(seq_along(cn), function(j) {
    g <- setNames(numeric(nrow(V)), rownames(V))
    g[cn[j]] <- 1 / s
    g[bn]    <- cc[j] / s^2
    v <- as.numeric(t(g) %*% V %*% g)
    data.frame(term = cn[j], c = cc[j], c_se = sqrt(V[cn[j], cn[j]]),
               longrun = cc[j] / s, longrun_se = sqrt(max(v, 0)))
  }))
  out$c_lo <- out$c - z * out$c_se;             out$c_hi <- out$c + z * out$c_se
  out$lr_lo <- out$longrun - z * out$longrun_se; out$lr_hi <- out$longrun + z * out$longrun_se
  se_s <- sqrt(sum(V[bn, bn]))
  attr(out, "shrinkage") <- c(s = s, se = se_s,
                              mult_lo = 1 / (s + z * se_s),
                              mult_hi = 1 / max(1e-4, s - z * se_s))
  out
}

# Dominant eigenvalue of the companion matrix - the honest persistence summary
# when the individual lag coefficients are weakly identified.
dominant_eigen <- function(mod) {
  p <- mod$p
  B <- matrix(0, p, p); B[1, ] <- mod$b
  if (p > 1) B[2:p, 1:(p - 1)] <- diag(p - 1)
  max(abs(eigen(B)$values))
}

# Drop-one LRTs with the ascent built in. Returns the table AND the (possibly
# improved) full model, because the ascent can change it.
drop_one <- function(full, X, U, ME, p = 2, su.fixed = 1, n_obs,
                     starts = c(0.001, 0.01, 0.05, 0.1, 0.25), quiet = FALSE,
                     se.min = 0) {
  U <- as.matrix(U)
  vars <- colnames(U)
  reduced <- vector("list", length(vars))

  for (i in seq_along(vars)) {
    reduced[[i]] <- fit_ms(X, U[, -i, drop = FALSE], ME, p = p, su.fixed = su.fixed,
                           starts = starts,
                           extra.starts = list(as.vector(full$c)[-i]),
                           b0.start = full$b0, b.start = full$b, se.min = se.min)
    if (!quiet) cat(sprintf("    drop %-9s reduced logLik = %10.4f  (se = %.4f)\n",
                            vars[i], reduced[[i]]$logLik, reduced[[i]]$se))
  }

  full <- ascend(full, reduced, X, U, ME, p = p, su.fixed = su.fixed,
                 quiet = quiet, se.min = se.min)

  dev <- 2 * (full$logLik - vapply(reduced, function(m) m$logLik, 0))
  assert_nonneg_dev(dev, vars)

  data.frame(var = vars,
             c = as.vector(full$c),
             longrun = longrun(full),
             reduced_logLik = vapply(reduced, function(m) m$logLik, 0),
             dev = dev,
             P = pchisq(pmax(dev, 0), df = 1, lower.tail = FALSE),
             row.names = NULL) -> tab
  list(table = tab, full = full, reduced = reduced)
}

# A negative deviance between nested models is arithmetically impossible at the
# maxima, so it means the optimiser failed. Stop rather than report a P-value.
# tol absorbs pure floating-point noise only.
assert_nonneg_dev <- function(dev, labels = NULL, tol = 1e-6) {
  bad <- which(dev < -tol)
  if (length(bad)) {
    lab <- if (is.null(labels)) bad else labels[bad]
    stop("negative likelihood-ratio deviance for: ",
         paste(sprintf("%s (%.4f)", lab, dev[bad]), collapse = ", "),
         "\n  The full model is not at its maximum. Widen the multi-start grid ",
         "or add ascent passes; do not report these tests.")
  }
  invisible(TRUE)
}

# =============================================================================
# 4b. The observation model: counting error plus a constant extra term
# =============================================================================

# ocfs_uncertainty.csv gives COUNTING error only - the binomial/Poisson noise in
# how many spores landed on the slide. Real observation error is larger, because
# a sample also carries dung-deposition patchiness, bioturbation, and age-model
# error, none of which counting statistics can see.
#
# Fitting su freely on top of ME = counting variance forces that extra error to
# be PROPORTIONAL to counting variance, which is the wrong shape: patchiness does
# not care how many Lycopodium grains were counted. On this dataset that
# misspecification is severe enough that the optimiser runs away to se -> 0,
# calling the entire series measurement noise (convergence code 10, se = 0.6% of
# sd(X)).
#
# The right shape is additive:
#
#     var_obs[t] = counting_var[t] + tau^2
#
# which is obtained by passing ME[t] = counting_var[t] + tau^2 with su.fixed = 1.
# tau is then one extra parameter, profiled out on a grid. tau = 0 recovers "the
# stated uncertainties are the whole story"; a large tau says counting error is a
# minor component.
me_with_tau <- function(var_obs, tau) var_obs + tau^2

# Profile the likelihood over tau. Returns the grid and the ML value. Profiling
# rather than joint optimisation because tau trades off against se, and a 1-D
# profile is far more reliable here than letting Nelder-Mead find it.
#
# Grid points are independent, so this parallelises perfectly. Set workers = 1
# for a serial run.
#
# Read the sum_b column as well as logLik. On this dataset se falls and sum(b)
# climbs towards 1 as tau rises: the likelihood has a ridge along which process
# noise is traded for observation noise, ending in a near-random-walk latent
# series with almost no innovation variance. That ridge is nearly flat, so tau is
# only weakly identified and the persistence estimate travels with it.
profile_tau <- function(X, U, var_obs, p = 2, se.min = 0,
                        taus = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2),
                        b0.start = NA, b.start = array(NA, dim = p),
                        starts = c(0.01, 0.05, 0.1, 0.25), verbose = TRUE,
                        workers = 1, root = getwd()) {

  one <- function(i) {
    ME <- me_with_tau(var_obs, taus[i])
    m <- try(fit_ms(X, U, ME, p = p, su.fixed = 1, starts = starts, se.min = se.min,
                    b0.start = b0.start, b.start = b.start), silent = TRUE)
    if (inherits(m, "try-error")) NULL else m
  }

  if (workers > 1 && requireNamespace("furrr", quietly = TRUE)) {
    Xw <- X; Uw <- U; vw <- var_obs; pw_ <- p; sw <- starts; smw <- se.min
    b0w <- b0.start; bw <- b.start; tw <- taus; rootw <- root
    future::plan(future::multisession, workers = workers)
    on.exit(future::plan(future::sequential), add = TRUE)
    fits <- furrr::future_map(seq_along(taus), function(i) {
      setwd(rootw)
      suppressMessages({
        library(forecast)
        source("kalman-smoother/TVARSS_12Aug26.r")
        source("kalman-smoother/ocfs_helpers.r")
      })
      ME <- me_with_tau(vw, tw[i])
      m <- try(fit_ms(Xw, Uw, ME, p = pw_, su.fixed = 1, starts = sw, se.min = smw,
                      b0.start = b0w, b.start = bw), silent = TRUE)
      if (inherits(m, "try-error")) NULL else m
    }, .options = furrr::furrr_options(seed = 1984))
    future::plan(future::sequential)
  } else {
    fits <- lapply(seq_along(taus), one)
  }

  res <- data.frame(
    tau    = taus,
    logLik = vapply(fits, function(m) if (is.null(m)) NA_real_ else m$logLik, 0),
    se     = vapply(fits, function(m) if (is.null(m)) NA_real_ else m$se, 0),
    sum_b  = vapply(fits, function(m) if (is.null(m)) NA_real_ else sum(m$b), 0)
  )
  if (verbose)
    for (i in seq_along(taus))
      cat(sprintf("  tau = %5.2f  logLik = %10.4f  se = %7.4f  sum(b) = %6.3f\n",
                  res$tau[i], res$logLik[i], res$se[i], res$sum_b[i]))

  best <- which.max(res$logLik)
  res$dlogLik <- res$logLik - res$logLik[best]
  # 95% profile-likelihood interval: within qchisq(.95, 1)/2 = 1.92 units
  inside <- res$tau[!is.na(res$dlogLik) & res$dlogLik > -qchisq(0.95, 1) / 2]
  list(profile = res, tau = res$tau[best], fit = fits[[best]],
       ci = if (length(inside)) range(inside) else c(NA, NA), fits = fits)
}

# =============================================================================
# 5b. Standalone likelihood and numerical standard errors
# =============================================================================

# TVARSS returns no Hessian, so model averaging has no SEs to work with. This is
# a standalone log-likelihood for the c-parameterisation (predictors in c, d = 0,
# sb0 = sb = 0), written independently of TVARSS_ml.
#
# It is exact rather than an approximation: with the b-block variances pinned at
# zero, the b/b0 part of the TVARSS state has zero variance for all time, so the
# B12/B13 blocks of the transition only ever multiply zeros and the (2p+1)-
# dimensional filter collapses to an ordinary p-dimensional linear Gaussian one
# in z = x - b0. The constant in the returned log-likelihood is TVARSS's own
# -((Tmax - p)/2) log(2 pi), so the two are directly comparable.
#
# ALWAYS check it against the fit before using it - check_loglik() does that.
#
# par is ordered as TVARSS's opt.par: (b0, b_1..b_p, se, c_1..c_nu).
ocfs_loglik <- function(par, X, U, ME, p, su) {
  X <- as.vector(X)
  b0 <- par[1]
  b  <- par[2:(p + 1)]
  se <- abs(par[p + 2])
  Tmax <- length(X)

  B <- matrix(0, p, p); B[1, ] <- b
  if (p > 1) B[2:p, 1:(p - 1)] <- diag(p - 1)
  Se <- matrix(0, p, p); Se[1, 1] <- se^2

  M <- diag(p * p) - kronecker(B, B)
  if (rcond(M) < 1e-12) return(-1e10)
  P <- matrix(solve(M) %*% as.vector(Se), p, p)

  if (is.null(U) || ncol(U) == 0) {
    uc <- rep(0, Tmax)
  } else {
    cc <- par[(p + 3):(p + 2 + ncol(U))]
    uc <- as.vector(as.matrix(U) %*% matrix(cc, ncol = 1))
  }

  a <- matrix(uc[1], p, 1)     # TVARSS's "stationary" initialisation
  sumlogF <- 0; sumvFv <- 0

  for (t in 1:Tmax) {
    if (t > 1) {
      a <- B %*% a
      a[1] <- a[1] + uc[t]
      P <- B %*% P %*% t(B) + Se
    }
    if (!is.na(X[t])) {
      F <- P[1, 1] + su^2 * ME[t]
      if (!is.finite(F) || F <= 0) return(-1e10)
      v <- X[t] - (a[1] + b0)
      K <- P[, 1, drop = FALSE] / F
      a <- a + K * v
      P <- P - K %*% P[1, , drop = FALSE]
      sumlogF <- sumlogF + log(F)
      sumvFv  <- sumvFv + v^2 / F
    }
  }
  -((Tmax - p) / 2) * log(2 * pi) - (sumlogF + sumvFv) / 2
}

check_loglik <- function(mod, X, U, ME, su, tol = 1e-6) {
  ll <- ocfs_loglik(mod$opt.par, X, U, ME, mod$p, su)
  ok <- abs(ll - mod$logLik) < tol
  if (!ok)
    stop(sprintf(paste("the standalone likelihood disagrees with TVARSS:",
                       "%.8f vs %.8f.\n  Do not use the numerical SEs until this",
                       "is resolved."), ll, mod$logLik))
  invisible(TRUE)
}

# Variance-covariance matrix from the numerical Hessian of the log-likelihood.
# Returns NULL (with a warning) rather than stopping if the Hessian is singular
# or gives negative variances, which happens when a coefficient sits against a
# flat direction of the likelihood.
param_vcov <- function(mod, X, U, ME, su, check = TRUE) {
  if (check) check_loglik(mod, X, U, ME, su)
  p <- mod$p
  nll <- function(par) -ocfs_loglik(par, X, U, ME, p, su)
  H <- try(numDeriv::hessian(nll, mod$opt.par), silent = TRUE)
  if (inherits(H, "try-error") || any(!is.finite(H))) {
    warning("param_vcov: could not evaluate the Hessian"); return(NULL)
  }
  V <- try(solve(H), silent = TRUE)
  if (inherits(V, "try-error")) {
    warning("param_vcov: singular Hessian"); return(NULL)
  }
  if (any(diag(V) < 0)) {
    warning("param_vcov: negative variance on the diagonal - the optimum is ",
            "not a proper interior maximum in every direction")
  }
  nm <- c("b0", paste0("b", seq_len(p)), "se",
          if (!is.null(U) && ncol(U) > 0) colnames(as.matrix(U)))
  dimnames(V) <- list(nm, nm)
  V
}

# SEs of the c coefficients only.
c_se <- function(mod, X, U, ME, su) {
  V <- param_vcov(mod, X, U, ME, su)
  if (is.null(U) || ncol(as.matrix(U)) == 0) return(numeric(0))
  nm <- colnames(as.matrix(U))
  if (is.null(V)) return(setNames(rep(NA_real_, length(nm)), nm))
  setNames(sqrt(pmax(diag(V)[nm], NA_real_)), nm)
}

# =============================================================================
# 6. Caching
# =============================================================================

# Fits take 30-60 s, so cache them. Delete results/cache/ to force a refit.
cached <- function(key, expr, dir = "results/cache", refresh = FALSE) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(dir, paste0(key, ".rds"))
  if (!refresh && file.exists(f)) return(readRDS(f))
  val <- force(expr)
  saveRDS(val, f)
  val
}
