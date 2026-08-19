## ocfs_analysis.r  --  Tier 2: the corrected single-model analysis
##
## Replaces kalman-smoother/ocfs_smoother.r, which is kept unchanged for
## reference. What is different, and why:
##
##   response      ocfs^0.25 instead of raw ocfs. Chosen by asking which scale
##                 makes the measurement error in ocfs_uncertainty.csv
##                 homoscedastic and symmetric; the optimum is 0.24-0.30.
##   observation   var_obs[t] = the measured per-bin counting variance from
##                 ocfs_uncertainty.csv, su fixed at 1, tau fixed at 0. The
##                 original used ME = 1 and su = 1, i.e. every sample equally
##                 and perfectly measured. Section 2 explains why tau is not
##                 estimated.
##   X[1]          removed. The series is trimmed to start at the first real
##                 observation instead of inventing a value of 100.
##   predictors    in c (process equation). Also fitted in d for comparison,
##                 because on this dataset the choice changes which predictors
##                 look important.
##   optimiser     multi-start plus lift-back ascent, with a hard assertion that
##                 no likelihood-ratio deviance comes out negative.
##   smoother      the corrected one, with an uncertainty band.
##
## Run:  Rscript kalman-smoother/ocfs_analysis.r
## Fits are cached in results/cache/; delete that folder to force a refit.

rm(list = ls())
suppressMessages({
  library(future); library(furrr)
  library(dplyr); library(readr); library(forecast)
})
source("kalman-smoother/TVARSS_12Aug26.r")
source("kalman-smoother/ocfs_helpers.r")

OUT <- "results/ocfs"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

p <- 2
PREDICTORS <- c("char_acc", "d18O", "heinrich", "mean_co2", "PrDens")
rule <- function(txt) cat("\n", strrep("=", 74), "\n", txt, "\n", strrep("=", 74), "\n", sep = "")

# =============================================================================
rule("1. Data")
# =============================================================================

all_composite <- read_csv("./data/all_composite.csv", show_col_types = FALSE) |> as.data.frame()
me  <- build_ocfs_me(all_composite, lambda = OCFS_LAMBDA)
dat <- build_ocfs_data(all_composite, me, predictors = PREDICTORS, lambda = OCFS_LAMBDA)

X <- dat$X; U <- dat$U; n_obs <- dat$n_obs

# Reject the degenerate zero-process-noise branch (see ocfs_helpers.r header).
SE_MIN <- SE_MIN_FRAC * sd(X, na.rm = TRUE)
cat(sprintf("se floor for admissible fits: %.4f (%.0f%% of sd(response) = %.3f)
",
            SE_MIN, 100 * SE_MIN_FRAC, sd(X, na.rm = TRUE)))
var_count <- dat$var_obs      # counting variance only

cat("\npredictor correlations:\n"); print(round(cor(U), 2))

# =============================================================================
rule("2. Observation model: how much error is there beyond counting error?")
# =============================================================================

# DECISION: tau is FIXED AT ZERO, so the observation variance is exactly the
# measured counting variance and su is not estimated at all.
#
# The reasoning, because this is a judgement call and not an obvious default:
#
#  * Estimating any observation-error parameter (su freely, or tau on a grid)
#    creates a trade-off against se, and that trade-off has a degenerate solution
#    at se = 0 in which the latent series becomes a deterministic AR(2) recursion
#    - flexible enough near the unit circle to trace almost any smooth curve, so
#    it wins on likelihood by 15-25 units. It is reached from only some starting
#    values, which made the tau profile non-reproducible (tau = 1.25 returned
#    logLik -347.87 on one run and -330.82 on the next) and produced a drop-one
#    deviance of -26.8 that halted an earlier version of this script.
#
#  * Fixing tau = 0 removes the trade-off rather than policing it. su is then
#    pinned entirely by external data - the whole point of having per-sample
#    quantiles - and se is the only variance parameter left. Empirically the model
#    is then well behaved: se ~ 0.31 of sd(response), sum(b) ~ 0.85, and every
#    drop-one deviance comes out positive.
#
#  * The cost is that tau = 0 understates observation error, since counting
#    statistics cannot see dung patchiness, bioturbation or age-model slop. But
#    that error is CONSERVATIVE for the question being asked: unmodelled
#    observation error is absorbed into se, making the latent series noisier and
#    predictor effects HARDER to detect. Allowing tau > 0 strengthens the omnibus
#    result (dev 16.3 at tau = 0 against 23.2 at tau = 1), so nothing is being
#    hidden by refusing to fit it.
#
# The profile is still computed and reported, because "AICc prefers extra
# observation error" is a real feature of these data and belongs in the write-up -
# it just should not be chased into an unidentified region. Section 7b then checks
# whether the predictor conclusions actually move along the profile.
TAU <- 0
SU  <- 1
ME  <- me_with_tau(var_count, TAU)
K_EXTRA <- 0          # no fitted observation-error parameter
aicc2 <- function(mod, n = n_obs) {
  k <- mod$npar + K_EXTRA
  if (n - k - 1 <= 0) return(NA_real_)
  -2 * mod$logLik + 2 * k + 2 * k * (k + 1) / (n - k - 1)
}
cat(sprintf("\ntau fixed at %.2f: observation variance IS the measured counting variance\n", TAU))
cat(sprintf("ME range %s (su fixed at 1)\n",
            paste(round(range(ME[dat$obs]), 3), collapse = " to ")))
cat(sprintf("mean observation SD %.3f vs sd(response) %.3f\n",
            sqrt(mean(var_count[dat$obs])), sd(X, na.rm = TRUE)))

# --- the profile, reported but not used to choose the model -------------------
WORKERS <- max(1, min(8, future::availableCores() - 2))
tp0 <- cached("tau_profile_null",
              profile_tau(X, NULL, var_count, p = p, se.min = SE_MIN, verbose = FALSE,
                          workers = WORKERS, root = getwd()))
tpF <- cached("tau_profile_full",
              profile_tau(X, U, var_count, p = p, se.min = SE_MIN, verbose = FALSE,
                          workers = WORKERS, root = getwd()))
prof <- tpF$profile
prof$k <- p + 2 + ncol(U) + (prof$tau > 0)
prof$AICc <- -2 * prof$logLik + 2 * prof$k +
             2 * prof$k * (prof$k + 1) / (n_obs - prof$k - 1)
prof$dAICc <- prof$AICc - min(prof$AICc, na.rm = TRUE)
prof$se_frac <- prof$se / sd(X, na.rm = TRUE)
prof$ok <- ifelse(prof$se >= SE_MIN, "yes", "NO-degenerate")
cat("\nfor reference only - profile over tau on the global model:\n")
print(prof[, c("tau", "logLik", "se", "se_frac", "sum_b", "AICc", "dAICc", "ok")],
      digits = 4, row.names = FALSE)
cat(sprintf("global-model ML tau = %.2f [%.2f, %.2f];  null-model ML tau = %.2f [%.2f, %.2f]\n",
            tpF$tau, tpF$ci[1], tpF$ci[2], tp0$tau, tp0$ci[1], tp0$ci[2]))
cat("Read the se_frac and sum_b columns: as tau rises, process noise vanishes and\n")
cat("persistence approaches a random walk. That is the unidentified region, and it\n")
cat("is why tau is not estimated. AICc preferring it is not evidence worth chasing.\n")

# =============================================================================
rule("3. Null and full models, multi-started")
# =============================================================================

mod0 <- cached("t2_mod0", fit_ms(X, NULL, ME, p = p, su.fixed = SU, se.min = SE_MIN))
cat(sprintf("null: logLik = %.4f  b0 = %.4f  b = (%s)  se = %.4f  sum(b) = %.4f\n",
            mod0$logLik, mod0$b0, paste(round(mod0$b, 3), collapse = ", "),
            mod0$se, sum(mod0$b)))
cat(sprintf("      dominant eigenvalue = %.4f\n", dominant_eigen(mod0)))

modF <- cached("t2_modF", fit_ms(X, U, ME, p = p, su.fixed = SU, se.min = SE_MIN,
                                b0.start = mod0$b0, b.start = mod0$b, quiet = FALSE))
cat("\nfull model, 5 predictors in c:\n"); print(modF)
check_loglik(modF, X, U, ME, SU)
cat("standalone likelihood agrees with TVARSS.\n")

# =============================================================================
rule("4. GATE: does the omnibus test survive at lambda = 0.25 with c?")
# =============================================================================

dev <- 2 * (modF$logLik - mod0$logLik)
assert_nonneg_dev(dev, "omnibus")
Pom <- pchisq(dev, df = ncol(U), lower.tail = FALSE)
cat(sprintf("omnibus LRT: dev = %.4f, df = %d, P = %.5g\n", dev, ncol(U), Pom))
cat(sprintf("AICc: null = %.2f  full = %.2f  difference = %+.2f   (n_obs = %d)\n",
            aicc2(mod0), aicc2(modF), aicc2(modF) - aicc2(mod0), n_obs))
cat(if (Pom <= 0.05) "\nGate passed.\n" else
    "\n*** GATE NOT PASSED - stop and re-examine before Tier 3.\n")

# =============================================================================
rule("5. c or d? The choice changes the answer, so fit both")
# =============================================================================
# c pushes the latent process, so its effect accumulates through the
# autoregression; d is an instantaneous offset on the observation. They are not
# jointly identifiable with 83 points, so they are fitted separately and compared
# by AICc. On this dataset the ranking of predictors differs between them, which
# is a substantive finding rather than a nuisance: a pulse-like covariate fits a
# process shock, a slow ramp fits a level offset.

fit_d <- function(X, U, ME, p, su.fixed, d.start, b0.start, b.start) {
  TVARSS(X = X, p = p, ME = ME, U = as.matrix(U), Tsamplefract = .9,
         show.fig = FALSE, annealing = FALSE, initial.points = "stationary",
         sb0.fixed = 0, sb.fixed = matrix(0, 1, p), su.fixed = su.fixed,
         c.fixed = rep(0, ncol(U)), d.fixed = rep(NA, ncol(U)),
         d.start = d.start, b0.start = b0.start, b.start = b.start)
}
modD <- cached("t2_modD", {
  best <- NULL
  for (s in c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5)) {
    m <- try(fit_d(X, U, ME, p, SU, rep(s, ncol(U)), mod0$b0, mod0$b), silent = TRUE)
    if (!se_admissible(m, SE_MIN)) next
    if (!inherits(m, "try-error") && (is.null(best) || m$logLik > best$logLik)) best <- m
  }
  best
})
cmp <- data.frame(
  parameterisation = c("null", "predictors in c", "predictors in d"),
  k      = c(mod0$npar, modF$npar, modD$npar) + K_EXTRA,
  logLik = c(mod0$logLik, modF$logLik, modD$logLik),
  AICc   = c(aicc2(mod0), aicc2(modF), aicc2(modD))
)
cmp$dAICc <- cmp$AICc - min(cmp$AICc)
print(cmp, digits = 5, row.names = FALSE)
cat("\ncoefficients side by side (not directly comparable in magnitude: c is a\n")
cat("per-step process input, d is a level shift):\n")
print(data.frame(term = colnames(U), c = as.vector(modF$c),
                 c_longrun = longrun(modF), d = as.vector(modD$d)),
      digits = 4, row.names = FALSE)
BEST_PAR <- if (aicc2(modF) <= aicc2(modD)) "c" else "d"
cat("\nbetter supported here:", BEST_PAR, "\n")

# =============================================================================
rule("6. Drop-one tests, with ascent")
# =============================================================================

d1 <- cached("t2_drop_one", drop_one(modF, X, U, ME, p = p, su.fixed = SU, n_obs = n_obs, se.min = SE_MIN))
modF <- d1$full
cat("\n")
print(d1$table[order(d1$table$P), ], digits = 4, row.names = FALSE)
cat(sprintf("\nfull model after ascent: logLik = %.4f  AICc = %.2f\n",
            modF$logLik, aicc2(modF)))
write.csv(d1$table, file.path(OUT, "drop_one_lambda025_c.csv"), row.names = FALSE)

k <- which(colnames(U) %in% c("d18O", "mean_co2"))
mj <- cached("t2_drop_clim", fit_ms(X, U[, -k, drop = FALSE], ME, p = p, su.fixed = SU, se.min = SE_MIN,
                                   b0.start = modF$b0, b.start = modF$b))
devj <- 2 * (modF$logLik - mj$logLik)
assert_nonneg_dev(devj, "d18O+mean_co2")
cat(sprintf("joint drop of d18O + mean_co2 (r = %.2f): dev = %.4f, df = 2, P = %.4g\n",
            cor(U[, "d18O"], U[, "mean_co2"]), devj,
            pchisq(devj, df = 2, lower.tail = FALSE)))

# =============================================================================
rule("7. Standard errors, and how to talk about effect sizes here")
# =============================================================================

V <- param_vcov(modF, X, U, ME, SU)
lr <- longrun_ci(modF, V)
sh <- attr(lr, "shrinkage")
cat(sprintf("sum(b) = %.4f, so 1 - sum(b) = %.4f (SE %.4f)\n",
            sum(modF$b), sh[["s"]], sh[["se"]]))
cat(sprintf("long-run multiplier 1/(1 - sum(b)) = %.2f, 95%% range %.1f to %.1f\n",
            1 / sh[["s"]], sh[["mult_lo"]], sh[["mult_hi"]]))
cat(sprintf("individual b SEs: %s, but SE of their SUM is only %.4f - the lag\n",
            paste(sprintf("%.3f", sqrt(diag(V)[paste0("b", seq_len(p))])), collapse = ", "),
            sh[["se"]]))
cat("coefficients trade off almost perfectly, so only the persistence is\n")
cat("identified, never the individual lags.\n")

NEAR_RW <- sum(modF$b) > 0.95
if (NEAR_RW) {
  cat("\nsum(b) is within 0.05 of 1, so the latent series is effectively a random\n")
  cat("walk. A random walk has no stationary mean, so 'the long-run effect of a\n")
  cat("sustained change' is not a defined quantity - the series never saturates.\n")
  cat("Report c as the effect per 200-yr bin, which IS estimable, and say that it\n")
  cat("accumulates for as long as the covariate persists. Do not divide by\n")
  cat("1 - sum(b); the long-run columns below are printed only to show how\n")
  cat("unstable that division is.\n")
} else {
  cat("\nsum(b) is far enough from 1 that long-run effects are meaningful, but they\n")
  cat("still carry the multiplier's uncertainty - use lr_lo/lr_hi, never a bare\n")
  cat("point estimate.\n")
}
cat("\n")
print(lr[, c("term", "c", "c_se", "c_lo", "c_hi", "longrun", "lr_lo", "lr_hi")],
      digits = 3, row.names = FALSE)
write.csv(lr, file.path(OUT, "coefficients_with_ci.csv"), row.names = FALSE)

# =============================================================================
rule("7b. Do the predictor conclusions survive the tau ridge?")
# =============================================================================
# tau is only weakly identified, and persistence moves with it. That only matters
# if the predictor inference moves too. Refit the drop-one tests at the ends and
# middle of the profile interval and compare. If the ranking holds, the weak
# identification of tau is a caveat about persistence, not about the predictors.

tau_check <- unique(c(0, 0.5, 1.0))
tau_check <- sort(tau_check[tau_check >= 0])
sens <- list()
for (tt in tau_check) {
  MEt <- me_with_tau(var_count, tt)
  key <- sprintf("t2_dropone_tau%03d", round(tt * 100))
  m0t <- cached(paste0("t2_null_tau", round(tt * 100)),
                fit_ms(X, NULL, MEt, p = p, su.fixed = 1, se.min = SE_MIN))
  mFt <- cached(paste0("t2_full_tau", round(tt * 100)),
                fit_ms(X, U, MEt, p = p, su.fixed = 1, se.min = SE_MIN,
                       b0.start = m0t$b0, b.start = m0t$b))
  dt  <- cached(key, drop_one(mFt, X, U, MEt, p = p, su.fixed = 1, se.min = SE_MIN,
                              n_obs = n_obs, quiet = TRUE))
  dv  <- 2 * (dt$full$logLik - m0t$logLik)
  sens[[sprintf("%.2f", tt)]] <- list(tab = dt$table, full = dt$full,
                                      omnibus = dv, sum_b = sum(dt$full$b),
                                      se = dt$full$se)
  cat(sprintf("\ntau = %.2f : omnibus dev = %.3f (P = %.4g), sum(b) = %.3f, se = %.4f\n",
              tt, dv, pchisq(max(dv, 0), ncol(U), lower.tail = FALSE),
              sum(dt$full$b), dt$full$se))
  print(dt$table[order(dt$table$P), c("var", "c", "dev", "P")],
        digits = 3, row.names = FALSE)
}

cat("\n---- P-values side by side across tau ----\n")
pmat <- sapply(sens, function(s) setNames(s$tab$P, s$tab$var))
print(round(pmat, 4))
cat("\nrank order of predictors at each tau (1 = strongest):\n")
print(apply(pmat, 2, rank))
stable <- all(apply(apply(pmat, 2, rank), 1, function(r) diff(range(r)) <= 1))
cat(if (stable)
      "\nRanking is stable to within one place across the tau interval: the weak\nidentification of tau is a caveat about persistence, not about the predictors.\n"
    else
      "\nWARNING: the predictor ranking MOVES along the tau ridge. Report the range\nacross tau rather than a single set of P-values.\n")
saveRDS(sens, file.path(OUT, "tau_sensitivity.rds"))

# =============================================================================
rule("8. Kalman smoother with uncertainty band")
# =============================================================================

fit <- TVARSS_KalmanSmoother(modF)
if (anyNA(fit$X.smoothed)) stop("smoother returned NA - is TVARSS_12Aug26.r sourced?")
if (is.null(fit$X.smoothed.se)) stop("no X.smoothed.se - old smoother in use")

sm <- data.frame(bins = dat$bins, age = dat$age,
                 ocfs = all_composite$ocfs[match(dat$bins, all_composite$bins)],
                 X = as.vector(fit$X), smooth = as.vector(fit$X.smoothed),
                 se = as.vector(fit$X.smoothed.se))
sm$lower <- sm$smooth - 1.96 * sm$se
sm$upper <- sm$smooth + 1.96 * sm$se
sm$ocfs_smooth <- ipw(sm$smooth); sm$ocfs_lower <- ipw(sm$lower); sm$ocfs_upper <- ipw(sm$upper)
if (any(sm$ocfs_smooth < 0, na.rm = TRUE)) stop("negative back-transformed concentration")

gap <- is.na(sm$X)
cat(sprintf("%d bins, all non-NA. mean SE: %.3f at observations, %.3f in gaps (sd = %.3f)\n",
            nrow(sm), mean(sm$se[!gap]), mean(sm$se[gap]), sd(X, na.rm = TRUE)))
cat(sprintf("variance resolved in the gaps: %.0f%%\n",
            100 * (1 - (mean(sm$se[gap]) / sd(X, na.rm = TRUE))^2)))
write.csv(sm, file.path(OUT, "smoothed_lambda025.csv"), row.names = FALSE)

# =============================================================================
rule("9. Figures")
# =============================================================================

png(file.path(OUT, "ocfs_smoothed.png"), width = 1600, height = 1000, res = 140)
op <- par(no.readonly = TRUE); par(mfrow = c(2, 1), mar = c(4, 4.6, 2.2, 1))
ageka <- approx(which(!is.na(sm$age)), sm$age[!is.na(sm$age)],
                xout = seq_len(nrow(sm)), rule = 2)$y / 1000
plot(ageka, sm$smooth, type = "n", xlim = rev(range(ageka)),
     ylim = range(c(sm$lower, sm$upper)), xlab = "Age (ka BP)",
     ylab = expression(ocfs^0.25), main = "Smoothed OCFS, analysis scale, 95% band")
polygon(c(ageka, rev(ageka)), c(sm$lower, rev(sm$upper)), col = "#1f5f8b22", border = NA)
lines(ageka, sm$smooth, col = "#1f5f8b", lwd = 2)
points(ageka, sm$X, pch = 21, bg = "white", cex = .6)
plot(ageka, sm$ocfs_smooth, type = "n", xlim = rev(range(ageka)),
     ylim = c(0, max(sm$ocfs_upper)), xlab = "Age (ka BP)",
     ylab = "OCFS concentration", main = "Back-transformed to concentration units")
polygon(c(ageka, rev(ageka)), c(sm$ocfs_lower, rev(sm$ocfs_upper)), col = "#18776622", border = NA)
lines(ageka, sm$ocfs_smooth, col = "#187766", lwd = 2)
points(ageka, sm$ocfs, pch = 21, bg = "white", cex = .6)
par(op); dev.off()

png(file.path(OUT, "tau_profile.png"), width = 1100, height = 750, res = 140)
par(mar = c(4.2, 4.4, 2.2, 1))
plot(tpF$profile$tau, tpF$profile$logLik, type = "b", pch = 19, col = "#1f5f8b",
     xlab = expression(tau), ylab = "log-likelihood",
     main = "Extra observation error beyond counting statistics")
abline(h = max(tpF$profile$logLik) - qchisq(0.95, 1) / 2, lty = 2, col = "grey50")
abline(v = TAU, lty = 3, col = "#a8322d")
legend("bottomright", bty = "n", cex = .85, lty = c(2, 3), col = c("grey50", "#a8322d"),
       legend = c("95% profile-likelihood cutoff", expression(paste("ML ", tau))))
dev.off()
cat("figures written to", OUT, "\n")

saveRDS(list(dat = dat, me = me, ME = ME, SU = SU, TAU = TAU, K_EXTRA = K_EXTRA,
             tau_profile = tpF, tau_profile_null = tp0, mod0 = mod0, modF = modF, modD = modD,
             drop_one = d1, vcov = V, longrun = lr, smoothed = sm,
             n_obs = n_obs, best_par = BEST_PAR),
        file.path(OUT, "tier2_fits.rds"))

rule("Done")
cat("outputs in", OUT, "\n")
