## ocfs_analysis.r  --  Tier 2: the corrected single-model analysis
##
## Replaces kalman-smoother/ocfs_smoother.r, which is kept unchanged for
## reference. What is different, and why:
##
##   response      ocfs^0.25 instead of raw ocfs. Chosen by asking which scale
##                 makes the measurement error in ocfs_uncertainty.csv
##                 homoscedastic and symmetric; the optimum is 0.24-0.30.
##   observation   var_obs[t] = counting_var[t] + tau^2, with counting_var from
##                 the stored quantiles and tau profiled. The original used
##                 ME = 1 and su = 1, i.e. every sample equally and perfectly
##                 measured.
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
var_count <- dat$var_obs      # counting variance only

cat("\npredictor correlations:\n"); print(round(cor(U), 2))

# =============================================================================
rule("2. Observation model: how much error is there beyond counting error?")
# =============================================================================
# ocfs_uncertainty.csv describes counting statistics only. Dung-deposition
# patchiness, bioturbation and age-model error are all invisible to it, so the
# real observation variance is counting_var[t] + tau^2 for some tau > 0. tau is
# profiled rather than jointly optimised, because it trades off against se and a
# 1-D profile is much more reliable than letting Nelder-Mead find it.
#
# Fitting su freely on top of ME = counting_var instead (i.e. assuming the extra
# error is PROPORTIONAL to counting error) fails badly here: se collapses to
# 0.6% of sd(X) and optim returns convergence code 10. That misspecification is
# what motivated the additive form.

tp <- cached("tau_profile_null", profile_tau(X, NULL, var_count, p = p, verbose = TRUE))
cat("\nprofile over tau, fitted on the null model:\n")
print(tp$profile, digits = 5, row.names = FALSE)
TAU <- tp$tau
cat(sprintf("\nML tau = %.2f, 95%% profile-likelihood CI [%.2f, %.2f]\n",
            TAU, tp$ci[1], tp$ci[2]))
cat(sprintf("counting-only observation SD averages %.3f; total is now %.3f\n",
            sqrt(mean(var_count[dat$obs])), sqrt(mean(var_count[dat$obs] + TAU^2))))
cat(sprintf("counting error is %.0f%% of the total observation variance\n",
            100 * mean(var_count[dat$obs]) / mean(var_count[dat$obs] + TAU^2)))
if (TAU == 0) cat("tau = 0: the stated uncertainties account for everything.\n")
if (!is.na(tp$ci[1]) && tp$ci[1] > 0)
  cat("tau is bounded away from 0, so the stated uncertainties are understated.\n")

# From here on: observation variance is exactly counting_var[t] + tau^2.
ME <- me_with_tau(var_count, TAU)
SU <- 1
cat(sprintf("\nME range %s (su fixed at 1, so ME IS the observation variance)\n",
            paste(round(range(ME[dat$obs]), 3), collapse = " to ")))
# tau costs one parameter that TVARSS does not know about; add it to every AICc
K_EXTRA <- if (TAU > 0) 1 else 0
aicc2 <- function(mod) {
  k <- mod$npar + K_EXTRA
  -2 * mod$logLik + 2 * k + 2 * k * (k + 1) / (n_obs - k - 1)
}

# =============================================================================
rule("3. Null and full models, multi-started")
# =============================================================================

mod0 <- cached("t2_mod0", fit_ms(X, NULL, ME, p = p, su.fixed = SU))
cat(sprintf("null: logLik = %.4f  b0 = %.4f  b = (%s)  se = %.4f  sum(b) = %.4f\n",
            mod0$logLik, mod0$b0, paste(round(mod0$b, 3), collapse = ", "),
            mod0$se, sum(mod0$b)))
cat(sprintf("      dominant eigenvalue = %.4f\n", dominant_eigen(mod0)))

modF <- cached("t2_modF", fit_ms(X, U, ME, p = p, su.fixed = SU,
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

d1 <- cached("t2_drop_one", drop_one(modF, X, U, ME, p = p, su.fixed = SU, n_obs = n_obs))
modF <- d1$full
cat("\n")
print(d1$table[order(d1$table$P), ], digits = 4, row.names = FALSE)
cat(sprintf("\nfull model after ascent: logLik = %.4f  AICc = %.2f\n",
            modF$logLik, aicc2(modF)))
write.csv(d1$table, file.path(OUT, "drop_one_lambda025_c.csv"), row.names = FALSE)

k <- which(colnames(U) %in% c("d18O", "mean_co2"))
mj <- cached("t2_drop_clim", fit_ms(X, U[, -k, drop = FALSE], ME, p = p, su.fixed = SU,
                                   b0.start = modF$b0, b.start = modF$b))
devj <- 2 * (modF$logLik - mj$logLik)
assert_nonneg_dev(devj, "d18O+mean_co2")
cat(sprintf("joint drop of d18O + mean_co2 (r = %.2f): dev = %.4f, df = 2, P = %.4g\n",
            cor(U[, "d18O"], U[, "mean_co2"]), devj,
            pchisq(devj, df = 2, lower.tail = FALSE)))

# =============================================================================
rule("7. Standard errors, and why the long-run effects are the weak part")
# =============================================================================

V <- param_vcov(modF, X, U, ME, SU)
lr <- longrun_ci(modF, V)
sh <- attr(lr, "shrinkage")
cat(sprintf("sum(b) = %.4f, so 1 - sum(b) = %.4f (SE %.4f)\n",
            sum(modF$b), sh[["s"]], sh[["se"]]))
cat(sprintf("the long-run multiplier 1/(1 - sum(b)) is %.2f, but its 95%% range is %.1f to %.1f\n",
            1 / sh[["s"]], sh[["mult_lo"]], sh[["mult_hi"]]))
cat("Individual b coefficients are almost unidentified (they trade off against\n")
cat("each other); only their sum is well determined. Quote c, not longrun.\n\n")
print(lr[, c("term", "c", "c_se", "c_lo", "c_hi", "longrun", "lr_lo", "lr_hi")],
      digits = 3, row.names = FALSE)
write.csv(lr, file.path(OUT, "coefficients_with_ci.csv"), row.names = FALSE)

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
plot(tp$profile$tau, tp$profile$logLik, type = "b", pch = 19, col = "#1f5f8b",
     xlab = expression(tau), ylab = "log-likelihood",
     main = "Extra observation error beyond counting statistics")
abline(h = max(tp$profile$logLik) - qchisq(0.95, 1) / 2, lty = 2, col = "grey50")
abline(v = TAU, lty = 3, col = "#a8322d")
legend("bottomright", bty = "n", cex = .85, lty = c(2, 3), col = c("grey50", "#a8322d"),
       legend = c("95% profile-likelihood cutoff", expression(paste("ML ", tau))))
dev.off()
cat("figures written to", OUT, "\n")

saveRDS(list(dat = dat, me = me, ME = ME, SU = SU, TAU = TAU, K_EXTRA = K_EXTRA,
             tau_profile = tp, mod0 = mod0, modF = modF, modD = modD,
             drop_one = d1, vcov = V, longrun = lr, smoothed = sm,
             n_obs = n_obs, best_par = BEST_PAR),
        file.path(OUT, "tier2_fits.rds"))

rule("Done")
cat("outputs in", OUT, "\n")
