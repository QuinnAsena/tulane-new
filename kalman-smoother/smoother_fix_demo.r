## smoother_fix_demo.r
##
## Self-contained demonstration that TVARSS_KalmanSmoother() in TVARSS_11Feb25.r
## returns incorrect smoothed states when the time series has missing values, and
## that the version in TVARSS_12Aug26.r fixes it.
##
## Run with:  Rscript kalman-smoother/smoother_fix_demo.r
## from the repository root. Takes about a minute. Prints a table and either
## "ALL CHECKS PASSED" or a failed assertion.
##
## The argument is made three ways so it does not rest on trusting any one
## implementation:
##
##   (a) An independent Kalman filter + Durbin-Koopman smoother is written from
##       scratch below (dk_smoother). It is legitimate as a reference because
##       when sb0 = sb = 0 the TVARSS state space reduces to an ordinary linear
##       Gaussian model in z = x - b0. We first prove it agrees with TVARSS's own
##       log-likelihood to machine precision, which establishes that it is
##       filtering the same model with the same parameters.
##
##   (b) Because the data are simulated, the true latent series is known, so both
##       smoothers can be scored against the truth rather than against each other.
##
##   (c) A smoother must beat the raw observations. That is the whole point of
##       running one, and it is a criterion no implementation can argue with.

rm(list = ls())
set.seed(1)

old_file <- "kalman-smoother/TVARSS_11Feb25.r"
new_file <- "kalman-smoother/TVARSS_12Aug26.r"
if (!file.exists(old_file)) stop("run this from the repository root")

old <- new.env(); new <- new.env()
suppressMessages({
  sys.source(old_file, envir = old)
  sys.source(new_file, envir = new)
})

# ---------------------------------------------------------------------------
# (a) Independent reference: Kalman filter + Durbin-Koopman smoother.
#     Valid only for sb0 = sb = 0, which is what this demo fits.
# ---------------------------------------------------------------------------
dk_smoother <- function(X, U, p, b0, b, se, su, dd, ME) {
  Tmax <- length(X)
  B <- matrix(0, p, p); B[1, ] <- b
  if (p > 1) B[2:p, 1:(p - 1)] <- diag(p - 1)
  Se <- matrix(0, p, p); Se[1, 1] <- se^2
  Z <- matrix(0, 1, p); Z[1, 1] <- 1
  Su <- su^2

  # stationary prior for z = x - b0
  P <- matrix(solve(diag(p * p) - kronecker(B, B)) %*% as.vector(Se), p, p)
  a <- matrix(0, p, 1)
  off <- if (is.null(dd)) rep(0, Tmax) else as.vector(U %*% matrix(dd, ncol = 1))

  apred <- matrix(NA, Tmax, p); Ppred <- array(NA, c(Tmax, p, p))
  vv <- numeric(Tmax); iF <- numeric(Tmax); L <- array(0, c(Tmax, p, p))
  sumlogF <- 0; sumvFv <- 0; nobs <- 0

  for (t in 1:Tmax) {
    apred[t, ] <- a; Ppred[t, , ] <- P
    if (!is.na(X[t])) {
      F <- as.numeric(Z %*% P %*% t(Z) + Su * ME[t])
      iF[t] <- 1 / F
      vv[t] <- X[t] - as.numeric(Z %*% a + b0 + off[t])
      K <- P %*% t(Z) / F
      au <- a + K * vv[t]
      Pu <- P - K %*% Z %*% P
      L[t, , ] <- B %*% (diag(p) - K %*% Z)
      sumlogF <- sumlogF + log(F); sumvFv <- sumvFv + vv[t]^2 / F; nobs <- nobs + 1
    } else {
      au <- a; Pu <- P; L[t, , ] <- B      # no update: gain is zero, L is the transition
    }
    a <- B %*% au
    P <- B %*% Pu %*% t(B) + Se
  }

  r <- matrix(0, p, 1); N <- matrix(0, p, p)
  ahat <- matrix(NA, Tmax, p); V <- numeric(Tmax)
  for (t in Tmax:1) {                       # every t, observed or not
    r <- t(Z) * iF[t] * vv[t] + t(L[t, , ]) %*% r
    N <- t(Z) %*% (iF[t] * Z) + t(L[t, , ]) %*% N %*% L[t, , ]
    ahat[t, ] <- apred[t, ] + Ppred[t, , ] %*% r
    V[t] <- (Ppred[t, , ] - Ppred[t, , ] %*% N %*% Ppred[t, , ])[1, 1]
  }

  list(smoothed = ahat[, 1] + b0 + off,
       se = sqrt(pmax(V, 0)),
       LL = sumlogF + sumvFv,
       nobs = nobs)
}

# TVARSS reports logLik = -((Tmax - p)/2) * log(2*pi) - LL/2, so invert that to
# recover the LL = sum(log F) + sum(v'F^-1 v) that the reference computes.
tvarss_LL <- function(mod) -2 * (mod$logLik + ((nrow(mod$X) - mod$p) / 2) * log(2 * pi))

# ---------------------------------------------------------------------------
# (b) Simulate an AR(2) with observation error. Truth is known.
# ---------------------------------------------------------------------------
Tmax <- 200; p <- 2
b_true <- c(0.6, 0.2); b0_true <- 10; se_true <- 1; su_true <- 1

burn <- 50
z <- numeric(Tmax + burn)
for (t in 3:(Tmax + burn)) {
  z[t] <- b_true[1] * z[t - 1] + b_true[2] * z[t - 2] + rnorm(1, sd = se_true)
}
z <- z[(burn + 1):(Tmax + burn)]
truth <- z + b0_true
Xfull <- matrix(truth + rnorm(Tmax, sd = su_true), ncol = 1)
ME <- rep(1, Tmax)

fit_it <- function(X) {
  new$TVARSS(X = X, p = p, ME = ME, Tsamplefract = .9, show.fig = FALSE,
             annealing = FALSE, sb0.fixed = 0, sb.fixed = matrix(0, 1, p),
             su.fixed = su_true)
}
rmse <- function(a, b) sqrt(mean((a - b)^2, na.rm = TRUE))

results <- list()
checks  <- character(0)
pass    <- TRUE
check <- function(label, ok) {
  pass <<- pass && isTRUE(ok)
  checks <<- c(checks, sprintf("  [%s] %s", if (isTRUE(ok)) "ok" else "FAIL", label))
}

# ---- Test A: fully observed -------------------------------------------------
cat("Test A: fully observed series (T = 200)\n")
fitA <- fit_it(Xfull)
cat(sprintf("  fitted b0=%.3f b=(%.3f, %.3f) se=%.3f su=%.3f\n",
            fitA$b0, fitA$b[1], fitA$b[2], fitA$se, fitA$su))

refA <- dk_smoother(as.vector(Xfull), NULL, p, fitA$b0, fitA$b, fitA$se, fitA$su, NULL, ME)
check("reference filter reproduces TVARSS log-likelihood",
      abs(refA$LL - tvarss_LL(fitA)) < 1e-6)

smA_old <- old$TVARSS_KalmanSmoother(fitA)
smA_new <- new$TVARSS_KalmanSmoother(fitA)

results$A <- c(data = rmse(as.vector(Xfull), truth),
               old  = rmse(as.vector(smA_old$X.smoothed), truth),
               new  = rmse(as.vector(smA_new$X.smoothed), truth),
               ref  = rmse(refA$smoothed, truth))

check("fixed smoother matches the independent reference (complete data)",
      max(abs(as.vector(smA_new$X.smoothed) - refA$smoothed)) < 1e-8)
check("fixed smoother beats the raw data (complete)",
      results$A[["new"]] < results$A[["data"]])
check("fixed smoother beats the old smoother (complete)",
      results$A[["new"]] < results$A[["old"]])

# ---- Test B: 70% of observations blanked ------------------------------------
cat("Test B: same series with 70% of observations set to NA\n")
Xmiss <- Xfull
drop <- sort(sample(2:Tmax, floor(0.7 * Tmax)))
Xmiss[drop] <- NA
keep <- setdiff(seq_len(Tmax), drop)

fitB <- fit_it(Xmiss)
cat(sprintf("  n observed = %d of %d;  fitted b0=%.3f b=(%.3f, %.3f) se=%.3f\n",
            sum(!is.na(Xmiss)), Tmax, fitB$b0, fitB$b[1], fitB$b[2], fitB$se))

refB <- dk_smoother(as.vector(Xmiss), NULL, p, fitB$b0, fitB$b, fitB$se, fitB$su, NULL, ME)
check("reference filter reproduces TVARSS log-likelihood (gappy)",
      abs(refB$LL - tvarss_LL(fitB)) < 1e-6)

smB_old <- old$TVARSS_KalmanSmoother(fitB)
smB_new <- new$TVARSS_KalmanSmoother(fitB)

results$B_obs <- c(data = rmse(as.vector(Xfull)[keep], truth[keep]),
                   old  = rmse(as.vector(smB_old$X.smoothed)[keep], truth[keep]),
                   new  = rmse(as.vector(smB_new$X.smoothed)[keep], truth[keep]),
                   ref  = rmse(refB$smoothed[keep], truth[keep]))
results$B_all <- c(data = NA,
                   old  = NA,
                   new  = rmse(as.vector(smB_new$X.smoothed), truth),
                   ref  = rmse(refB$smoothed, truth))

check("old smoother returns nothing in the gaps",
      sum(!is.na(smB_old$X.smoothed)) == sum(!is.na(Xmiss)))
check("fixed smoother returns a value at every time point",
      all(!is.na(smB_new$X.smoothed)))
check("fixed smoother matches the independent reference (gappy)",
      max(abs(as.vector(smB_new$X.smoothed) - refB$smoothed)) < 1e-8)
check("fixed smoother beats the raw data at observed times (gappy)",
      results$B_obs[["new"]] < results$B_obs[["data"]])

# The old smoother is not merely worse, it is close to useless here: compare how
# much of the achievable RMSE reduction each one actually recovers. Note this must
# be measured at observed times only, since the old smoother returns nothing else.
recovered <- (results$B_obs[["data"]] - results$B_obs[["old"]]) /
             (results$B_obs[["data"]] - results$B_obs[["new"]])
cat(sprintf("  of the RMSE reduction a correct smoother achieves, 11Feb25 recovers %.1f%%\n",
            100 * recovered))
check("old smoother recovers less than 10% of the achievable improvement",
      recovered < 0.10)
check("smoothed SE is returned and is larger in gaps than at observations",
      mean(smB_new$X.smoothed.se[drop]) > mean(smB_new$X.smoothed.se[keep]))

# ---- Test C: su -> 0 means the smoother must pass through the data ----------
cat("Test C: near-zero observation error, so the state must interpolate the data\n")
fitC <- new$TVARSS(X = Xmiss, p = p, ME = ME, Tsamplefract = .9, show.fig = FALSE,
                   annealing = FALSE, sb0.fixed = 0, sb.fixed = matrix(0, 1, p),
                   su.fixed = 1e-4)
smC_old <- old$TVARSS_KalmanSmoother(fitC)
smC_new <- new$TVARSS_KalmanSmoother(fitC)
dev_old <- max(abs(as.vector(smC_old$X.smoothed)[keep] - as.vector(Xmiss)[keep]))
dev_new <- max(abs(as.vector(smC_new$X.smoothed)[keep] - as.vector(Xmiss)[keep]))
cat(sprintf("  max |smoothed - data| at observed times:  old = %.4g   fixed = %.4g\n",
            dev_old, dev_new))
check("with su -> 0 the fixed smoother passes through the observations",
      dev_new < 1e-3)
check("with su -> 0 the old smoother does not", dev_old > 1e-3)

# ---- Test D: the filter and log-likelihood are untouched --------------------
cat("Test D: the fix must not change any fit\n")
fitB_old_engine <- old$TVARSS(X = Xmiss, p = p, ME = ME, Tsamplefract = .9,
                              show.fig = FALSE, annealing = FALSE, sb0.fixed = 0,
                              sb.fixed = matrix(0, 1, p), su.fixed = su_true)
check("TVARSS() gives an identical fit from both files",
      isTRUE(all.equal(fitB_old_engine$par.full, fitB$par.full)) &&
        abs(fitB_old_engine$logLik - fitB$logLik) < 1e-10)
check("both smoothers report the same log-likelihood (filter unchanged)",
      abs(smB_old$logLik - smB_new$logLik) < 1e-10)

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
tab <- rbind(
  "fully observed"                 = results$A[c("data", "old", "new", "ref")],
  "70% missing, at observed times" = results$B_obs[c("data", "old", "new", "ref")],
  "70% missing, all times"         = results$B_all[c("data", "old", "new", "ref")]
)
colnames(tab) <- c("raw data", "11Feb25", "12Aug26", "independent ref")

cat("\n")
cat("RMSE against the known latent series (lower is better)\n")
cat("------------------------------------------------------------------------\n")
print(round(tab, 4), na.print = "-")
cat("------------------------------------------------------------------------\n")
cat("'-' in the 11Feb25 column of the last row: no values are returned in gaps.\n")

cat("\nChecks\n")
cat(paste(checks, collapse = "\n"), "\n")
cat("\n")
if (pass) {
  cat("ALL CHECKS PASSED\n")
} else {
  stop("one or more checks FAILED - see above")
}

# ---- optional figure -------------------------------------------------------
if (!interactive()) {
  dir.create("reports/figs", recursive = TRUE, showWarnings = FALSE)
  png("reports/figs/smoother_fix_demo.png", width = 1500, height = 700, res = 130)
}
op <- par(no.readonly = TRUE)
par(mfrow = c(1, 1), mar = c(4, 4.2, 2.2, 1))
tt <- 1:Tmax
plot(tt, truth, type = "n", xlab = "Time", ylab = "X",
     ylim = range(c(truth, smB_new$X.smoothed, smB_old$X.smoothed), na.rm = TRUE),
     main = "70% missing: fixed smoother (blue) vs TVARSS_11Feb25 (red x)")
polygon(c(tt, rev(tt)),
        c(smB_new$X.smoothed - 1.96 * smB_new$X.smoothed.se,
          rev(smB_new$X.smoothed + 1.96 * smB_new$X.smoothed.se)),
        col = "#1f5f8b22", border = NA)
lines(tt, truth, col = "grey55", lwd = 2)
lines(tt, smB_new$X.smoothed, col = "#1f5f8b", lwd = 2)
points(tt, smB_old$X.smoothed, col = "#a8322d", pch = 4, cex = .8)
points(tt[keep], as.vector(Xmiss)[keep], pch = 21, bg = "white", cex = .7)
legend("topleft", bty = "n", cex = .85,
       legend = c("true latent series", "observations", "12Aug26 smoother +/- 95%", "11Feb25 smoother"),
       col = c("grey55", "black", "#1f5f8b", "#a8322d"),
       pch = c(NA, 21, NA, 4), lty = c(1, NA, 1, NA), lwd = c(2, NA, 2, NA), pt.bg = "white")
par(op)
if (!interactive()) {
  dev.off()
  cat("figure written to reports/figs/smoother_fix_demo.png\n")
}
