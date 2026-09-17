## ocfs_smoother.r -- OCFS drivers, TVARSS
##
## STEP 1: data prep, plus one null model and one full model.
## Step 2 will replace the single full fit with an expand.grid over covariate
## scenarios. Tony's example run() loop that used to sit at the bottom of this
## file is superseded by that; it is still in git history if needed.

library(dplyr)
library(readr)
library(future)
library(furrr)
library(forecast)

# TVARSS_12Aug26.r is TVARSS_11Feb25.r with the Kalman SMOOTHER corrected; the
# filter and the log-likelihood are byte-identical, so no fit changes. The old
# smoother returned values only at observed times and got them wrong on gappy
# series - see smoother_fix_demo.r, which reproduces the bug and the fix.
source("kalman-smoother/TVARSS_12Aug26.r")

all_composite <- read_csv("./data/all_composite.csv", show_col_types = FALSE)

# all_composite is stored sorted bins 309 -> 1, i.e. oldest first, so row order
# already runs forward in time. If that ever changes, the predictors would be
# silently misaligned with the response and every coefficient would be wrong -
# hence a hard stop rather than a re-sort.
stopifnot(all(diff(all_composite$bins) < 0))
plot(all_composite$ocfs)

# Response ----------------------------------------------------------------

# Fourth root, not sqrt and not log. Chosen from the per-sample uncertainties in
# ocfs_uncertainty.csv by asking which exponent makes the MEASUREMENT ERROR
# homoscedastic and symmetric: the optimum is 0.24-0.30. sqrt (0.5) leaves the
# error growing with the value (r = 0.60), log overshoots so it shrinks with the
# value (r = -0.77). 0^0.25 = 0, so the seven zero-count bins need no offset.
d <- all_composite |> mutate(ocfs_t = ocfs^0.25)

# Clip to the first real observation. The previous version set X[1] <- 100
# because TVARSS's initial-updating block has no NA check and breaks on a
# missing first value - but with su small that invented number is then treated
# as near-exact data anchoring the series at 62.7 ka.
d <- d[min(which(!is.na(d$ocfs_t))):nrow(d), ]

X <- matrix(d$ocfs_t, ncol = 1)
plot(X)

# Predictors --------------------------------------------------------------

# char_acc and mean_co2 have internal gaps (18 and 14 bins here) which are
# interpolated because TVARSS requires U complete at every time step.
# heinrich is a 0/1 stadial indicator, so it is left unscaled - its coefficient
# then reads per stadial rather than per SD.
U <- d |>
  select(char_acc, d18O, heinrich, mean_co2, PrDens) |>
  mutate(across(c(char_acc, d18O, mean_co2), ~ as.numeric(forecast::na.interp(.))),
         across(-heinrich, ~ as.numeric(scale(.)))) |>
  as.matrix()

stopifnot(nrow(U) == nrow(X), !anyNA(U))


# Model settings ----------------------------------------------------------

p <- 2                    # lags, so 200 and 400 yr

# Measurement error is a single number for this first pass, pending a
# conversation with the analysts about how the csv's uncertainties were derived.
# su = 1 happens to be close to right on this scale: the mean per-bin
# observation SD implied by the csv is 1.013. (On the raw concentration scale,
# where sd(ocfs) = 5867, su = 1 was absurd.) Fixing su rather than estimating it
# also avoids a degenerate solution where se collapses to 0 and the latent
# series becomes a deterministic AR(2) that traces the data.
# To move to per-bin weighting later, replace this one line with a vector.
ME <- rep(1, nrow(X))
su.fixed <- 1

# sb0 = sb = 0 pins the autoregression coefficients, so this is an ordinary
# AR(2) state-space model. The "time-varying" part of TVARSS is switched off;
# freeing it cost 5-12 AICc when tested.
sb0.fixed <- 0
sb.fixed <- matrix(0, 1, p)

n_obs <- sum(!is.na(X))   # 83 observations in 293 bins; use n_obs for AICc


# Fitting -----------------------------------------------------------------

# Predictors enter through c (the process equation), so their effect is filtered
# through the autoregression and accumulates while the covariate persists.
# d.fixed = 0 keeps them out of the observation equation; c and d are not
# jointly identifiable with 83 points.
fit_ocfs <- function(U, c.start, b0.start = NA, b.start = rep(NA, p)) {
  args <- list(X = X, p = p, ME = ME, su.fixed = su.fixed, Tsamplefract = .9,
               show.fig = FALSE, annealing = FALSE, initial.points = "stationary",
               sb0.fixed = sb0.fixed, sb.fixed = sb.fixed,
               b0.start = b0.start, b.start = b.start)
  if (!is.null(U)) {
    U <- as.matrix(U)
    # c.start is normally one number used for every covariate, but a full vector
    # can be passed to start from another model's solution.
    args <- c(args, list(U = U,
                         c.fixed = rep(NA, ncol(U)),
                         d.fixed = rep(0, ncol(U)),
                         c.start = if (length(c.start) == ncol(U)) c.start
                                   else rep(c.start[1], ncol(U))))
  }
  do.call(TVARSS, args)
}

# One fit is never enough here. TVARSS optimises with Nelder-Mead, which reports
# convergence = 0 when its simplex stops moving, not when it has found the
# optimum - different starting values land on log-likelihoods several units
# apart and can flip coefficient signs. So try a grid and keep the best.
#
# Every fit is kept, not just the winner: the spread across starts is the
# evidence that a single fit cannot be trusted, so discarding it would discard
# the diagnostic. Inspect it with attr(mod, "starts") for the summary table, or
# attr(mod, "fits") for the fitted objects themselves.
#
# `extra` takes a list of full c.start vectors, for starting from a smaller
# model's solution with zeros for the new terms. That guarantees the larger model
# can always reach at least the smaller one's likelihood, which is what stops
# nested comparisons coming out negative.
best_fit <- function(U, starts = c(0.01, 0.05, 0.1, 0.25, 0.5), extra = list(), ...) {
  if (is.null(U)) starts <- NA          # nothing in c to start, so one fit only
  grid <- c(as.list(starts), extra)
  fits <- lapply(grid, function(s) try(fit_ocfs(U, s, ...), silent = TRUE))
  ok <- !vapply(fits, inherits, logical(1), "try-error")
  stopifnot(any(ok))
  fits <- fits[ok]

  summ <- data.frame(
    start       = vapply(grid[ok], function(s) s[1], numeric(1)),
    logLik      = vapply(fits, function(m) m$logLik, numeric(1)),
    se          = vapply(fits, function(m) m$se, numeric(1)),
    sum_b       = vapply(fits, function(m) sum(m$b), numeric(1)),
    convergence = vapply(fits, function(m) m$opt.convergence, numeric(1))
  )
  best <- fits[[which.max(summ$logLik)]]
  attr(best, "starts") <- summ
  attr(best, "fits")   <- fits
  best
}

mod0 <- best_fit(NULL)                                        # no covariates
mod1 <- best_fit(U, b0.start = mod0$b0, b.start = mod0$b)     # all five

attr(mod1, "starts")   # how much did the starting value matter?


# Result ------------------------------------------------------------------

# A negative deviance between nested models is arithmetically impossible at the
# maxima, so it means the optimiser failed on one of them. pchisq() of a
# negative deviance returns 1, which is how the previous version of this script
# produced a table of P = 1 and read it as "no effect". Stop instead.
dev <- 2 * (mod1$logLik - mod0$logLik)
stopifnot(dev >= -1e-6)

mod0
mod1
c(dev = dev, df = ncol(U), P = pchisq(dev, df = ncol(U), lower.tail = FALSE))

# Verification: se should be a sensible fraction of the response variation. If
# it comes back near zero the degenerate branch has appeared and the fit is not
# describing a stochastic process.
c(sd_X = sd(X, na.rm = TRUE), se = mod1$se, se_frac = mod1$se / sd(X, na.rm = TRUE),
  sum_b = sum(mod1$b), n_obs = n_obs, n_bins = nrow(X))


# Covariate scenarios -----------------------------------------------------

# Covariates grouped into hypothesis blocks, switched on and off together.
# d18O and mean_co2 share a block because they correlate at 0.89: fitted
# separately each looks weak while the pair is informative, so testing them
# individually would understate the climate signal.
blocks <- list(
  CLIM   = c("d18O", "mean_co2"),
  ABRUPT = "heinrich",
  FIRE   = "char_acc",
  HUMAN  = "PrDens",
  TIME   = "time"
)

# TIME is a plain linear trend, and it exists so that HUMAN has to earn its
# keep. PrDens is exactly zero before ~13 ka (bins 62-309, so only 18 of the 83
# observations lie in its support), which means on its own it cannot be told
# apart from "the late record" - any late decline in OCFS would be credited to
# human population by default. An indicator for bins <= 61 is the harsher
# version of this test and can be swapped in here later.
U_all <- cbind(U, time = as.numeric(scale(seq_len(nrow(X)))))

# Every on/off combination of the five blocks. Row 1 is all FALSE, the null.
scen <- expand.grid(CLIM   = c(FALSE, TRUE), ABRUPT = c(FALSE, TRUE),
                    FIRE   = c(FALSE, TRUE), HUMAN  = c(FALSE, TRUE),
                    TIME   = c(FALSE, TRUE))

# Each block appears in exactly half the models, which is what makes the summed
# weights at the bottom comparable between blocks.
stopifnot(colSums(scen[names(blocks)]) == nrow(scen) / 2)

scen$label <- apply(scen[names(blocks)], 1, \(on)
  if (any(on)) paste(names(blocks)[on], collapse = "+") else "null")

# Warm-started from the null fit, and multi-started within each scenario.
#
# 32 scenarios x 5 starts is ~160 TVARSS fits at roughly 30 s each, so this is
# the slow step: about 85 min sequentially, ~15 min on 6 workers. The result is
# cached because the downstream table gets edited far more often than the fits
# change - delete the file to force a refit. It lives under results/cache/,
# which is gitignored.
fits_file <- "results/cache/ocfs_scenario_fits.rds"

if (file.exists(fits_file)) {
  fits <- readRDS(fits_file)
} else {
  plan(multisession, workers = 6)
  # Progress is not streamed back from workers, so this runs quietly.
  fits <- future_map(seq_len(nrow(scen)), \(i) {
    cols <- unlist(blocks[unlist(scen[i, names(blocks)])], use.names = FALSE)
    Ui <- if (length(cols)) U_all[, cols, drop = FALSE] else NULL
    best_fit(Ui, b0.start = mod0$b0, b.start = mod0$b)
  }, .options = furrr_options(seed = 1984))
  plan(sequential)

  dir.create(dirname(fits_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(fits, fits_file)
}
names(fits) <- scen$label

# Did any scenario land on a meaningfully different optimum depending on where
# it started? If the spread is wide for some scenario, widen `starts` rather
# than trusting its winner.
spread <- vapply(fits, \(m) diff(range(attr(m, "starts")$logLik)), numeric(1))
sort(spread, decreasing = TRUE)[1:5]


# Model comparison --------------------------------------------------------

scen$k      <- vapply(fits, \(m) m$npar,   numeric(1))
scen$logLik <- vapply(fits, \(m) m$logLik, numeric(1))

# AICc on the 83 OBSERVATIONS, not the 293 rows of the series. The correction
# term is about 2.5 here, so it is not negligible. TVARSS's own $AIC uses a
# constant based on the series length, so it is only comparable within one
# response definition - don't reuse it across transforms or clippings.
scen$AICc  <- with(scen, -2 * logLik + 2 * k + 2 * k * (k + 1) / (n_obs - k - 1))
scen$dAICc <- scen$AICc - min(scen$AICc)
scen$w     <- exp(-scen$dAICc / 2) / sum(exp(-scen$dAICc / 2))

# Every scenario nests the null, so each one also gets a likelihood-ratio test
# against it. A negative deviance would mean the optimiser failed, not that the
# covariates are useless - stop rather than report P = 1.
i0 <- which(scen$label == "null")
scen$dev <- 2 * (scen$logLik - scen$logLik[i0])
scen$df  <- scen$k - scen$k[i0]
scen$P   <- pchisq(scen$dev, scen$df, lower.tail = FALSE)
scen$P[i0] <- NA          # the null against itself is not a test
stopifnot(scen$dev >= -1e-6)

results <- scen[order(scen$AICc),
                c("label", "k", "logLik", "AICc", "dAICc", "w", "dev", "df", "P")]
head(results, 12)

# Summed Akaike weight per block: the weight of evidence that a block belongs in
# the model, pooled over every model containing it. Comparable between blocks
# only because the design above is balanced. ~0.5 is the no-information value.
importance <- vapply(names(blocks), \(b) sum(scen$w[scen[[b]]]), numeric(1))
sort(importance, decreasing = TRUE)

# TVARSS returns c as an unnamed matrix, so take the names from the U it stored.
best <- fits[[results$label[1]]]
coefs <- data.frame(term = colnames(best$U), c = as.vector(best$c))


# Drop-one variable importance --------------------------------------------

# Every model above is compared with the NULL. This instead compares each model
# with itself minus one block, which is the usual variable-importance test and
# what Tony's template does (it tests `mod` against `mod1`, the main effects).
#
# No refitting is needed. The grid is every subset of the five blocks, so for any
# model and any block inside it, the model with that block removed is already in
# the grid - each test is a lookup.
ll <- setNames(scen$logLik, scen$label)
kk <- setNames(scen$k, scen$label)

drop1 <- do.call(rbind, lapply(seq_len(nrow(scen)), \(i) {
  on <- unlist(scen[i, names(blocks)])
  if (!any(on)) return(NULL)                      # nothing to drop from the null
  do.call(rbind, lapply(names(blocks)[on], \(b) {
    keep <- on; keep[b] <- FALSE
    red <- if (any(keep)) paste(names(blocks)[keep], collapse = "+") else "null"
    data.frame(model = scen$label[i], dropped = b, reduced = red,
               dev = unname(2 * (ll[scen$label[i]] - ll[red])),
               df  = unname(kk[scen$label[i]] - kk[red]))
  }))
}))
drop1$P <- pchisq(pmax(drop1$dev, 0), drop1$df, lower.tail = FALSE)

# Same reasoning as the omnibus guard: a negative deviance between nested models
# means the optimiser failed on the larger one. A warning rather than a stop,
# because values of order 1e-3 are just numerical noise - but anything larger
# means that row's test cannot be reported.
if (any(drop1$dev < -1e-6)) {
  warning("negative drop-one deviance(s) - the larger model is not at its maximum")
  print(drop1[drop1$dev < -1e-6, ], row.names = FALSE, digits = 4)
}

# CLIM is a 2-df test because d18O and mean_co2 enter and leave together. That is
# deliberate: they correlate at 0.89, and splitting them makes both look weak
# while the pair is strong.
drop1[drop1$model == results$label[1], ]


# Interactions within the best model ---------------------------------------

# Are the covariate effects additive? Six pairwise interactions among the four
# variables in the best-supported model, built as raw products the way Tony's
# template does (cbind(U1, inter = a*b)), and each tested against the
# MAIN-EFFECTS model rather than against the null.
best_blocks <- names(blocks)[unlist(scen[scen$label == results$label[1], names(blocks)])]
best_cols   <- unlist(blocks[best_blocks], use.names = FALSE)
U_best      <- U_all[, best_cols, drop = FALSE]

prs <- combn(best_cols, 2, simplify = FALSE)
INT <- sapply(prs, \(pr) U_all[, pr[1]] * U_all[, pr[2]])
colnames(INT) <- vapply(prs, paste, character(1), collapse = ":")

# READ THIS BEFORE INTERPRETING ANY INTERACTION COEFFICIENT.
# These products are near-duplicates of their own main effects. PrDens was
# standardised, so its 232 structural zeros became a constant -0.229, which makes
# PrDens x d18O very nearly a rescaled copy of d18O across 79% of the record.
#   The likelihood-ratio tests and AICc below are still valid - adding the raw
#   product spans the same model space however it is parameterised, so the
#   maximised likelihood and the df are unaffected.
#   The individual c coefficients are NOT interpretable for the terms near 1
#   below, and the main effects in those models will shift a long way from the
#   base. Check the `spread` column: a wide spread across starting values means
#   the optimiser is wandering along a near-flat ridge.
round(apply(INT, 2, \(z) max(abs(cor(z, U_best)))), 3)

# Note on the three X:time terms: they are the parametric version of "does this
# driver's effect change through time", which is what TVARSS's sb tests
# non-parametrically - and that was tested and rejected (freeing sb0 cost 5.09
# AICc, sb 6.92, both 12.09). A significant X:time interaction would be in
# tension with that and is worth investigating, not reporting at face value.

int_specs <- c(list(none = character(0)),
               setNames(as.list(colnames(INT)), colnames(INT)),
               list(all = colnames(INT)))

int_U <- \(cols) if (length(cols)) cbind(U_best, INT[, cols, drop = FALSE]) else U_best

int_file <- "results/cache/ocfs_interaction_fits.rds"
if (file.exists(int_file)) {
  int_fits <- readRDS(int_file)
} else {
  # Fit the base first, then start every interaction model FROM the base solution
  # with the new coefficient(s) at zero, as well as from the usual grid. That
  # starting point reproduces the base exactly, so an interaction model can never
  # come out worse than the model it nests - which is what happened on the first
  # attempt here (mean_co2:PrDens landed 0.0099 short) and what leaves the
  # all-six model wandering, since its log-likelihood spanned 16 units across
  # plain starting values.
  ibase <- best_fit(U_best, b0.start = mod0$b0, b.start = mod0$b)
  plan(multisession, workers = 6)
  int_fits <- future_map(int_specs[-1], \(cols)
    best_fit(int_U(cols),
             extra = list(c(as.vector(ibase$c), rep(0, length(cols)))),
             b0.start = ibase$b0, b.start = ibase$b),
    .options = furrr_options(seed = 1984))
  plan(sequential)
  int_fits <- c(list(none = ibase), int_fits)
  dir.create(dirname(int_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(int_fits, int_file)
}

int_res <- data.frame(
  interaction = names(int_fits),
  k      = vapply(int_fits, \(m) m$npar, integer(1)),
  logLik = vapply(int_fits, \(m) m$logLik, numeric(1)),
  spread = vapply(int_fits, \(m) diff(range(attr(m, "starts")$logLik)), numeric(1))
)
int_res$AICc  <- with(int_res, -2*logLik + 2*k + 2*k*(k+1)/(n_obs - k - 1))
int_res$dAICc <- int_res$AICc - min(int_res$AICc)
ib <- which(int_res$interaction == "none")
int_res$dev <- 2 * (int_res$logLik - int_res$logLik[ib])
int_res$df  <- int_res$k - int_res$k[ib]
int_res$P   <- pchisq(pmax(int_res$dev, 0), int_res$df, lower.tail = FALSE)
int_res$P[ib] <- NA

# The base is the same model as the grid's best, on the same data, so its AICc
# must match. If it does not, the refit found a different optimum and nothing in
# this table can be compared with the grid.
stopifnot(abs(int_res$AICc[ib] - results$AICc[1]) < 1e-6)
stopifnot(int_res$dev >= -1e-6)

int_res[order(int_res$AICc), ]


# Smoothed series ---------------------------------------------------------

# The corrected smoother, applied to the best-supported model. This is what the
# fix was for: an interpolated OCFS series across the 210 bins with no
# observation, each with an honest standard error. The old smoother returned
# nothing at those bins at all.
sm <- TVARSS_KalmanSmoother(best)
stopifnot(!anyNA(sm$X.smoothed),          # a value at every bin, not just the 83
          !is.null(sm$X.smoothed.se))     # NULL here means the old smoother

smoothed <- data.frame(
  bin      = d$bins,
  age      = d$age,
  ocfs     = d$ocfs,                      # raw concentration, NA where unobserved
  obs      = as.vector(sm$X),             # the response actually fed in
  smooth   = as.vector(sm$X.smoothed),
  se       = as.vector(sm$X.smoothed.se)
)
smoothed$lower <- smoothed$smooth - 1.96 * smoothed$se
smoothed$upper <- smoothed$smooth + 1.96 * smoothed$se

# Back to concentration units. The interval is transformed, not the SE, because
# ^4 is non-linear - so the band is asymmetric on this scale, correctly.
smoothed$ocfs_smooth <- pmax(smoothed$smooth, 0)^4
smoothed$ocfs_lower  <- pmax(smoothed$lower,  0)^4
smoothed$ocfs_upper  <- pmax(smoothed$upper,  0)^4
stopifnot(smoothed$ocfs_smooth >= 0)

# How much does the smoother actually know where there is no data? The SE in the
# gaps against the sd of the response is the honest answer.
gap <- is.na(smoothed$obs)
c(se_at_obs = mean(smoothed$se[!gap]), se_in_gaps = mean(smoothed$se[gap]),
  sd_X = sd(X, na.rm = TRUE),
  var_resolved_in_gaps = 1 - (mean(smoothed$se[gap]) / sd(X, na.rm = TRUE))^2)


# Save --------------------------------------------------------------------

out <- "results/ocfs"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.csv(results, file.path(out, "ocfs_model_ranking.csv"), row.names = FALSE)
write.csv(data.frame(block = names(importance), summed_weight = importance),
          file.path(out, "ocfs_importance.csv"), row.names = FALSE)
write.csv(data.frame(label = names(spread), logLik_spread = spread),
          file.path(out, "ocfs_start_spread.csv"), row.names = FALSE)
write.csv(coefs, file.path(out, "ocfs_top_coefficients.csv"), row.names = FALSE)
write.csv(drop1, file.path(out, "ocfs_drop_one.csv"), row.names = FALSE)
write.csv(int_res, file.path(out, "ocfs_interactions.csv"), row.names = FALSE)
write.csv(smoothed, file.path(out, "ocfs_smoothed.csv"), row.names = FALSE)

png(file.path(out, "ocfs_smoothed.png"), width = 1600, height = 1000, res = 140)
op <- par(no.readonly = TRUE); par(mfrow = c(2, 1), mar = c(4, 4.6, 2.2, 1))
ka <- approx(which(!is.na(smoothed$age)), smoothed$age[!is.na(smoothed$age)],
             xout = seq_len(nrow(smoothed)), rule = 2)$y / 1000
plot(ka, smoothed$smooth, type = "n", xlim = rev(range(ka)),
     ylim = range(c(smoothed$lower, smoothed$upper)),
     xlab = "Age (ka BP)", ylab = expression(ocfs^0.25),
     main = paste0("Smoothed OCFS, 95% band  [", results$label[1], "]"))
polygon(c(ka, rev(ka)), c(smoothed$lower, rev(smoothed$upper)),
        col = "#1f5f8b22", border = NA)
lines(ka, smoothed$smooth, col = "#1f5f8b", lwd = 2)
points(ka, smoothed$obs, pch = 21, bg = "white", cex = .6)
plot(ka, smoothed$ocfs_smooth, type = "n", xlim = rev(range(ka)),
     ylim = c(0, max(smoothed$ocfs_upper)),
     xlab = "Age (ka BP)", ylab = "OCFS concentration",
     main = "Back-transformed to concentration units")
polygon(c(ka, rev(ka)), c(smoothed$ocfs_lower, rev(smoothed$ocfs_upper)),
        col = "#18776622", border = NA)
lines(ka, smoothed$ocfs_smooth, col = "#187766", lwd = 2)
points(ka, smoothed$ocfs, pch = 21, bg = "white", cex = .6)
par(op); dev.off()
