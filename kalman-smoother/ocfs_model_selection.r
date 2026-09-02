## ocfs_model_selection.r  --  Tier 3: study design and multi-model inference
##
## Requires results/ocfs/tier2_fits.rds, so run ocfs_analysis.r first.
##
## Design. Five a priori hypothesis blocks, giving 2^5 = 32 models. Each block
## appears in exactly 16 of the 32, which is what makes summed Akaike weights
## comparable between blocks (Burnham & Anderson 2002, sec. 6.9.6 - summed
## weights are only interpretable when the model set is balanced).
##
##   CLIM    d18O + mean_co2   kept together: they correlate at r = 0.88 and are
##                             not separately identifiable. Splitting them would
##                             make both look weak while the pair is strong.
##   ABRUPT  heinrich          0/1 stadial indicator
##   FIRE    char_acc          charcoal accumulation
##   HUMAN   PrDens            human population density (Kelly 2025 SPD)
##   TIME    a late-record term
##
## Why TIME exists. PrDens is exactly zero for bins 62-309 and nonzero only in
## the contiguous block 1-61 (1.1-13.1 ka); only 18 of the 83 OCFS observations
## fall in its support. So on its own it cannot be told apart from "the last 13
## kyr", and a significant PrDens might be nothing more than the terminal OCFS
## decline. TIME is the competitor that makes HUMAN interpretable: PrDens earns
## weight only where it outperforms a generic late-record term.
##
## The set is fitted twice, once per TIME definition:
##   time_lin  a linear trend in time (standardised row index, forward in time)
##   time_hol  an indicator for bins 1-61, i.e. exactly PrDens's own support
## time_hol is the harsher test. If HUMAN keeps its weight against it, the claim
## is about the magnitude of human population rather than merely about being in
## the late record.
##
## Run:  Rscript kalman-smoother/ocfs_model_selection.r
## Fits are cached per model in results/cache/; the script is restartable.

rm(list = ls())
suppressMessages({
  library(readr); library(forecast); library(future); library(furrr)
})
source("kalman-smoother/TVARSS_12Aug26.r")
source("kalman-smoother/ocfs_helpers.r")

ROOT <- normalizePath(getwd(), winslash = "/")
OUT  <- "results/ocfs"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

RUN_SENSITIVITY <- TRUE
WORKERS <- max(1, min(10, future::availableCores() - 2))

rule <- function(txt) cat("\n", strrep("=", 74), "\n", txt, "\n", strrep("=", 74), "\n", sep = "")

t2 <- readRDS(file.path(OUT, "tier2_fits.rds"))
dat <- t2$dat; n_obs <- t2$n_obs; p <- t2$modF$p
X <- dat$X
# observation variance is counting_var[t] + tau^2 with su fixed at 1, as settled
# in ocfs_analysis.r section 2
ME <- t2$ME; SU <- t2$SU; TAU <- t2$TAU; K_EXTRA <- t2$K_EXTRA
cat("observation model: tau =", TAU, " su =", SU, " extra parameters =", K_EXTRA, "\n")

# tau is a parameter TVARSS does not count, so every AICc here adds K_EXTRA. It
# is the same constant for every model in the set, so rankings are unaffected -
# but the absolute values must match the ones in ocfs_analysis.r.
aicc2 <- function(mod, n = n_obs) {
  k <- mod$npar + K_EXTRA
  if (n - k - 1 <= 0) return(NA_real_)
  -2 * mod$logLik + 2 * k + 2 * k * (k + 1) / (n - k - 1)
}

# ---- predictor matrix, with both TIME variants appended ----------------------
Ubase <- dat$U
time_lin <- as.numeric(scale(seq_len(nrow(X))))       # rows run forward in time
time_hol <- as.numeric(dat$bins <= 61)                # PrDens's exact support
cat("time_hol is 1 in", sum(time_hol), "of", nrow(X), "bins;",
    sum(time_hol == 1 & dat$obs), "of", n_obs, "observations\n")
cat("cor(time_lin, PrDens) =", round(cor(time_lin, Ubase[, "PrDens"]), 3),
    "  cor(time_hol, PrDens) =", round(cor(time_hol, Ubase[, "PrDens"]), 3), "\n")

BLOCKS <- list(CLIM = c("d18O", "mean_co2"), ABRUPT = "heinrich",
               FIRE = "char_acc", HUMAN = "PrDens", TIME = "time")

make_U <- function(blocks, time_var) {
  cols <- unlist(BLOCKS[blocks], use.names = FALSE)
  if (!length(cols)) return(NULL)
  M <- matrix(nrow = nrow(X), ncol = 0)
  for (cl in cols) {
    v <- if (cl == "time") time_var else Ubase[, cl]
    M <- cbind(M, v)
  }
  colnames(M) <- cols
  M
}

# ---- parallel fitting with per-model caching --------------------------------
# Everything sent to a worker is plain data: the design matrices are built here
# in the parent, so no closures or environments are serialised. Each model is
# cached separately, which makes the whole script restartable.
run_set <- function(spec_list, time_var, time_name,
                    starts = c(0.001, 0.01, 0.05, 0.1, 0.25)) {
  keys <- vapply(spec_list, function(b)
    paste0("ms_", time_name, "_", if (length(b)) paste(b, collapse = "-") else "null"), "")
  cache_files <- file.path("results/cache", paste0(keys, ".rds"))
  dir.create("results/cache", recursive = TRUE, showWarnings = FALSE)
  have <- file.exists(cache_files)
  cat(sprintf("  %d models, %d already cached, fitting %d on %d worker(s)\n",
              length(keys), sum(have), sum(!have), WORKERS))

  if (any(!have)) {
    todo  <- which(!have)
    Ulist <- lapply(spec_list[todo], make_U, time_var = time_var)
    Xw <- X; MEw <- ME; pw_ <- p; SUw <- SU
    b0s <- t2$mod0$b0; bs <- t2$mod0$b; rootw <- ROOT; startsw <- starts

    plan(multisession, workers = WORKERS)
    on.exit(plan(sequential), add = TRUE)
    res <- future_map(seq_along(todo), function(j) {
      setwd(rootw)
      suppressMessages({
        library(forecast)
        source("kalman-smoother/TVARSS_12Aug26.r")
        source("kalman-smoother/ocfs_helpers.r")
      })
      out <- try(fit_ms(Xw, Ulist[[j]], MEw, p = pw_, su.fixed = SUw,
                        starts = startsw, b0.start = b0s, b.start = bs),
                 silent = TRUE)
      if (inherits(out, "try-error")) NULL else out
    }, .options = furrr_options(seed = 1984))
    plan(sequential)

    for (j in seq_along(todo))
      if (!is.null(res[[j]])) saveRDS(res[[j]], cache_files[todo[j]])
  }

  fits <- lapply(cache_files, function(f) if (file.exists(f)) readRDS(f) else NULL)
  names(fits) <- keys
  fits
}

# =============================================================================
rule("1. HUMAN vs TIME head-to-head (run first: it can change the headline)")
# =============================================================================

for (tn in c("time_lin", "time_hol")) {
  tv <- get(tn)
  cat("\n--", tn, "--\n")
  specs <- list(character(0), "HUMAN", "TIME", c("HUMAN", "TIME"))
  fits <- run_set(specs, tv, paste0("h2h_", tn))
  lab <- c("null", "HUMAN only", "TIME only", "HUMAN + TIME")
  tabh <- data.frame(model = lab,
                     k = vapply(fits, function(m) if (is.null(m)) NA_integer_ else m$npar, 1L),
                     logLik = vapply(fits, function(m) if (is.null(m)) NA_real_ else m$logLik, 0),
                     AICc = vapply(fits, function(m) if (is.null(m)) NA_real_ else aicc2(m), 0))
  tabh$dAICc <- tabh$AICc - min(tabh$AICc, na.rm = TRUE)
  print(tabh, digits = 5, row.names = FALSE)
  # is PrDens worth anything once TIME is in the model?
  if (!is.null(fits[[3]]) && !is.null(fits[[4]])) {
    dv <- 2 * (fits[[4]]$logLik - fits[[3]]$logLik)
    cat(sprintf("  adding HUMAN to TIME: dev = %.4f, df = 1, P = %.4g\n",
                dv, pchisq(max(dv, 0), 1, lower.tail = FALSE)))
    if (dv < -1e-6) cat("  (negative deviance: optimiser trouble, treat as ~0)\n")
  }
  if (!is.null(fits[[2]]) && !is.null(fits[[4]])) {
    dv <- 2 * (fits[[4]]$logLik - fits[[2]]$logLik)
    cat(sprintf("  adding TIME to HUMAN: dev = %.4f, df = 1, P = %.4g\n",
                dv, pchisq(max(dv, 0), 1, lower.tail = FALSE)))
  }
}

# =============================================================================
rule("2. The full balanced set: 32 models per TIME definition")
# =============================================================================

bnames <- names(BLOCKS)
all_specs <- lapply(0:(2^length(bnames) - 1), function(m)
  bnames[as.logical(bitwAnd(m, 2^(seq_along(bnames) - 1)))])
stopifnot(length(all_specs) == 32)
# balance check: every block must appear in exactly half the models
appear <- vapply(bnames, function(b) sum(vapply(all_specs, function(s) b %in% s, TRUE)), 1L)
cat("block appearances (must all be 16):", paste(names(appear), appear, sep = "=", collapse = ", "), "\n")
stopifnot(all(appear == 16))

sets <- list()
for (tn in c("time_lin", "time_hol")) {
  cat("\n--", tn, "--\n")
  sets[[tn]] <- run_set(all_specs, get(tn), tn)
}

# ---- AICc table, weights, evidence ratios ----------------------------------
build_table <- function(fits, specs) {
  ok <- !vapply(fits, is.null, TRUE)
  tab <- data.frame(
    model  = vapply(specs, function(s) if (length(s)) paste(s, collapse = "+") else "null", ""),
    nblock = vapply(specs, length, 1L),
    k      = vapply(fits, function(m) if (is.null(m)) NA_integer_ else m$npar, 1L),
    logLik = vapply(fits, function(m) if (is.null(m)) NA_real_ else m$logLik, 0),
    stringsAsFactors = FALSE
  )
  tab$AICc  <- vapply(fits, function(m) if (is.null(m)) NA_real_ else aicc2(m), 0)
  tab$dAICc <- tab$AICc - min(tab$AICc, na.rm = TRUE)
  tab$w     <- exp(-tab$dAICc / 2); tab$w <- tab$w / sum(tab$w, na.rm = TRUE)
  tab$evid  <- max(tab$w, na.rm = TRUE) / tab$w
  tab$failed <- !ok
  tab[order(tab$dAICc), ]
}

tabs <- lapply(names(sets), function(tn) build_table(sets[[tn]], all_specs))
names(tabs) <- names(sets)

for (tn in names(tabs)) {
  cat("\n---- model ranking,", tn, "(top 12 of 32) ----\n")
  tt <- tabs[[tn]]
  if (any(tt$failed)) cat("  NOTE:", sum(tt$failed), "model(s) failed to fit and are excluded\n")
  print(head(tt[, c("model", "k", "logLik", "AICc", "dAICc", "w", "evid")], 12),
        digits = 4, row.names = FALSE)
  cat("  models within 2 AICc:", sum(tt$dAICc <= 2, na.rm = TRUE),
      " within 4:", sum(tt$dAICc <= 4, na.rm = TRUE), "\n")
  write.csv(tt, file.path(OUT, paste0("model_ranking_", tn, ".csv")), row.names = FALSE)
}

# ---- summed Akaike weights per block ---------------------------------------
cat("\n---- variable importance: summed Akaike weights ----\n")
imp <- data.frame(block = bnames)
for (tn in names(tabs)) {
  tt <- tabs[[tn]]
  w_by_model <- setNames(tt$w, tt$model)
  imp[[tn]] <- vapply(bnames, function(b) {
    inc <- vapply(all_specs, function(s) b %in% s, TRUE)
    lbl <- vapply(all_specs, function(s) if (length(s)) paste(s, collapse = "+") else "null", "")
    sum(w_by_model[lbl[inc]], na.rm = TRUE)
  }, 0)
}
print(imp, digits = 3, row.names = FALSE)
cat("(balanced set, so these are comparable between blocks; ~0.5 means no support)\n")
write.csv(imp, file.path(OUT, "variable_importance.csv"), row.names = FALSE)

# =============================================================================
rule("3. Model-averaged coefficients with unconditional SEs")
# =============================================================================
# Burnham & Anderson eq. 4.9: the unconditional variance of a model-averaged
# estimate adds the between-model spread to the within-model variance,
#   var_uncond = sum_i w_i ( var_i + (theta_i - theta_bar)^2 ).
# Two averages are reported because they answer different questions:
#   conditional - averaged only over models containing the term, i.e. "given the
#                 term is in the model, how big is it"
#   full        - averaged over all 32 with zero substituted where the term is
#                 absent, a shrinkage estimate that folds in selection
#                 uncertainty. This is the one to quote for effect sizes.

model_average <- function(fits, specs, tab, time_var) {
  lbl <- vapply(specs, function(s) if (length(s)) paste(s, collapse = "+") else "null", "")
  w <- setNames(tab$w, tab$model)[lbl]
  terms_all <- unlist(BLOCKS, use.names = FALSE)
  est <- se <- matrix(NA_real_, length(specs), length(terms_all),
                      dimnames = list(lbl, terms_all))

  for (i in seq_along(fits)) {
    m <- fits[[i]]
    if (is.null(m) || is.null(m$c)) next
    U <- make_U(specs[[i]], time_var)
    s <- tryCatch(c_se(m, X, U, ME, SU), error = function(e) {
      warning("SE failed for ", lbl[i], ": ", conditionMessage(e)); NULL })
    cn <- colnames(U)
    est[i, cn] <- as.vector(m$c)
    if (!is.null(s)) se[i, cn] <- s[cn]
  }

  out <- do.call(rbind, lapply(terms_all, function(v) {
    inc <- !is.na(est[, v])
    if (!any(inc)) return(NULL)
    wi <- w[inc]; th <- est[inc, v]; vi <- se[inc, v]^2
    vi[is.na(vi)] <- 0                      # a missing SE contributes spread only
    wc <- wi / sum(wi)
    th_cond <- sum(wc * th)
    se_cond <- sqrt(sum(wc * (vi + (th - th_cond)^2)))
    th_full <- sum(wi * th)                 # zero where the term is absent
    se_full <- sqrt(sum(wi * (vi + (th - th_full)^2)) +
                    (1 - sum(wi)) * th_full^2)
    data.frame(term = v, n_models = sum(inc), sum_w = sum(wi),
               cond = th_cond, cond_se = se_cond,
               full = th_full, full_se = se_full,
               n_missing_se = sum(is.na(se[inc, v])))
  }))
  out$z_full <- out$full / out$full_se
  out
}

avg <- list()
for (tn in names(tabs)) {
  cat("\n---- model-averaged c coefficients,", tn, "----\n")
  avg[[tn]] <- model_average(sets[[tn]], all_specs, tabs[[tn]], get(tn))
  print(avg[[tn]], digits = 3, row.names = FALSE)
  if (any(avg[[tn]]$n_missing_se > 0))
    cat("  NOTE: some models gave no usable Hessian; their SE contributes 0 and\n",
        "  the unconditional SE is therefore slightly optimistic for those terms.\n")
  write.csv(avg[[tn]], file.path(OUT, paste0("model_averaged_", tn, ".csv")), row.names = FALSE)
}

# =============================================================================
rule("4. Top model: ascent, drop-one, and a test of time-varying coefficients")
# =============================================================================

tn_best <- names(tabs)[which.min(vapply(tabs, function(t) min(t$AICc, na.rm = TRUE), 0))]
cat("TIME definition with the better-supported set:", tn_best, "\n")
tt <- tabs[[tn_best]]
best_lbl <- tt$model[1]
best_spec <- all_specs[[which(vapply(all_specs, function(s)
  identical(if (length(s)) paste(s, collapse = "+") else "null", best_lbl), TRUE))]]
cat("top model:", best_lbl, "\n")

Ubest <- make_U(best_spec, get(tn_best))
mbest <- sets[[tn_best]][[paste0("ms_", tn_best, "_",
          if (length(best_spec)) paste(best_spec, collapse = "-") else "null")]]
if (is.null(mbest)) stop("could not retrieve the top model from the cache")

if (!is.null(Ubest)) {
  d1 <- cached(paste0("top_drop_one_", tn_best),
               drop_one(mbest, X, Ubest, ME, p = p, su.fixed = SU, n_obs = n_obs))
  mbest <- d1$full
  cat("\ndrop-one tests within the top model:\n")
  print(d1$table, digits = 4, row.names = FALSE)
  se_b <- c_se(mbest, X, Ubest, ME, SU)
  cat("\nWald SEs from the numerical Hessian, with long-run effects:\n")
  print(data.frame(term = colnames(Ubest), c = as.vector(mbest$c), se = as.vector(se_b),
                   longrun = longrun(mbest),
                   longrun_se = as.vector(se_b) / (1 - sum(mbest$b))),
        digits = 4, row.names = FALSE)
  cat(sprintf("\ndominant eigenvalue = %.4f;  1 - sum(b) = %.4f\n",
              dominant_eigen(mbest), 1 - sum(mbest$b)))
}

cat("\n-- is there support for time-varying coefficients? --\n")
tv_tests <- list(
  "sb0 free"      = list(sb0.fixed = NA, sb.fixed = matrix(0, 1, p)),
  "sb free"       = list(sb0.fixed = 0,  sb.fixed = matrix(NA, 1, p)),
  "sb0 + sb free" = list(sb0.fixed = NA, sb.fixed = matrix(NA, 1, p))
)
tvtab <- data.frame(variant = c("fixed (as fitted)", names(tv_tests)),
                    logLik = NA_real_, AICc = NA_real_, sb0 = NA_real_)
tvtab$logLik[1] <- mbest$logLik; tvtab$AICc[1] <- aicc2(mbest); tvtab$sb0[1] <- 0
for (i in seq_along(tv_tests)) {
  a <- tv_tests[[i]]
  m <- cached(paste0("tv_", gsub("[^a-z0-9]", "", names(tv_tests)[i]), "_", tn_best),
    try(fit_ocfs(X, Ubest, ME, p = p, su.fixed = SU,
                 c.start = if (is.null(Ubest)) NULL else as.vector(mbest$c),
                 b0.start = mbest$b0, b.start = mbest$b,
                 sb0.fixed = a$sb0.fixed, sb.fixed = a$sb.fixed), silent = TRUE))
  if (!inherits(m, "try-error")) {
    tvtab$logLik[i + 1] <- m$logLik; tvtab$AICc[i + 1] <- aicc2(m); tvtab$sb0[i + 1] <- m$sb0
  }
}
tvtab$dAICc <- tvtab$AICc - min(tvtab$AICc, na.rm = TRUE)
print(tvtab, digits = 5, row.names = FALSE)
cat("A negative deviance against the fixed model means the optimiser could not\n")
cat("beat the constrained fit even though it is nested - read that as no support.\n")

# =============================================================================
rule("5. Sensitivity")
# =============================================================================

if (RUN_SENSITIVITY) {
  all_composite <- read_csv("./data/all_composite.csv", show_col_types = FALSE) |> as.data.frame()

  cat("\n-- (a) transform exponent: is the conclusion knife-edge on lambda = 0.25? --\n")
  # lambda is varied with tau held at 0 so the comparison isolates lambda alone.
  # Deviances are then directly comparable with the primary analysis.
  for (lam in c(0.20, 0.25, 0.30)) {
    me_l  <- build_ocfs_me(all_composite, lambda = lam, verbose = FALSE)
    dat_l <- build_ocfs_data(all_composite, me_l, lambda = lam, verbose = FALSE)
    # tau stays at 0, matching the primary analysis: this tests lambda alone.
    ME_l <- me_with_tau(dat_l$var_obs, 0)
    Ul <- make_U(best_spec, if (tn_best == "time_lin")
                   as.numeric(scale(seq_len(nrow(dat_l$X)))) else as.numeric(dat_l$bins <= 61))
    m0 <- cached(sprintf("sens_lam%03d_null", lam * 100),
                 fit_ms(dat_l$X, NULL, ME_l, p = p, su.fixed = 1))
    mF <- cached(sprintf("sens_lam%03d_top", lam * 100),
                 fit_ms(dat_l$X, Ul, ME_l, p = p, su.fixed = 1,
                        b0.start = m0$b0, b.start = m0$b))
    dv <- 2 * (mF$logLik - m0$logLik)
    cat(sprintf("  lambda = %.2f (tau = %.2f): dev = %7.3f, df = %d, P = %-9.4g  c = %s\n",
                lam, 0, dv, ncol(Ul),
                pchisq(max(dv, 0), ncol(Ul), lower.tail = FALSE),
                paste(sprintf("%+.3f", as.vector(mF$c)), collapse = " ")))
  }

  cat("\n-- (b) autoregressive order --\n")
  for (pp in 1:3) {
    m0 <- cached(sprintf("sens_p%d_null", pp), fit_ms(X, NULL, ME, p = pp, su.fixed = SU))
    mF <- cached(sprintf("sens_p%d_top", pp),
                 fit_ms(X, Ubest, ME, p = pp, su.fixed = SU, b0.start = m0$b0, b.start = m0$b))
    cat(sprintf("  p = %d: logLik = %9.4f  AICc = %8.2f  dominant eigen = %.4f  c = %s\n",
                pp, mF$logLik, aicc2(mF), dominant_eigen(mF),
                paste(sprintf("%+.3f", as.vector(mF$c)), collapse = " ")))
  }

  cat("
-- (c) without the Holocene (bins >= 54, as in tulane.R) --
")
  # HUMAN is dropped from this test BY DESIGN, not for numerical convenience.
  # PrDens is identically zero before ~13 ka, so on a pre-Holocene subset it has
  # almost no support left and cannot be estimated: an earlier version left it in
  # and the fit returned c = -126.9 with a deviance of 1.5e42. Asking whether
  # human population predicts OCFS in a window that predates the population is
  # not a meaningful question, so the test is restricted to the blocks that exist
  # throughout the record.
  keep <- dat$bins >= 54
  spec_h <- intersect(best_spec, c("CLIM", "ABRUPT", "FIRE",
                                   if (tn_best == "time_lin") "TIME"))
  dropped <- setdiff(best_spec, spec_h)
  if (length(dropped))
    cat("  excluded (undefined or unsupported outside the Holocene):",
        paste(dropped, collapse = ", "), "
")
  Xh  <- X[keep, , drop = FALSE]
  MEh <- ME[keep]
  Ubh <- if (length(spec_h)) make_U(spec_h, time_lin)[keep, , drop = FALSE] else NULL
  sdh <- sd(Xh, na.rm = TRUE)
  cat("  rows:", nrow(Xh), " observations:", sum(!is.na(Xh)),
      " predictors:", if (is.null(Ubh)) 0 else ncol(Ubh), "
")
  if (!is.null(Ubh) && ncol(Ubh) > 0) {
    obs_h <- !is.na(Xh[, 1])
    vr <- apply(Ubh[obs_h, , drop = FALSE], 2, sd)
    if (any(vr < 1e-8)) {
      cat("  dropping zero-variance column(s):",
          paste(colnames(Ubh)[vr < 1e-8], collapse = ", "), "
")
      Ubh <- Ubh[, vr >= 1e-8, drop = FALSE]
    }
    m0h <- cached("sens_woholo_null",
                  fit_ms(Xh, NULL, MEh, p = p, su.fixed = SU,
                         se.min = SE_MIN_FRAC * sdh))
    mFh <- cached("sens_woholo_top",
                  fit_ms(Xh, Ubh, MEh, p = p, su.fixed = SU,
                         se.min = SE_MIN_FRAC * sdh,
                         b0.start = m0h$b0, b.start = m0h$b))
    dv <- 2 * (mFh$logLik - m0h$logLik)
    bad <- !is.finite(dv) || dv < -1e-6 || dv > 500 ||
           max(abs(as.vector(mFh$c))) > 10 * sdh
    if (bad) {
      cat("  FIT FAILED - implausible estimates, not reporting. Check for a",
          " covariate that is near-degenerate on this subset.
")
    } else {
      cat(sprintf("  dev = %.3f, df = %d, P = %.4g   sum(b) = %.3f  se = %.4f
",
                  dv, ncol(Ubh), pchisq(dv, ncol(Ubh), lower.tail = FALSE),
                  sum(mFh$b), mFh$se))
      print(data.frame(term = colnames(Ubh), c = as.vector(mFh$c),
                       longrun = longrun(mFh)), digits = 4, row.names = FALSE)
    }
  } else cat("  nothing left to test after exclusions
")

  cat("\n-- (d) dropping the na.interp-filled predictor bins --\n")
  filled_bins <- unique(unlist(dat$filled))
  keep2 <- !(dat$bins %in% filled_bins)
  Xd <- X; Xd[!keep2] <- NA          # blank the response there rather than reindexing time
  cat("  filled bins:", length(filled_bins), "; observations lost:",
      sum(!is.na(X) & !keep2), "of", n_obs, "\n")
  m0d <- cached("sens_nofill_null", fit_ms(Xd, NULL, ME, p = p, su.fixed = SU))
  mFd <- cached("sens_nofill_top",
                fit_ms(Xd, Ubest, ME, p = p, su.fixed = SU, b0.start = m0d$b0, b.start = m0d$b))
  dv <- 2 * (mFd$logLik - m0d$logLik)
  cat(sprintf("  dev = %.3f, df = %d, P = %.4g   c = %s\n", dv, ncol(Ubest),
              pchisq(max(dv, 0), ncol(Ubest), lower.tail = FALSE),
              paste(sprintf("%+.3f", as.vector(mFd$c)), collapse = " ")))
}

# =============================================================================
rule("Done")
# =============================================================================
saveRDS(list(sets = sets, tabs = tabs, imp = imp, avg = avg, best = mbest,
             best_spec = best_spec, tn_best = tn_best),
        file.path(OUT, "tier3_selection.rds"))
cat("outputs in", OUT, "\n")
cat("\nReporting checklist:\n")
cat("  * quote the model-averaged 'full' estimates with their unconditional SEs\n")
cat("  * quote summed weights for variable importance, noting the set is balanced\n")
cat("  * state how many models fall within 2 AICc of the top one\n")
cat("  * if HUMAN loses to TIME in part 1, say so explicitly rather than quoting\n")
cat("    PrDens from a model set that never made it compete\n")
