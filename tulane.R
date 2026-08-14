if (!require("pacman")) install.packages("pacman", repos="http://cran.r-project.org")
pacman::p_load(multinomialTS, tidyverse, ggtext, patchwork, forecast,
               future, furrr, viridis, viridisLite, RColorBrewer)

packageVersion("multinomialTS")

# Data hanling ------------------------------------------------------------

# Function to re-fit the model if necessary
refit_func <- function(mod, n_refit = 10) {

  Tsample <- mod$Tsample
  Y <- mod$Y
  X <- mod$X
  p <- ncol(X) + 1 # Number of independent variables plus intercept
  n <- ncol(Y)

  ss_seq_list <- vector(mode = "list", length = n_refit)
  ss_seq_list[[1]] <- mod
  idx <- 1

  while (idx < length(ss_seq_list)) {
    
    B0.start <- ss_seq_list[[idx]]$B0
    B.start <- ss_seq_list[[idx]]$B

    sigma.start <- ss_seq_list[[idx]]$sigma

    V.fixed = matrix(NA, n, n) # Covariance matrix of environmental variation in process eq
    V.fixed[1] = 1

    V.start <- ss_seq_list[[idx]]$V

    B.fixed <- matrix(NA, ncol(X), n)
    B.fixed[,1] <- 0
    B0.fixed = matrix(c(0, rep(NA, n - 1)), nrow = 1, ncol = n)

    C.start <- ss_seq_list[[idx]]$C
    C.fixed <- C.start
    C.fixed[C.fixed != 0] <- NA
    print(C.start)
    print(C.fixed)

    ssm <- multinomialTS::mnTS(Y = Y, X = X, Tsample = Tsample, B0.start = B0.start, B.start = B.start,
                               C.start = C.start, C.fixed = C.fixed, B0.fixed = B0.fixed,
                               V.fixed = V.fixed, V.start = V.start,
                               B.fixed = B.fixed, dispersion.fixed = 1, maxit.optim = 1e+07)
    print(ssm$AIC)
    idx <- idx + 1
    print(idx)
    ss_seq_list[[idx]] <- ssm

  }
  return(ss_seq_list)
}

# necessary as global variable
bin_width = 200

# Download and save Grimm 06 core from Neotoma
# grimm_06_pollen <- get_sites(siteid = 2570) |> 
#   get_datasets() |>
#   neotoma2::filter(datasettype == "pollen" & !is.na(age_range_young)) |> 
#   get_downloads() |> 
#   samples()
# 
# saveRDS(grimm_06_pollen, "./data/data_230814/grimm_06_pollen.rds")

# Load pollen data
grimm_06_pollen <- readRDS("./data/tula94_pollen_files/grimm_06_pollen.rds")
nora_pollen_wide <- read_csv("./data/tula94_pollen_files/tula94_pollen_wide_count.csv")

# Tulane 2020 core
chron <- read_csv("./data/TULA20_age-depth_files/TULA20_compsiteCore_w-ages.csv")
core_20_char <- read_csv("./data/TULA20_age-depth_files/TULA20_CHAR_w-ages.csv")
# core_20_spore <- read_csv("./data/TULA20_age-depth_files/TULA20_CFS_w-ages.csv")
core_20_spore <- read_csv("./data/TULA20_age-depth_files/ocfs_uncertainty.csv")
# Isotope cores
co2 <- read_csv("./data/isotope/co2_merged.csv")
oxy18 <- read_csv("./data/isotope/ngrip-d18o-50yr.csv")

chron <- chron |> 
  select(depth_composite, median_age, quant_5perc, quant_95perc)
colnames(chron) <- c("depth", "cov_age", "quant_5perc", "quant_95perc")

core_20_char <- core_20_char |> 
  select(depth_composite, accu_rate)
colnames(core_20_char) <- c("depth", "char_acc")

core_20_spore <- core_20_spore |> 
  select(OCFS_conc, depth_core)
colnames(core_20_spore) <- c("ocfs", "depth")

humans <- read_csv("data/humans-kelly2025/kelly2025_recalcSPD_HUC3.csv") |>
  select(-c(`...1`)) |>
  arrange(desc(calBP))

# Handling state-variables
# Variations in tree cover in North America since the LGM
# graminoids/grass
# Herbs/Forbs
herbs <- c("Apiaceae.*|Ambrosia.*|Artemisia.*|Asteraceae.*|
            Tubuliflorae.*|Brassicaceae.*|Caryophyllaceae.*|
            Chenopodiaceae.*|Amaranthaceae.*|Dryas.*|Ephedra.*|
            Eriogonum.*|Euphorbiaceae.*|Oxyria.*|Ranunculaceae.*|
            Sarcobatus.*|Saxifragaceae.*")

grass <- c("Cyperaceae.*|Poaceae.*")

herb_forb <- grimm_06_pollen |> 
  dplyr::filter(datasetid == 19620,
                stringr::str_detect(variablename, herbs)) |> 
  mutate(variablename = "Herbs") |> 
  group_by(siteid, sitename,
           sampleid, variablename, units, age,
           agetype, depth, datasetid,
           long, lat) |>
  summarise(value = sum(value), .groups='keep')

gram_grass <- grimm_06_pollen |> 
  dplyr::filter(datasetid == 19620,
                stringr::str_detect(variablename, grass)) |> 
  mutate(variablename = "Grass") |> 
  group_by(siteid, sitename,
           sampleid, variablename, units, age,
           agetype, depth, datasetid,
           long, lat) |>
  summarise(value = sum(value), .groups='keep')

grimm_06_tree <- grimm_06_pollen |> 
  dplyr::filter(datasetid == 19620,
                ecologicalgroup %in% c("UPHE", "TRSH"),
                stringr::str_detect(variablename, c("Pinus.*|Quercus.*"))
  ) |>
  mutate(variablename = replace(variablename, stringr::str_detect(variablename, "Pinus.*"), "Pinus"),
         variablename = replace(variablename, stringr::str_detect(variablename, "Picea.*"), "Picea")) |>
  group_by(siteid, sitename,
           sampleid, variablename, units, age,
           agetype, depth, datasetid,
           long, lat) |>
  summarise(value = sum(value), .groups='keep')

other_spp <- grimm_06_pollen |> 
  dplyr::filter(datasetid == 19620,
                ecologicalgroup %in% c("UPHE", "TRSH"),
                !stringr::str_detect(variablename, herbs),
                !stringr::str_detect(variablename, grass),
                !stringr::str_detect(variablename, c("Pinus.*|Quercus.*"))) |>
  mutate(variablename = "other") |> 
  group_by(siteid, sitename,
           sampleid, variablename, units, age,
           agetype, depth, datasetid,
           long, lat) |>
  summarise(value = sum(value), .groups='keep')

grimm_06_groups <- bind_rows(other_spp, gram_grass, herb_forb, grimm_06_tree)

pollen_wide <- grimm_06_groups |>
  pivot_wider(id_cols = c(depth, age), names_from = variablename, values_from = value) |> 
  mutate(across(c(other, Grass, Herbs, Pinus, Quercus), ~ replace_na(., 0))) |> 
  rename("grimm_age" = "age") |> 
  rename("grimm_depth" = "depth") |> 
  arrange(grimm_age) |>
  ungroup() |> 
  mutate(age = nora_pollen_wide$new_age) |> 
  arrange(desc(age))

# write_csv(pollen_wide, "./data/pollen_wide.csv")

# nora_poll <-  read_csv("./data/tula94_pollen_files/tula94_pollen.csv")
# 
# nora_poll <- nora_poll |>
#   pivot_wider(id_cols = c(depth, age, new_age), names_from = variablename, values_from = value) |> 
#   select(depth, age, new_age) |> 
#   rename("nora_new_age" = "new_age") |> 
#   rename("nora_old_age" = "age") |> 
#   rename("nora_depth" = "depth") |> 
#   arrange(desc(nora_depth))
# 
# check_ages <- full_join(pollen_wide, nora_poll, by = c("grimm_depth" = "nora_depth"))

# Plot counts
# Yes, they counted >500 quercus in that spike... crazy
pollen_wide |>
  pivot_longer(-c(age, grimm_depth, grimm_age)) |> 
  ggplot(aes(x = age, y = value)) +
    geom_area(fill = "grey20") +
    geom_segment(aes(x = age, xend = age,
                     y = 0, yend = value), colour = "grey30", linewidth = 0.6) +
    scale_x_reverse(breaks = scales::breaks_pretty(n = 6)) +
    coord_flip() +
    labs(y = "Pollen counts", x = "Time (ybp)") +
    facet_wrap(~name, nrow = 1) +
    theme_minimal() +
    theme(
      text = element_text(size = 9),
    )

bins <- cut(pollen_wide$age,
            breaks = seq(from = min(pollen_wide$age), 
                         to = max(pollen_wide$age + bin_width), 
                         by = bin_width),
            include.lowest = T, labels = F)

pollen_wide_binned <- bind_cols(bins = bins, pollen_wide) |> 
  group_by(bins) |> 
  summarise(
    age = mean(age, na.rm = T),
    other = sum(other, na.rm = T),
    Grass = sum(Grass, na.rm = T),
    Herbs = sum(Herbs, na.rm = T),
    Pinus = sum(Pinus, na.rm = T),
    Quercus = sum(Quercus, na.rm = T)) |> 
  arrange(desc(bins))

dim(pollen_wide)
dim(pollen_wide_binned)

# Handling covariates
composite_covariate_join <- chron |> 
  full_join(core_20_char, by = "depth") |> 
  full_join(core_20_spore, by = "depth") |> 
  mutate(heinrich = NA,
         heinrich = case_when(cov_age <= 62400 & cov_age >= 59700 ~ 1,
                              cov_age <= 49300 & cov_age >= 47600 ~ 1,
                              cov_age <= 40200 & cov_age >= 38300 ~ 1,
                              cov_age <= 31300 & cov_age >= 30000 ~ 1,
                              cov_age <= 24700 & cov_age >= 23400 ~ 1,
                              cov_age <= 18300 & cov_age >= 15100 ~ 1,
                              cov_age <= 12900 & cov_age >= 11600 ~ 1,
                              .default = NA)) |>
  arrange(desc(cov_age))
# Binary variables left as NA for summarise later
dim(composite_covariate_join)
# write_csv(composite_covariate_join, "./data/composite_covariate_join.csv")

composite_covariate_join |> 
  select(-c(quant_5perc, quant_95perc)) |> 
  pivot_longer(-c(depth, cov_age)) |> 
  ggplot(aes(x = cov_age, y = value)) +
  geom_point() +
  geom_line() +
  scale_x_reverse(breaks = scales::breaks_pretty(n = 6)) +
  coord_flip() +
  labs(y = "Values", x = "Time (ybp)") +
  facet_wrap(~name, nrow = 1, scales = "free") +
  theme_minimal() +
  theme(
    text = element_text(size = 9),
  )


cov_bins <- cut(composite_covariate_join$cov_age,
                breaks = seq(from = min(pollen_wide_binned$age), 
                             to = max(pollen_wide_binned$age + bin_width), 
                             by = bin_width),
                include.lowest = T, labels = F)


composite_covariate_join_bin <- bind_cols(bins = cov_bins, composite_covariate_join) |> 
  drop_na(bins) |> # we lose one ocfs observation at 64720
  group_by(bins) |> 
  summarise(
    cov_age = mean(cov_age, na.rm = T),
    char_acc = mean(char_acc, na.rm = T),
    heinrich = mean(heinrich, na.rm = T),
    ocfs = mean(ocfs, na.rm = T)) |> 
  mutate(heinrich = ifelse(is.nan(heinrich), 0, heinrich)) |>
  arrange(desc(cov_age))


# CO2
co2 <- co2[nrow(co2):1, , drop = F]

co2_bins <- cut(co2$age_calBP, 
                breaks = seq(from = min(pollen_wide_binned$age), 
                             to = max(pollen_wide_binned$age + bin_width), 
                             by = bin_width),
                include.lowest = T, labels = F)

co2_mean <- cbind(bins = co2_bins, co2) |> 
  drop_na(bins) |> 
  select(bins, age_calBP, CO2_blank_gravity_corrected) |> 
  group_by(bins) |> 
  summarise(mean_age_co2 = mean(age_calBP),
            mean_co2 = mean(CO2_blank_gravity_corrected)) |> 
  arrange(desc(mean_age_co2))

# Oxy18
oxy18_mean <- oxy18 |> 
  rename(age = Age) |> 
  group_by(age) |> 
  summarise(d18O = mean(d18O)) |> 
  arrange(desc(age))

oxy_bins <- cut(oxy18_mean$age, breaks = seq(from = min(pollen_wide_binned$age), 
                                             to = max(pollen_wide_binned$age + bin_width), 
                                             by = bin_width),
                include.lowest = T, labels = F)

oxy18_mean <- cbind(bins = oxy_bins, oxy18_mean) |> 
  drop_na(bins) |> 
  group_by(bins) |> 
  summarise(mean_age_oxy = mean(age),
            d18O = mean(d18O)) |> 
  arrange(desc(mean_age_oxy))

###

human_bins <- cut(humans$calBP,
            breaks = seq(from = min(pollen_wide$age), 
                         to = max(pollen_wide$age + bin_width), 
                         by = bin_width),
            include.lowest = T, labels = F)


humans_binned <- bind_cols(human_bins = human_bins, humans) |> 
  group_by(human_bins) |> 
  summarise(
    calBP = mean(calBP, na.rm = T),
    PrDens = mean(PrDens, na.rm = T)) |> 
  arrange(desc(human_bins))

###

all_composite <- pollen_wide_binned |>
  full_join(composite_covariate_join_bin, by = "bins") |>
  full_join(co2_mean, by = "bins") |>
  full_join(oxy18_mean, by = "bins") |>
  left_join(humans_binned |> select(-c(calBP)), by = c("bins" = "human_bins")) |>
  mutate(PrDens = ifelse(is.na(PrDens), 0, PrDens)) |>
  arrange(desc(bins))

# write_csv(all_composite, "./data/all_composite.csv")

# all_composite <- pollen_wide_binned |>
#   full_join(composite_covariate_join_bin, by = "bins") |>
#   full_join(co2_mean, by = "bins") |>
#   full_join(oxy18_mean, by = "bins") |>
#   arrange(desc(bins))

# all_composite |> select(bins, age, ocfs) |> 
#   arrange(bins) |> 
#   print(n = 50)
# 
# all_composite |> select(bins, age, ocfs) |> 
#   print(n = 50)
# 
# all_composite |> select(bins, age, ocfs) |> 
#   arrange(bins) |> 
#   mutate(across(c(ocfs), forecast::na.interp)) |> 
#   print(n = 20) |> 
#   tail(n = 20)
# 
# all_composite |> select(bins, age, ocfs) |> 
#   arrange(desc(bins)) |> 
#   mutate(across(c(ocfs), forecast::na.interp)) |> 
#   print(n = 20) |> 
#   tail(n = 20)


all_composite |> 
  select(bins, char_acc, ocfs, d18O, heinrich, mean_co2, PrDens) |> 
  pivot_longer(-c(bins)) |> 
  ggplot(aes(x = bins, y = value)) +
  geom_point() +
  geom_line() +
  facet_wrap(~name, scales = "free", ncol = 1)

all_composite |> 
  select(bins, char_acc, ocfs, d18O, heinrich, mean_co2, PrDens) |> 
  mutate(across(c(char_acc, ocfs, d18O, mean_co2), forecast::na.interp)) |> 
  pivot_longer(-c(bins)) |> 
  ggplot(aes(x = bins, y = value)) +
  geom_point() +
  geom_line() +
  facet_wrap(~name, scales = "free", ncol = 1)

all_composite |>
  select(bins, char_acc, ocfs, d18O, heinrich, mean_co2, PrDens) |> 
  mutate(across(c(char_acc, ocfs, d18O, mean_co2), forecast::na.interp)) |> 
  mutate(across(-c(heinrich, bins), ~ as.numeric(scale(.)))) |> 
  pivot_longer(-c(bins)) |> 
  ggplot(aes(x = bins, y = value)) +
  geom_point() +
  geom_line() +
  facet_wrap(~name, scales = "free", ncol = 1)

# Try sqrt transform?
plot(all_composite$ocfs)
plot(log(all_composite$ocfs))
plot(sqrt(all_composite$ocfs))
plot(scale(forecast::na.interp(all_composite$ocfs)))

all_composite |>
  select(bins, ocfs, age) |> 
  mutate(across(c(ocfs), forecast::na.interp)) |>
  arrange(age) |> 
  ggplot(aes(x = bins, y = ocfs)) +
    geom_point() +
    geom_line()

x <- sqrt(all_composite$ocfs)
x <- forecast::na.interp(x)
plot(scale(x))

# all_composite |> 
#   select(bins, char_acc, ocfs, d18O, humans, heinrich, mean_co2) |> 
#   mutate(ocfs = log(ocfs),
#   ocfs = ifelse(ocfs == -Inf, 0, ocfs)) |>
#   mutate(across(c(char_acc, ocfs, d18O, mean_co2), forecast::na.interp)) |>
#   pivot_longer(-c(bins)) |> 
#   ggplot(aes(x = bins, y = value)) +
#   geom_point() +
#   geom_line() +
#   facet_wrap(~name, scales = "free", ncol = 1)

# X_co2 <- data.frame(y = co2$CO2_blank_gravity_corrected, x = co2$age_calBP)
# X_co2_gam <- mgcv::gam(y ~ s(x, bs = "bs", k = nrow(X_co2)), method = "REML", data =  X_co2)
# pred_co2 <- predict(X_co2_gam, newdata = data.frame(x = oxy_co2$mean_age_oxy[which(is.na(oxy_co2$mean_co2))] ))
# 
# oxy_co2$mean_co2[which(is.na(oxy_co2$mean_co2))] <- pred_co2

# X_char <- data.frame(y = core_20_char$char_acc, x = core_20_char$age)
# X_char_gam <- mgcv::gam(y ~ s(x, bs = "bs", k = nrow(X_char)), method = "REML", data =  X_char[nrow(X_char):1, , drop = F])
# pred <- predict(X_char_gam, newdata = data.frame(x = char_pollen_ocfs_bin$age[which(is.na(char_pollen_ocfs_bin$char_acc))] ))
# char_pollen_ocfs_bin$char_acc[which(is.na(char_pollen_ocfs_bin$char_acc))] <- pred

# pred_ofc <- predict(X_ofc_gam, newdata = data.frame(x = char_pollen_ocfs_bin$age[which(is.nan(char_pollen_ocfs_bin$ocfs))] ))
# ofc_pred <- as_tibble(char_pollen_ocfs_bin)
# # ofc_pred$ocfs <- forecast::na.interp(ofc_pred$ocfs)
# ofc_pred$ocfs[is.na(ofc_pred$ocfs)] <- pred_ofc

# Model scenarios ---------------------------------------------------------

# Scenario design: humans (PrDens in/out) x Holocene (in/out) = 4 datasets.
# Each dataset is fit twice, without and with the species interactions in the
# C matrix (Pinus <-> Quercus), so 8 models in total. Every model is refitted,
# bootstrapped, and written to its own folder under ./results_scenarios.

# Covariates common to every scenario. PrDens is added for the "humans" runs.
base_covariates <- c("char_acc", "ocfs", "d18O", "heinrich", "mean_co2")

# Bin 54 is ~11,732 ybp. Everything younger is dropped for the "without
# Holocene" scenarios.
holocene_bin_cutoff <- 54

# Covariates that are interpolated before scaling (heinrich is binary)
interp_covariates <- c("char_acc", "ocfs", "d18O", "mean_co2")

# Note on humans x without-Holocene: the Kelly 2025 SPD only reaches 13,075 ybp,
# so in that scenario PrDens is non-zero in 8 of 256 bins. The model will run,
# but treat the PrDens coefficient there as near-unidentifiable.

# set up Y
make_Y <- function(dat) {
  dat |>
    select(bins, other, Grass, Herbs, Pinus, Quercus) |>
    arrange(desc(bins)) |>
    select(-bins) |>
    as.matrix()
}

# set up X
make_X <- function(dat, covariates) {
  dat |>
    select(bins, all_of(covariates)) |>
    arrange(desc(bins)) |>
    mutate(across(any_of("ocfs"), sqrt)) |>
    mutate(across(any_of(interp_covariates), forecast::na.interp),
           across(!any_of(c("heinrich", "bins")), ~ as.numeric(scale(.)))) |>
    select(-bins) |>
    as.matrix()
}

# Fit and refit the two models (no interactions, with interactions) for one
# scenario and save them. Sequential by design: the interaction model starts
# from the no-interaction fit, and refit_func chains each refit onto the
# previous one, so there is nothing here to run in parallel.
# refit_index = 1 is the original fit (what the previous version of this
# script bootstrapped); raise it to carry a later refit forward instead.
# interactions: list of C-matrix off-diagonal pairs to free. c(4, 5) is
# Pinus <-> Quercus given the Y column order above. Pass list() for a purely
# diagonal C in the "with interactions" model.
fit_scenario <- function(dat, covariates, out_dir,
                         n_refit = 5, refit_index = 1,
                         interactions = list(c(4, 5))) {

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  message("\n===== fitting ", out_dir, " =====")
  message("covariates: ", paste(covariates, collapse = ", "),
          " | bins: ", nrow(dat))

  Y <- make_Y(dat)
  X <- make_X(dat, covariates)
  Tsample <- which(rowSums(Y) != 0)

  stopifnot(!anyNA(X), nrow(X) == nrow(Y))
  matplot(X, type = 'l', main = basename(out_dir))

  p <- ncol(X) + 1 # Number of independent variables plus intercept
  n <- ncol(Y)

  # GLMM for starting values
  V.fixed = diag(n) # Covariance matrix of environmental variation in process eq

  B.fixed <- matrix(c(rep(0, p), rep(NA, (n - 1) * p)), p, n)
  B.start <- matrix(c(rep(0, p), rep(.01, (n - 1) * p)), p, n)

  glmm_mod <- multinomialTS::mnGLMM(Y = Y[Tsample, ],
                                    X = X[Tsample, , drop = F],
                                    B.start = B.start, B.fixed = B.fixed,
                                    V.fixed = V.fixed)
  print(summary(glmm_mod))

  B0.start <- glmm_mod$B[1, , drop = F]
  B.start <- glmm_mod$B[2:p, , drop = F]

  sigma.start <- glmm_mod$sigma

  V.fixed = matrix(NA, n, n) # Covariance matrix of environmental variation in process eq
  V.fixed[1] = 1

  V.start <- glmm_mod$V

  B.fixed <- matrix(NA, ncol(X), n)
  B.fixed[, 1] <- 0
  B0.fixed = matrix(c(0, rep(NA, n - 1)), nrow = 1, ncol = n)

  # Set-up C without interactions
  C.start.diag = .5 * diag(n)
  C.fixed.diag <- C.start.diag
  C.fixed.diag[C.fixed.diag != 0] <- NA

  # Set-up C with interactions
  C.start.diag.int = .5 * diag(n)
  for (pair in interactions) {
    C.start.diag.int[pair[1], pair[2]] = C.start.diag.int[pair[2], pair[1]] = .5
  }
  C.fixed.diag.int <- C.start.diag.int
  C.fixed.diag.int[C.fixed.diag.int != 0] <- NA

  # Model with no interactions
  message("fitting mnTS, no interactions")
  start_time <- Sys.time()
  mnTS_mod <- mnTS(Y = Y[Tsample, ],
                   X = X, Tsample = Tsample,
                   B0.start = B0.start, B0.fixed = B0.fixed,
                   B.start = B.start, B.fixed = B.fixed,
                   C.start = C.start.diag, C.fixed = C.fixed.diag,
                   V.start = V.start, V.fixed = V.fixed,
                   dispersion.fixed = 1, maxit.optim = 1e+6)
  message("  ", format(Sys.time() - start_time))

  mnTS_mod_refit <- refit_func(mnTS_mod, n_refit)
  saveRDS(mnTS_mod_refit, file.path(out_dir, "mnTS_mod_refit.rds"))

  # Model with interactions, started from the no-interaction fit
  message("fitting mnTS, with interactions")
  start_time <- Sys.time()
  mnTS_mod_int <- mnTS(Y = Y[Tsample, ],
                       X = X, Tsample = Tsample,
                       B0.start = mnTS_mod$B0, B0.fixed = B0.fixed,
                       B.start = mnTS_mod$B, B.fixed = B.fixed,
                       C.start = C.start.diag.int, C.fixed = C.fixed.diag.int,
                       V.start = mnTS_mod$V, V.fixed = V.fixed,
                       dispersion.fixed = 1, maxit.optim = 1e+6)
  message("  ", format(Sys.time() - start_time))

  mnTS_mod_int_refit <- refit_func(mnTS_mod_int, n_refit)
  saveRDS(mnTS_mod_int_refit, file.path(out_dir, "mnTS_mod_int_refit.rds"))

  # Original fit plus each refit, to check the optimiser settled
  aic_tbl <- bind_rows(
    tibble(model = "mnTS_mod",
           refit = seq_along(mnTS_mod_refit),
           logLik = map_dbl(mnTS_mod_refit, "logLik"),
           AIC = map_dbl(mnTS_mod_refit, "AIC"),
           opt.convergence = map_dbl(mnTS_mod_refit, "opt.convergence")),
    tibble(model = "mnTS_mod_int",
           refit = seq_along(mnTS_mod_int_refit),
           logLik = map_dbl(mnTS_mod_int_refit, "logLik"),
           AIC = map_dbl(mnTS_mod_int_refit, "AIC"),
           opt.convergence = map_dbl(mnTS_mod_int_refit, "opt.convergence"))
  )
  write_csv(aic_tbl, file.path(out_dir, "model_aic.csv"))
  print(aic_tbl)

  # Coefficient tables of the models that get bootstrapped
  walk2(list(mnTS_mod_refit[[refit_index]], mnTS_mod_int_refit[[refit_index]]),
        c("mnTS_mod", "mnTS_mod_int"),
        \(mod, nm) {
          as_tibble(coef(mod), rownames = "cov") |>
            write_csv(file.path(out_dir, paste0(nm, "_coef.csv")))
        })

  invisible(list(mnTS_mod_refit = mnTS_mod_refit,
                 mnTS_mod_int_refit = mnTS_mod_int_refit,
                 Y = Y, X = X, Tsample = Tsample,
                 covariates = covariates,
                 out_dir = out_dir))
}

# Read a scenario's fits back off disk, so bootstrapping (or re-bootstrapping
# with more replicates) can be run in a session that did not do the fitting.
load_scenario_fits <- function(out_dir) {
  list(mnTS_mod_refit = readRDS(file.path(out_dir, "mnTS_mod_refit.rds")),
       mnTS_mod_int_refit = readRDS(file.path(out_dir, "mnTS_mod_int_refit.rds")),
       out_dir = out_dir)
}

# Bootstrap one scenario. This is the expensive step, so this is where the
# parallelism goes.
#
# boot_reps is the total number of replicates per model. It is split into
# n_chunks jobs of boot_reps / n_chunks replicates, and both models (no
# interactions, with interactions) are queued together, so 2 * n_chunks jobs
# share the workers. More, smaller chunks balance the load better across
# workers; fewer, larger chunks carry less per-job overhead. The chunks are
# kept as separate list elements named mnTS_mod<i> / mnTS_mod_int<i> and are
# pooled by the plotting code below.
#
# workers = NULL leaves whatever future plan is already set alone; otherwise a
# multisession plan is set here and the previous plan restored on exit.
bootstrap_scenario <- function(fits, out_dir = fits$out_dir,
                               boot_reps = 1000, n_chunks = 10, workers = 10,
                               refit_index = 1, keep_all_mods = FALSE,
                               seed = 1984) {

  stopifnot(boot_reps %% n_chunks == 0)
  reps_per_chunk <- boot_reps %/% n_chunks

  if (!is.null(workers)) {
    oplan <- future::plan(strategy = multisession, workers = workers)
    on.exit(future::plan(oplan), add = TRUE)
  }

  mods <- c(setNames(rep(list(fits$mnTS_mod_refit[[refit_index]]), n_chunks),
                     paste0("mnTS_mod", seq_len(n_chunks))),
            setNames(rep(list(fits$mnTS_mod_int_refit[[refit_index]]), n_chunks),
                     paste0("mnTS_mod_int", seq_len(n_chunks))))

  message("\n===== bootstrapping ", out_dir, " =====")
  message(length(mods), " jobs x ", reps_per_chunk, " reps on ",
          future::nbrOfWorkers(), " workers (", boot_reps, " reps per model)")
  start_time <- Sys.time()

  res <- furrr::future_map(mods, multinomialTS::bootstrap,
                           reps = reps_per_chunk,
                           .options = furrr_options(seed = seed))

  # all_mods holds every bootstrap refit, ~200 MB per scenario. The summaries
  # and plots below only need coef.table.boot and all_mods_pars, so drop it
  # unless keep_all_mods = TRUE.
  if (!keep_all_mods) {
    res <- map(res, ~ .x[c("coef.table.boot", "all_mods_pars")])
  }
  saveRDS(res, file.path(out_dir, "bootstraps.rds"))
  message("  ", format(Sys.time() - start_time))

  invisible(res)
}

# The four scenarios ------------------------------------------------------

scenarios <- tibble::tribble(
  ~scenario,             ~humans, ~holocene,
  "humans_holocene",     TRUE,    TRUE,
  "humans_woholo",       TRUE,    FALSE,
  "nohumans_holocene",   FALSE,   TRUE,
  "nohumans_woholo",     FALSE,   FALSE
)

results_root <- "./results_scenarios"
dir.create(results_root, showWarnings = FALSE)
write_csv(scenarios, file.path(results_root, "scenarios.csv"))

scenario_dirs <- set_names(file.path(results_root, scenarios$scenario),
                           scenarios$scenario)

scenario_data <- map(scenarios$scenario, \(s) {
  scen <- scenarios[scenarios$scenario == s, ]
  list(
    dat = if (scen$holocene) {
      all_composite
    } else {
      all_composite |> filter(bins >= holocene_bin_cutoff)
    },
    covariates = if (scen$humans) c(base_covariates, "PrDens") else base_covariates,
    out_dir = scenario_dirs[[s]]
  )
}) |> set_names(scenarios$scenario)

scenario_fits <- vector(mode = "list", length = nrow(scenarios))
names(scenario_fits) <- scenarios$scenario
scenario_boot <- scenario_fits

## Stage 1: fitting, scenarios in sequence, single core --------------------

# A single mnTS fit takes roughly an hour on this data, and n_refit = 5 means
# ten fits per scenario, so stage 1 is a long job. Two knobs:
#
# refit_existing = FALSE picks up any scenario already written to disk instead
# of refitting it, so an interrupted run can be resumed.
#
# fit_workers > 1 runs the scenarios' fits concurrently, one core each. It
# does not conflict with stage 2, which is sequential across scenarios and
# uses the cores there. Worker messages are not streamed, so progress is only
# visible in the saved files until each scenario finishes.
refit_existing <- FALSE
fit_workers <- 1

fit_one <- function(d) {
  saved <- file.path(d$out_dir, c("mnTS_mod_refit.rds", "mnTS_mod_int_refit.rds"))
  if (!refit_existing && all(file.exists(saved))) {
    message("using saved fits in ", d$out_dir)
    return(load_scenario_fits(d$out_dir))
  }
  fit_scenario(dat = d$dat, covariates = d$covariates, out_dir = d$out_dir)
}

start_time_all <- Sys.time()
if (fit_workers > 1) {
  oplan <- future::plan(strategy = multisession,
                        workers = min(fit_workers, nrow(scenarios)))
  scenario_fits <- furrr::future_map(scenario_data, fit_one,
                                     .options = furrr_options(seed = 1984))
  future::plan(oplan)
} else {
  for (scen in scenarios$scenario) {
    scenario_fits[[scen]] <- fit_one(scenario_data[[scen]])
  }
}
Sys.time() - start_time_all

## Stage 2: bootstrapping, scenarios in sequence, replicates in parallel ---

# Fitting done in an earlier session? Pick the models back up with:
# scenario_fits <- map(scenario_dirs, load_scenario_fits)

start_time_all <- Sys.time()
for (scen in scenarios$scenario) {
  scenario_boot[[scen]] <- bootstrap_scenario(
    scenario_fits[[scen]],
    boot_reps = 1000,  # per model, matching the old 5 x 200
    n_chunks = 10,     # 20 jobs in total (2 models x 10 chunks)
    workers = 10
  )
}
Sys.time() - start_time_all

# Plotting ----------------------------------------------------------------

# Wald coefficients of the eight fitted models (4 scenarios x 2 C matrices)
ssms <- unlist(map(scenario_fits, \(f) {
  list(mnTS_mod = f$mnTS_mod_refit[[1]],
       mnTS_mod_int = f$mnTS_mod_int_refit[[1]])
}), recursive = FALSE)

wald <- lapply(ssms, \(hyp) {
  wald <- coef(hyp)
  as_tibble(wald, rownames = "cov")
})
wald_bind <- bind_rows(wald, .id = "model") |>
  separate_wider_delim(cols = model, delim = ".", names = c("scenario", "hyp")) |>
  mutate(sig = P < 0.05,
         hyp = forcats::fct(hyp),
         scenario = fct(scenario, levels = scenarios$scenario))

wald_bind_x <- wald_bind |>
  filter(!grepl("sp.|^Grass|^Herbs|^Pinus|^Quercus", cov))


B_plot <- ggplot(wald_bind_x, aes(x = hyp, y = Coef., colour = as_factor(sig),
                                  shape = scenario)) +
  geom_point(position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = Coef. - se, ymax = Coef. + se),
                position = position_dodge(width = 0.6)) +
  geom_hline(yintercept = 0) +
  scale_color_manual(name = "Significance", labels = c("> 0.05", "< 0.05"),
                     values = c("#BF0606", "#5ab4ac")) +
  scale_shape_manual(name = "Scenario", values = c(17, 19, 2, 1)) +
  labs(x = "Taxa", y = "Coefficient") +
  facet_wrap(~cov) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90),
    legend.position = "bottom"
  )

B_plot

wald_bind_c <- wald_bind |>
  filter(grepl("sp.", cov))

C_plot <- ggplot(wald_bind_c, aes(x = as_factor(cov), y = Coef., colour = as_factor(sig))) +
  geom_point() +
  geom_errorbar(aes(ymin = Coef. - se, ymax = Coef. + se)) +
  geom_hline(yintercept = 0) +
  scale_color_manual(name = "Significance", labels = c("> 0.05", "< 0.05"),
                     values = c("#BF0606", "#5ab4ac")) +
  labs(x = "Taxa", y = "Coefficient") +
  facet_grid(scenario ~ hyp) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90),
    legend.position = "bottom"
  )

C_plot


# ggplot(composite_join_pollen |> filter(name == "Pinus"), aes(x = age, y = value)) +
#   geom_area(colour = "grey90") +
#   geom_col() +
#   geom_vline(xintercept = c(62400,
#                             49300,
#                             40200,
#                             31300,
#                             24700,
#                             18300,
#                             12900,
#                             59700,
#                             47600,
#                             36800,
#                             30000,
#                             23400,
#                             15100,
#                             11000), colour = "red") +
#   scale_x_reverse(breaks = scales::breaks_pretty(n = 6)) +
#   # coord_flip() +
#   # ylim(0, 0.5) +
#   labs(y = "Pollen counts", x = "Time (ybp)") +
#   # facet_wrap(~variablename,
#   #            nrow = 1) +
#   theme_minimal() +
#   theme(
#     text = element_text(size = 10),
#   )

## Bootstrap plotting -----------------------------------------------------
# Runs from the saved .rds files, so this section works in a fresh session
# without refitting the models.

results_root <- "./results_scenarios"

scenario_labels <- c(
  humans_holocene   = "Humans, with Holocene",
  humans_woholo     = "Humans, without Holocene",
  nohumans_holocene = "No humans, with Holocene",
  nohumans_woholo   = "No humans, without Holocene"
)

hyp_labels <- c(
  mnTS_mod = "Without species interaction",
  mnTS_mod_int = "With species interaction"
)

X_names_list <- c(
  heinrich ="A: Heinrich events",
  d18O = "B: &delta;<sup>18</sup>O",
  mean_co2 ="C: CO<sub>2</sub>",
  char_acc ="D: Charcoal accumulation",
  ocfs ="E: Fungal spores",
  PrDens ="F: Human population density"
)

boot_res <- map(set_names(names(scenario_labels)),
                ~ readRDS(file.path(results_root, .x, "bootstraps.rds")))

# check on convergence
imap(boot_res, \(res, scen) {
  tibble(scenario = scen,
         hyp = str_remove(names(res), pattern = '[[:digit:]]+'),
         n_fail_converged = map_dbl(res, ~ sum(.x$all_mods_pars[, "opt.convergence"] != 0)),
         n_total = map_dbl(res, ~ nrow(.x$all_mods_pars)))
}) |>
  bind_rows() |>
  print(n = Inf)

mods_boot <- imap(boot_res, \(res, scen) {
  map(res, ~ {
    as_tibble(.x$all_mods_pars) |>
      pivot_longer(-c(logLik, opt.convergence))
  }) |>
    bind_rows(.id = "hyp") |>
    mutate(hyp = str_remove(hyp, pattern = '[[:digit:]]+'))
}) |>
  bind_rows(.id = "scenario") |>
  mutate(scenario = fct(scenario, levels = names(scenario_labels)))

mods_boot_68 <- mods_boot |>
#  filter(opt.convergence == 0) |>
  group_by(scenario, hyp, name) |>
  summarise(boot_mean = mean(value),
            boot_sd = sd(value),
            upper_68 = quantile(value, probs = 0.84),
            lower_68 = quantile(value, probs = 0.16),
            .groups = "drop") |>
  mutate(t_scores = boot_mean / boot_sd,
         p_vals = 2 * pnorm(q = abs(t_scores), lower.tail = F),
         sig = p_vals < 0.05)

mods_boot_table <- mods_boot_68 |>
  mutate(name = str_replace_all(name,
    pattern = "y1|y2|y3|y4|y5",
    replacement = function(x) case_when(
      x == "y1" ~ "Other",
      x == "y2" ~ "Grass",
      x == "y3" ~ "Herbs",
      x == "y4" ~ "Pinus",
      x == "y5" ~ "Quercus",
      TRUE ~ x  # Keep other values unchanged
    ))) |>
    filter(!str_detect(name, "v."))

write_csv(mods_boot_table, file.path(results_root, "boot_coef_allmods.csv"))

# same table split back out into each scenario's folder
walk(names(scenario_labels), \(scen) {
  mods_boot_table |>
    filter(scenario == scen) |>
    write_csv(file.path(results_root, scen, "boot_coef.csv"))
})

mods_boot_68_B <- mods_boot_68 |>
  filter(grepl(paste(names(X_names_list), collapse = "|"), name)) |>
  separate_wider_delim(cols = name, delim = ".", names = c("cov", "name")) |>
  mutate(name = str_replace_all(name,
    pattern = "y2|y3|y4|y5",
    replacement = function(x) case_when(
      x == "y2" ~ "Grass",
      x == "y3" ~ "Herbs",
      x == "y4" ~ "_Pinus_",
      x == "y5" ~ "_Quercus_",
      TRUE ~ x  # Keep other values unchanged
    )),
        name = fct(name, levels = c("Grass", "Herbs", "_Pinus_", "_Quercus_")),
        cov = fct(cov, levels = c("heinrich", "d18O", "mean_co2", "char_acc",
                                  "ocfs", "PrDens")))

# One scenario per plot
plot_boot_B <- function(dat, title = NULL) {
  ggplot(dat, aes(x = name, y = boot_mean, colour = as_factor(sig))) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.2) +
    geom_errorbar(aes(ymin = lower_68, ymax = upper_68),
                  width = .4, alpha = 0.5) +
    scale_color_manual(name = "Significance", labels = c("> 0.05", "< 0.05"),
                       values = c("#202020", "#d80000")) +
    labs(x = "Taxa", y = "MultinomialTS coefficient estimate", title = title) +
    facet_wrap(~ cov, labeller = as_labeller(X_names_list)) +
    theme_bw() +
    theme(
      strip.text = element_markdown(size = 9),
      strip.background = element_rect(fill = NA),
      legend.position = "inside",
      legend.position.inside = c(.09, .92),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.background = element_rect(fill = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x.bottom = element_markdown(size = 8),
      axis.title = element_text(size = 10),
      panel.spacing.x = unit(0, "lines"),
      panel.spacing.y = unit(0, "lines")
    )
}

dir.create("./images", showWarnings = FALSE)

# One B plot per scenario x C-matrix variant, saved to ./images
boot_plots <- tidyr::expand_grid(scenario = names(scenario_labels),
                                 hyp = names(hyp_labels)) |>
  pmap(\(scenario, hyp) {
    p <- mods_boot_68_B |>
      filter(scenario == .env$scenario, hyp == .env$hyp) |>
      mutate(cov = fct_drop(cov)) |>
      plot_boot_B(title = paste0(scenario_labels[[scenario]], ", ",
                                 tolower(hyp_labels[[hyp]])))

    stem <- paste0("boot_plot_", scenario, "_",
                   ifelse(hyp == "mnTS_mod_int", "int", "noint"))
    ggsave(file.path("images", paste0(stem, ".svg")), p,
           height = 15, width = 14, units = "cm", device = svg)
    ggsave(file.path("images", paste0(stem, ".png")), p,
           height = 15, width = 14, units = "cm", device = png)
    p
  })

boot_plot_int <- boot_plots[[2]] # humans, with Holocene, with interactions
boot_plot_int

## Bootstrap plotting supp info --------------------------------------------
# All four scenarios overlaid, one figure per C-matrix variant

boot_plot_compare <- function(dat, title = NULL) {
  ggplot(dat, aes(x = name, y = boot_mean, colour = as_factor(sig),
                  shape = scenario)) +
    geom_point(position = position_dodge(width = 0.6), size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.2) +
    geom_errorbar(aes(ymin = lower_68, ymax = upper_68),
                  width = .4, alpha = 0.5, position = position_dodge(width = 0.6)) +
    scale_color_manual(name = "Significance", labels = c("> 0.05", "< 0.05"),
                       values = c("#202020", "#d80000")) +
    scale_shape_manual(name = NULL, labels = scenario_labels, values = c(17, 19, 2, 1),
                       drop = FALSE) +
    guides(shape = guide_legend(nrow = 2, byrow = TRUE)) +
    labs(x = "Taxa", y = "MultinomialTS coefficient estimate", title = title) +
    facet_wrap(~ cov, labeller = as_labeller(X_names_list)) +
    theme_bw() +
    theme(
      strip.text = element_markdown(size = 9),
      strip.background = element_rect(fill = NA),
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.background = element_rect(fill = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x.bottom = element_markdown(size = 8),
      axis.title = element_text(size = 10),
      panel.spacing.x=unit(0, "lines"),
      panel.spacing.y=unit(0, "lines")
    )
}

boot_plot_all_int <- mods_boot_68_B |>
  filter(hyp == "mnTS_mod_int") |>
  boot_plot_compare(title = "With species interaction")
boot_plot_all_int

boot_plot_all_noint <- mods_boot_68_B |>
  filter(hyp == "mnTS_mod") |>
  boot_plot_compare(title = "Without species interaction")
boot_plot_all_noint

# Everything on one figure: covariates down, C-matrix variant across
boot_plot <- mods_boot_68_B |>
  boot_plot_compare() +
  facet_grid(cov ~ hyp,
             labeller = labeller(cov = as_labeller(X_names_list),
                                 hyp = as_labeller(hyp_labels)))
boot_plot

walk2(list(boot_plot_all_int, boot_plot_all_noint),
      c("boot_plot_all_int", "boot_plot_all_noint"),
      \(p, nm) {
        ggsave(file.path("images", paste0(nm, ".svg")), p,
               height = 19, width = 27, units = "cm", device = svg)
        ggsave(file.path("images", paste0(nm, ".png")), p,
               height = 19, width = 27, units = "cm", device = png)
      })

ggsave(
  "./images/boot_plot_all.svg",
  boot_plot,
  height = 32,
  width = 27,
  units = "cm",
  device = svg)

ggsave(
  "./images/boot_plot_all.png",
  boot_plot,
  height = 32,
  width = 27,
  units = "cm",
  device = png)

## Bootstrap C plotting supp info ------------------------------------------

mods_boot_68_C <- mods_boot_68 |>
  filter(grepl("sp.", name)) |>
  mutate(name = str_replace_all(name,
    pattern = "y1|y2|y3|y4|y5",
    replacement = function(x) case_when(
      x == "y1" ~ "Other",
      x == "y2" ~ "Grass",
      x == "y3" ~ "Herbs",
      x == "y4" ~ "_Pinus_",
      x == "y5" ~ "_Quercus_",
      TRUE ~ x
    )),
   name = str_remove(name, "sp."),
   name = fct(name, levels = c("Other.Other", "Grass.Grass", "Herbs.Herbs",
                                "_Pinus_._Pinus_", "_Quercus_._Quercus_",
                                "_Pinus_._Quercus_", "_Quercus_._Pinus_")))


mods_boot_68_C_plot <- ggplot(mods_boot_68_C |> filter(hyp == "mnTS_mod_int"),
                              aes(x = name, y = boot_mean, colour = as_factor(sig), shape = scenario)) +
  geom_point(position = position_dodge(width = 0.6), size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.2) +
  geom_errorbar(aes(ymin = lower_68, ymax = upper_68),
                    width = .4, alpha = 0.5, position = position_dodge(width = 0.6)) +
  scale_color_manual(name = "Significance", labels = c("> 0.05", "< 0.05"),
                     values = c("#202020", "#d80000")) +
  scale_shape_manual(name = NULL, labels = scenario_labels, values = c(17, 19, 2, 1)) +
  guides(shape = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(x = "Taxa", y = "MultinomialTS coefficient estimate") +
  # facet_wrap(~ hyp) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.background = element_rect(fill = NA),
    axis.text.x.bottom = element_markdown(size = 8, angle = 45, hjust = 1),
    axis.title = element_text(size = 10))

mods_boot_68_C_plot

ggsave(
  "./images/boot_plot_C_SI.svg",
  mods_boot_68_C_plot,
  height = 15,
  width = 14,
  units = "cm",
  device = svg)


ggsave(
  "./images/boot_plot_C_SI.png",
  mods_boot_68_C_plot,
  height = 15,
  width = 14,
  units = "cm",
  device = png)



# Covariate and pollen plots ----------------------------------------------

X_names_list2 <- c(
  heinrich ="Heinrich",
  d18O = "&delta;<sup>18</sup>O",
  mean_co2 ="CO<sub>2</sub>",
  char_acc ="Charcoal",
  ocfs ="Spores"
)

cov_plot <- all_composite |>
  select(bins, char_acc, ocfs, d18O, heinrich, mean_co2) |> 
  mutate(across(c(char_acc, ocfs, d18O, mean_co2), forecast::na.interp)) |> 
  pivot_longer(-c(bins)) |> 
  ggplot(aes(x = bins, y = value)) +
  geom_point(size = 0.4) +
  geom_line(linewidth = 0.3) +
  scale_x_reverse() +
  scale_y_continuous(breaks = scales::breaks_pretty(2)) +
  coord_flip() +
  facet_wrap(~name, scales = "free", nrow = 1,labeller = as_labeller(X_names_list2)) +
  theme_minimal() +
  theme(
    text = element_text(size = 10),
    strip.text = element_markdown(size = 9),
    axis.text.y = element_blank(),
    axis.title.y = element_blank())


all_composite_plot <- all_composite |>
  select(age, bins, other, Grass, Herbs, Pinus, Quercus) |>
  pivot_longer(-c(age, bins)) |>
  mutate(name = str_replace_all(name, 
                                pattern = "Pinus|Quercus", 
                                replacement = function(x) case_when(
                                  x == "Pinus" ~ "_Pinus_",
                                  x == "Quercus" ~ "_Quercus_",
                                  TRUE ~ x
                                )),
         name = fct(name, levels = c("other", "Grass", "Herbs", "_Pinus_", "_Quercus_")))

pol_count_plot <- ggplot(all_composite_plot, aes(x = age, y = value)) +
  geom_area(colour = "grey90") +
  geom_segment(data = all_composite_plot,
               aes(x = age, xend = age,
                   y = 0, yend = value), colour = "grey30", linewidth = 0.6) +
  scale_x_reverse() +
  coord_flip() +
  ylim(0, 400) +
  labs(y = "Pollen counts", x = "Time (ybp)") +
  facet_wrap(~name,
             nrow = 1) +
  theme_minimal() +
  theme(
    text = element_text(size = 10),
    strip.text = element_markdown()
  )

tulane_esa <- pol_count_plot + cov_plot + plot_layout(widths = c(0.6, 0.4))
ggsave(
  filename = "./images/tulane_esa.svg",
  plot = tulane_esa,
  device = "svg",
  width = 10.67,
  height = 6,
  units = "in"
)


###

# Scenario used for the headline (ESA) figures below
esa_scenario <- "humans_holocene"

boot_C_plot_int <- mods_boot_68_C |>
  filter(
    scenario == esa_scenario,
    hyp %in% c("mnTS_mod_int"),
    name %in% c("_Pinus_._Quercus_", "_Quercus_._Pinus_")
  ) %>%
  mutate(panel = "Interaction coefficients") %>%  # fake facet label
  ggplot(aes(x = name, y = boot_mean)) +
  geom_point(position = position_dodge(width = 0.6), size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.2) +
  geom_errorbar(
    aes(ymin = lower_68, ymax = upper_68),
    width = .4, alpha = 0.5, position = position_dodge(width = 0.6)
  ) +
  labs(x = "Taxa", y = NULL) +
  facet_wrap(~panel) +  # adds strip title
  theme_bw() +
  theme(
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.background = element_rect(fill = NA),
    axis.text.x.bottom = element_markdown(size = 8, angle = 45, hjust = 1),
    axis.title = element_text(size = 10),
    strip.background = element_rect(fill = NA),
  )

X_names_list3 <- c(
  heinrich ="Heinrich events",
  d18O = "&delta;<sup>18</sup>O",
  mean_co2 ="CO<sub>2</sub>",
  char_acc ="Charcoal accumulation",
  ocfs ="Fungal spores"
)
boot_plot_int2 <- ggplot(mods_boot_68_B |> filter(scenario == esa_scenario,
                                                  hyp == "mnTS_mod_int"),
                        aes(x = name, y = boot_mean, colour = as_factor(sig))) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.2) +
  geom_errorbar(aes(ymin = lower_68, ymax = upper_68),
                width = .4, alpha = 0.5) +
  scale_color_manual(name = "Significance", labels = c("> 0.05", "< 0.05"),
                     values = c("#202020", "#d80000")) +
  labs(x = "Taxa", y = "MultinomialTS coefficient estimate") +
  facet_wrap(~ cov, labeller = as_labeller(X_names_list3)) +
  theme_bw() +
  theme(
    strip.text = element_markdown(size = 9),
    strip.background = element_rect(fill = NA),
    legend.position = "inside",
    legend.position.inside = c(.09, .92),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.background = element_rect(fill = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x.bottom = element_markdown(size = 8),
    axis.title = element_text(size = 10),
    panel.spacing.x=unit(0, "lines"),
    panel.spacing.y=unit(0, "lines")
  )
boot_plot_int2

tulane_B_C_ESA <- boot_plot_int2 + boot_C_plot_int + plot_layout(widths = c(0.8, 0.2))

ggsave(
  filename = "./images/tulane_B_C_ESA.svg",
  plot = tulane_B_C_ESA,
  device = "svg",
  width = 10,
  height = 6,
  units = "in"
)
