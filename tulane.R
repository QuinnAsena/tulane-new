if (!require("pacman")) install.packages("pacman", repos="http://cran.r-project.org")
pacman::p_load(multinomialTS, tidyverse, ggtext, patchwork, forecast,
               future, furrr, viridis, viridisLite, RColorBrewer)


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
core_20_spore <- read_csv("./data/TULA20_age-depth_files/TULA20_CFS_w-ages.csv")
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
  select(depth_composite, `OCFS Concentration`)
colnames(core_20_spore) <- c("depth" , "ocfs")

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
                              cov_age <= 40200 & cov_age >= 36800 ~ 1,
                              cov_age <= 31300 & cov_age >= 30000 ~ 1,
                              cov_age <= 24700 & cov_age >= 23400 ~ 1,
                              cov_age <= 18300 & cov_age >= 15100 ~ 1,
                              cov_age <= 12900 & cov_age >= 11000 ~ 1,
                              .default = NA),
         humans = NA,
         humans = case_when(cov_age  <= 14500 ~ 1, .default = NA)) |> 
  arrange(desc(cov_age))
# Binary variables left as NA for summarise later
dim(composite_covariate_join)

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
    humans = mean(humans, na.rm = T),
    ocfs = mean(ocfs, na.rm = T)) |> 
  mutate(heinrich = ifelse(is.nan(heinrich), 0, heinrich),
         humans = ifelse(is.nan(humans), 0, humans)) |>
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

all_composite <- pollen_wide_binned |> 
  full_join(composite_covariate_join_bin, by = "bins") |> 
  full_join(co2_mean, by = "bins") |> 
  full_join(oxy18_mean, by = "bins") |> 
  arrange(desc(bins))

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
  select(bins, char_acc, ocfs, d18O, humans, heinrich, mean_co2) |> 
  pivot_longer(-c(bins)) |> 
  ggplot(aes(x = bins, y = value)) +
  geom_point() +
  geom_line() +
  facet_wrap(~name, scales = "free", ncol = 1)

all_composite |> 
  select(bins, char_acc, ocfs, d18O, humans, heinrich, mean_co2) |> 
  mutate(across(c(char_acc, ocfs, d18O, mean_co2), forecast::na.interp)) |> 
  pivot_longer(-c(bins)) |> 
  ggplot(aes(x = bins, y = value)) +
  geom_point() +
  geom_line() +
  facet_wrap(~name, scales = "free", ncol = 1)

all_composite |>
  select(bins, char_acc, ocfs, d18O, humans, heinrich, mean_co2) |> 
  mutate(across(c(char_acc, ocfs, d18O, mean_co2), forecast::na.interp)) |> 
  mutate(across(-c(humans, heinrich, bins), ~ as.numeric(scale(.)))) |> 
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

# Without interactions ----------------------------------------------------

# set up Y
Y <- all_composite |>
  select(bins, other, Grass, Herbs, Pinus, Quercus) |>
  arrange(desc(bins)) |> 
  select(-bins) |> 
  as.matrix()

Tsample <- which(rowSums(Y) != 0)

# set up X
X <- all_composite |>
  select(bins, char_acc, ocfs, d18O, humans, heinrich, mean_co2) |> 
  arrange(desc(bins)) |> 
  mutate(ocfs = sqrt(ocfs)) |> 
  mutate(across(c(char_acc, ocfs, d18O, mean_co2), forecast::na.interp),
         across(-c(humans, heinrich, bins), ~ as.numeric(scale(.)))) |>
  select(-bins) |> 
  as.matrix()

which(is.na(X))

matplot(X, type = 'l')

p <- ncol(X) + 1 # Number of independent variables plus intercept
n <- ncol(Y)

V.fixed = diag(n) # Covariance matrix of environmental variation in process eq
# V.fixed = matrix(NA, n, n)
# V.fixed[1] = 1

B.fixed <- matrix(c(rep(0,p),rep(NA, (n - 1) * p)), p, n)
B.start <- matrix(c(rep(0,p),rep(.01, (n - 1) * p)), p, n)

glmm_mod <- multinomialTS::mnGLMM(Y = Y[which(rowSums(Y) != 0),],
                                  X = X[which(rowSums(Y) != 0), ,drop = F],
                                  B.start = B.start, B.fixed = B.fixed,
                                  V.fixed = V.fixed)
summary(glmm_mod)

B0.start <- glmm_mod$B[1, , drop = F]
B.start <- glmm_mod$B[2:p, , drop = F]

sigma.start <- glmm_mod$sigma

V.fixed = matrix(NA, n, n) # Covariance matrix of environmental variation in process eq
V.fixed[1] = 1

V.start <- glmm_mod$V
# V.start <- diag(diag(V.start))

B.fixed <- matrix(NA, ncol(X), n)
B.fixed[,1] <- 0
B0.fixed = matrix(c(0, rep(NA, n - 1)), nrow = 1, ncol = n)

# Set-up C without interactions
C.start.diag = .5 * diag(n)
C.fixed.diag <- C.start.diag
C.fixed.diag[C.fixed.diag != 0] <- NA


# Model with no interactions
start_time <- Sys.time()
mnTS_mod <- mnTS(Y = Y[Tsample, ],
                 X = X, Tsample = Tsample,
                 B0.start = B0.start, B0.fixed = B0.fixed,
                 B.start = B.start, B.fixed = B.fixed,
                 C.start = C.start.diag, C.fixed = C.fixed.diag,
                 V.start = V.start, V.fixed = V.fixed,
                 dispersion.fixed = 1, maxit.optim = 1e+6)

end_time <- Sys.time()
end_time - start_time

mnTS_mod_refit <- refit_func(mnTS_mod, 5)
lapply(mnTS_mod_refit, coef)

# With interactions -------------------------------------------------------

C.start.diag.int = .5 * diag(n)
C.start.diag.int[4, 5] = C.start.diag.int[5, 4] = .5

C.fixed.diag.int <- C.start.diag.int
C.fixed.diag.int[C.fixed.diag.int != 0] <- NA


start_time <- Sys.time()
mnTS_mod_int <- mnTS(Y = Y[Tsample, ],
                     X = X, Tsample = Tsample,
                     B0.start = mnTS_mod$B0, B0.fixed = B0.fixed,
                     B.start = mnTS_mod$B, B.fixed = B.fixed,
                     C.start = C.start.diag.int, C.fixed = C.fixed.diag.int,
                     V.start = mnTS_mod$V, V.fixed = V.fixed,
                     dispersion.fixed = 1, maxit.optim = 1e+6)

end_time <- Sys.time()
end_time - start_time

mnTS_mod_int_refit <- refit_func(mnTS_mod_int, 5)
lapply(mnTS_mod_int_refit, coef)


# bootstrapping -----------------------------------------------------------

start_time <- Sys.time()

future::plan(strategy = multisession, workers = 10)
mods <- c(setNames(rep(list(mnTS_mod_refit[[5]]), 5),
                    paste0("mnTS_mod", 1:5)),
          setNames(rep(list(mnTS_mod_int_refit[[3]]), 5),
                     paste0("mnTS_mod_int", 1:5)))

res <- furrr::future_map(mods, multinomialTS::boot, rep = 200,
                         .options = furrr_options(seed = 1984))
saveRDS(res, "./results/bootstraps_1000.rds")
end_time <- Sys.time()
end_time - start_time

# Plotting ----------------------------------------------------------------

ssms <- list(mnTS_mod = mnTS_mod_refit[[5]],
             mnTS_mod_int = mnTS_mod_int_refit[[3]])

wald <- lapply(ssms, \(hyp) {
  wald <- coef(hyp)
  as_tibble(wald, rownames = "cov")
})
wald_bind <- bind_rows(wald, .id = "hyp") |> 
  mutate(sig = P < 0.05,
         hyp = forcats::fct(hyp))

wald_bind_x <- wald_bind |> 
  filter(!grepl("sp.|^Grass|^Herbs|^Pinus|^Quercus", cov))


B_plot <- ggplot(wald_bind_x, aes(x = hyp, y = Coef., colour = as_factor(sig))) +
  geom_point() +
  geom_errorbar(aes(ymin = Coef. - se, ymax = Coef. + se)) +
  geom_hline(yintercept = 0) +
  scale_color_manual(name = "Significance", labels = c("> 0.05", "< 0.05"),
                     values = c("#BF0606", "#5ab4ac")) +
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
  facet_wrap(~hyp) +
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
res <- readRDS("./results/bootstraps_1000.rds")

X_names_list <- c(
  heinrich ="A: Heinrich events",
  d18O = "B: &delta;<sup>18</sup>O",
  mean_co2 ="C: CO<sub>2</sub>",
  char_acc ="D: Charcoal accumulation",
  ocfs ="E: Fungal spores",
  humans = "F: Human presence"
)

mods_boot <- map(res, ~ {
  as_tibble(.x[[2]]) |> 
  pivot_longer(-c(logLik, opt.convergence))
}) |> 
  bind_rows(.id = "hyp") |> 
  mutate(hyp = str_remove(hyp, pattern = '[[:digit:]]+'))

mods_boot_68 <- mods_boot |> 
#  filter(opt.convergence == 0) |> 
  group_by(hyp, name) |> 
  summarise(boot_mean = mean(value),
            boot_sd = sd(value),
            upper_68 = quantile(value, probs = 0.84),
            lower_68 = quantile(value, probs = 0.16)) |> 
  mutate(t_scores = boot_mean / boot_sd,
         p_vals = 2 * pnorm(q = abs(t_scores), lower.tail = F),
         sig = p_vals < 0.05)
#

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
        cov = fct(cov, levels = c("heinrich", "d18O", "mean_co2", "char_acc", "ocfs", "humans")))


boot_plot_int <- ggplot(mods_boot_68_B |> filter(hyp == "mnTS_mod_int"),
                        aes(x = name, y = boot_mean, colour = as_factor(sig))) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.2) +
  geom_errorbar(aes(ymin = lower_68, ymax = upper_68),
                    width = .4, alpha = 0.5) +
  scale_color_manual(name = "Significance", labels = c("> 0.05", "< 0.05"),
                     values = c("#202020", "#d80000")) +
  labs(x = "Taxa", y = "MultinomialTS coefficient estimate") +
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
    axis.text = element_markdown(size = 8),
    axis.title = element_text(size = 10),
    panel.spacing.x=unit(0, "lines"),
    panel.spacing.y=unit(0, "lines") 
  )
boot_plot_int


ggsave(
  "./results/boot_plot_int.svg",
  boot_plot_int,
  height = 15,
  width = 14,
  units = "cm",
  device = svg)

ggsave(
  "./results/boot_plot_int.png",
  boot_plot_int,
  height = 15,
  width = 14,
  units = "cm",
  device = png)

# Without holocene --------------------------------------------------------
all_composite_woholo <- all_composite |>
  filter(bins >= 54) |>
  arrange(desc(bins))

# Without interactions ----------------------------------------------------

# set up Y
Y <- all_composite_woholo |>
  select(bins, other, Grass, Herbs, Pinus, Quercus) |>
  arrange(desc(bins)) |> 
  select(-bins) |> 
  as.matrix()

Tsample <- which(rowSums(Y) != 0)

# set up X
X <- all_composite_woholo |>
  select(bins, char_acc, ocfs, d18O, humans, heinrich, mean_co2) |> 
  arrange(desc(bins)) |> 
  mutate(ocfs = sqrt(ocfs)) |> 
  mutate(across(c(char_acc, ocfs, d18O, mean_co2), forecast::na.interp),
         across(-c(humans, heinrich, bins), ~ as.numeric(scale(.)))) |>
  select(-bins) |> 
  as.matrix()

which(is.na(X))

matplot(X, type = 'l')

p <- ncol(X) + 1 # Number of independent variables plus intercept
n <- ncol(Y)

V.fixed = diag(n) # Covariance matrix of environmental variation in process eq
# V.fixed = matrix(NA, n, n)
# V.fixed[1] = 1

B.fixed <- matrix(c(rep(0,p),rep(NA, (n - 1) * p)), p, n)
B.start <- matrix(c(rep(0,p),rep(.01, (n - 1) * p)), p, n)

glmm_mod <- multinomialTS::mnGLMM(Y = Y[which(rowSums(Y) != 0),],
                                  X = X[which(rowSums(Y) != 0), ,drop = F],
                                  B.start = B.start, B.fixed = B.fixed,
                                  V.fixed = V.fixed)
summary(glmm_mod)

B0.start <- glmm_mod$B[1, , drop = F]
B.start <- glmm_mod$B[2:p, , drop = F]

sigma.start <- glmm_mod$sigma

V.fixed = matrix(NA, n, n) # Covariance matrix of environmental variation in process eq
V.fixed[1] = 1

V.start <- glmm_mod$V
# V.start <- diag(diag(V.start))

B.fixed <- matrix(NA, ncol(X), n)
B.fixed[,1] <- 0
B0.fixed = matrix(c(0, rep(NA, n - 1)), nrow = 1, ncol = n)

# Set-up C without interactions
C.start.diag = .5 * diag(n)
C.fixed.diag <- C.start.diag
C.fixed.diag[C.fixed.diag != 0] <- NA


# Model with no interactions
start_time <- Sys.time()
mnTS_mod <- mnTS(Y = Y[Tsample, ],
                 X = X, Tsample = Tsample,
                 B0.start = B0.start, B0.fixed = B0.fixed,
                 B.start = B.start, B.fixed = B.fixed,
                 C.start = C.start.diag, C.fixed = C.fixed.diag,
                 V.start = V.start, V.fixed = V.fixed,
                 dispersion.fixed = 1, maxit.optim = 1e+6)

end_time <- Sys.time()
end_time - start_time

mnTS_mod_refit <- refit_func(mnTS_mod, 5)
lapply(mnTS_mod_refit, coef)

# With interactions -------------------------------------------------------

C.start.diag.int = .5 * diag(n)
C.start.diag.int[4, 5] = C.start.diag.int[5, 4] = .5

C.fixed.diag.int <- C.start.diag.int
C.fixed.diag.int[C.fixed.diag.int != 0] <- NA


start_time <- Sys.time()
mnTS_mod_int <- mnTS(Y = Y[Tsample, ],
                     X = X, Tsample = Tsample,
                     B0.start = mnTS_mod$B0, B0.fixed = B0.fixed,
                     B.start = mnTS_mod$B, B.fixed = B.fixed,
                     C.start = C.start.diag.int, C.fixed = C.fixed.diag.int,
                     V.start = mnTS_mod$V, V.fixed = V.fixed,
                     dispersion.fixed = 1, maxit.optim = 1e+6)

end_time <- Sys.time()
end_time - start_time

mnTS_mod_int_refit <- refit_func(mnTS_mod_int, 5)
lapply(mnTS_mod_int_refit, coef)


# bootstrapping -----------------------------------------------------------

start_time <- Sys.time()

future::plan(strategy = multisession, workers = 10)
mods <- c(setNames(rep(list(mnTS_mod_refit[[5]]), 5),
                    paste0("mnTS_mod", 1:5)),
          setNames(rep(list(mnTS_mod_int_refit[[3]]), 5),
                     paste0("mnTS_mod_int", 1:5)))

res <- furrr::future_map(mods, multinomialTS::boot, rep = 200,
                         .options = furrr_options(seed = 1984))
saveRDS(res, "./results/woholo_bootstraps_1000.rds")
end_time <- Sys.time()
end_time - start_time


## Bootstrap plotting supp info --------------------------------------------

res_woholo <- readRDS("./results/woholo_bootstraps_1000.rds")
res <- readRDS("./results/bootstraps_1000.rds")

# check on convergence
lapply(res_woholo, \(f) {
  x <- f$all_mods_pars
  sum(x[ ,colnames(x) %in% "opt.convergence"])
})

names(res_woholo) <- paste0(names(res_woholo), "_woholo")
res <- c(res, res_woholo)

mods_boot <- map(res, ~ {
  as_tibble(.x[[2]]) |> 
  pivot_longer(-c(logLik, opt.convergence))
}) |> 
  bind_rows(.id = "hyp") |> 
  mutate(hyp = str_remove(hyp, pattern = '[[:digit:]]+'))

X_names_list <- c(
  heinrich ="A: Heinrich events",
  d18O = "B: &delta;<sup>18</sup>O",
  mean_co2 ="C: CO<sub>2</sub>",
  char_acc ="D: Charcoal accumulation",
  ocfs ="E: Fungal spores",
  humans = "F: Human presence"
)

mods_boot <- map(res, ~ {
  as_tibble(.x[[2]]) |> 
  pivot_longer(-c(logLik, opt.convergence))
}) |> 
  bind_rows(.id = "hyp") |> 
  mutate(hyp = str_remove(hyp, pattern = '[[:digit:]]+'))

mods_boot_68 <- mods_boot |> 
#  filter(opt.convergence == 0) |> 
  group_by(hyp, name) |> 
  summarise(boot_mean = mean(value),
            boot_sd = sd(value),
            upper_68 = quantile(value, probs = 0.84),
            lower_68 = quantile(value, probs = 0.16)) |> 
  mutate(t_scores = boot_mean / boot_sd,
         p_vals = 2 * pnorm(q = abs(t_scores), lower.tail = F),
         sig = p_vals < 0.05)
#

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
        cov = fct(cov, levels = c("heinrich", "d18O", "mean_co2", "char_acc", "ocfs", "humans")))


boot_plot <- ggplot(mods_boot_68_B, aes(x = name, y = boot_mean, colour = as_factor(sig), shape = as_factor(hyp))) +
  geom_point(position = position_dodge(width = 0.6), size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.2) +
  geom_errorbar(aes(ymin = lower_68, ymax = upper_68),
                    width = .4, alpha = 0.5, position = position_dodge(width = 0.6)) +
  scale_color_manual(name = "Significance", labels = c("> 0.05", "< 0.05"),
                     values = c("#202020", "#d80000")) +
  scale_shape_manual(name = NULL, labels = c(
    "With species interaction", " Without species interaction",
    "With species interaction, without Holocene", " Without species interaction, without Holocene"),
                     values = c(17, 19, 2, 1)) +
  guides(shape = guide_legend(nrow = 2,byrow = TRUE)) +
  labs(x = "Taxa", y = "MultinomialTS coefficient estimate") +
  facet_wrap(~ cov, labeller = as_labeller(X_names_list)) +
  theme_bw() +
  theme(
    strip.text = element_markdown(size = 9),
    strip.background = element_rect(fill = NA),
    legend.position = "bottom",
    # legend.position = "inside",
    # legend.position.inside = c(.09, .92),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.background = element_rect(fill = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_markdown(size = 8),
    axis.title = element_text(size = 10),
    panel.spacing.x=unit(0, "lines"),
    panel.spacing.y=unit(0, "lines")
  )
boot_plot


ggsave(
  "./results/boot_plot_all.svg",
  boot_plot,
  height = 19,
  width = 27,
  units = "cm",
  device = svg)

ggsave(
  "./results/boot_plot_all.png",
  boot_plot,
  height = 19,
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


mods_boot_68_C_plot <- ggplot(mods_boot_68_C |> filter(hyp %in% c("mnTS_mod_int", "mnTS_mod_int_woholo")), aes(x = name, y = boot_mean, colour = as_factor(sig), shape = as_factor(hyp))) +
  geom_point(position = position_dodge(width = 0.6), size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.2) +
  geom_errorbar(aes(ymin = lower_68, ymax = upper_68),
                    width = .4, alpha = 0.5, position = position_dodge(width = 0.6)) +
  scale_color_manual(name = "Significance", labels = c("> 0.05", "< 0.05"),
                     values = c("#202020", "#d80000")) +
  scale_shape_manual(name = NULL, labels = c("With Holocene", " Without Holocene"), values = c(17, 2)) +
  labs(x = "Taxa", y = "MultinomialTS coefficient estimate") +
  # facet_wrap(~ hyp) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.background = element_rect(fill = NA),
    axis.text = element_markdown(size = 8, angle = 45, hjust = 1),
    axis.title = element_text(size = 10))


ggsave(
  "./results/boot_plot_C_SI.svg",
  mods_boot_68_C_plot,
  height = 15,
  width = 14,
  units = "cm",
  device = svg)


ggsave(
  "./results/boot_plot_C_SI.png",
  mods_boot_68_C_plot,
  height = 15,
  width = 14,
  units = "cm",
  device = png)
