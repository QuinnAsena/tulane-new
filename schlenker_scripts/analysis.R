# Lake Tulane Combined Age-depth Analysis
# Using LOI tie points with the Grimm 1994 core (downloaded from Neotoma)

rm(list = ls())

# Libraries
library(neotoma2)
library(Bchron)
library(tidyverse)
library(tidytext)
library(ggtext)
library(ggrepel)
library(readxl)
library(topicmodels)
library(riojaPlot)
library(rioja)
library(tapas)
library(gridExtra)
library(bcp)
library(ggtext)
library(RColorBrewer)
# library(devtools)
# devtools::install_github("wfinsinger/tapas")

dir.create("./output")
dir.create("./figures")
dir.create("./figures/topic_plots")
dir.create("./TULA20_output")
dir.create("./TULA94_output")

heinrichs <- data.frame(heinrichs = c(0:6), 
                        agemax = c(12900, 18300, 24700, 31300, 40200, 49300, 62400), 
                        agemin = c(11600, 15100, 23400, 30000, 38300, 47600, 59700),
                        name = c("**H0** <br> 12.9-11.6", "**H1** <br> 18.3-15.1", "**H2** <br> 24.7-23.4",
                                 "**H3** <br> 31.1-30.0", "**H4** <br> 40.2-38.3", "**H6** <br> 49.3-47.6", 
                                 "**H7** <br> 62.4-59.7"))

# 1. Age-Depth Models -----------
## TULA94 Age-depth Model -----------
# # Get information for Lake Tulane - siteID = 2570
tula_site <- get_sites(2570) #Download site information
tula_datasets <- get_datasets(tula_site)
tula_downloads <- get_downloads(tula_site)
tula_chron <- chronologies(tula_downloads) #Download chronology information

tula_chron %>% as.data.frame() #View chronologies, we are using the controls from the 2006 analysis but we need uncalibrated dates which is 11349
 
controls <- neotoma2::chroncontrols(tula_downloads) %>%
   dplyr::filter(chronologyid == 11349) %>%
   arrange(depth)
controls

#Remove outliers per Grimm 2006 and old Heinrich dates, then add back new heinrich dates
controls_outrm <- controls %>% dplyr::filter(depth != 665.0 & depth != 1021.0 & depth != 1165.0 & depth != 1341.0 & depth != 1365.0 &
                                              depth != 1377.0 & depth != 1393.0 & depth != 1401.0 & depth != 1676.0) %>%
  add_row(siteid = NA, chronologyid = NA, depth = 1369, thickness = 1, agelimitolder = 47650, chroncontrolid = 1,
          agelimityounger = 47550, chroncontrolage = 47600, chroncontroltype = "Heinrich stadials") %>%
  add_row(siteid = NA, chronologyid = NA, depth = 1605, thickness = 1, agelimitolder = 59750, chroncontrolid = 2,
          agelimityounger = 59650, chroncontrolage = 59700, chroncontroltype = "Heinrich stadials")
write.csv(controls_outrm, "./TULA94_output/TULA94_chroncontrols.csv")

## re-run Grimm age-depth model using Bchron
predictDepths <- c(0:1800)

# Confirmed that chroncontrolage and the agelimityounger equal to radiocarbon error

tula94_newChron <- Bchron::Bchronology(ages = controls_outrm$chroncontrolage,
                                ageSds = abs(controls_outrm$agelimityounger -
                                               controls_outrm$chroncontrolage),
                                calCurves = c("normal", rep("intcal20", 38), "normal", "normal"),
                                positionThicknesses = controls_outrm$thickness,
                                positions = controls_outrm$depth,
                                allowOutside = TRUE,
                                ids = controls_outrm$chroncontrolid)

# Predict ages at each depth for which we have samples.
tula94_chronPredicts <- predict(tula94_newChron, predictDepths)
tula94_dates <- data.frame(depth = predictDepths, median_age = apply(tula94_chronPredicts, 2, median), quant_5perc = apply(tula94_chronPredicts, 2, quantile, probs = c(0.05)),
                           quant_95perc = apply(tula94_chronPredicts, 2, quantile, probs = c(0.95)))

write.csv(tula94_dates, "./TULA94_output/TULA94_newAges.csv")
tula94_agedepthmodel_new <- list(tula94_newChron = tula94_newChron, tula94_dates = tula94_dates)
save(tula94_agedepthmodel_new, file = "./TULA94_output/tula94_agedepthmodel_new.RData")
load("./TULA94_output/tula94_agedepthmodel_new.RData")
tula94_newChron <- tula94_agedepthmodel_new[["tula94_newChron"]]

tula94_chronPlot <- 
  plot(tula94_newChron) + coord_cartesian(xlim = c(0, 65000), ylim = c(0,1650)) + 
    scale_x_continuous(breaks = c(seq(0, 65000,5000)), labels = c(seq(0, 65,5))) +
    scale_y_continuous(breaks = c(seq(0, 1650, 200)), labels = c(seq(0, 16,2))) +
    xlab("Age (yr. BP x 1000)") +
    ylab("Depth (m)")+
    coord_flip()
ggsave("./TULA94_output/TULA94_chronPlot.jpeg", tula94_chronPlot, width = 3.5, height = 5, units = "in")

### Update TULA94 Pollen and updated age-depth model -----------
#datasets(tula_datasets)
get_publications(datasetid = 19620) # Check to make sure this is the correct pollen record from the Grimm 2006 paper
tula94_pollen_dataset <- tula_datasets %>% neotoma2::filter(datasetid == 19620)
tula94_pollen_dl <- get_downloads(tula94_pollen_dataset)
tula94_pollen <- samples(tula94_pollen_dl)
load(file = "./TULA94_output/tula94_agedepthmodel_new.RData")
tula94_dates <- tula94_agedepthmodel_new[["tula94_dates"]]
tula94_datestoadd <- data.frame(depth = tula94_dates$depth, new_age = tula94_dates$median_age)
tula94_pollen <- tula94_pollen %>% inner_join(tula94_datestoadd, by = c("depth" = "depth"))
tula94_agecomp <- unique(tula94_pollen[,c("age", "new_age")])
write.csv(tula94_agecomp, "./TULA94_output/tula94_agecompare.csv")


# Just using exactly the pollen counts and taxonomy from Grimm - only updating the ages

tula94_pollen_wide_prop <- tula94_pollen %>% toWide(ecologicalgroup = c("TRSH", "UPHE"),
         unit = c("NISP"),
         elementtypes = c("pollen"),
         groupby = "new_age",
         operation = "prop") %>%
  arrange(new_age)

tula94_pollen_wide_count <- tula94_pollen %>% toWide(ecologicalgroup = c("TRSH", "UPHE"),
                                               unit = c("NISP"),
                                               elementtypes = c("pollen"),
                                               groupby = "new_age",
                                               operation = "sum") %>%
  arrange(new_age)

save(tula94_pollen, file = "./TULA94_output/tula94_pollen.RData")
write.csv(tula94_pollen, file = "./TULA94_output/tula94_pollen.csv")
save(tula94_pollen_wide_prop, file = "./TULA94_output/tula94_pollen_wide_prop.RData")
write.csv(tula94_pollen_wide_prop, file = "./TULA94_output/tula94_pollen_wide_prop.csv")
save(tula94_pollen_wide_count, file = "./TULA94_output/tula94_pollen_wide_count.RData")
write.csv(tula94_pollen_wide_count, file = "./TULA94_output/tula94_pollen_wide_count.csv")

## TULA20 Age-depth model --------------
### Read in raw data -----------
sheet_names <- excel_sheets("./data_raw/TULA20_RawData.xlsx")
TULA20_data <- list()
for (i in 1:length(sheet_names)) {
  st <- sheet_names[i]
  temp <- read_xlsx("./TULA20_RawData.xlsx", sheet = st, col_names = TRUE)
  TULA20_data[[i]] <- temp
}

names(TULA20_data) <- sheet_names

rm(i, st, temp, sheet_names)

## TULA20 Age-depth model 
### add ages to the LOI tie points from the TULA94 age depth model
tula20_loi_tiepoints <- TULA20_data[["LOI_TiePoints"]] %>% left_join(tula94_dates, by = c("Depth_Grimm" = "depth"))
tula20_loi_tiepoints$error <- (tula20_loi_tiepoints$quant_95perc - tula20_loi_tiepoints$quant_5perc)/2 #calculate error
tula20_loi_tiepoints$curve <- "normal" #specify what curve to use

tula20_otherdates <- TULA20_data[["Other_Dates"]] #load in the other dates


tula20_controls <- data.frame(labID = c(tula20_loi_tiepoints$Event, tula20_otherdates$LabID),
                              depth = c(tula20_loi_tiepoints$TULA20_Depth_Composite, tula20_otherdates$depth_comp),
                              age = round(c(tula20_loi_tiepoints$median_age, tula20_otherdates$age), 0),
                              error = round(c(tula20_loi_tiepoints$error, tula20_otherdates$error), 0),
                              curve = c(tula20_loi_tiepoints$curve, tula20_otherdates$curve)) %>% arrange(depth) %>%
                    add_row(labID = "coretop", depth = 0, age = -70, error = 15, curve = "normal") %>%
                    arrange(depth)
  
tula20_controls$thickness <- 1
tula20_controls$ageid <- 1:nrow(tula20_controls)

write.csv(tula20_controls , "./TULA20_output/TULA20_controls.csv")

tula20_newChron <- Bchron::Bchronology(ages =tula20_controls$age,
                                       ageSds = tula20_controls$error,
                                       calCurves = tula20_controls$curve,
                                       positionThicknesses = tula20_controls$thickness,
                                       positions = tula20_controls$depth,
                                       allowOutside = TRUE,
                                       ids = tula20_controls$ageid)

# Predict ages at each depth for which we have samples. 
tula20_chronPredicts <- predict(tula20_newChron, predictDepths)
tula20_dates <- data.frame(depth = predictDepths, median_age = apply(tula20_chronPredicts, 2, median), quant_5perc = apply(tula20_chronPredicts, 2, quantile, probs = c(0.05)), 
                           quant_95perc = apply(tula20_chronPredicts, 2, quantile, probs = c(0.95)))
tula20_chronPlot <- plot(tula20_newChron)
ggsave("./TULA20_output/TULA20_chronPlot.jpeg", tula20_chronPlot, height = 6, width = 8, units = "in", dpi = 300)

tula20_composite <- TULA20_data[["Composite_Core"]] %>%
  left_join(tula20_dates, by = c("depth_composite" = "depth")) %>%
  mutate(accu_rate = median_age - lag(median_age, default = median_age[1]))

tula20_agedepthmodel <- list(tula20_newChron = tula20_newChron, tula20_chronPredicts = tula20_chronPredicts, tula20_composite = tula20_composite)

save(tula20_agedepthmodel, file = "./TULA20_output/tula20_agedepthmodel.RData")
write.csv(tula20_composite , "./TULA20_output/TULA20_compsiteCore_w-ages.csv")

load("./TULA20_output/tula20_agedepthmodel.RData")
tula20_newChron <- tula20_agedepthmodel[["tula20_newChron"]]
tula20_chronPlot <- plot(tula20_newChron) + coord_cartesian(xlim = c(0, 65000), ylim = c(0,1650)) + 
  scale_x_continuous(breaks = c(seq(0, 65000,5000)), labels = c(seq(0, 65,5))) +
  scale_y_continuous(breaks = c(seq(0, 1650, 200)), labels = c(seq(0, 16,2))) +
  xlab("Age (yr. BP x 1000)") +
  ylab("Depth (m)")+
  coord_flip()
ggsave("./TULA20_output/TULA20_chronPlot.jpeg", tula20_chronPlot, height = 5, width = 4, units = "in", dpi = 300)

### Update TULA20 ages in LOI, CHAR, and CFS ----------
load(file = "./TULA20_output/tula20_agedepthmodel.RData")
tula20_composite <- tula20_agedepthmodel[["tula20_composite"]]

tula20_loi <- TULA20_data[["TULA_LOI"]] %>% inner_join(tula20_composite, by = c("core" = "core", "depth_core" = "totalHolDepth"))
tula20_char <- TULA20_data[["TULA_CHAR"]] %>% inner_join(tula20_composite, by = c("core" = "core", "depth_core" = "totalHolDepth"))
tula20_cfs <- TULA20_data[["TULA_CFS"]] %>% inner_join(tula20_composite, by = c("Core" = "core", "depth_core" = "totalHolDepth"))

tula20_char$CHAR <- tula20_char$charCount / tula20_char$accu_rate

write.csv(tula20_loi , "./TULA20_output/TULA20_LOI_w-ages.csv")
write.csv(tula20_cfs , "./TULA20_output/TULA20_CFS_w-ages.csv")
write.csv(tula20_char , "./TULA20_output/TULA20_CHAR_w-ages.csv")
#write.csv(tula20_char500 , "./TULA20_output/TULA20_CHAR500_w-ages.csv")

save(tula20_loi, file = "./TULA20_output/TULA20_LOI_w-ages.RData")
save(tula20_cfs, file = "./TULA20_output/TULA20_CFS_w-ages.RData")
save(tula20_char, file = "./TULA20_output/TULA20_CHAR_w-ages.RData")


# 2. Topic Modelling ----------
rm(list = ls()) #Remove files from age-depth models that are no longer needed
## Aggregate pollen into genera 
load(file = "./TULA94_output/tula94_pollen.RData")

tula94_pollen_ag_count <- tula94_pollen %>% 
  mutate(variablename = replace(variablename, stringr::str_detect(variablename, "Pinus*"), "Pinus") ) %>%
  toWide(ecologicalgroup = c("TRSH", "UPHE"),
         unit = c("NISP"),
         elementtypes = c("pollen"),
         groupby = "new_age",
         operation = "sum") %>%
  pivot_longer(-new_age) %>% 
  mutate(value = as.integer(value)) %>% #need to convert to integer for topic models
  pivot_wider(names_from = name, values_from = value) %>%
  arrange(new_age)

tula94_pollen_ag_prop <- tula94_pollen %>% 
  mutate(variablename = replace(variablename, stringr::str_detect(variablename, "Pinus*"), "Pinus") ) %>%
  toWide(ecologicalgroup = c("TRSH", "UPHE"),
         unit = c("NISP"),
         elementtypes = c("pollen"),
         groupby = "new_age",
         operation = "prop") %>%
  arrange(new_age)



#create defaults for multiple numbers of topics 
nlist <- 3:6 #set number of topics
reptimes <- length(nlist) #set number or replications for below
seeds <- read.csv(file = "./data_raw/CTMseeds.csv")
seeds <- seeds[,2]

#Used controls suggested for CTM from Grün and Hornik 2011
myctrl <-
  list(
    estimate.beta = TRUE,
    verbose = 0,
    prefix = tempfile(),
    save = 0,
    keep = 0,
    seed = seeds,
    nstart = 20L, 
    best = TRUE,
    var = list(iter.max = 500, tol = 10 ^ -6),
    em = list(iter.max = 1000, tol = 10 ^ -4),
    initialize = "random",
    cg = list(iter.max = 500, tol = 10 ^ -5)
  )

#make repeated lists of pollen counts and controls to input into mapply 
poll_only <- dplyr::select(tula94_pollen_ag_count, -new_age)
poll_reps <- do.call("list", replicate(reptimes, as.matrix(poll_only), simplify = FALSE))
ctrls <- do.call("list", replicate(reptimes, myctrl, simplify = FALSE))
ages <- tula94_pollen_ag_count$new_age

# Run the CTM topic model
ctm_mods <- mapply(CTM, k=nlist, x=poll_reps, control = ctrls) #CTM run with controls and using 20 seeds with the best model selected

save(ctm_mods, file = "./output/ctm_mods.RData")

# Evaluate the CTM models
aicsctm <- do.call(rbind,lapply(ctm_mods, AIC))
bicsctm <- do.call(rbind,lapply(ctm_mods, BIC))
eval_ctm <- data.frame(k = nlist, aic = aicsctm, bic = bicsctm)
eval_ctm_plot <- ggplot(eval_ctm) + geom_point(aes(k, aicsctm), col = "blue") + geom_point(aes(k, bicsctm), col = "orange") + xlim(c(3,8)) + xlab("Number of Topics") + ylab("Model Fit (AIC and BIC)") + theme_minimal() 
eval_ctm_plot
ggsave("./figures/topic_plots/eval_ctm.jpeg", eval_ctm_plot)

#function to extract terms and topics and make basic plots for each number of groups
termtopic_func <- function(tm_result, age_vector = ages){
  postr <- lapply(tm_result, posterior)
  
  ## Plot "terms" 
  terms_dfs <- lapply(postr, function(x){
    x <- x[["terms"]]
    x <- data.frame(x)
    x$ngroups <- c(rep(nrow(x), nrow(x)))
    x$topic <- factor(1:nrow(x))
    pivot_longer(x, c(-ngroups, -topic), names_to = "taxa", values_to = "beta")
  })
  
  term_plot <- lapply(terms_dfs, function(x){
    x<- x %>%
      group_by(ngroups, topic) %>%
      top_n(8, beta) %>%
      ungroup() %>%
      mutate(taxa = reorder_within(taxa, beta, topic))
    ggplot(x, aes(x = reorder(taxa, beta), beta, fill = topic)) + geom_col() +
      facet_wrap(~topic, scales = "free", drop = TRUE, ncol = 2) +
      scale_fill_brewer(palette = "Set3") +
      scale_x_reordered() +
      coord_flip()
  })
  
  ## Plot "topics"
  topic_dfs <- lapply(postr, function(x){
    x <- x[["topics"]]
    x <- data.frame(x)
    colnames(x) <- factor(1:ncol(x))
    x$sample <- 1:nrow(x)
    x$age <- age_vector
    x<- pivot_longer(x, c(-sample, -age), names_to = "topic", values_to = "gamma")
    x
  })
  
  
  topic_plot <- lapply(topic_dfs, function(x){
    ggplot(x, aes(x = age, y = gamma, fill = topic)) + geom_area()+
      scale_fill_brewer(palette = "Set3") +
      scale_x_reverse()+
      theme_minimal()
  })
  list(terms_dfs, term_plot, topic_dfs, topic_plot)
}


ctm_termtopic <- termtopic_func(ctm_mods)
ctm_termtopic[[4]][[2]]
ctm_termtopic[[2]][[2]]

# Save all graphs for sup figs
# save term graphs
for (i in 1:length(nlist)) {
  n <- i + (min(nlist)-1)
  p <- ctm_termtopic[[2]][[i]]
  jpeg(paste0("./figures/topic_plots/terms_N", n, ".jpeg"), height = 4, width = 6, units = "in", res = 300)
  print(p)
  dev.off()
}
# save topic graphs
for (i in 1:length(nlist)) {
  n <- i + (min(nlist)-1)
  p <- ctm_termtopic[[4]][[i]]
  jpeg(paste0("./figures/topic_plots/topics_N", n, ".jpeg"), height = 3, width = 6, units = "in", res = 300)
  print(p)
  dev.off()
}

saveRDS(ctm_termtopic, file = "./output/ctm_termtopic.RDS")
#ctm_termtopic <- readRDS("./output/ctm_termtopic.RDS")

rm(p, n, i, nlist, reptimes, ctrls, aicsctm, bicsctm, myctrl, poll_only, poll_reps, ctm_mods)
#Will use 4 topics moving forward with main analysis
ctm_termtopic <- readRDS(file = "./output/ctm_termtopic.RDS")
ctm_topics4 <- ctm_termtopic[[3]][[2]]
ctm_topics4 <- ctm_topics4 %>%
  pivot_wider(names_from = topic, values_from = gamma)
ctm_terms4 <- ctm_termtopic[[1]][[2]]
write.csv(ctm_topics4, "./output/ctm_topics4.csv")
write.csv(ctm_terms4, "./output/ctm_terms4.csv")

# 3. Charcoal ------------
load(file = "./TULA20_output/TULA20_CHAR_w-ages.RData")
tapas_input <- data.frame(
  CmTop	= tula20_char$depth_composite,
  CmBot	= tula20_char$depth_composite + 1,
  AgeTop = tula20_char$median_age,
  AgeBot = tula20_char$median_age + tula20_char$accu_rate,
  Volume = 1,
  CH = tula20_char$charCount)

peaks_full <- peak_detection(tapas_input, thresh_type = "local")

tapas_in_7a <- dplyr::filter(tapas_input, AgeTop < 10000)
tapas_in_4a <- dplyr::filter(tapas_input, AgeTop > 10000)

peaks7a <- peak_detection(tapas_in_7a, thresh_type = "local")
peaks4a <- peak_detection(tapas_in_4a, thresh_type = "local")

peaks_analysis <- list("peaks7a" = peaks7a, "peaks4a" = peaks4a)
saveRDS(peaks_analysis, "./output/CHAR_peaks_analysis.RDS")
peaks_analysis <- readRDS("./output/CHAR_peaks_analysis.RDS")
peaks7a <- peaks_analysis[["peaks7a"]]
peaks4a <- peaks_analysis[["peaks4a"]]

# Plot.Anomalies(series = peaks7a , plot.crosses = TRUE, plot.x = TRUE, plot.neg = FALSE)
# Plot_ReturnIntervals(series = peaks7a, plot.x = TRUE, plot.neg = FALSE)
# Plot.Anomalies(series = peaks4a, plot.crosses = TRUE, plot.x = TRUE, plot.neg = FALSE)
# Plot_ReturnIntervals(series = peaks4a, plot.x = TRUE, plot.neg = FALSE)

ri_tbl_7a <- data.frame(age = peaks7a[["thresh"]]$peaks.pos.ages, RI = peaks7a[["thresh"]]$RI_pos, time = "Holocene")
ri_tbl_4a <- data.frame(age = peaks4a[["thresh"]]$peaks.pos.ages, RI = peaks4a[["thresh"]]$RI_pos, time = "Pleistocene")
ri_tbl <- rbind(ri_tbl_7a, ri_tbl_4a)
write.csv(ri_tbl, "./output/ri_tbl.csv")

char_tbl_7a <- data.frame(age = peaks7a[["int"]][["series.int"]]$age, CHAR = peaks7a[["int"]][["series.int"]]$CHAR, time = "Holocene")
char_tbl_4a <- data.frame(age = peaks4a[["int"]][["series.int"]]$age, CHAR = peaks4a[["int"]][["series.int"]]$CHAR, time = "Pleistocene")
char_tbl <- rbind(char_tbl_7a, char_tbl_4a)
write.csv(char_tbl, "./output/char_tbl.csv")

char_peaks_7a <- data.frame(peaks = peaks7a[["thresh"]]$peaks.pos.ages, loc = 45)
char_peaks_4a <- data.frame(peaks = peaks4a[["thresh"]]$peaks.pos.ages, loc = 45)
char_peaks <- rbind(char_peaks_7a, char_peaks_4a)
write.csv(char_peaks, "./output/char_peaks.csv")


sni_tbl_7a <- data.frame(age = peaks7a[["thresh"]]$ages.thresh, SNI = peaks7a[["thresh"]][["SNI_pos"]]$SNI_raw)
sni_tbl_4a <- data.frame(age = peaks4a[["thresh"]]$ages.thresh, SNI = peaks4a[["thresh"]][["SNI_pos"]]$SNI_raw)
sni_tbl <- rbind(sni_tbl_7a, sni_tbl_4a)
write.csv(sni_tbl, "./output/sni_tbl.csv")

png(file = "./figures/supp_CHAR_fullFig.png", width = 6, height = 6, units = "in", res = 600)
layout(mat = matrix(c(1,3,4,2,3,4), ncol = 2))
Plot.Anomalies(series = peaks4a, plot.crosses = TRUE, plot.x = TRUE, plot.neg = FALSE)
Plot.Anomalies(series = peaks7a, plot.crosses = TRUE, plot.x = TRUE, plot.neg = FALSE)
base::plot(x = ri_tbl$age, y = ri_tbl$RI, type = "l",  xlim = rev(range(ri_tbl$age)), ylab = ("Return Interval (years)"), xlab = "Age (calendar years BP)", frame = FALSE)
base::plot(x = sni_tbl$age, y = sni_tbl$SNI, type = "l", ylim = c(0, 10), xlim = rev(range(sni_tbl$age)), ylab = ("Signal-to-Noise Ratio"), xlab = "Age (calendar years BP)", frame = FALSE)
abline(h = 3, col = "red")
dev.off()

char_plot <- ggplot(char_tbl) + 
  geom_col(aes(age, CHAR), col = "black") + 
  geom_point(data = char_peaks, mapping = aes(x = peaks, y = loc), col = "red", shape = 3) + 
  coord_cartesian(ylim=c(0, 55), expand = c(0)) + 
  scale_x_reverse(breaks = seq(0,65000, 10000)) + 
  ylab("Charcoal Accumulation Rate <br> (CHAR)") +
  theme_minimal()+
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.title.y = element_markdown(size = 10),
        axis.line = element_line(colour = "grey"))
char_plot 

ri_plot <- 
  ri_tbl %>% ggplot(aes(age, RI)) + geom_line() + 
  scale_x_reverse(breaks = seq(0,65000, 10000)) + 
  scale_y_continuous(breaks = c(0, 500, 1000, 1500))+
  coord_cartesian(expand = c(0))+
  ylab("Return Interval <br> (years)") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.title.y = element_markdown(size = 10),
        axis.line = element_line(colour = "grey"))

sni_plot <- 
  sni_tbl %>% ggplot(aes(age, SNI)) + geom_line() + 
  scale_x_reverse(breaks = seq(0,65000, 10000)) + 
  geom_hline(yintercept = 3, col = "red") +
  xlab("Age (calendar years BP)") +
  ylab("Signal-to-Noise <br> Ratio") +
  scale_y_continuous(breaks = c(0,2,4,6,8,10), expand = c(0,0))+
  coord_cartesian(ylim=c(0, 10), expand = c(0))+
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.title.y = element_markdown(size = 10),
        axis.line = element_line(colour = "grey"))

p1 <- ggplotGrob(char_plot)
p2 <- ggplotGrob(ri_plot)
p3 <- ggplotGrob(sni_plot)

g1 <- rbind(p1,p2,p3, size = "first")
g1$widths <- unit.pmax(p1$widths, p2$widths, p3$widths)
grid.newpage()

png(file = "./figures/supp_CHAR_fullFigx2.png", width = 6, height = 6, units = "in", res = 600)
grid.arrange(g1)
dev.off()

mean(ri_tbl_7a$RI, na.rm = TRUE)
mean(ri_tbl_4a$RI, na.rm = TRUE)
mean(ri_tbl$RI, na.rm = TRUE)

mean(char_tbl$CHAR)
mean(char_tbl$CHAR[which(char_tbl$age>=11700)])
mean(char_tbl$CHAR[which(char_tbl$age<11700)])

max(char_tbl$CHAR)
max(char_tbl$CHAR[which(char_tbl$age>=11700)])
max(char_tbl$CHAR[which(char_tbl$age<11700)])

mean(sni_tbl$SNI)
length(which(sni_tbl$SNI >3))
length(sni_tbl$SNI)
length(which(sni_tbl$SNI >3))/length(sni_tbl$SNI)

## Charcoal tests to figure out what is happening with the last 10k
### Split record in two


# 4. CFS --------------
load(file = "./TULA20_output/TULA20_CFS_w-ages.RData")

ocfs <- tula20_cfs%>% select("depth_core", "median_age", "Delitschia", "Sporormiella", "Podospora", "Cercophora", "Sordaria", "Lycopodium", "Tracer_Added") %>%
  mutate(OCFS_total = (Delitschia + Sporormiella + Podospora + Cercophora + Sordaria),
    OCFS_conc = ((Tracer_Added * OCFS_total) /(Lycopodium * 1)),
    total_counted = (OCFS_total + Lycopodium))

ComputeCI <- function(samplesize, tracerfound, traceradded, n = 10^6) {
  
  prop <- tracerfound / samplesize
  max <- 1 / (prop * (1 - prop)) * traceradded * (1 / prop - 1)
  step <- ceiling(max / 500)
  spore <- seq(0, max, step)
  
  # point estimate of spore density
  estimspore <- ceiling(traceradded * samplesize / tracerfound) - traceradded
  
  # compute the maximum likelihood
  # use this to scale (normalize) the likelihood function
  maxlik <- dhyper(x = tracerfound,
                   m = traceradded,
                   n = estimspore,
                   k = samplesize)
  
  # compute the (normalized) likelihood for each spore density
  ylik <- dhyper(x = tracerfound,
                 m = traceradded,
                 n = spore,
                 k = samplesize) / maxlik
  
  # the interpolated function for the rejection sampling
  flik <- splinefun(spore, ylik)
  
  # generate random bivariate data
  bestestim  <- traceradded * (1 / prop - 1)
  ymax <- floor(1 / (prop * (1 - prop)) * bestestim)
  dat <- matrix(c(runif(n, 0, ymax), runif(n, 0, 1.005)), ncol = 2)
  
  accept  <- c()
  # rejection sampling
  for (i in 1:n){
    xx <- dat[i, 1]
    yy <- dat[i, 2]
    if (yy <= flik(xx)) {accept <- c(accept, xx)}
  }
  
  # compute the quantile for the confidence interval
  quant  <- quantile(accept, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975))
  
  names(bestestim)  <- "best estimate"
  result <- c(bestestim, quant)
  
  return(result)
}

ocfs_no0s <- filter(ocfs, OCFS_total != 0)

ll <- nrow(ocfs_no0s) # number of data
out <- data.frame("best_estimate"= numeric(), "Q2.5" = numeric(), "Q5" = numeric(), 
                  "Q25"= numeric(), "Q50"= numeric(), "Q75"= numeric(), "Q95"= numeric(), "Q97.5"= numeric())
for (i in 1:ll){
  ss <- ocfs_no0s$total_counted[i]
  tf <- ocfs_no0s$Lycopodium[i]
  ta <- ocfs_no0s$Tracer_Added[i]
  out[i,] <- ComputeCI(ss, tf, ta, n = 10^3)
}
out$age <- ocfs_no0s$median_age

ocfs_uncertainty <- left_join(ocfs, out, by = c("median_age" = "age"))
saveRDS(ocfs_uncertainty, "./output/ocfs_uncertainty.RDS" )
write.csv(ocfs_uncertainty, "./output/ocfs_uncertainty.csv")

ocfs_long <-tula20_cfs %>% 
  select("median_age", "Delitschia", "Sporormiella", "Podospora", "Cercophora", "Sordaria", "Lycopodium", "Tracer_Added", "accu_rate") %>%
  pivot_longer(cols = -c(median_age, Tracer_Added, Lycopodium, accu_rate)) %>%
  mutate(OCFS_conc = ((Tracer_Added * value) /(Lycopodium * 1)), 
         OCFS_accumuation = OCFS_conc/accu_rate) 
saveRDS(ocfs_long, "./output/ocfs_long.RDS" )

ocfs_uncertainty <- readRDS("./output/ocfs_uncertainty.RDS" )
ocfs_long <- readRDS("./output/ocfs_long.RDS" )

ocfs_taxa_conc <- ocfs_long %>% 
  ggplot(aes(median_age, OCFS_conc, colour = name)) +
  geom_line() +
  scale_x_reverse() +
  ylab("Spore Concentration <br> (spores/cm<sup>3</sup>)")+
  xlab("Age (calendar year BP)")+
  theme_minimal() +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_markdown(),
        panel.grid = element_blank())
ggsave("./figures/supp_OCFS_taxa_conc.jpeg")

ocfs_taxa_accu <- ocfs_long %>% 
  ggplot(aes(median_age, OCFS_accumuation, colour = name)) +
  geom_line() +
  scale_x_reverse() +
  ylab("Spore Accumulation Rate <br> (spores/cm<sup>2</sup>/year)")+
  xlab("Age (calendar year BP)")+
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(), axis.title.y = element_markdown(),
        panel.grid = element_blank())
ggsave("./figures/supp_OCFS_taxa_accu.jpeg")

p1 <- ggplotGrob(ocfs_taxa_conc)
p2 <- ggplotGrob(ocfs_taxa_accu)

g4 <- rbind(p1,p2, size = "first")
g4$widths <- unit.pmax(p1$widths, p2$widths)
grid.newpage()

jpeg("./figures/supp_OCFS_TaxaPlot.jpeg", width = 7, height = 5, units = "in", res = 300)
grid.arrange(g4, ncol =1)
dev.off()

ocfs_hienrichs <- ocfs_uncertainty %>% select(median_age, OCFS_conc) %>%
  mutate(hienrichs = case_when(
    median_age > 11600 & median_age < 12900 ~ "H0",
    median_age > 15100 & median_age < 18300 ~ "H1",
    median_age > 23400 & median_age < 24700 ~ "H2",
    median_age > 30000 & median_age < 31100 ~ "H3",
    median_age > 38300 & median_age < 40200 ~ "H4",
    median_age > 47600 & median_age < 49300 ~ "H5",
    median_age > 59700 & median_age < 62400 ~ "H6",
    .default = "NO"))

mean(ocfs_hienrichs$OCFS_conc[which(ocfs_hienrichs$hienrichs %in% c("H1", "H2", "H3", "H4"))])

mean(ocfs_hienrichs$OCFS_conc[which(ocfs_hienrichs$hienrichs == "NO" & 
                                      ocfs_hienrichs$median_age < 38300 &
                                      ocfs_hienrichs$median_age > 18300 )])

# 5. BCP ----------
### the package bcp is not archived, downloaded latest version 
# devtools::install_version("bcp",version="4.0.3")
pinus_bcp <- bcp(tula94_pollen_ag_prop$Pinus, burnin = 100, mcmc = 1000, p0 = 0.12)
quercus_bcp <- bcp(tula94_pollen_ag_prop$Quercus, burnin = 100, mcmc = 1000, p0 = 0.12)
ocfs_bcp <- bcp(ocfs_uncertainty$OCFS_conc, burnin = 100, mcmc = 1000, p0 = 0.12)
char_bcp <-  bcp(char_tbl$CHAR, burnin = 100, mcmc = 1000, p0 = 0.12)

bcp_tbl_pollen <- data.frame(tula94_pollen_ag_prop$new_age, 
                      pinus_bcp$posterior.mean, 
                      pinus_bcp$posterior.prob, 
                      quercus_bcp$posterior.mean, 
                      quercus_bcp$posterior.prob)
colnames(bcp_tbl_pollen) <- c("age", "Pinus_mean", "Pinus_prob", "Quercus_mean", "Quercus_prob")
bcp_tbl_ocfs <- data.frame(ocfs_uncertainty$median_age, 
                      ocfs_bcp$posterior.mean, 
                      ocfs_bcp$posterior.prob)
colnames(bcp_tbl_ocfs) <- c("age", "OCFS_mean", "OCFS_prob")
write.csv(bcp_tbl_pollen, "./output/bcp_tbl_pollen.csv")
write.csv(bcp_tbl_ocfs, "./output/bcp_tbl_ocfs.csv")
# bcp_tbl <- read.csv("./output/bcp_tbl_pollen.csv")
# pinus_bcp_dates <- bcp_tbl %>% select(age, Pinus_prob) %>% filter(Pinus_prob >= 0.5)

jpeg("./figures/supp_bcp_pinus.jpeg")
plot(pinus_bcp, main = "Pinus")
dev.off()

jpeg("./figures/supp_bcp_quercus.jpeg")
plot(quercus_bcp, main = "Quercus")
dev.off()

jpeg("./figures/supp_bcp_ocfs.jpeg")
plot(ocfs_bcp, main = "OCFS")
dev.off()

## Sup-fig combined
pinus_bcp_plot <- ggplot() +
  geom_area(data = tula94_pollen_ag_prop, aes(x = ages, y = Pinus), col = "#6D8B3D", fill = "#6D8B3D") +
  geom_line(data = bcp_tbl_pollen, aes(x = ages, y = Pinus_prob*0.81), linetype = "dashed") +
  scale_y_continuous(name = "*Pinus* Pollen (%)", breaks = c(0,.2,.4,.6,.8), labels = c(0,20,40,60,80), 
                     sec.axis = sec_axis( trans=~./0.81, name="Posterior Probability", breaks = c(0,.5,1)), expand = c(0, 0)) +
  scale_x_reverse(breaks = seq(0,65000, 10000), expand = c(0, .01)) + 
  coord_cartesian(ylim=c(0, .85), xlim=c(65000, 0)) + 
  theme_minimal()+
  theme(panel.grid = element_blank(), axis.line = element_line(color = "grey"),
        axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_markdown(size = 10)) +
  geom_rect(data = heinrichs, aes(xmax=agemax, xmin=agemin, ymin=-Inf, ymax=Inf), alpha = 0.25, fill='black')

quercus_bcp_plot <- ggplot() +
  geom_area(data = tula94_pollen_ag_prop, aes(x = ages, y = Quercus), col = "#EC7825", fill = "#EC7825") +
  geom_line(data = bcp_tbl_pollen, aes(x = ages, y = Quercus_prob*0.57), linetype = "dashed") +
  scale_y_continuous(name = "*Quercus* Pollen (%)", breaks = c(0,.2,.4,.6), labels = c(0,20,40,60), 
                     sec.axis = sec_axis( trans=~./0.57, name="Posterior Probability", breaks = c(0,.5,1)), expand = c(0, 0)) +
  scale_x_reverse(breaks = seq(0,65000, 10000), expand = c(0, .01)) + 
  coord_cartesian(ylim=c(0, .6), xlim=c(65000, 0)) + 
  theme_minimal()+
  theme(panel.grid = element_blank(), axis.line = element_line(color = "grey"),
        axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_markdown(size = 10)) +
  geom_rect(data = heinrichs, aes(xmax=agemax, xmin=agemin, ymin=-Inf, ymax=Inf), alpha = 0.25, fill='black')

ocfs_bcp_plot <- ggplot() +
  geom_area(data = ocfs_uncertainty, aes(x = median_age, y = OCFS_conc/1000), col = "#F2696C", fill = "#F2696C") +
  geom_line(data = bcp_tbl_ocfs, aes(x = age, y = OCFS_prob*39.010), linetype = "dashed") +
  scale_y_continuous(name = "OCFS <br> (1000xspores/cm<sup>3</sup>)", breaks = seq(0,40,10), 
                     sec.axis = sec_axis( trans=~./39.010, name="Posterior Probability", breaks = c(0,.5,1)), expand = c(0, 0)) +
  scale_x_reverse(breaks = seq(0,65000, 10000), labels = seq(0,65000, 10000), expand = c(0, .01)) + 
  coord_cartesian(ylim=c(0, 40), xlim=c(65000, 0)) + 
  xlab("Age (calendar years BP)")+
  theme_minimal()+
  theme(panel.grid = element_blank(), axis.line = element_line(color = "grey"), 
        axis.title.y = element_markdown(size = 10)) +
  geom_rect(data = heinrichs, aes(xmax=agemax, xmin=agemin, ymin=-Inf, ymax=Inf), alpha = 0.25, fill='black')
  
p1 <- ggplotGrob(pinus_bcp_plot)
p2 <- ggplotGrob(quercus_bcp_plot)
p3 <- ggplotGrob(ocfs_bcp_plot)

g1 <- rbind(p1,p2,p3, size = "first")
g1$widths <- unit.pmax(p1$widths, p2$widths, p3$widths)
grid.newpage()

png(file = "./figures/supp_BCP_full.png", width = 8, height = 6, units = "in", res = 600)
grid.arrange(g1)
dev.off()


# 6. Pollen ----------------
poll_plot <- tula94_pollen_ag_prop %>% 
  select(c("new_age", "Pinus",  "Quercus",  "Ulmus", "Carya", "Myrica", "Cupressaceae undiff.",
           "Ericaceae", "Poaceae", "Artemisia", "Iva frutescens-type",
           "Asteroideae undiff.",  "Amaranthaceae undiff.", "Cyperaceae")) %>%
  pivot_longer(cols = -new_age) %>%
  mutate(name_fact = factor(name, levels = c("Pinus",  "Quercus",  "Ulmus", "Carya", "Myrica", "Cupressaceae undiff.", 
                                             "Ericaceae","Poaceae", "Artemisia", "Iva frutescens-type",
                                             "Asteroideae undiff.",  "Amaranthaceae undiff.", "Cyperaceae"),
                            labels = c("*Pinus*",  "*Quercus*",  "*Ulmus*", "*Carya*", "*Myrica*", "Cupressaceae", 
                                       "Ericaceae", "Poaceae", "*Artemisia*", "*Iva*",
                                       "Asteroideae",  "Amaranthaceae", "Cyperaceae")
                            ),
         ecol_group = factor(name, levels = c("Pinus",  "Quercus",  "Ulmus", "Carya", "Myrica", "Cupressaceae undiff.", 
                                             "Ericaceae","Poaceae", "Artemisia", "Iva frutescens-type",
                                             "Asteroideae undiff.",  "Amaranthaceae undiff.", "Cyperaceae"),
                            labels = c("Pine",  "TRSH",  "TRSH", "TRSH", "TRSH", "TRSH", 
                                       "TRSH", "UPHE", "UPHE", "UPHE",
                                       "UPHE",  "UPHE", "AQVP")
         )) %>%
  ggplot() + 
  geom_area(aes(x = new_age, y = round(value*100, 0), fill = ecol_group)) + 
  geom_vline(xintercept = pinus_bcp_dates$age, linetype = "dashed", linewidth = 0.25, color = "red") +
  geom_rect(data = heinrichs, aes(xmax=agemax, xmin=agemin, ymin=-Inf, ymax=Inf), alpha = 0.3, fill='black') +
  scale_fill_discrete(type = c("#015245", "#e3c027", "#b07eed", "#83c0f2")) +
  scale_x_reverse(expand = c(0,0, 0, 0), limits = c(62750, -1), breaks = seq(from = 60000, to = 0, by = -10000))+
  scale_y_continuous(labels = scales::number_format(accuracy = 1))+
  coord_flip()+
  facet_wrap(~name_fact, scales = "free_x", nrow = 1) +
  xlab("Age (Calendar Years BP)") +
  theme_minimal()+
  theme(strip.text.x = element_markdown(angle=70), 
        axis.text.x = element_text(size = 5, angle = 315, margin=margin(-2,0,0,0)),
        axis.title.x = element_blank(), 
        panel.grid = element_blank(), panel.spacing.x = unit(3, "mm"), 
        axis.ticks.x = element_line(color = "black", linewidth = 0.1),
        axis.line = element_line(color = "black", linewidth = 0.1),
        legend.position = "none")
poll_plot 

ggsave("./figures/pollen_diagram.pdf", poll_plot, width = 7.5, height = 4)

# 7. Other Data ----------
humans <- read.csv("./data_raw/humanpop.csv")
O18 <- read.csv("./data_raw/ngrip-d18o-50yr.csv")
CO2 <- read.csv("./data_raw/co2_merged.csv")
sealevel <- read.csv("./data_raw/sea_level.csv")

# Manuscript Figure ------
plot_margins <- c(-3, 5, 0, 5)
h_plot <- geom_rect(data = heinrichs, aes(xmax=agemax, xmin=agemin, ymin=-Inf, ymax=Inf), fill='lightgrey')
combo_theme <- theme(panel.border = element_rect(colour = "#00000000", fill = "#00000000"),
                     panel.background = element_rect(colour = "#00000000", fill = "#00000000"),
                      panel.grid = element_blank(), 
                     axis.title = element_markdown(size = 8))

terms_plot <- ctm_terms4 %>% select(-ngroups) %>% 
  group_by(topic) %>%
  top_n(8, beta) %>%
  ungroup() %>%
  mutate(topic_fact = factor(topic, levels = c(3, 1, 4, 2), labels = c("Pine Forest", "Oak-Grass <br> Savanna", "Oak-Pine <br> Forest", "Grass-herbs"))) %>%
  mutate(
         taxa = replace(taxa, stringr::str_detect(taxa, "Ambrosia*"), "*Ambrosia*"),
         taxa = replace(taxa, stringr::str_detect(taxa, "Quercus*"), "<i>Quercus</i>"),
         taxa = replace(taxa, stringr::str_detect(taxa, "Iva*"), "<i>Iva</i>"),
         taxa = replace(taxa, stringr::str_detect(taxa, "Fagus*"), "<i>Fagus</i>"),
         taxa = replace(taxa, stringr::str_detect(taxa, "Fraxinus*"), "<i>Fraxinus</i>"),
         taxa = replace(taxa, stringr::str_detect(taxa, "Juglans*"), "<i>Juglans</i>"),
         taxa = replace(taxa, stringr::str_detect(taxa, "Ulmus*"), "<i>Ulmus</i>"),
         taxa = replace(taxa, stringr::str_detect(taxa, "Polygonum*"), "<i>Polygonum</i>"),
         taxa = replace(taxa, stringr::str_detect(taxa, "Carya*"), "<i>Carya</i>"),
         taxa = replace(taxa, stringr::str_detect(taxa, "Pinus*"), "<i>Pinus</i>"),
         taxa = replace(taxa, stringr::str_detect(taxa, "Myrica*"), "<i>Myrica</i>"), 
         taxa = replace(taxa, stringr::str_detect(taxa, "Amaranthaceae*"), "Amaranth."),
         taxa = replace(taxa, stringr::str_detect(taxa, "Asteroideae*"), "Astero."),
         taxa = replace(taxa, stringr::str_detect(taxa, "Cupressaceae*"), "Cupress."),
         taxa = replace(taxa, stringr::str_detect(taxa, "Cyperaceae*"), "Cyper."),
         taxa = reorder_within(taxa, beta, topic)) %>%
  ggplot(aes(x = reorder(taxa, beta), beta, fill = topic_fact)) + 
  geom_col() +
  scale_fill_discrete(type = c("#048772", "#81CDC3", "#A76214", "#E0C37F")) +
  scale_x_reordered(limits=rev) + 
  scale_y_continuous(expand = c(0,0, 0, 0), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0.00", "0.25", "0.50", "0.75", "1.00")) +
  facet_wrap(~topic_fact, scales = "free_x", drop = TRUE, ncol = 4) +
  ylab("Beta") +
  theme_minimal() + 
  theme(axis.text.x = element_markdown(angle = 315, hjust = 0, size = 5),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 8),
        strip.background = element_blank(), 
        strip.text = element_markdown(size = 7), 
        legend.position = "none", 
        panel.grid = element_blank(),
        plot.margin = margin(0, 30, 5, 10))
terms_plot

ctm_plot <- ctm_topics4 %>%
  pivot_longer(c( -sample, -age), names_to = "topic", values_to = "gamma") %>%
  mutate(topic_fact = factor(topic, levels = c(2, 4, 1, 3), labels = c("Mixed Forest???", "Oak-Hickory <br> Forb Forest", "Oak-Grass Savanna", "Pine Forest"))) %>%
  ggplot()+
  geom_area(aes(age, gamma, fill = topic_fact)) +
  scale_fill_discrete(type = c("#E0C37F", "#A76214", "#81CDC3", "#048772")) +
  scale_x_reverse(expand = c(0,0), limits = c(65000, 0), )+
  scale_y_continuous(expand = c(0.01,0,0,0))+
  ylab("Gamma") +
  theme_minimal()+
  theme(axis.text.y=element_text(size=6), axis.text.x=element_blank(), axis.title.x = element_blank(), 
        legend.position = "none", plot.margin = margin(plot_margins)) + combo_theme
ctm_plot <- ctm_plot +  geom_rect(data = heinrichs, aes(xmax=agemax, xmin=agemin, ymin=-Inf, ymax=Inf), alpha = 0.25, fill='black')
#ctm_plot

cfs_plot <- ocfs_uncertainty %>% 
  ggplot() +
  h_plot +
  geom_line(aes(x = median_age, y = OCFS_conc)) +
  geom_ribbon(aes(x = median_age, ymin = Q2.5, ymax = Q97.5), fill = "darkgreen", alpha = 0.5) +
  scale_x_reverse(expand = c(0,0), limits = c(65000, 0))+
  scale_y_continuous(expand = c(0.1,0,0,0), limits = c(0, 40000), position = "right")+
  ylab("CFS") +
  theme_minimal()+
  theme(axis.text.y=element_text(size=6), 
        axis.text.x=element_blank(), 
        axis.title.x = element_blank(),
        legend.position = "none", plot.margin = margin(plot_margins)) + combo_theme
 
#cfs_plot 
char_peaks$loc <- 80

char_plot <- ggplot(char_tbl) + 
  h_plot +
  geom_col(aes(age, CHAR), col = "black", width = 0.5) +
  # geom_line(aes(age, CHAR), col = "black") +
  # geom_area(aes(age, CHAR), fill = "black") + 
  geom_point(data = char_peaks, mapping = aes(x = peaks, y = loc), col = "red", shape = 3, size = 0.5) + 
  coord_cartesian(ylim=c(0, 120)) + 
  ylab("CHAR") +
  scale_x_reverse(expand = c(0,0), limits = c(65000, 0))+
  scale_y_continuous(expand = c(0,0,0,0), limits = c(0, 90), breaks = c(0, 20, 40, 60, 80))+
  theme_minimal()+
  theme(axis.text.y=element_text(size=6), 
        axis.text.x=element_blank(), 
        axis.title.x = element_blank(),
        axis.line.x = element_line(color = "black"),
        legend.position = "none", plot.margin = margin(plot_margins)) + combo_theme

char_plot 

o18_plot <- ggplot(O18) + 
  h_plot +
  geom_line(aes(Age, d18O), col = "black", linewidth = 0.05) +
  scale_x_reverse(expand = c(0,0), limits = c(65000, 0))+
  scale_y_continuous(expand = c(.2,0,0,0), limits = c(-45, -33), breaks = c(-44, -40, -36), position = "right")+
  labs(y = "&delta;<sup>18</sup>O") +
  theme_minimal()+
  theme(axis.text.y=element_markdown(size=6), 
        axis.title.y = ggtext::element_markdown(),
        axis.text.x=element_blank(), 
        axis.title.x = element_blank(), 
        legend.position = "none", plot.margin = margin(plot_margins)) + combo_theme

#o18_plot

sl_plot <- sealevel %>%
  mutate(age = t.ka.*-1000) %>%
  ggplot() + 
  h_plot +
  geom_line(aes(x = age, y = SL), col = "black") +
  #geom_richtext(data = heinrichs, mapping = aes(x = agemean, y = -15, label = name), 
                #size = 1.25, colour = "black", fill = "#00000000", label.colour = "#00000000") +
  scale_x_reverse(expand = c(0,0), limits = c(65000, 0))+
  scale_y_continuous(expand = c(0.1,0,0,1), limits = c(-125, 0))+
  ylab("Sea Level (m)") +
  theme_minimal()+
  theme(axis.text.y=element_text(size=6), 
        axis.text.x=element_blank(), 
        axis.title.x = element_blank(), 
        legend.position = "none", plot.margin = margin(plot_margins)) + combo_theme

#sl_plot

co2_plot <- ggplot(CO2) + 
  h_plot +
  geom_line(aes(x = age_calBP, y = CO2_blank_gravity_corrected), col = "black", size = 0.3) +
  scale_x_reverse(expand = c(0,0), limits = c(65000, 0))+
  scale_y_continuous(expand = c(0.1,0,0,2), limits = c(175, 275), position = "right")+
  ylab("CO2") +
  theme_minimal()+
  theme(axis.text.y=element_text(size=6), 
        axis.text.x=element_blank(), 
        axis.title.x = element_blank(), 
        legend.position = "none", plot.margin = margin(plot_margins)) + combo_theme

#co2_plot

human_annos <- data.frame(age = c(14550, 22000), loc = 0.17, loc2 = c(0.21, 0.14), name = c("Page-Ladson", "White Sands"))
human_plot <- ggplot(humans) + 
  h_plot +
  geom_line(aes(x = Date, y = NA.Human), col = "black") +
  geom_point(data = human_annos, mapping = aes(x = age, y = loc), col = "darkred", fill = "darkred", shape = 25) + 
  geom_text(data = human_annos, mapping = aes(x = age+5000, y = loc2, label = name), size = 2) +
  scale_x_reverse(expand = c(0,0), limits = c(65000, 0), 
                  breaks = seq(from = 65000, to = 0, by = -5000), 
                  minor_breaks = seq(from = 65000, to = 0, by = -1000),
                  guide = guide_axis(minor.ticks = TRUE))+
  scale_y_continuous(expand = c(0.1,0.01,0,0.01), limits = c(0, 0.25))+
  ylab("Humans")+
  xlab("Age (calendar years BP)") +
  #theme_minimal()+
  theme(axis.text.x = element_text(size = 5), axis.text.y=element_text(size=6), 
        legend.position = "none", plot.margin = margin(plot_margins), axis.ticks = element_line(linewidth = 0.2)) + combo_theme

#human_plot 

## Combine plot -----------
p1 <- ggplotGrob(terms_plot)
p2 <- ggplotGrob(ctm_plot)
p3 <- ggplotGrob(cfs_plot)
p4 <- ggplotGrob(char_plot)
p5 <- ggplotGrob(o18_plot) 
p6 <- ggplotGrob(sl_plot) 
p7 <- ggplotGrob(co2_plot) 
p8 <- ggplotGrob(human_plot) 

g4 <- rbind(p2,p3,p4,p5,p6,p7,p8, size = "first")
g4$widths <- unit.pmax(p2$widths, p3$widths, p4$widths, p5$widths, p6$widths, p7$widths, p8$widths)
grid.newpage()

layout_in <- matrix(data = c(1, 2, 2, 2), ncol = 1)

pdf("./figures/fig2_combo_plot.pdf", width = 4, height = 6)
grid.arrange(p1, g4, layout_matrix = layout_in, padding = unit(0, "line"))
dev.off()  

jpeg("./figures/fig2_combo_plot.jpeg", width = 4, height = 6, units = "in", res = 300)
grid.arrange(p1, g4, layout_matrix = layout_in, padding = unit(0, "line"))
dev.off()  

  
  
