library(dplyr)
library(readr)
library(forecast)
source("kalman-smoother/TVARSS_11Feb25.r")

# Read-in wrangled data from tulane.R
all_composite <- read_csv("./data/all_composite.csv")

# necessary as global variable
bin_width = 200

# Predictors
# Use Nora's topic model output instead of pollen counts 
topic <- read_csv("data/topic_model/ctm_topics4.csv") |>
  select(-c(`...1`)) |>
  rename(topic1 = `1`, topic2 = `2`, topic3 = `3`, topic4 = `4`) |>
  arrange(desc(age))

# Bin topics using original pollen data age
topic_bins <- cut(topic$age,
            breaks = seq(from = min(all_composite$age, na.rm = T), 
                         to = max((all_composite$age + bin_width), na.rm = T), 
                         by = bin_width),
            include.lowest = T, labels = F)

topic_binned <- bind_cols(topic_bins = topic_bins, topic) |> 
  group_by(topic_bins) |> 
  summarise(
    age = mean(age, na.rm = T),
    across(starts_with("topic"), \(x) mean(x, na.rm = TRUE))) |> 
  arrange(desc(topic_bins))

# Join bins, ages are almost identical but not exact so joining by bins instead of age
all_composite |>
  left_join(topic_binned |> select(-c(age)), by = c("bins" = "topic_bins")) |>
  select(-c(other, Grass, Herbs, Pinus, Quercus)) |>
  arrange(desc(bins))

# set up X - ocfs as response variable
X <- all_composite |>
  select(bins, ocfs) |>
  arrange(desc(bins)) |> 
  select(-bins) |> 
  as.matrix()

# First point in X needs a value so setting to 100
X[1] <- 100

# X_df <- data.frame(y = all_composite$ocfs, x = all_composite$bins)
# X_df_gam <- mgcv::gam(y ~ s(x, bs = "bs", k = nrow(X_df)), method = "REML", data =  X_df)
# pred_df <- predict(X_df_gam, newdata = data.frame(x = X_df$x[which(is.na(X_df$y))] ))
# X_df$y[which(is.na(X_df$y))] <- pred_df
# plot(X_df$x, X_df$y, type = "l", xlab = "bins", ylab = "ocfs")

# set up U - predictors
U <- all_composite |>
  select(bins, char_acc, d18O, heinrich, mean_co2, PrDens) |> 
  mutate(across(c(char_acc, d18O, mean_co2), forecast::na.interp),
         across(-c(heinrich, bins), ~ as.numeric(scale(.)))) |>
  arrange(desc(bins)) |> 
  select(-bins) |> 
  as.matrix()


# Set up parameters -------------------------------------------------------

# number of lags in the autocorrelation function
p <- 2
# standard error of the measurement error.
# Setting su.fixed = 1 assumes that the measurement error has SD = 1.
# Setting su.fixed = NA will get the fitting to estimate su.
ME <- rep(1, nrow(U))
su.fixed <- 1

# -----------------------------------------------------------------
# Analyses: I set this up so you have to pick and choose the
# X and U variables. I like this better than doing all combinations,
# because you can go through each model separately to see 
# what is going on.
# -----------------------------------------------------------------

# It is possible to have the autoregression coefficients change
# through time by setting sb0.fixed = NA and/or sb.fixed = matrix(NA,1,p).
# Setting them to zero means the auoregression coeffients don't change.
sb0.fixed <- 0
sb.fixed <- matrix(0,1,p)

mod0 <- TVARSS(X = X, p = p, ME = ME, Tsamplefract = .9, show.fig = F, annealing = F,
               sb0.fixed = sb0.fixed,
               sb.fixed = sb.fixed,
               su.fixed = su.fixed)
mod0

mod <- TVARSS(X = X, p = p, ME = ME, U = U, Tsamplefract = .9, show.fig = F, annealing = F,
                  sb0.fixed = sb0.fixed,
                  sb.fixed = sb.fixed,
                  su.fixed = su.fixed,
                  c.fixed = rep(0, ncol(U)),               
                  d.fixed = rep(NA, ncol(U)),
                  d.start = rep(.01, ncol(U)),
                  b0.start = mod0$b0,
                  b.start = mod0$b)
mod
dev <- 2*(mod$logLik - mod0$logLik)
df <- ncol(U)
c(dev = dev, df = df, P = pchisq(dev, df = df, lower.tail = F))

P_indep_vars_TVARSS(mod)



#---------- Below is Tony's code to run as a loop, not used for now ----------#

run <- function(p){
	
	w <- list()
	count <- 0
	for(i.var in c("Hemlock_d13c", "Beech_d13c")){
		count <- count + 1
		print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
		print(i.var)
		print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
	
		X <- d[,i.var]
		
		###############################################
		print("tests of LL and temp")
		
		pick.list <- list(c("LL"), c("temp"), c("LL", "temp"))
		for(i in 1:3){
			
			U <- U.list[,pick.list[[i]], drop = F]
		
			mod0 <- TVARSS(X = X, p = p, ME = ME, Tsamplefract = .9, show.fig = F, annealing = F,
			               sb0.fixed = sb0.fixed,
			               sb.fixed = sb.fixed,
			               su.fixed = su.fixed)
			mod0
#browser()
			mod <- TVARSS(X = X, p = p, ME = ME, U = U, Tsamplefract = .9, show.fig = F, annealing = F,
			                  sb0.fixed = sb0.fixed,
			                  sb.fixed = sb.fixed,
			                  su.fixed = su.fixed,
			                  c.fixed = rep(0, ncol(U)),               
			                  d.fixed = rep(NA, ncol(U)),
			                  b0.start = mod0$b0,
			                  b.start = mod0$b)
			show(pick.list[[i]])
			show(mod)
			
			dev <- 2*(mod$logLik - mod0$logLik)
			df <- ncol(U)
			show(c(dev = dev, df = df, P = pchisq(dev, df = df, lower.tail = F)))
		}
		
		###############################################
		print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
		print("tests of interactions")
		
		U1 <- U.list[,c("LL","temp"), drop = F]
		U <- cbind(U1, inter = U.list[,"LL_temp"])
		
		mod0 <- TVARSS(X = X, p = p, ME = ME, Tsamplefract = .9, show.fig = F, annealing = F,
		               sb0.fixed = sb0.fixed,
		               sb.fixed = sb.fixed,
		               su.fixed = su.fixed)
		               
		mod1 <- TVARSS(X = X, p = p, ME = ME, U = U1, Tsamplefract = .9, show.fig = F, annealing = F,
		                  sb0.fixed = sb0.fixed,
		                  sb.fixed = sb.fixed,
		                  su.fixed = su.fixed,
		                  d.fixed = rep(NA, ncol(U1)),
		                  c.fixed = rep(0, ncol(U1)),               
		                  b0.start = mod0$b0,
		                  b.start = mod0$b)
		                  
		mod <- TVARSS(X = X, p = p, ME = ME, U = U, Tsamplefract = .9, show.fig = F, annealing = F,
		                  sb0.fixed = sb0.fixed,
		                  sb.fixed = sb.fixed,
		                  su.fixed = su.fixed,
		                  c.fixed = rep(0, ncol(U)),               
		                  d.fixed = rep(NA, ncol(U)),
		                  b0.start = mod0$b0,
		                  b.start = mod0$b)		
		show(mod)
		
		# the first test is for all variables in U together
		dev <- 2*(mod$logLik - mod0$logLik)
		df <- ncol(U)
		show(c(dev = dev, df = df, P = pchisq(dev, df = df, lower.tail = F)))
		
		# the second test is for the variables in U that aren't in U1
		dev <- 2*(mod$logLik - mod1$logLik)
		df <- ncol(U) - ncol(U1)
		show(c(dev = dev, df = df, P = pchisq(dev, df = df, lower.tail = F)))
		
		
		###############################################
		print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
		print("tests of changes in time")
		
		# pick whichever pair you want to test; the code is set up to test the significance of the
		# variable(s) added to U1 to produce U.
		for(i in c("LL","temp")){
			if(i == "LL"){
				U1 <- U.list[,c("LL", "Age5000")]
				U <- cbind(U1, inter = U.list[,"LL_Age5000"])
			}else{
				U1 <- U.list[,c("temp", "Age5000")]
				U <- cbind(U1, inter = U.list[,"temp_Age5000"])
			}
			mod0 <- TVARSS(X = X, p = p, ME = ME, Tsamplefract = .9, show.fig = F, annealing = F,
			               sb0.fixed = sb0.fixed,
			               sb.fixed = sb.fixed,
			               su.fixed = su.fixed)
			               
			mod1 <- TVARSS(X = X, p = p, ME = ME, U = U1, Tsamplefract = .9, show.fig = F, annealing = F,
			                  sb0.fixed = sb0.fixed,
			                  sb.fixed = sb.fixed,
			                  su.fixed = su.fixed,
			                  c.fixed = rep(0, ncol(U1)),               
			                  d.fixed = rep(NA, ncol(U1)),
			                  b0.start = mod0$b0,
			                  b.start = mod0$b)		
			                  
			mod <- TVARSS(X = X, p = p, ME = ME, U = U, Tsamplefract = .9, show.fig = F, annealing = F,
			                  sb0.fixed = sb0.fixed,
			                  sb.fixed = sb.fixed,
			                  su.fixed = su.fixed,
			                  c.fixed = rep(0, ncol(U)),               
			                  d.fixed = rep(NA, ncol(U)),
			                  b0.start = mod0$b0,
			                  b.start = mod0$b)		
			show(i)
			show(mod)
			
			# the first test is for all variables in U together
			dev <- 2*(mod$logLik - mod0$logLik)
			df <- ncol(U)
			show(c(dev = dev, df = df, P = pchisq(dev, df = df, lower.tail = F)))
			
			# the second test is for the variables in U that are not in U1
			dev <- 2*(mod$logLik - mod1$logLik)
			df <- ncol(U) - ncol(U1)
			show(c(dev = dev, df = df, P = pchisq(dev, df = df, lower.tail = F)))
		}
		
		###############################################
		print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
		print("tests of interactions with a change point")
		
		U1 <- U.list[,c("LL","temp", "Age5000", "LL_temp")]
		U2 <- U.list[,c("LL","temp", "Age5000", "LL_temp", "LL_Age5000", "temp_Age5000")]
		U <- cbind(U2, inter_change = U.list[,"LL_temp_Age5000"])
		
		mod0 <- TVARSS(X = X, p = p, ME = ME, Tsamplefract = .9, show.fig = F, annealing = F,
		               sb0.fixed = sb0.fixed,
		               sb.fixed = sb.fixed,
		               su.fixed = su.fixed)
		               
		mod1 <- TVARSS(X = X, p = p, ME = ME, U = U1, Tsamplefract = .9, show.fig = F, annealing = F,
		                  sb0.fixed = sb0.fixed,
		                  sb.fixed = sb.fixed,
		                  su.fixed = su.fixed,
		                  c.fixed = rep(0, ncol(U1)),               
		                  d.fixed = rep(NA, ncol(U1)),
		                  b0.start = mod0$b0,
		                  b.start = mod0$b)		

		mod2 <- TVARSS(X = X, p = p, ME = ME, U = U2, Tsamplefract = .9, show.fig = F, annealing = F,
		                  sb0.fixed = sb0.fixed,
		                  sb.fixed = sb.fixed,
		                  su.fixed = su.fixed,
		                  c.fixed = rep(0, ncol(U2)),               
		                  d.fixed = rep(NA, ncol(U2)),
		                  b0.start = mod0$b0,
		                  b.start = mod0$b)		
		                  
		mod <- TVARSS(X = X, p = p, ME = ME, U = U, Tsamplefract = .9, show.fig = F, annealing = F,
		                  sb0.fixed = sb0.fixed,
		                  sb.fixed = sb.fixed,
		                  su.fixed = su.fixed,
		                  c.fixed = rep(0, ncol(U)),               
		                  d.fixed = rep(NA, ncol(U)),
		                  b0.start = mod0$b0,
		                  b.start = mod0$b)		
		show(mod)
		
		# test for all variables in the full mod
		dev <- 2*(mod$logLik - mod0$logLik)
		df <- ncol(U)
		show(c(dev = dev, df = df, P = pchisq(dev, df = df, lower.tail = F)))
		
		# test for all variables in U2 that are not in U1
		dev <- 2*(mod2$logLik - mod1$logLik)
		df <- ncol(U2) - ncol(U1)
		show(c(dev = dev, df = df, P = pchisq(dev, df = df, lower.tail = F)))
		
		# test for the variables in U that are not in U2
		dev <- 2*(mod$logLik - mod2$logLik)
		df <- ncol(U) - ncol(U2)
		show(c(dev = dev, df = df, P = pchisq(dev, df = df, lower.tail = F)))
		
		show(c(mod$logLik, mod2$logLik, mod1$logLik, mod0$logLik))
		
		###############################################
		# plot the Kalman-filtered and smoothed data
		fit <- TVARSS_KalmanSmoother(mod)
		
		# blue = Kalman filter, red = Kalman smoother
		pdf(paste0(i.var," figure.pdf"), height = 4, width = 6)
			par(mfrow = c(1,1), mai = c(.8,.8,.1,.1))
			plot(X ~ time, data = fit, xlab = "Time", ylab="data (blk), filter (blu), smooth (red), pred (gr)")
			lines(X ~ time, data = fit, lwd=2)
			lines(X.filtered ~ time, data = fit, col = "blue", lwd=2)
			lines(X.smoothed ~ time, data = fit, col = "red", lwd=2)
			lines(X.mean ~ time, data = fit, col = "green", lwd=2)
			points(X ~ time, data = fit, lty=2)
			points(X.filtered ~ time, data = fit, typ = "l", col = "blue")
			points(X.smoothed ~ time, data = fit, col = "red")
			points(X.mean ~ time, data = fit, col = "green")
		dev.off()
		
		w[[count]] <- data.frame(var = i.var, bin = fit$time, X = fit$X, X.filtered = fit$X.filtered, X.smoothed = fit$X.smoothed)
	}
	return(w)
}
# run analyses
p <- 3
sb.fixed <- matrix(0,1,p)
output <- run(p)
write.csv(output, "isotope Kalmen smoother output p = 3.csv", row.names = F)

p <- 2
sb.fixed <- matrix(0,1,p)
output <- run(p)
write.csv(output, "isotope Kalmen smoother output p = 2.csv", row.names = F)
