## TVARSS.R  ---  12 Aug 2026 revision
##
## Fork of TVARSS_11Feb25.r by Anthony R. Ives. The ONLY function changed is
## TVARSS_KalmanSmoother(); TVARSS(), summary/print/plot.TVARSS, simulateTS()
## and P_indep_vars_TVARSS() are byte-identical to the 11Feb25 original. The
## Kalman FILTER and the log-likelihood are unchanged, so every fit and every
## AIC from the original file is reproduced exactly.
##
## Why: the Kalman smoother gave wrong answers on series with missing
## observations. On a 73%-missing spore series it returned values only at the
## observed times, and those values had RMSE 4147 against data it should have
## reproduced to four decimals, including negative values for a strictly
## positive concentration. See kalman-smoother/smoother_fix_demo.r for a
## self-contained reproduction on simulated data.
##
## Five changes, each marked CHANGE n at the relevant line:
##
##   1. The Kalman gain used to build L[t] is now formed from the PREDICTED
##      covariance P(t|t-1). Previously PP had already been overwritten with
##      the updated covariance on the preceding line, so L[t] was wrong even
##      with complete data (RMSE 0.738 vs 0.653 for a correct smoother on a
##      simulated fully-observed AR(2)).
##
##   2. v.array, invF.array and L.array are now filled at MISSING times too
##      (v = 0, invF = 0, L = the transition matrix). Previously they kept the
##      zero matrix they were initialised to, which deleted the transition
##      from the backward pass.
##
##   3. The backward recursion runs over every t, not only observed t. This is
##      what makes interpolation across gaps possible at all: the smoothed
##      series used to be NA wherever the data were NA.
##
##   4. L[t] uses the transition mapping t -> t+1, built from the state AFTER
##      updating at t. Previously it used the matrix mapping t-1 -> t. Inert
##      while sb0 = sb = 0 (the b block has zero variance so only the constant
##      B sub-block matters) but wrong as soon as those are freed.
##
##   5. New: the N recursion, returning X.smoothed.se, the SE of the smoothed
##      latent state. The return list is otherwise unchanged, so this file is a
##      drop-in replacement for source().
##
## Note: with initial.points = "observed" the first p-1 rows of the smoothed
## series remain NA, as before; with the default "stationary" every row is
## returned.


## TVARSS.R 

## Time-Varying Autoregressive State Space model
## Copyright 2024 Anthony R. Ives

# This code performs a Kalman filter and Kalman smoother on univariate AR(p) time series in which the autoregression coefficients are allowed to vary as a random walk. The function simulateTS() simulates time series with the option of having breakpoints when the dynamics change. Missing data points in the fitted time series X are replaced with NA, although values of U are required for all data points. 

# Process equation where x[t] is the p x 1 vector of predictions of (X[t], X[t-1], ..., X[t-p]) from x[t-1]

	# x[t] <- b0[t-1]
			 # + b1[t-1] * (x[t-1] - b0[t-1])
			 # + b2[t-1] * (x[t-2] - b0[t-1])
			 # + ... + bp[t-1] * (x[t-p] - b0[t-1]) 
			 # + U[t,] %*% c + rnorm(sd = se)
	
	# b0[t] <- b0[t-1] + rnorm(sd = sb0)
	# bi[t] <- bi[t-1] + rnorm(sd = sbp)

# Updating equation where X[t] is the observation at time t

	# X[t] <- x[t] + U[t,] %*% d + rnorm(sd = su)
	
# Initial estimates of x are either taken from the stationary distribution (x[1]) or from the observed values X[1:p] with options initial.points = "stationary" versus "observed". 

# If real values are given to par.fixed, these values are not estimated; setting any par.fixed = NA makes these parameters estimated.

# Initial regression parameters are obtained from arima(order = c(p,0,0)) applied to the time series X[1:floor(Tsamplefract * Tmax)]. Setting Tsamplefract < 1 is useful if autoregression coefficients are expected to vary (i.e., sb0 and/or sb > 0).

# OPtimization is performed by optim() with the default method "BFGS"
	
# The function simulateTS(B, B0, C, D, X0, U, se, su=0, break.times) simulates the time series in segments defined by break.times. For  given segment j, the rows of B, B0, C, and D correspond to the paramters b, b0, c, and d.

require("GenSA")

TVARSS <- function(X, p = 1, ME = NULL, U = NULL, initial.points = "stationary", b0.start = NA, b.start = array(NA, dim = p), se.start = NA, su.start = 0.01, sb0.start = .05, sb.start = array(0.05, dim = p), c.start = if(is.null(U)) NULL else array(0, dim=dim(U)[2]), d.start = if(is.null(U)) NULL else array(0, dim=dim(U)[2]), b0.fixed = NA, b.fixed = array(NA, dim = p), se.fixed = NA, su.fixed = NA, sb0.fixed = NA, sb.fixed = array(NA, dim = p), c.fixed = if(is.null(U)) NULL else array(NA, dim=dim(U)[2]), d.fixed = if(is.null(U)) NULL else array(0, dim=dim(U)[2]), Tsamplefract = .9, annealing = F, show.count = 10^10, show.fig = T, method = "Nelder-Mead", maxit.BFGS = 10^4, maxit.SANN = 10^2, optim.BFGS.control = NULL, optim.SANN.control = NULL) {
	
	require(GenSA)
	
	####################################################
	# Begin TVARSS.ml
	TVARSS_ml <- function(par, X, ME, U, p, par.fixed, initial.points, show.count = 10^10, show.fig = T, LL.only = T) {

		Tmax <- dim(X)[1]

		par.full <- par.fixed
		par.full[is.na(par.fixed)] <- par

		b0 <- par.full[1]
		b <- par.full[2:(p+1)]
		se <- par.full[p+2]
		su <- par.full[p+3]
		sb <- par.full[(p+4):(p+4+p)]
		if(!is.null(U)){
			nu <- dim(U)[2]
			cc <- matrix(par.full[(p+4+p+1):(p+4+p+nu)], ncol=1)
			dd <- matrix(par.full[(p+4+p+nu+1):(p+4+p+nu+nu)], ncol=1)
		}		

		B0 <- b0		
		B <- as.matrix(rbind(array(0, dim=p), diag(1, p - 1, p)))
		B[1, ] <- b
		
		Se <- diag(0, p)
		Se[1, 1] <- se^2
		
		Su <- su^2	
		Sb <- diag(sb^2)
		
		S <- as.matrix(rbind(cbind(Se, matrix(0, p, p+1)), cbind(matrix(0, p+1, p), Sb)))
		
		Z <- matrix(0, 1, 2*p+1)
		Z[1,1] <- 1

		if(initial.points == "stationary"){
			if(is.null(U)){
				x <- matrix(B0, p, 1)
			}else{	
				x <- matrix(B0 + U[1,] %*% cc, p, 1)
			}
			
			if(rcond(diag(p*p)-kronecker(B,B)) < 1e-12){
				show(paste0("rcond for initial point = ",rcond(diag(p*p)-kronecker(B,B))))
				return(1e10)
			}
			PP <- solve(diag(p*p)-kronecker(B,B)) %*% matrix(Se, nrow = p*p)		
			PP <- matrix(PP, p, p)
			PP <- as.matrix(rbind(cbind(PP, matrix(0, p, p+1)), cbind(matrix(0, p+1, p), Sb)))
			t.start <- 1
		}
		if(initial.points == "observed"){
			if(is.null(U)){
				x <- X[p:1]
			}else{	
				x <- X[p:1] - U[p:1,,drop=F] %*% dd
			}
			
			PP <- Se		
			PP <- as.matrix(rbind(cbind(PP, matrix(0, p, p+1)), cbind(matrix(0, p+1, p), Sb)))
			t.start <- p
		}
		
		# INITIAL UPDATING
		FF <- Z %*% PP %*% t(Z) + Su * ME[t.start]
		invF <- 1/FF
		
		y <- matrix(c(x, B0, B[1,]), ncol=1)	
		if(is.null(U)){
			v <- X[t.start] - Z %*% y
		}else{	
			v <- X[t.start] - (Z %*% y + U[t.start,,drop = F] %*% dd)
		}
		y <- y + PP %*% t(Z) %*% invF %*% v
		PP <- PP - PP %*% t(Z) %*% invF %*% Z %*% PP
		
		x <- y[1:p]
		B0 <- y[p+1]
		B[1,] <- y[(p+2):length(y)]
		
		# START logLik
		logFt <- 0
		vFv <- 0
		if(FF > 0) {
			logFt <- logFt + log(FF)
		}else{
			logFt <- 10^10
		}	
		vFv <- vFv + t(v) %*% invF %*% v

		# information used for plotting
		if(counter %% show.count == 0){
			eiglist <- matrix(NA, Tmax, 1)
			eigB <- eigen(B)$values
			eiglist[p:1] <- max(abs(eigB))
			
			Xlist <- matrix(NA, nrow=Tmax, ncol=1)	
			meanX <- matrix(NA, nrow=Tmax, ncol=1)
			if(is.null(U)){
				Xlist[1:p] <- x[p:1]
				b0cclist <- matrix(NA, nrow=Tmax, ncol=1)
				b0cclist[p:1] <- B0
			}else{	
				Xlist[1:p] <- x[p:1] + U[p:1,,drop=F] %*% dd
				b0cclist <- matrix(NA, nrow=Tmax, ncol=1)
				b0cclist[p:1] <- B0 + U[p:1,,drop = F] %*% cc
				meanX[p:1] <- b0cclist[p:1] + U[p:1,,drop=F] %*% dd

			}
			
			b0list <- matrix(NA, nrow=Tmax, ncol=1)
			b0list[p:1] <- B0
			
			
			blist <- matrix(NA, nrow=Tmax, ncol=p)
			blist[p:1] <- sum(b)
			
			PPlist <- matrix(NA, nrow=Tmax, ncol=dim(PP)[1]^2)
			sdX <- sd(X)
		}
		
		# ITERATIONS
		for(t in (t.start+1):Tmax) {
				
			# PREDICTION EQUATIONS
			
			B12 <- as.matrix(rbind(1 - sum(B[1,]), matrix(0, p - 1, 1)))
			B13 <- as.matrix(rbind(t(x) - B0, matrix(0, p-1, p)))
			
			BB <- as.matrix(rbind(cbind(B, B12, B13), cbind(matrix(0, 1, p), 1, matrix(0, 1, p)), cbind(matrix(0, p, p+1), diag(p))))
			PP <- BB %*% PP %*% t(BB) + S
			x <- B0 + B %*% (x-B0)
			if(!is.null(U)){
				x[1] <- x[1] + U[t,] %*% cc
			}
						
			# UPDATING EQUATIONS
			if(!any(is.na(X[t]))){
				FF <- Z %*% PP %*% t(Z) + Su * ME[t]
				invF <- 1/FF
				
				y <- matrix(c(x, B0, B[1,]), ncol=1)	
				if(is.null(U)){
					v <- X[t] - Z %*% y
				}else{	
					v <- X[t] - (Z %*% y + U[t,,drop = F] %*% dd)
				}
				y <- y + PP %*% t(Z) %*% invF %*% v
				PP <- PP - PP %*% t(Z) %*% invF %*% Z %*% PP
				
				x <- y[1:p]
				B0 <- y[p+1]
				B[1,] <- y[(p+2):length(y)]
				
				# TERMS OF LIKELIHOOD FUNCTION
				if(FF > 0) {
					logFt <- logFt + log(FF)
				}else{
					logFt <- 10^10
				}
				
				vFv <- vFv + t(v) %*% invF %*% v
				
				if(counter %% show.count == 0){
					eigB <- eigen(B)$values
					eiglist[t] <- max(abs(eigB))				
					if(is.null(U)){
						Xlist[t] <- x[1]
						b0cclist[t] <- B0
						meanX[t] <- b0cclist[t]
					}else{	
						Xlist[t] <- x[1] + as.numeric(U[t,,drop=F] %*% dd)
						b0cclist[t] <- B0 + U[t,,drop=F] %*% cc
						meanX[t] <- b0cclist[t] + U[t,,drop=F] %*% dd
					}
					PPlist[t,] <- matrix(PP, nrow=1)
					b0list[t] <- B0
					blist[t,] <- B[1,]
				}
			}
		}
		
		LL <- logFt + vFv
		if(is.complex(LL)) LL <- 10^10

		if((counter %% show.count == 0) & (show.fig == T)){
			npar <- length(par)
			logLik <- -((Tmax - p)/2) * log(2*pi) - LL/2
			
			par.op <- par(no.readonly = TRUE)
			par(mfrow=c(2,1), mar=c(4, 4, 2, .5), cex.axis=1, cex.lab=1, cex.main=1)
			
			plot(1:Tmax, X, typ="p", xlab="Time", ylab="X (o), X.fitted (b), and b0.fitted (g)", main=paste('logLik=', .001*round(1000*logLik)))
			lines(1:Tmax, Xlist, col = "blue")
			if(any(is.na(X))) points(1:Tmax, Xlist, col="blue", pch = 20, cex = .5)
			
			if(is.null(U)){
				lty <- 1				
			}else{
				lty <- 1 + any(cc != 0)				
			}		
			lines(1:Tmax, meanX, col = "green", lty = lty)
			
			main.text <- NULL
			for(i in 0:p) main.text <- paste(main.text, "   sb[", i, "]=", abs(.001*round(1000*sb[i+1])), sep="")
			
			plot(1:Tmax, eiglist, typ="l", col="red", xlab="Time", ylab=expression(lambda), ylim=c(0,1.3), main=eval(expression(main.text)))
			if(any(is.na(X))) points(1:Tmax, eiglist, col="red", pch = 20, cex = .5)
			lines(c(0, Tmax), c(1, 1), col="black")
			
			par(par.op)

		}		
		
		counter <<- counter + 1

		if(LL.only){
			return(LL)
		}else{
			return(list(LL = LL, 
				X.fitted = Xlist,
				X.mean = meanX,
				b0.fitted = b0list,
				b0cc.fitted = b0cclist,
				b.fitted = blist,
				PP.fitted = PPlist,
				eigen.fitted = eiglist))
		}		
		
	}			
	# End TVARSS_ml
	####################################################

	if (var(X, na.rm = TRUE) == 0) {
		stop("The response (dependent variable) has no variation.")
	}
	
	X <- as.matrix(X)
	if(!is.null(U)) {
		U <- as.matrix(U)
		nu <- ncol(U)
	}
	
	if (is.null(ME)) {
		ME <- matrix(1, dim(X))
	}else{
		ME <- as.matrix(ME)
		if (length(ME) != length(X)) {
			stop("The measurement error matrix ME should have the same length as X.")
		}
	}

	Tmax <- dim(X)[1]	
	Tsample <- floor(Tsamplefract * Tmax)
	
	if(length(sb.start) != p) stop("The length of sb.start must be p.")

	ar.init <- arima(X, xreg = U, order = c(p, 0, 0))
	b.init <- ar.init$coef[1:p]
	b0.init <- ar.init$coef[p+1]
	s2 <- ar.init$sigma2
	se.init <- (s2/2)^.5
	su.init <- (s2/2)^.5
	sb0.init <- sb0.start
	sb.init <- sb.start
	if(!is.null(U)) {
		if(any(is.na(c.fixed))){
			#c.init <- ar.init$coef[(p+2):(p+1+nu)]/(1-sum(b.init))
			c.init <- ar.init$coef[(p+2):(p+1+nu)]
			d.init <- 0 * ar.init$coef[(p+2):(p+1+nu)]
		}else{
			c.init <- 0 * ar.init$coef[(p+2):(p+1+nu)]
			d.init <- ar.init$coef[(p+2):(p+1+nu)]
		}
	}else{
		c.init <- NULL
		d.init <- NULL
	}
	par.init <- c(b0.init, b.init, se.init, su.init, sb0.init, sb.init, c.init, d.init) 
	par.start <- c(b0.start, b.start, se.start, su.start, sb0.start, sb.start, c.start, d.start)
	par.fixed <- c(b0.fixed, b.fixed, se.fixed, su.fixed, sb0.fixed, sb.fixed, c.fixed, d.fixed)

	# set up variables for fitting
	par.full <- par.init
	par.full[!is.na(par.start)] <- par.start[!is.na(par.start)]
	
	par.full[!is.na(par.fixed)] <- par.fixed[!is.na(par.fixed)]	
	par <- par.full[is.na(par.fixed)]

	if(is.null(optim.BFGS.control)) optim.BFGS.control = list(maxit = maxit.BFGS)
	if(is.null(optim.SANN.control)) optim.SANN.control = list(temp = 0.1, tmax = 10, maxit = maxit.SANN)
		
	counter <- 1

	if(annealing == T){
		if(!is.null(U)) {
			par.upper <- c(100, array(2, dim=p), 10, 10, 10, array(10, dim=p), array(100, dim=2*nu))
			par.lower <- c(-100, array(-2, dim=p), 0, 0, 0, array(0, dim=p), array(-100, dim=2*nu))
		}else{
			par.upper <- c(100, array(2, dim=p), 10, 10, 10, array(10, dim=p))
			par.lower <- c(-100, array(-2, dim=p), 0, 0, 0, array(0, dim=p))
		}
			
		par.upper <- par.upper[is.na(par.fixed)]
		par.lower <- par.lower[is.na(par.fixed)]
		
		optSANN <- GenSA(fn = TVARSS_ml, par = par, lower = par.lower, upper = par.upper, X = X, U = U, ME = ME, p = p, initial.points = initial.points, par.fixed = par.fixed, show.count = show.count, show.fig = show.fig, control=list(smooth = F, maxit = maxit.SANN))	
		par <- optSANN$par
	}

	opt <- optim(fn = TVARSS_ml, par = par, X = X, U = U, ME = ME, p = p, initial.points = initial.points, par.fixed = par.fixed, show.count = show.count, show.fig = show.fig, method = method, control = optim.BFGS.control)
	
	if(opt$convergence != 0) cat("/nconvergence failed")

	# retrieve final fitted values
	fit <- TVARSS_ml(opt$par, X = X, U = U, ME = ME, p = p, initial.points = initial.points, par.fixed = par.fixed, show.count = 1, show.fig = show.fig, LL.only = F)

	par.full <- par.fixed
	par.full[is.na(par.fixed)] <- opt$par

	b0 <- par.full[1]
	b <- par.full[2:(p+1)]
	se <- abs(par.full[p+2])
	su <- abs(par.full[p+3])
	sb0 <- abs(par.full[p+4])
	sb <- abs(par.full[(p+5):(p+4+p)])
	
	if(any(is.na(c.fixed))){
		cc <- matrix(par.full[(p+4+p+1):(p+4+p+nu)], ncol=1)
	}else{
		cc <- NULL
	}		
	if(any(is.na(d.fixed))){
		dd <- matrix(par.full[(p+4+p+nu+1):(p+4+p+nu+nu)], ncol=1)
	}else{
		dd <- NULL
	}		

	LL <- opt$value
	npar <- length(par)
	logLik <- -((Tmax - p)/2) * log(2*pi) - LL/2
	AIC <- -2*logLik + 2*npar;

	results <- list(X = X, p = p, ME = ME, U = U, se = se, su = su, sb0 = sb0, sb = sb, c = cc, d = dd, b0 = b0, b = b, logLik = logLik, AIC = AIC, npar = npar, initial.points = initial.points, X.fitted = fit$X.fitted, X.mean = fit$X.mean, b0.fitted = fit$b0.fitted, b0cc.fitted = fit$b0cc.fitted, b.fitted = fit$b.fitted, PP.fitted = fit$PP.fitted, eigen.fitted = fit$eigen.fitted, b0.start = b0.start, b.start = b.start, se.start = se.start, su.start = su.start, sb.start = sb.start, c.start = c.start, d.start = d.start, b0.fixed = b0.fixed, b.fixed = b.fixed, se.fixed = se.fixed, su.fixed = su.fixed, sb0.fixed = sb0.fixed, sb.fixed = sb.fixed, c.fixed = c.fixed, d.fixed = d.fixed, opt.par = opt$par, par.full = par.full, par.fixed = par.fixed, Tsamplefract = Tsamplefract, annealing = annealing, method = method, optim.BFGS.control = optim.BFGS.control, optim.SANN.control = optim.SANN.control, opt.convergence = opt$convergence)
	class(results) <- "TVARSS"

	return(results)
}

######################################################
######################################################
# summary.TVARSS
######################################################
######################################################

summary.TVARSS <- function(x, ...) {

	cat("\nCall: TVARSS with p = ", x$p, "\n")

	cat("\nlogLik =", x$logLik)
	cat(",  AIC = ", x$AIC, " [df = ",x$npar, "]\n", sep="")

	cat("\nAutoregression coefficients\n")
	cat("\tb[0] = ", x$b0, "\n", sep="")
	for(i in 1:length(x$b))
		cat("\tb[",i,"] = ", x$b[i], "\n", sep="")

	cat("\nVariation in autoregression coefficients (as standard deviations)\n")
	cat("\tsb[0] = ", x$sb0, "\n", sep="")
	for(i in 1:length(x$sb))
		cat("\tsb[",i,"] = ", x$sb[i], "\n", sep="")

	cat("\nProcess and Measurement errors (as standard deviations)")
	cat("\n\tse =", x$se)
	if(!is.null(x$su)) cat("\n\tsu =", x$su)
	
	if(!is.null(x$c)){
		cat("\n\nCoefficients of c for U\n")
		if(all(is.null(colnames(x$U)))){
			for(i in 1:length(x$c)) cat("\tc[",i-1,"] = ", x$c[i], "\n", sep="")
		}else{	
			for(i in 1:length(x$c)) cat("\t",colnames(x$U)[i], " = ", x$c[i], "\n", sep="")
		}
	}
	
	if(!is.null(x$d)){
		cat("\n\nCoefficients of d for U\n")
		if(all(is.null(colnames(x$U)))){
			for(i in 1:length(x$d)) cat("\td[",i-1,"] = ", x$d[i], "\n", sep="")
		}else{	
			for(i in 1:length(x$d)) cat("\t",colnames(x$U)[i], " = ", x$d[i], "\n", sep="")
		}
	}

	cat("\n")
}
print.TVARSS <- summary.TVARSS

######################################################
######################################################
# plot.TVARSS
######################################################
######################################################

plot.TVARSS <- function(x, ...) {

	####################################################
	# Begin TVARSS.fig
	TVARSS.fig <- function(par, X, ME, U, initial.points, p, par.fixed) {

		Tmax <- dim(X)[1]

		par.full <- par.fixed
		par.full[is.na(par.fixed)] <- par

		b0 <- par.full[1]
		b <- par.full[2:(p+1)]
		se <- par.full[p+2]
		su <- par.full[p+3]
		sb <- par.full[(p+4):(p+4+p)]
		if(!is.null(U)){
			nu <- dim(U)[2]
			cc <- matrix(par.full[(p+4+p+1):(p+4+p+nu)], ncol=1)
			dd <- matrix(par.full[(p+4+p+nu+1):(p+4+p+nu+nu)], ncol=1)
		}		

		B0 <- b0		
		B <- as.matrix(rbind(array(0, dim=p), diag(1, p - 1, p)))
		B[1, ] <- b
		
		Se <- diag(0, p)
		Se[1, 1] <- se^2
		
		Su <- su^2	
		Sb <- diag(sb^2)
		
		S <- as.matrix(rbind(cbind(Se, matrix(0, p, p+1)), cbind(matrix(0, p+1, p), Sb)))
		
		Z <- matrix(0, 1, 2*p+1)
		Z[1,1] <- 1

		if(initial.points == "stationary"){
			if(is.null(U)){
				x <- matrix(B0, p, 1)
			}else{	
				x <- matrix(B0 + U[1,] %*% cc, p, 1)
			}
			
			if(rcond(diag(p*p)-kronecker(B,B)) < 1e-12){
				show(paste0("rcond for initial point = ",rcond(diag(p*p)-kronecker(B,B))))
				return(1e10)
			}
			PP <- solve(diag(p*p)-kronecker(B,B)) %*% matrix(Se, nrow = p*p)		
			PP <- matrix(PP, p, p)
			PP <- as.matrix(rbind(cbind(PP, matrix(0, p, p+1)), cbind(matrix(0, p+1, p), Sb)))
			t.start <- 1
		}
		if(initial.points == "observed"){
			if(is.null(U)){
				x <- X[p:1]
			}else{	
				x <- X[p:1] - U[p:1,,drop=F] %*% dd
			}
			
			PP <- Se		
			PP <- as.matrix(rbind(cbind(PP, matrix(0, p, p+1)), cbind(matrix(0, p+1, p), Sb)))
			t.start <- p
		}
		
		# INITIAL UPDATING
		FF <- Z %*% PP %*% t(Z) + Su * ME[t.start]
		invF <- 1/FF
		
		y <- matrix(c(x, B0, B[1,]), ncol=1)	
		if(is.null(U)){
			v <- X[t.start] - Z %*% y
		}else{	
			v <- X[t.start] - (Z %*% y + U[t.start,,drop = F] %*% dd)
		}
		y <- y + PP %*% t(Z) %*% invF %*% v
		PP <- PP - PP %*% t(Z) %*% invF %*% Z %*% PP
		
		x <- y[1:p]
		B0 <- y[p+1]
		B[1,] <- y[(p+2):length(y)]
		
		# START logLik
		logFt <- 0
		vFv <- 0
		if(FF > 0) {
			logFt <- logFt + log(FF)
		}else{
			logFt <- 10^10
		}	
		vFv <- vFv + t(v) %*% invF %*% v

		# information used for plotting
		eiglist <- matrix(NA, Tmax, 1)
		eigB <- eigen(B)$values
		eiglist[p:1] <- max(abs(eigB))
		
		Xlist <- matrix(NA, nrow=Tmax, ncol=1)	
		meanX <- matrix(NA, nrow=Tmax, ncol=1)
		if(is.null(U)){
			Xlist[1:p] <- x[p:1]
			b0cclist <- matrix(NA, nrow=Tmax, ncol=1)
			b0cclist[p:1] <- B0
		}else{	
			Xlist[1:p] <- x[p:1] + U[p:1,,drop=F] %*% dd
			b0cclist <- matrix(NA, nrow=Tmax, ncol=1)
			b0cclist[p:1] <- B0 + U[p:1,,drop = F] %*% cc
			meanX[p:1] <- b0cclist[p:1] + U[p:1,,drop=F] %*% dd

		}
		
		b0list <- matrix(NA, nrow=Tmax, ncol=1)
		b0list[p:1] <- B0
		
		
		blist <- matrix(NA, nrow=Tmax, ncol=p)
		blist[p:1] <- sum(b)
		
		PPlist <- matrix(NA, nrow=Tmax, ncol=dim(PP)[1]^2)
		sdX <- sd(X)
		
		# ITERATIONS
		for(t in (t.start+1):Tmax) {
				
			# PREDICTION EQUATIONS
			
			B12 <- as.matrix(rbind(1 - sum(B[1,]), matrix(0, p - 1, 1)))
			B13 <- as.matrix(rbind(t(x) - B0, matrix(0, p-1, p)))
			
			BB <- as.matrix(rbind(cbind(B, B12, B13), cbind(matrix(0, 1, p), 1, matrix(0, 1, p)), cbind(matrix(0, p, p+1), diag(p))))
			PP <- BB %*% PP %*% t(BB) + S
			x <- B0 + B %*% (x-B0)
			if(!is.null(U)){
				x[1] <- x[1] + U[t,] %*% cc
			}
						
			# UPDATING EQUATIONS
			if(!any(is.na(X[t]))){
				FF <- Z %*% PP %*% t(Z) + Su * ME[t]
				invF <- 1/FF
				
				y <- matrix(c(x, B0, B[1,]), ncol=1)	
				if(is.null(U)){
					v <- X[t] - Z %*% y
				}else{	
					v <- X[t] - (Z %*% y + U[t,,drop = F] %*% dd)
				}
				y <- y + PP %*% t(Z) %*% invF %*% v
				PP <- PP - PP %*% t(Z) %*% invF %*% Z %*% PP
				
				x <- y[1:p]
				B0 <- y[p+1]
				B[1,] <- y[(p+2):length(y)]
				
				# TERMS OF LIKELIHOOD FUNCTION
				if(FF > 0) {
					logFt <- logFt + log(FF)
				}else{
					logFt <- 10^10
				}
				
				vFv <- vFv + t(v) %*% invF %*% v
				
				eigB <- eigen(B)$values
				eiglist[t] <- max(abs(eigB))				
				if(is.null(U)){
					Xlist[t] <- x[1]
					b0cclist[t] <- B0
					meanX[t] <- b0cclist[t]
				}else{	
					Xlist[t] <- x[1] + as.numeric(U[t,,drop=F] %*% dd)
					b0cclist[t] <- B0 + U[t,,drop=F] %*% cc
					meanX[t] <- b0cclist[t] + U[t,,drop=F] %*% dd
				}
				PPlist[t,] <- matrix(PP, nrow=1)
				b0list[t] <- B0
				blist[t,] <- B[1,]
			}
		}
		
		LL <- logFt + vFv

		npar <- length(par)
		logLik <- -((Tmax - p)/2) * log(2*pi) - LL/2
		
		par.op <- par(no.readonly = TRUE)
		par(mfrow=c(2,1), mar=c(4, 4, 2, .5), cex.axis=1, cex.lab=1, cex.main=1)
		
		plot(1:Tmax, X, typ="p", xlab="Time", ylab="X (o), X.fitted (b), and b0.fitted (g)", main=paste('logLik=', .001*round(1000*logLik)))
		lines(1:Tmax, Xlist, col = "blue")
		if(any(is.na(X))) points(1:Tmax, Xlist, col="blue", pch = 20, cex = .5)
		
		if(is.null(U)){
			lty <- 1				
		}else{
			lty <- 1 + any(cc != 0)				
		}		
		lines(1:Tmax, meanX, col = "green", lty = lty)
		
		main.text <- NULL
		for(i in 0:p) main.text <- paste(main.text, "   sb[", i, "]=", abs(.001*round(1000*sb[i+1])), sep="")
		
		plot(1:Tmax, eiglist, typ="l", col="red", xlab="Time", ylab=expression(lambda), ylim=c(0,1.3), main=eval(expression(main.text)))
		if(any(is.na(X))) points(1:Tmax, eiglist, col="red", pch = 20, cex = .5)
		lines(c(0, Tmax), c(1, 1), col="black")
		
		par(par.op)
	}		

	# End TVARSS.fig
	####################################################

	# plot for par = opt.par
	par.op <- par(no.readonly = TRUE)		
	par(mfrow=c(2,1), mar=c(4, 4, 2, .5), cex.axis=1, cex.lab=1, cex.main=1)

	TVARSS.fig(par=x$opt.par, X=x$X, ME = x$ME, U=x$U, initial.points=x$initial.points, p=x$p, par.fixed=x$par.fixed)
	
	par(par.op)
}

TVARSS_KalmanSmoother <- function(mod, show.fig = F){

	require(MASS)

	# ####################################################
	# # Begin TVARSSsmoother_ml
	TVARSSsmoother_ml <- function(par, X, ME, U, p, par.fixed, initial.points) {

		Tmax <- dim(X)[1]
		m <- 2*p + 1

		par.full <- par.fixed
		par.full[is.na(par.fixed)] <- par

		b0 <- par.full[1]
		b <- par.full[2:(p+1)]
		se <- par.full[p+2]
		su <- par.full[p+3]
		sb <- par.full[(p+4):(p+4+p)]
		if(!is.null(U)){
			nu <- dim(U)[2]
			cc <- matrix(par.full[(p+4+p+1):(p+4+p+nu)], ncol=1)
			dd <- matrix(par.full[(p+4+p+nu+1):(p+4+p+nu+nu)], ncol=1)
		}

		B0 <- b0
		B <- as.matrix(rbind(array(0, dim=p), diag(1, p - 1, p)))
		B[1, ] <- b

		Se <- diag(0, p)
		Se[1, 1] <- se^2

		Su <- su^2
		Sb <- diag(sb^2)

		S <- as.matrix(rbind(cbind(Se, matrix(0, p, p+1)), cbind(matrix(0, p+1, p), Sb)))

		Z <- matrix(0, 1, m)
		Z[1,1] <- 1

		# CHANGE 4: the transition matrix is now built by a helper so that the SAME
		# matrix is used for L[t] and for the prediction step into t+1. BB maps the
		# state at t to the state at t+1, so it must be built from the state AFTER
		# the updating step at t.
		make.BB <- function(x, B0, B){
			B12 <- as.matrix(rbind(1 - sum(B[1,]), matrix(0, p - 1, 1)))
			B13 <- as.matrix(rbind(t(x) - B0, matrix(0, p-1, p)))
			as.matrix(rbind(cbind(B, B12, B13),
			                cbind(matrix(0, 1, p), 1, matrix(0, 1, p)),
			                cbind(matrix(0, p, p+1), diag(p))))
		}

		# All four arrays hold PREDICTED (prior) quantities at t, which is what the
		# Durbin-Koopman backward recursion requires.
		y.array <- array(NA, dim=c(Tmax, m))
		PP.array <- array(0, dim=c(Tmax, m, m))
		v.array <- array(0, dim=Tmax)
		invF.array <- array(0, dim=Tmax)
		L.array <- array(0, dim=c(Tmax, m, m))

		if(initial.points == "stationary"){
			if(is.null(U)){
				x <- matrix(B0, p, 1)
			}else{
				x <- matrix(B0 + U[1,] %*% cc, p, 1)
			}

			if(rcond(diag(p*p)-kronecker(B,B)) < 1e-12){
				stop("TVARSS_KalmanSmoother: the fitted B is too close to non-stationary for the stationary initial condition.")
			}
			PP <- solve(diag(p*p)-kronecker(B,B)) %*% matrix(Se, nrow = p*p)
			PP <- matrix(PP, p, p)
			PP <- as.matrix(rbind(cbind(PP, matrix(0, p, p+1)), cbind(matrix(0, p+1, p), Sb)))
			t.start <- 1
		}
		if(initial.points == "observed"){
			if(is.null(U)){
				x <- X[p:1]
			}else{
				x <- X[p:1] - U[p:1,,drop=F] %*% dd
			}

			PP <- Se
			PP <- as.matrix(rbind(cbind(PP, matrix(0, p, p+1)), cbind(matrix(0, p+1, p), Sb)))
			t.start <- p
		}

		# START logLik
		logFt <- 0
		vFv <- 0
		BB <- NULL

		# ITERATIONS
		# The filter is unchanged from TVARSS_11Feb25.r; only the bookkeeping needed
		# by the smoother has been corrected. The loop now starts at t.start so that
		# the initial time point is handled by the same code as every other one.
		for(t in t.start:Tmax) {

			# PREDICTION EQUATIONS
			# At t.start the prior is the initial condition set above, so there is
			# nothing to predict from.
			if(t > t.start){
				PP <- BB %*% PP %*% t(BB) + S
				x <- B0 + B %*% (x-B0)
				if(!is.null(U)){
					x[1] <- x[1] + U[t,] %*% cc
				}
			}

			y.array[t,] <- c(t(x), B0, B[1,])
			PP.array[t,,] <- PP

			# UPDATING EQUATIONS
			if(!any(is.na(X[t]))){
				FF <- Z %*% PP %*% t(Z) + Su * ME[t]
				invF <- 1/FF

				# CHANGE 1: the Kalman gain is formed from the PREDICTED covariance.
				# The previous version reused PP after it had been overwritten with
				# the updated covariance on the line below, which made L[t] wrong
				# even when no data were missing.
				KK <- PP %*% t(Z) %*% invF

				y <- matrix(c(x, B0, B[1,]), ncol=1)
				if(is.null(U)){
					v <- X[t] - Z %*% y
				}else{
					v <- X[t] - (Z %*% y + U[t,,drop = F] %*% dd)
				}
				y <- y + KK %*% v
				PP <- PP - KK %*% Z %*% PP

				x <- y[1:p]
				B0 <- y[p+1]
				B[1,] <- y[(p+2):length(y)]

				v.array[t] <- v
				invF.array[t] <- invF

				# TERMS OF LIKELIHOOD FUNCTION
				if(FF > 0) {
					logFt <- logFt + log(FF)
				}else{
					logFt <- 10^10
				}

				vFv <- vFv + t(v) %*% invF %*% v
			}else{
				# CHANGE 2: at a missing observation there is no update, so the gain
				# is zero and L[t] reduces to the transition matrix. The previous
				# version left L.array[t,,] as the zero matrix it was initialised
				# to, which silently deleted the transition from the backward pass.
				KK <- matrix(0, m, 1)
				v.array[t] <- 0
				invF.array[t] <- 0
			}

			BB <- make.BB(x, B0, B)
			L.array[t,,] <- BB %*% (diag(m) - KK %*% Z)
		}

		LL <- logFt + vFv
		logLik <- -((Tmax - p)/2) * log(2*pi) - LL/2

		return(list(logLik = logLik,
				Z = Z,
				y.array = y.array,
				v.array = v.array,
				PP.array = PP.array,
				invF.array = invF.array,
				L.array = L.array,
				t.start = t.start))

	}
	# End TVARSSsmoother_ml
	####################################################


	Tmax <- dim(mod$X)[1]
	p <- mod$p
	m <- 2*p + 1

	est <- TVARSSsmoother_ml(par = mod$opt.par, mod$X, mod$ME, mod$U, mod$p, mod$par.fixed, initial.points = mod$initial.points)

	Z <- est$Z
	y <- est$y.array
	v <- est$v.array
	PP <- est$PP.array
	invF <- est$invF.array
	L <- est$L.array
	t.start <- est$t.start

	# CHANGE 3: the backward recursion now runs over EVERY time point, not only the
	# observed ones. At a missing observation the first term vanishes (invF = 0) but
	# the accumulator R must still be carried through the transition by t(L[t,,]),
	# otherwise information from later observations reaches earlier ones undecayed.
	# This is what produced the large errors on series with gaps, and it is why the
	# smoothed series used to be NA wherever the data were NA.
	#
	# CHANGE 5: the companion recursion in N gives the smoothed state covariance, so
	# the smoothed series can be plotted with an uncertainty band. In long gaps this
	# SE approaches the prior SD, which is the honest statement that the smoother
	# knows nothing there.
	ys <- array(NA, dim = dim(y))
	Vs <- array(NA, dim = Tmax)
	R <- matrix(0, m, 1)
	N <- matrix(0, m, m)
	for(t in Tmax:t.start){
		R <- t(Z) %*% invF[t] %*% v[t] + t(L[t,,]) %*% R
		N <- t(Z) %*% (invF[t] * Z) + t(L[t,,]) %*% N %*% L[t,,]
		ys[t,] <- y[t,] + PP[t,,] %*% R
		Vs[t] <- (PP[t,,] - PP[t,,] %*% N %*% PP[t,,])[1,1]
	}

	if(!is.null(mod$d)){
		xs <- ys[,1] + mod$U %*% matrix(mod$d, ncol = 1)
	}else{
		xs <- ys[,1]
	}

	# SE of the smoothed latent state, conditional on the fitted parameters. It does
	# not include uncertainty in b, c, d, se or su.
	xs.se <- sqrt(pmax(Vs, 0))

  time <- 1:Tmax

  if(show.fig){
    par.op <- par(no.readonly = TRUE)
    par(mfrow = c(1,1))
    plot(mod$X ~ time, typ = "p",
         ylim = range(c(mod$X, xs - 1.96*xs.se, xs + 1.96*xs.se), na.rm = TRUE))
    polygon(c(time, rev(time)),
            c(xs - 1.96*xs.se, rev(xs + 1.96*xs.se)),
            col = "#0000ff22", border = NA)
    lines(mod$X.fitted ~ time, col = "lightblue", lwd = 2)
    lines(xs ~ time, col = "blue", lwd = 2)
    points(mod$X ~ time)
    par(par.op)
  }

  return(list(X = mod$X, U = mod$U, X.filtered = mod$X.fitted, X.smoothed = xs, X.smoothed.se = xs.se, X.mean = mod$X.mean, Y.filtered = y, Y.smoothed = ys, time = time, logLik = est$logLik))

}

simulateTS <- function(B, B0, C, D, X0, U, se, su=0, break.times){

	p <- ncol(B)
	Tmax <- max(break.times)
	for(i.segment in 1:length(break.times)){		
		if(i.segment == 1){
			if(is.null(U)){
				x <- X0
			}else{	
				x <- X0 - U[p:1,,drop=F] %*% D[i.segment,,drop = F] 
			}
	
			X <- X0
			meanX <- X0
			t <- p
		}
		while(t < break.times[i.segment]){
			t <- t + 1
			BB <- diag(p+1)[-(p+1),-1]
			BB[1,] <- B[i.segment, , drop = F]
			
			x <- B0[i.segment] + BB %*% (x-B0[i.segment])
			x[1] <- x[1] + U[t,,drop = F] %*% C[i.segment,,drop = F] + rnorm(1, sd = se)

			X <- c(X, x[1] + B0[i.segment] + U[t,,drop = F] %*% D[i.segment,,drop = F] + rnorm(1, sd = su))
			meanX <- c(meanX, (B0[i.segment] + U[t,,drop = F] %*% C[i.segment,,drop = F])/(1-sum(B[i.segment, , drop = F])) + U[t,,drop = F]%*% D[i.segment,,drop = F])
		}
	}
	return(data.frame(time = 1:t, X = X, meanX = meanX))
}

# simple code for producing a list of P-values from LRTs
P_indep_vars_TVARSS <- function(mod){
	out <- data.frame(var = colnames(mod$U))
	if(ncol(mod$U) == 1){
		mod0 <- TVARSS(X = mod$X, ME = mod$ME, p = mod$p, Tsamplefract = .9, show.fig = T, annealing = F,
			initial.points = "stationary",
			su.fixed = mod$su.fixed,
			sb0.fixed = mod$sb0.fixed,
			sb.fixed = mod$sb.fixed,
			b.start = mod$b,
			b0.start = mod$b0	
			)
		dev <- 2*(mod$logLik - mod0$logLik)
		P <- pchisq(dev, df = 1, lower.tail = F)
		if(!is.null(mod$c)) out$val[1] <- mod$c[1]
		if(!is.null(mod$d)) out$val[1] <- mod$d[1]
		out$dev[1] <- dev
		out$P[1] <- P	
			
	}else{
		
		for(i.col in 1:ncol(mod$U)){
			mod0 <- TVARSS(X = mod$X, U = mod$U[,-i.col], ME = mod$ME, p = mod$p, Tsamplefract = .9, show.fig = T, annealing = F,
				initial.points = "stationary",
				su.fixed = mod$su.fixed,
				c.fixed = mod$c.fixed[-i.col],
				d.fixed = mod$d.fixed[-i.col],
				sb0.fixed = mod$sb0.fixed,
				sb.fixed = mod$sb.fixed,
				c.start = if(any(!is.null(mod$c))) mod$c[-i.col] else rep(0,ncol(mod$U) - 1),
				d.start = if(any(!is.null(mod$d))) mod$d[-i.col] else rep(0,ncol(mod$U) - 1),
				b.start = mod$b,
				b0.start = mod$b0	
				)
			dev <- 2*(mod$logLik - mod0$logLik)
			P <- pchisq(dev, df = 1, lower.tail = F)
			if(!is.null(mod$c)) out$val[i.col] <- mod$c[i.col]
			if(!is.null(mod$d)) out$val[i.col] <- mod$d[i.col]
			out$dev[i.col] <- dev
			out$P[i.col] <- P		
		}
	}
	df <- out[,-1]
	rownames(df) <- out[,1]

	return(df)
}
