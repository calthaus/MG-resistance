# Modeling the spread of macrolide-resistance in Mycoplasma genitalium (MG)
# Christian L. Althaus & Dominique Cadosch, 2 October 2019

# Load libraries
library(deSolve)
library(bbmle)
library(mvtnorm)
library(plotrix)

set.seed(71309)

# Create data sets with MG resistance levels from published studies
france <- data.frame(time =  c(2003:2010 + 0.5, # Chrisment et al. (2012, J Antimicrob Chemother)
							   2011:2012 + 0.5, # Touati et al. (2014, J Clin Microbiol)
							   2013:2014 + 0.5, # Le Roy et al. (2016, Emerg Infect Dis)
							   2016 + 0.5), # Le Roy et al. (2017, J Clin Microbiol)
					 size = c(1, 10, 6, 10, 15, 13, 21, 39,
					 		  69, 65,
					 		  112, 109,
					 		  72),
					 pos = c(0, 0, 0, 1, 2, 2, 3, 5,
					 		 10, 9,
					 		 19, 19,
					 		 6))

denmark <- data.frame(time = c(2007:2010 + 0.5), # Salado-Rasmussen & Jensen (2014, Clin Infect Dis)
					  size = c(11, 226, 378, 454),
					  pos = c(3, 81, 135, 191))

sweden <- data.frame(time = c(2006:2013 + 0.5), # Anagrius et al. (2013, PLOS ONE)
					 size = c(18, 53, 58, 81, 98, 100, 71, 114),
					 pos = c(0, 0, 1, 5, 14, 21, 8, 10))

sets <- list(France = france, Denmark = denmark, Sweden = sweden)
countries <- names(sets)
for(i in countries) {
	li <- numeric()
	ui <- numeric()
	for(j in 1:length(sets[[i]]$time)) {
		ci <- binom.test(sets[[i]]$pos[j], sets[[i]]$size[j])
		li[j] <- ci$conf.int[1]
		ui[j] <- ci$conf.int[2]
	}
	sets[[i]] <- cbind(sets[[i]], li, ui)
}

rm(france, denmark, sweden, li, ui, ci, i, j)

# MG transmission model
model <- function(t, x, parms) {
	with(as.list(c(parms, x)),{
		dS   <- - beta*S*(I_S + I_A + I_T) + gamma*(I_S + I_A + I_T) + tau*(1 - mu)*I_S		# Susceptible individuals
		dI_S <- beta*S*I_S - gamma*I_S - tau*I_S											# Individuals with sensitive strain
		dI_A <- tau*mu*I_S - gamma*I_A														# Individuals with acquired resistance
		dI_T <- beta*S*(I_A + I_T) - gamma*I_T												# Individuals with transmitted resistance
		der <- c(dS, dI_S, dI_A, dI_T)
		list(der)
	})
}

# Log-likelihood function
nll <- function(beta, gamma, mu, tau, intro, country) {
	pars <- c(beta = beta, gamma = gamma, mu = mu, tau = exp(tau))
	set <- sets[[countries[country]]]
	set <- set[order(set$time),]
	intro <- min_time + plogis(intro)*(min(set$time[set$pos > 0]) - min_time)
	times <- sort(c(0, set$time - intro))
	times <- times[times >= 0]
	if(length(times) > 1) {
		simulation <- as.data.frame(ode(init, times, model, parms = pars))
		p <- (simulation$I_A + simulation$I_T)/(simulation$I_S + simulation$I_A + simulation$I_T)
	} else p <- 0
	p <- c(rep(0, 20), p)
	p <- tail(p, length(set$time))
	ll <- sum(dbinom(set$pos, set$size, p, log = TRUE))
	return(-ll)
}

# Fit transmission model to data from France, Denmark and Sweden
min_time <- 1990
max_time <- 2025
init <- c(S = 0.98, I_S = 0.02, I_A = 0, I_T = 0)
fixed <- c(beta = 0.8/(1 - 0.02), gamma = 0.8, mu = 0.12)
free <- c(tau = log(0.4), intro = qlogis(0.5))
fit <- list()
for(i in 1:3) fit[[i]] <- mle2(nll, start = as.list(free), fixed = as.list(c(fixed, country = i)), method = "Nelder-Mead")

# Bootstrap sampling to calculate confidence intervals of model simulations and parameter estimates
n_sim <- 1e4
timepoints <- length(seq(min_time, max_time, 0.1))
sims_resistant <- array(NA, c(3, n_sim, timepoints))
sims_denovo <- array(NA, c(3, n_sim, timepoints))
interval_resistant <- list()
interval_denovo <- list()
interval_estimates <- list()
for(i in 1:3) {
	# Parameter sampling
	m <- coef(fit[[i]], exclude.fixed = TRUE)
	sigma <- vcov(fit[[i]])
	sim_coef <- data.frame(rmvnorm(n_sim, mean = m, sigma = sigma))
	sim_coef$tau <- exp(sim_coef$tau)
	sim_coef$intro <- plogis(sim_coef$intro)
	set <- sets[[i]]	
	sim_coef$intro <- min_time + sim_coef$intro*(min(set$time[set$pos > 0]) - min_time)
	for(j in 1:n_sim) {
		pars <- c(unlist(sim_coef[j, ]), fixed)
		times <- seq(round(pars["intro"], 1), max_time, 0.1)
		simulation <- as.data.frame(ode(init, times, model, parms = pars))
		p_resistant <- c(rep(0, timepoints), (simulation$I_A + simulation$I_T)/(simulation$I_S + simulation$I_A + simulation$I_T))
		p_denovo <- c(rep(1, timepoints), simulation$I_A/(simulation$I_A + simulation$I_T))
		p_denovo <- ifelse(is.nan(p_denovo), 1, p_denovo)
		sims_resistant[i, j, ] <- tail(p_resistant, timepoints)
		sims_denovo[i, j, ] <- tail(p_denovo, timepoints)
	}
	interval_resistant[[i]] <- apply(sims_resistant[i, , ], MAR = 2, FUN = quantile, probs = c(0.025, 0.975))
	interval_denovo[[i]] <- apply(sims_denovo[i, , ], MAR = 2, FUN = quantile, probs = c(0.025, 0.975), na.rm = TRUE)
	interval_estimates[[i]] <- apply(sim_coef, MAR = 2, FUN = quantile, probs = c(0.025, 0.975))
}

# Print parameter estimates
for(i in 1:3) {
	set <- sets[[countries[i]]]
	pars <- coef(fit[[i]])
	pars["tau"] <- exp(pars["tau"])
	intro <- min_time + plogis(pars["intro"])*(min(set$time[set$pos > 0]) - min_time)
	times <- seq(intro, max_time, 0.1)
	simulation <- as.data.frame(ode(init, times, model, parms = pars))
	
	print(pars["tau"])
	print(intro)
	print(interval_estimates[i])
	print(tail((simulation$I_A + simulation$I_T)/(simulation$I_S + simulation$I_A + simulation$I_T), 1))
	print(interval_resistant[[i]][, dim(interval_resistant[[i]])[2]])
}

# Plot model fits (Figure 3)
par(mfrow = c(3, 2))
for(i in 1:3) {
	# Proportion resistance MG
	plot(NA,
		xlim = c(min_time, max_time), ylim = c(0, 1),
		xlab = NA, ylab="Proportion resistant MG", main = countries[i], frame = FALSE)
	plotCI(sets[[i]]$time, sets[[i]]$pos/sets[[i]]$size, ui = sets[[i]]$ui, li = sets[[i]]$li,
		pch = 16,
		add = TRUE)
	polygon(c(seq(min_time, max_time, 0.1), rev(seq(min_time, max_time, 0.1))), c(interval_resistant[[i]][1, ], rev(interval_resistant[[i]][2, ])),
			col = rgb(1, 0, 0, alpha = 0.2), border = NA)
	pars <- coef(fit[[i]])
	pars["tau"] <- exp(pars["tau"])
	set <- sets[[countries[i]]]	
	intro <- min_time + plogis(pars["intro"])*(min(set$time[set$pos > 0]) - min_time)
	times <- seq(intro, max_time, 0.1)
	simulation <- as.data.frame(ode(init, times, model, parms = pars))
	lines(simulation$time, (simulation$I_A + simulation$I_T)/(simulation$I_S + simulation$I_A + simulation$I_T), col = "red")
	pars["mu"] <- 0.01
	simulation <- as.data.frame(ode(init, times, model, parms = pars))
	lines(simulation$time, (simulation$I_A + simulation$I_T)/(simulation$I_S + simulation$I_A + simulation$I_T), lty = 2, col = "red")
	pars["mu"] <- 0.001
	simulation <- as.data.frame(ode(init, times, model, parms = pars))
	lines(simulation$time, (simulation$I_A + simulation$I_T)/(simulation$I_S + simulation$I_A + simulation$I_T), lty = 3, col = "red")
	legend("topleft", c(expression(mu~"= 12%"), expression(mu~"= 1%"), expression(mu~"= 0.1%")), lty = 1:3, col = "red", bty = "n")

	# Proportion de novo resistance
	plot(NA,
		 xlim = c(min_time, max_time), ylim = c(0, 1),
		 xlab = NA, ylab="Proportion de novo resistance", main = countries[i], frame = FALSE)
	polygon(c(seq(min_time, max_time, 0.1), rev(seq(min_time, max_time, 0.1))), c(interval_denovo[[i]][1, ], rev(interval_denovo[[i]][2, ])),
			col = rgb(0, 0, 1, alpha = 0.2), border = NA)
	lines(simulation$time, simulation$I_A/(simulation$I_A + simulation$I_T), col = "blue")
}

# Plot resistance growth rate (Figure 4)
par(mfrow = c(1, 1))
cols <- c("red", "blue", "darkgreen")
plot(NA,
	 xlim = c(0, 1), ylim = c(0, 0.5),
	 xlab = "Proportion of resistant infections", ylab = expression(paste(Relative~growth~rate~of~resistant~infections~(y^-1))), frame = FALSE)
for(i in 1:3) {
	pars <- coef(fit[[i]])
	tau <- exp(pars["tau"])
	mu <- pars["mu"]
	p <- seq(0.01, 1, 0.01)
	lines(p, tau*(1 + mu*(1-p)/p), col = cols[i], lty = i)
	abline(h = tau, col = "gray", lty = i)
}
legend("topright", inset = 0.1, countries, lty = 1:3, col = cols, bty = "n")
