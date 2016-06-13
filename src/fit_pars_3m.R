## This script is a part of the supporting source code for paper:             ##
## Hu et. al. Computational identification of protein binding sites on        ##
##            RNA using high-throughput RNA structure-probing data            ##
## It is free for non-commercial academic use.                                ##

## Compare three Models

## requre some functions in the mseq package
library(methods) ## a bug in Rscript
source("../src/mseq.R")

## replace the corresponding function in the mseq package
expData <- function(oriData, llen, rlen) {
	cat("Expanding the surrounding sequences...\n")
	seq <- as.character(oriData$seq)
	run_i <- 0
	num_i <- sum(oriData$tag >= 0)
	index <- numeric(num_i)
	sseq <- character(num_i * (llen + rlen))
	count <- numeric(num_i)
	for (i in 1 : length(seq)) {
		if (oriData$tag[i] >= 0) {
			run_i <- run_i + 1
			index[run_i] <- oriData$index[i]
			sseq[((run_i - 1) * (llen + rlen) + 1) : (run_i * (llen + rlen))] <-
				seq[(i - llen) : (i + rlen - 1)] ## check for negative range
			count[run_i] <- oriData$count[i]
		}
	}
	sseq <- factor(sseq)
	cat("set of characters = ", levels(sseq), "\n")
	sseq <- matrix(sseq, ncol = llen + rlen, byrow = TRUE)
	data <- data.frame(index = index, count = count, sseq)
	cname <- character(2 + llen + rlen)
	cname[1] <- "index"
	cname[2] <- "count"
	for (i in 3 : (2 + llen + rlen)) {
		j <- i - 3 - llen
		if (j < 0) {
			cname[i] <- paste("pM", -j, sep = '')
		} else {
			cname[i] <- paste("p", j, sep = '')
		}
	}
	colnames(data) <- cname
	cat("number of genes =", length(unique(data$index)), "\n")
	cat("length of surrounding sequences =", dim(data)[2] - 2, "\n")
	cat("number of counts (positions) =", dim(data)[1], "\n")
	cat("total number of reads =", sum(data$count), "\n")
	return(data)
}

## replace the corresponding function in the mseq package
glmPred <- function(model.glm, newdata) {
	#m <- length(characters)
	oriCoef <- model.glm$coefficients[-1]
	p <- ncol(newdata) - 2 ## remove 'index' and 'count' column
	m <- length(oriCoef) / p + 1 ## m levels can produce m-1 coefs
	coef <- rep(0, p * m)
	for (i in 1 : p) {
		coef[(i - 1) * m + 2 : m] <- oriCoef[(i - 1) * (m-1) + 1 : (m-1)]
	}
	## get the predicted log preferences
	log_pref <- rep(0, nrow(newdata))
	for (i in 1 : p) {
		## Skip 'index' and 'count', so + 2
		log_pref <- log_pref + coef[(i - 1) * m + as.numeric(newdata[, i + 2])]
	}
	return(log_pref)
}

## replace the corresponding function in the mseq package
getPredCount <- function(data, model.glm) {
	uniq_cat <- unique(data$index)
	pred_pref <- exp(glmPred(model.glm, data))
	pred_count <- rep(0, length(data$count))
	for (k in uniq_cat) {
		pred_count[data$index == k] <- sum(data$count[data$index == k]) / 
			sum(pred_pref[data$index == k]) * pred_pref[data$index == k]
	}
	return(pred_count)
}

getPredCountWeight <- function(data, model.glm, rho) {
	uniq_cat <- unique(data$index)
	pred_pref <- exp(glmPred(model.glm, data))
	pred_count <- rep(0, length(data$count))
	for (k in uniq_cat) {
		pred_count[data$index == k] <- 
		    sum(rho[data$index == k] * data$count[data$index == k]) / sum(rho[data$index == k]) *
		    sum(data$index == k) * pred_pref[data$index == k] / sum(pred_pref[data$index == k])
	}
	return(pred_count)
}

##-----------------------
##-- Poisson Linear Model
train.PoissonLinear <- function(data, thrd=1e-4, max_iter=80) {
	cat("Train: Poisson Linear\n")

	## Step 1
	log_expr <- setOriOffset(data)
	for_dev <- 0
	cat("Iter", "mu", "alpha", "dev", "\n", sep="\t")
	for (k in 1:max_iter) {
		## Step 2
		model.glm <- glm(count ~ . - index, data = data, family = poisson(link = "log"), offset = log_expr)

		## Step 3
		log_expr <- updateOffset(data, model.glm$fitted.values, log_expr)

		## Step 4
		now_dev <- getDev(getPredCount(data, model.glm), data$count)
		cat(k, round(mean(exp(log_expr)),3), round(model.glm$coef[1],3), now_dev, "\n", sep = "\t")
		if (k > 1 && abs((now_dev - for_dev) / for_dev) < thrd) { break }
		for_dev <- now_dev
	}
	null_dev <- getDev(getNullCount(data), data$count)
	cat("R-squared =", 1 - now_dev/null_dev,"\n")

	return(list(model.glm, dev= now_dev, dev0= null_dev))
}

cv.PoissonLinear <- function(data, fold=5, thrd=1e-4, max_iter=80, seed=2012) {
	cat("Cross validation: Poisson Linear\n")
	train_index <- CVTrainIndexGene(data, fold, seed)
	
	dev <- rep(0, fold)
	null_dev <- rep(0, fold)
	
	for (j in 1 : fold) {
		train_data <- data[train_index[j, ],]
		test_data <- data[!train_index[j, ],]
		cat("Fold ", j, "of", fold, nrow(train_data), "->", nrow(test_data), "\n")
		## Train
		outcomes <- train.PoissonLinear(train_data, thrd, max_iter)
		## Predict
		dev[j] <- getDev(getPredCount(test_data, outcomes[[1]]), test_data$count)
		null_dev[j] <- getDev(getNullCount(test_data), test_data$count)
	}
	
	cat("deviance =", dev, "\n")
	cat("null deviance =", null_dev, "\n")
	
	R_squared <- 1 - sum(dev) / sum(null_dev)
	cat("R_squared =", R_squared, "\n")
	
	return(R_squared)
}

meanWeight <- function(data, rho) {
	mu <- rep(0, length(data$index))
	for(i in unique(data$index)) {
		mu[data$index == i] <- weighted.mean(
			data$count[data$index == i],
			rho[data$index == i])
	}
	return(mu)
}

getDev <- function(pred_count, real_count) {
	dev <- sum(real_count[real_count != 0] * log(real_count[real_count != 0] / pred_count[real_count != 0]))
	dev <- dev - sum(real_count) + sum(pred_count)
	dev <- 2 * dev
	return(dev)
}

getDevMix <- function(tau0, lambda0, tau1, lambda1, real_count) {
	lambda0 <- lambda0[real_count != 0]
	lambda1 <- lambda1[real_count != 0]
	x <- real_count[real_count != 0]
	a <- log(tau0) + (-lambda0 + x * log(lambda0))
	b <- log(tau1) + (-lambda1 + x * log(lambda1))
	dev <- rep(0, length(x))
	for(i in 1:length(x)) {
		if(abs(a[i] - b[i]) > 100) {
			dev[i] <- max(a[i], b[i])
		}else{
			dev[i] <- a[i] + log(1 + exp(b[i]-a[i]))
		}
	}
	dev <- -2 * sum(dev + x - x * log(x))
	return(dev)
}


##----------------------------
##-- Mixture of Poisson Models 
train.MixturePoisson <- function(data, thrd=1e-4, max_iter=80, rep_seed=1) {
	for(i in 1:rep_seed) {
		if(i == 1)
			best <- train.MixturePoisson_seed(data, thrd=thrd, max_iter=max_iter, seed=1)
		new <- train.MixturePoisson_seed(data, thrd=thrd, max_iter=max_iter, seed=i)
		if(new$r2 > best$r2)
			best <- new
	}
	return(best)
}

train.MixturePoisson_seed <- function(data, thrd=1e-4, max_iter=80, seed=2012) {
	cat("Train: Mixture of Poisson with seed =", seed, "\n")

	## Step 1
	set.seed(seed)
	rho0 <- runif(nrow(data), 0.1, 0.9)
	rho1 <- 1 - rho0

	## Step 2
	tau0 <- mean(rho0)
	tau1 <- 1 - tau0
	
	mu0 <- meanWeight(data, rho0)
	mu1 <- meanWeight(data, rho1)

	for_dev <- 0
	cat("Iter", "tau0", "tau1", "mu0", "mu1", "dev", "\n", sep="\t")
	for (k in 1:max_iter) {
		## Step 3
		log_pred_prob0 <- dpois(data$count, mu0, log=TRUE)
		log_pred_prob1 <- dpois(data$count, mu1, log=TRUE)

		rho0 <- tau0 / (tau0  + tau1 * exp(log_pred_prob1 - log_pred_prob0))
		rho1 <- 1 - rho0

		## Step 4
		mu0 <- meanWeight(data, rho0)
		mu1 <- meanWeight(data, rho1)

		## Step 5
		swap <- mu0 > mu1 ## mu0 should less than mu1, otherwise swap
		if(sum(swap) > 0) {
		cat("Swap the labels of", length(unique(data$index[swap])), "genes\n")
		tmp <- rho0[swap]; rho0[swap] <- rho1[swap]; rho1[swap] <- tmp
		tmp <- mu0[swap]; mu0[swap] <- mu1[swap]; mu1[swap] <- tmp
		}

		## Step 6
		tau0 <- mean(rho0)
		tau1 <- 1 - tau0

		## Step 7
		now_dev <- getDevMix(tau0, mu0, tau1, mu1, data$count)
		cat(k, round(tau0,3), round(tau1,3), round(mean(mu0),3), round(mean(mu1),3), now_dev, "\n", sep = "\t")
		if (k > 3 && abs((now_dev - for_dev) / for_dev) < thrd) { break }
		for_dev <- now_dev
	}
	null_dev <- getDev(getNullCount(data), data$count)
	r_squared = 1 - now_dev/null_dev
	cat("R-squared =", r_squared,"\n")

	return(list(r2= r_squared, dev= now_dev, dev0= null_dev, prob= rho1, pred= tau0*mu0 + tau1*mu1))
}

cv.MixturePoisson <- function(data, fold=5, thrd=1e-4, max_iter=80, seed=2012) {
	cat("Cross validation using Mixture of Poisson ...\n")
	train_index <- CVTrainIndexGene(data, fold, seed)
	dev <- rep(0, fold)
	null_dev <- rep(0, fold)
	for (j in 1 : fold) {
		test_data <- data[!train_index[j, ],]
		cat("Fold ", j, "of", fold, "0 ->", nrow(test_data), "\n")
		## Predict
		outcomes <- train.MixturePoisson(test_data, thrd, max_iter, rep_seed=1)
		dev[j] <- outcomes$dev
		null_dev[j] <- outcomes$dev0
	}
	
	cat("deviance =", dev, "\n")
	cat("null deviance =", null_dev, "\n")
	
	R_squared <- 1 - sum(dev) / sum(null_dev)
	cat("R_squared =", R_squared, "\n")
	
	return(R_squared)
}


##------------------------------------
##-- Mixture of Poisson Linear Models
train.MixturePoissonLinear <- function(data, thrd=1e-4, max_iter=80, rep_seed=1) {
	for(i in 1:rep_seed) {
		if(i == 1)
			best <- train.MixturePoissonLinear_seed(data, thrd=thrd, max_iter=max_iter, seed=1)
		new <- train.MixturePoissonLinear_seed(data, thrd=thrd, max_iter=max_iter, seed=i)
		if(new$r2 > best$r2)
			best <- new
	}
	return(best)
}

train.MixturePoissonLinear_seed <- function(data, thrd=1e-4, max_iter=80, seed=2012) {
	cat("Train: Mixture of Poisson Linear with seed =", seed, "\n")

	## Step 1
	set.seed(seed)
	rho0 <- runif(nrow(data), 0.1, 0.9)
	rho1 <- 1 - rho0

	## Step 2
	tau0 <- mean(rho0)
	tau1 <- 1 - tau0

	log_expr0 <- setOffsetWeight(data, rho0)
	log_expr1 <- setOffsetWeight(data, rho1)

	for_dev <- 0
	cat("Iter", "tau0", "tau1", "mu0", "mu1", "alpha0", "alpha1", "dev", "\n", sep="\t")
	for (k in 1:max_iter) {
		## Step 3
		model.glm0 <- glm(count ~ . - index, data = data, family = poisson(link = "log"), offset = log_expr0, weights = rho0)
		model.glm1 <- glm(count ~ . - index, data = data, family = poisson(link = "log"), offset = log_expr1, weights = rho1)

		## Step 4
		log_pref0 <- log_expr0 + model.glm0$coef[1] + glmPred(model.glm0, data)
		log_pref1 <- log_expr1 + model.glm1$coef[1] + glmPred(model.glm1, data)

		log_pred_prob0 <- dpois(data$count, exp(log_pref0), log=TRUE)
		log_pred_prob1 <- dpois(data$count, exp(log_pref1), log=TRUE)

		rho0 <- tau0 / (tau0 + tau1 * exp(log_pred_prob1 - log_pred_prob0))
		rho1 <- 1 - rho0

		## Step 5
		log_expr0 <- updateOffsetWeight(data, log_pref0, log_expr0, rho0)
		log_expr1 <- updateOffsetWeight(data, log_pref1, log_expr1, rho1)

		## Step 6
		swap <- log_expr0 > log_expr1 ## log_expr0 should less than log_expr1, otherwise swap
		if(sum(swap) > 0) {
		cat("Swap the labels of", length(unique(data$index[swap])), "genes\n")
		tmp <- rho0[swap]; rho0[swap] <- rho1[swap]; rho1[swap] <- tmp
		tmp <- log_expr0[swap]; log_expr0[swap] <- log_expr1[swap]; log_expr1[swap] <- tmp
		}

		## Step 7
		tau0 <- mean(rho0)
		tau1 <- 1 - tau0

		## Step 8
		now_dev <- getDevMix(tau0, getPredCountWeight(data, model.glm0, rho0), tau1, getPredCountWeight(data, model.glm1, rho1), data$count)
		cat(k, round(tau0,3), round(tau1,3), round(mean(exp(log_expr0)),2), round(mean(exp(log_expr1)),2), round(model.glm0$coef[1],3), round(model.glm1$coef[1],3), now_dev, "\n", sep = "\t")
		if (k > 1 && abs((now_dev - for_dev) / for_dev) < thrd) { break }
		for_dev <- now_dev
	}

	null_dev <- getDev(getNullCount(data), data$count)
	r_squared <- 1 - now_dev/null_dev
	cat("R-squared =", r_squared,"\n")

	return(list(model.glm0, model.glm1, tau0, r2=r_squared, dev= now_dev, dev0= null_dev))
}

predict.MixturePoissonLinear <- function(model, data, thrd=1e-4, max_iter=80, seed=2012) {
	cat("Predict: Mixture of Poisson Linear\n")

	model.glm0 <- model[[1]]
	model.glm1 <- model[[2]]
	tau0 <- model[[3]]

	## Step 1
	set.seed(seed)
	rho0 <- runif(nrow(data), 0.1, 0.9)
	rho1 <- 1 - rho0

	## Step 2
	#tau0 <- mean(rho0)
	tau1 <- 1 - tau0

	log_expr0 <- setOffsetWeight(data, rho0)
	log_expr1 <- setOffsetWeight(data, rho1)

	for_dev <- 0
	cat("Iter", "tau0", "tau1", "mu0", "mu1", "alpha0", "alpha1", "dev", "\n", sep="\t")
	for (k in 1:max_iter) {
		## Step 3
		#model.glm0 <- glm(count ~ . - index, data = data, family = poisson(link = "log"), offset = log_expr0, weights = rho0)
		#model.glm1 <- glm(count ~ . - index, data = data, family = poisson(link = "log"), offset = log_expr1, weights = rho1)

		## Step 4
		log_pref0 <- log_expr0 + model.glm0$coef[1] + glmPred(model.glm0, data)
		log_pref1 <- log_expr1 + model.glm1$coef[1] + glmPred(model.glm1, data)

		log_pred_prob0 <- dpois(data$count, exp(log_pref0), log=TRUE)
		log_pred_prob1 <- dpois(data$count, exp(log_pref1), log=TRUE)

		rho0 <- tau0 / (tau0 + tau1 * exp(log_pred_prob1 - log_pred_prob0))
		rho1 <- 1 - rho0

		## Step 5
		log_expr0 <- updateOffsetWeight(data, log_pref0, log_expr0, rho0)
		log_expr1 <- updateOffsetWeight(data, log_pref1, log_expr1, rho1)

		## Step 6
		swap <- log_expr0 > log_expr1 ## log_expr0 should less than log_expr1, otherwise swap
		if(sum(swap) > 0) {
		cat("Swap the labels of", length(unique(data$index[swap])), "genes\n")
		tmp <- rho0[swap]; rho0[swap] <- rho1[swap]; rho1[swap] <- tmp
		tmp <- log_expr0[swap]; log_expr0[swap] <- log_expr1[swap]; log_expr1[swap] <- tmp
		}

		## Step 7
		#tau0 <- mean(rho0)
		#tau1 <- 1 - tau0

		## Step 8
		now_dev <- getDevMix(tau0, getPredCountWeight(data, model.glm0, rho0), tau1, getPredCountWeight(data, model.glm1, rho1), data$count)
		cat(k, round(tau0,3), round(tau1,3), round(mean(exp(log_expr0)),2), round(mean(exp(log_expr1)),2), round(model.glm0$coef[1],3), round(model.glm1$coef[1],3), now_dev, "\n", sep = "\t")
		if (k > 1 && abs((now_dev - for_dev) / for_dev) < thrd) { break }
		for_dev <- now_dev
	}
	null_dev <- getDev(getNullCount(data), data$count)
	cat("R-squared =", 1 - now_dev/null_dev,"\n")

	pred_count <- rho0 * getPredCountWeight(data, model.glm0, rho0) + rho1 * getPredCountWeight(data, model.glm1, rho1)
	return(list(dev= now_dev, dev0= null_dev, prob= rho1, pred= pred_count))
}



cv.MixturePoissonLinear <- function(data, fold=5, thrd=1e-4, max_iter=80, seed=2012) {
	cat("Cross validation: Mixture of Poisson Linear\n")
	train_index <- CVTrainIndexGene(data, fold, seed)
	
	dev <- rep(0, fold)
	null_dev <- rep(0, fold)
	
	for (j in 1 : fold) {
		train_data <- data[train_index[j, ],]
		test_data <- data[!train_index[j, ],]
		cat("Fold ", j, "of", fold, nrow(train_data), "->", nrow(test_data), "\n")
		## Train
		model <- train.MixturePoissonLinear(train_data, thrd, max_iter)
		## Predict
		outcomes <- predict.MixturePoissonLinear(model, test_data, thrd, max_iter)
		dev[j] <- outcomes$dev
		null_dev[j] <- outcomes$dev0
	}
	
	cat("deviance =", dev, "\n")
	cat("null deviance =", null_dev, "\n")
	
	R_squared <- 1 - sum(dev) / sum(null_dev)
	cat("R_squared =", R_squared, "\n")
	
	return(R_squared)
}

setOffsetWeight <- function(data, rho) {
	uniq_cat <- unique(data$index)
	log_expr <- rep(0, length(data$index))
	
	for (k in uniq_cat) {
		log_expr[data$index == k] <- log(
			sum(rho[data$index ==k] * data$count[data$index == k]) /
			sum(rho[data$index == k]))
	}
	return(log_expr)
}

updateOffsetWeight <- function(data, log_pred_count, log_expr, rho) {
	uniq_cat <- unique(data$index)
	log_pref <- log_pred_count - log_expr
	
	for (k in uniq_cat) {
		log_expr[data$index == k] <- log(
			sum(rho[data$index == k] * data$count[data$index == k]) /
			sum(rho[data$index == k] * exp(log_pref[data$index == k])))
	}
	return(log_expr)
}

## generate a new random data set
generateData <- function(oldData, tau=0.7, seed=2011) {
	set.seed(seed)
	newData <- oldData ## a new copy
	for (i in unique(oldData$index)) {
		geneIndex <- oldData$index == i
		counts <- oldData$count[geneIndex]
		mu0 <- mean(counts[counts <= mean(counts)])
		mu1 <- 2 * mu0
		num0 <- round(sum(geneIndex) * tau)
		num1 <- sum(geneIndex) - num0
		newData$count[geneIndex] <- sample(c(rpois(num0, mu0), rpois(num1, mu1)))
	}
	return(newData)
}

## random the local sequence
randomSeq <- function(data, seed=2011) {
	set.seed(seed)
	index <- seq(1, nrow(data))
	index <- sample(index)
	data$seq <- data$seq[index]
	return(data)
}

##--------
##-- laboring function
fit_pars <- function(data_file, ratio= -1) {
	data_table <- read.csv(data_file)

	## remove elements tagged with 0, which is the sign of low expression genes
	data_table <- data_table[data_table$tag != 0, ]

	## get the maximum span size
	tag_list <- data_table$tag[data_table$index == data_table$index[1]]
	max_span <- 0
	for(one_tag in tag_list) {
		if(one_tag > 0) {
			break
		}
		max_span <- max_span + 1
	}
	cat("The maximum span size is", max_span, "\n")

	## generate a random data set from current data
	if(ratio >= 0) { ## use synthetic
		data_table <- generateData(data_table, ratio)
		## default span size for testing
		left_span <- 2
		right_span <- 2
		encoded <- expData(data_table, left_span, right_span)
		R_sq1 <- cv.PoissonLinear(encoded)
		R_sq2 <- cv.MixturePoisson(encoded)
		R_sq3 <- cv.MixturePoissonLinear(encoded)

		## save result
		write(paste(ratio, left_span, right_span, R_sq1, R_sq2, R_sq3, sep="\t"), file="log.txt", append=TRUE)
	}else{ ## use real
	## find the best region of surrounding sequences
	write(data_file, file="log.txt", append=TRUE)
	write("Ratio\tLeftSpan\tRightSpan\tPL\tMP\tMPL", file="log.txt", append=TRUE)
	for(left_span in seq(from= 1, to= max_span, by= 1)) {
		right_span <- left_span
#	for(right_span in seq(10,100,10)) {
		cat(left_span, right_span, "\n")

		#data_table <- randomSeq(data_table)
		encoded <- expData(data_table, left_span, right_span)

		R_sq1 <- cv.PoissonLinear(encoded)
		R_sq2 <- cv.MixturePoisson(encoded)
		R_sq3 <- cv.MixturePoissonLinear(encoded)

		## save result
		write(paste(ratio, left_span, right_span, R_sq1, R_sq2, R_sq3, sep="\t"), file="log.txt", append=TRUE)
	}
	}
}

### MAIN 
main <- function() 
{
	## use synthetic data 
	write("Synthetic data PARS_V1", file="log.txt", append=TRUE)
	write("Ratio\tLeftSpan\tRightSpan\tPL\tMP\tMPL", file="log.txt", append=TRUE)
	for(ratio in seq(0, 1, 0.1))
		fit_pars("../work/PARS_V1.csv", ratio)

	## use real data
	write("Real data PARS_V1", file="log.txt", append=TRUE)
	fit_pars("../work/PARS_V1.csv")
	write("Read data PARS_S1", file="log.txt", append=TRUE)
	fit_pars("../work/PARS_S1.csv")
}

## only run when there is no need to be reloaded 
if(length(commandArgs(TRUE)) == 0){main()}
