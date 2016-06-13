## This script is a part of the supporting source code for paper:             ##
## Hu et. al. Computational identification of protein binding sites on        ##
##            RNA using high-throughput RNA structure-probing data            ##
## It is free for non-commercial academic use.                                ##

## requre some functions in the mseq package
library(methods) ## a bug in Rscript
source("../src/fit_pars_3m.R")

## replace the corresponding function in the mseq package
expData <- function(oriData, choose, llen, rlen) {
	cat("Expanding the surrounding sequences...\n")
	seq <- as.character(oriData$seq)
	run_i <- 0
	num_i <- sum(choose)
	index <- numeric(num_i)
	sseq <- character(num_i * (llen + rlen))
	count <- numeric(num_i)
	for (i in 1 : length(seq)) {
		if (choose[i]) {
			run_i <- run_i + 1
			index[run_i] <- oriData$index[i]
			sseq[((run_i - 1) * (llen + rlen) + 1) : (run_i * (llen + rlen))] <- seq[(i - llen) : (i + rlen - 1)]
			count[run_i] <- rep(0, length(oriData$index[i])) ## skip the count data, processed other ways
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
	oriCoef[is.na(oriCoef)] <- 0 ## ignore NA coef
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

meanWeight <- function(data, rho) {
	mu <- rep(0, length(data$index))
	for(i in unique(data$index)) {
		mu[data$index == i] <- weighted.mean(data$count[data$index == i], rho[data$index == i])
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
		if(abs(a[i] - b[i]) > 100) { ## approximate
			dev[i] <- max(a[i], b[i])
		}else{
			dev[i] <- a[i] + log(1 + exp(b[i]-a[i]))
		}
	}
	dev <- -2 * sum(dev + x - x * log(x))
	return(dev)
}

#
#
train.MixturePoissonLinearCombine <- function(data, same_dir=TRUE, thrd=1e-4, max_iter=80, rep_seed=1) {
	for(i in 1:rep_seed) {
		if (i == 1)
			best <- train.MixturePoissonLinearCombine_seed(data, same_dir, thrd, max_iter, seed=1)
		new <- train.MixturePoissonLinearCombine_seed(data, same_dir, thrd, max_iter, seed=i)
		if(new$r2 > best$r2)
			best <- new
	}
	return(best)
}

train.MixturePoissonLinearCombine_seed <- function(data, same_dir=TRUE, thrd=1e-4, max_iter=80, seed=2012) {
	cat("Train: Mixture of Poisson Linear Combine with seed =", seed, "\n")

	## Step 1
	set.seed(seed)
	rho0 <- runif(nrow(data$v1), 0.1, 0.9)
	rho1 <- 1 - rho0

	## Step 2
	tau0 <- mean(rho0)
	tau1 <- 1 - tau0

	log_expr0v <- setOffsetWeight(data$v1, rho0)
	log_expr1v <- setOffsetWeight(data$v1, rho1)
	log_expr0s <- setOffsetWeight(data$s1, rho0)
	log_expr1s <- setOffsetWeight(data$s1, rho1)

	for_dev <- 0
	cat("Iter", "tau0", "tau1", "mu0v", "mu1v", "mu0s", "mu1s", "dev", "\n", sep="\t")
	for (k in 1 : max_iter) {
		## Step 3
		model.glm0v <- glm(count ~ .-index, data = data$v1, family = poisson(link = "log"), offset = log_expr0v, weights = rho0)
		model.glm1v <- glm(count ~ .-index, data = data$v1, family = poisson(link = "log"), offset = log_expr1v, weights = rho1)
		model.glm0s <- glm(count ~ .-index, data = data$s1, family = poisson(link = "log"), offset = log_expr0s, weights = rho0)
		model.glm1s <- glm(count ~ .-index, data = data$s1, family = poisson(link = "log"), offset = log_expr1s, weights = rho1)

		## Step 4
		log_pref0v <- log_expr0v + model.glm0v$coef[1] + glmPred(model.glm0v, data$v1)
		log_pref1v <- log_expr1v + model.glm1v$coef[1] + glmPred(model.glm1v, data$v1)
		log_pref0s <- log_expr0s + model.glm0s$coef[1] + glmPred(model.glm0s, data$s1)
		log_pref1s <- log_expr1s + model.glm1s$coef[1] + glmPred(model.glm1s, data$s1)

		log_pred_prob0v <- dpois(data$v1$count, exp(log_pref0v), log=TRUE)
		log_pred_prob1v <- dpois(data$v1$count, exp(log_pref1v), log=TRUE)
		log_pred_prob0s <- dpois(data$s1$count, exp(log_pref0s), log=TRUE)
		log_pred_prob1s <- dpois(data$s1$count, exp(log_pref1s), log=TRUE)

		rho0v <- 1 / (1 + (tau1/tau0) * exp(log_pred_prob1v - log_pred_prob0v))
		rho0s <- 1 / (1 + (tau1/tau0) * exp(log_pred_prob1s - log_pred_prob0s))
		rho1v <- 1 - rho0v
		rho1s <- 1 - rho0s

		## Step 5
		log_expr0v <- updateOffsetWeight(data$v1, log_pref0v, log_expr0v, rho0)
		log_expr1v <- updateOffsetWeight(data$v1, log_pref1v, log_expr1v, rho1)
		log_expr0s <- updateOffsetWeight(data$s1, log_pref0s, log_expr0s, rho0)
		log_expr1s <- updateOffsetWeight(data$s1, log_pref1s, log_expr1s, rho1)

		## Step 6
		swap <- log_expr0v > log_expr1v
		if(sum(swap) > 0) {
		cat("Swap the labels of", length(unique(data$v1$index[swap])), "genes in V1\n")
		tmp <- rho0v[swap]; rho0v[swap] <- rho1v[swap]; rho1v[swap] <- tmp
		}
		if(same_dir) {
		swap <- log_expr0s > log_expr1s
		}else{ ## opposite direction
		swap <- log_expr0s < log_expr1s
		}
		if(sum(swap) > 0) {
		cat("Swap the labels of", length(unique(data$s1$index[swap])), "genes in S1\n")
		tmp <- rho0s[swap]; rho0s[swap] <- rho1s[swap]; rho1s[swap] <- tmp
		}
		
		## Step 7
		rho0 <- (rho0v+rho0s)/2
		rho1 <- 1 - rho0

		## Step 8
		tau0 <- mean(rho0)
		tau1 <- 1 - tau0

		log_expr0v <- updateOffsetWeight(data$v1, log_pref0v, log_expr0v, rho0)
		log_expr1v <- updateOffsetWeight(data$v1, log_pref1v, log_expr1v, rho1)
		log_expr0s <- updateOffsetWeight(data$s1, log_pref0s, log_expr0s, rho0)
		log_expr1s <- updateOffsetWeight(data$s1, log_pref1s, log_expr1s, rho1)

		## Step 9
		now_dev <- 
		  getDevMix(tau0, getPredCountWeight(data$v1, model.glm0v, rho0), tau1, getPredCountWeight(data$v1, model.glm1v, rho1), data$v1$count) +
		  getDevMix(tau1, getPredCountWeight(data$s1, model.glm0s, rho0), tau0, getPredCountWeight(data$s1, model.glm1s, rho1), data$s1$count)
		cat(k, round(tau0,3), round(tau1,3), round(mean(exp(log_expr0v))), round(mean(exp(log_expr1v))), round(mean(exp(log_expr0s))), round(mean(exp(log_expr1s))), now_dev, "\n", sep = "\t")
		if (k > 3 && (for_dev - now_dev) / for_dev < thrd) { break }
		for_dev <- now_dev
	}
	null_dev <- getDev(getNullCount(data$v1), data$v1$count) + getDev(getNullCount(data$s1), data$s1$count)
	r_squared <- 1 - now_dev/null_dev
	cat("R-squared =", r_squared,"\n")
	return(list(model.glm0v, model.glm1v, model.glm0s, model.glm1s, tau0, r2= r_squared, dev= now_dev, dev0= null_dev))
}

predict.MixturePoissonLinearCombine <- function(model, data, same_dir=TRUE, thrd=1e-4, max_iter=80, seed=2012) {
	cat("Predict: Mixture of Poisson Linear Combine\n")
	model.glm0v <- model[[1]]
	model.glm1v <- model[[2]]
	model.glm0s <- model[[3]]
	model.glm1s <- model[[4]]
	tau0 <- model[[5]]

	## Step 1
	set.seed(seed)
	rho0 <- runif(nrow(data$v1), 0.1, 0.9)
	rho1 <- 1 - rho0

	## Step 2
	#tau0 <- mean(rho0)
	tau1 <- 1 - tau0

	log_expr0v <- setOffsetWeight(data$v1, rho0)
	log_expr1v <- setOffsetWeight(data$v1, rho1)
	log_expr0s <- setOffsetWeight(data$s1, rho0)
	log_expr1s <- setOffsetWeight(data$s1, rho1)

	for_dev <- 0
	cat("Iter", "tau0", "tau1", "mu0v", "mu1v", "mu0s", "mu1s", "dev", "\n", sep="\t")
	for (k in 1 : max_iter) {
		## Step 3
		#model.glm0v <- glm(count ~ .-index, data = data$v1, family = poisson(link = "log"), offset = log_expr0v, weights = rho0)
		#model.glm1v <- glm(count ~ .-index, data = data$v1, family = poisson(link = "log"), offset = log_expr1v, weights = rho1)
		#model.glm0s <- glm(count ~ .-index, data = data$s1, family = poisson(link = "log"), offset = log_expr0s, weights = rho0)
		#model.glm1s <- glm(count ~ .-index, data = data$s1, family = poisson(link = "log"), offset = log_expr1s, weights = rho1)

		## Step 4
		log_pref0v <- log_expr0v + model.glm0v$coef[1] + glmPred(model.glm0v, data$v1)
		log_pref1v <- log_expr1v + model.glm1v$coef[1] + glmPred(model.glm1v, data$v1)
		log_pref0s <- log_expr0s + model.glm0s$coef[1] + glmPred(model.glm0s, data$s1)
		log_pref1s <- log_expr1s + model.glm1s$coef[1] + glmPred(model.glm1s, data$s1)

		log_pred_prob0v <- dpois(data$v1$count, exp(log_pref0v), log=TRUE)
		log_pred_prob1v <- dpois(data$v1$count, exp(log_pref1v), log=TRUE)
		log_pred_prob0s <- dpois(data$s1$count, exp(log_pref0s), log=TRUE)
		log_pred_prob1s <- dpois(data$s1$count, exp(log_pref1s), log=TRUE)

		rho0v <- 1 / (1 + (tau1/tau0) * exp(log_pred_prob1v - log_pred_prob0v))
		rho0s <- 1 / (1 + (tau1/tau0) * exp(log_pred_prob1s - log_pred_prob0s))
		rho1v <- 1 - rho0v
		rho1s <- 1 - rho0s

		## Step 5
		log_expr0v <- updateOffsetWeight(data$v1, log_pref0v, log_expr0v, rho0)
		log_expr1v <- updateOffsetWeight(data$v1, log_pref1v, log_expr1v, rho1)
		log_expr0s <- updateOffsetWeight(data$s1, log_pref0s, log_expr0s, rho0)
		log_expr1s <- updateOffsetWeight(data$s1, log_pref1s, log_expr1s, rho1)

		## Step 6
		swap <- log_expr0v > log_expr1v
		if(sum(swap) > 0) {
		cat("Swap the labels of", length(unique(data$v1$index[swap])), "genes in V1\n")
		tmp <- rho0v[swap]; rho0v[swap] <- rho1v[swap]; rho1v[swap] <- tmp
		}
		if(same_dir) {
		swap <- log_expr0s > log_expr1s
		}else{ ## opposite direction
		swap <- log_expr0s < log_expr1s
		}
		if(sum(swap) > 0) {
		cat("Swap the labels of", length(unique(data$s1$index[swap])), "genes in S1\n")
		tmp <- rho0s[swap]; rho0s[swap] <- rho1s[swap]; rho1s[swap] <- tmp
		}
		
		## Step 7
		rho0 <- (rho0v+rho0s)/2
		rho1 <- 1 - rho0

		## Step 8
		#tau0 <- mean(rho0)
		#tau1 <- 1 - tau0

		log_expr0v <- updateOffsetWeight(data$v1, log_pref0v, log_expr0v, rho0)
		log_expr1v <- updateOffsetWeight(data$v1, log_pref1v, log_expr1v, rho1)
		log_expr0s <- updateOffsetWeight(data$s1, log_pref0s, log_expr0s, rho0)
		log_expr1s <- updateOffsetWeight(data$s1, log_pref1s, log_expr1s, rho1)

		## Step 9
		now_dev_v <- getDevMix(tau0, getPredCountWeight(data$v1, model.glm0v, rho0), tau1, getPredCountWeight(data$v1, model.glm1v, rho1), data$v1$count)
		now_dev_s <- getDevMix(tau1, getPredCountWeight(data$s1, model.glm0s, rho0), tau0, getPredCountWeight(data$s1, model.glm1s, rho1), data$s1$count)
		cat(k, round(tau0,3), round(tau1,3), round(mean(exp(log_expr0v))), round(mean(exp(log_expr1v))), round(mean(exp(log_expr0s))), round(mean(exp(log_expr1s))), now_dev_v + now_dev_s, "\n", sep = "\t")
		if (k > 3 && (for_dev - now_dev_v - now_dev_s) / for_dev < thrd) { break }
		for_dev <- now_dev_v + now_dev_s
	}
	null_dev_v <- getDev(getNullCount(data$v1), data$v1$count)
	null_dev_s <- getDev(getNullCount(data$s1), data$s1$count)
	cat("R-squared =", 1 - (now_dev_v + now_dev_s) / (null_dev_v + null_dev_s),"\n")

	predCountV <- rho0 * getPredCountWeight(data$v1, model.glm0v, rho0) + rho1 * getPredCountWeight(data$v1, model.glm1v, rho1)
	predCountS <- rho0 * getPredCountWeight(data$s1, model.glm0s, rho0) + rho1 * getPredCountWeight(data$s1, model.glm1s, rho1)
	return(list(ProbV1=rho1v, ProbS1=rho1s,
	            PredV1=predCountV, PredS1=predCountS,
				DevV1=now_dev_v, DevS1=now_dev_s,
				Dev0V1=null_dev_v, Dev0S1=null_dev_s))
}


cv.MixturePoissonLinearCombine <- function(data, fold=5, same_dir=TRUE, seed=2011, thrd=1e-4, max_iter=80) {
	cat("Cross validation: Mixture of Poisson Linear Combine\n")
	train_index <- CVTrainIndexGene(data$v1, fold, seed)
	dev_v1 <- rep(0, fold); dev_s1 <- rep(0, fold);
	null_dev_v1 <- rep(0, fold); null_dev_s1 <- rep(0, fold);
	for (j in 1 : fold) {
		train_data <- list(v1= data$v1[train_index[j, ],], s1= data$s1[train_index[j, ],])
		test_data <- list(v1= data$v1[!train_index[j, ],], s1= data$s1[!train_index[j, ],])
		cat("Fold ", j, "of", fold, nrow(train_data$v1), "->", nrow(test_data$v1), "\n")
		## Train
		model <- train.MixturePoissonLinearCombine(train_data, same_dir, thrd, max_iter)
		## Predict
		outcomes <- predict.MixturePoissonLinearCombine(model, test_data, same_dir, thrd, max_iter)
		dev_v1[j] = outcomes$DevV1
		dev_s1[j] = outcomes$DevS1
		null_dev_v1[j] = outcomes$Dev0V1
		null_dev_s1[j] = outcomes$Dev0S1
	}
	return(list(r2v = 1-sum(dev_v1)/sum(null_dev_v1),
	            r2s = 1-sum(dev_s1)/sum(null_dev_s1)))
}


##############################################################
## fitting seperately
train.MixturePoissonSeperate <- function(data, thrd=1e-4, max_iter=80, rep_seed=1) {
	v1_fitted <- train.MixturePoisson(data$v1, thrd, max_iter, rep_seed)
	s1_fitted <- train.MixturePoisson(data$s1, thrd, max_iter, rep_seed)
	return(list(ProbV1=v1_fitted$prob, ProbS1=s1_fitted$prob, 
				PredV1=v1_fitted$pred, PredS1=s1_fitted$pred, 
				DevV1=v1_fitted$dev, DevS1=s1_fitted$dev, 
				Dev0V1=v1_fitted$dev0, Dev0S1=s1_fitted$dev0))
}

train.MixturePoissonLinearSeperate <- function(data, thrd=1e-4, max_iter=80, rep_seed=1) {
	v1_model <- train.MixturePoissonLinear(data$v1, thrd, max_iter, rep_seed)
	s1_model <- train.MixturePoissonLinear(data$s1, thrd, max_iter, rep_seed)
	return(list(v1_model, s1_model))
}

predict.MixturePoissonLinearSeperate <- function(model, data, thrd=1e-4, max_iter=80) {
	v1_fitted <- predict.MixturePoissonLinear(model[[1]], data$v1, thrd, max_iter)
	s1_fitted <- predict.MixturePoissonLinear(model[[2]], data$s1, thrd, max_iter)
	return(list(ProbV1=v1_fitted$prob, ProbS1=s1_fitted$prob, 
	            PredV1=v1_fitted$pred, PredS1=s1_fitted$pred, 
				DevV1=v1_fitted$dev, DevS1=s1_fitted$dev,
				Dev0V1=v1_fitted$dev0, Dev0S1=s1_fitted$dev0))
}


## format data
expDataCombine <- function(data_table, choose, left_span, right_span, duplicate=TRUE) {
	encoded <- expData(data_table, choose, left_span, right_span)
	if (duplicate) { ## default option
		encoded_copy <- encoded ## a new copy
		encoded$count <- data_table$v1[choose]
		encoded_copy$count <- data_table$s1[choose]
		if(FALSE) { ## reverse the order to break the correlation
			idx <- sample(1:nrow(encoded_copy))
			encoded_copy <- encoded_copy[idx,]
			encoded_copy$index <- encoded$index
		}
		return(list(v1=encoded, s1=encoded_copy))
	} else {
		encoded$v1 <- data_table$v1[choose]
		encoded$s1 <- data_table$s1[choose]
		return(encoded)
	}
}

## generate a new random data set
generateData <- function(choose, index, count, tau=0.7, relation=TRUE, seed=2011) {
	newCount <- rep(0, length(count)) ## a new copy
	for (i in unique(index)) {
		geneIndex <- choose & index == i
		counts <- count[geneIndex]
		num0 <- round(sum(geneIndex) * tau)
		num1 <- sum(geneIndex) - num0
		mu0 <- mean(counts) * 0.8
		mu1 <- mean(counts) * 1.2
		set.seed(seed)
		randomIndex <- sample(c(1:sum(geneIndex)))
		if (relation) { ## positive
			mixCount <- c(rpois(num0, mu0), rpois(num1, mu1))
		}else{
			mixCount <- c(rpois(num0, mu1), rpois(num1, mu0))
		}
		newCount[geneIndex] <- mixCount[randomIndex]# + round(rnorm(length(mixCount), 10, 3))
	}
	newCount[newCount < 0] <- 0 ## trim
	return(newCount)
}

fit_pars <- function(action, combine, window, data_file) {
	cat(action, combine, window, data_file, "\n")
	data_table <- read.csv(data_file)
	log_file <- paste(data_file, ".log", sep="")
	
	## generate data set
	#data_table$v1 <- generateData(data_table$tag>=0, data_table$index, data_table$v1, 0.5, TRUE)
	#data_table$s1 <- generateData(data_table$tag>=0, data_table$index, data_table$s1, 0.5, TRUE) 

	########################################################################
	left_span <- window/2
	right_span <- window/2

	########################################################################
	## Train
	if (action == "Train")
	{
		encoded <- expDataCombine(data_table, data_table$tag >= 0, left_span, right_span)

		if (combine == "MixPoiSep")
			model <- train.MixturePoissonSeperate(encoded)
		if (combine == "MixPoiLin") 
			model <- train.MixturePoissonLinearSeperate(encoded)
		if (combine == "MPLComSam")
			model <- train.MixturePoissonLinearCombine(encoded, TRUE)
		if (combine == "MPLComOpp")
			model <- train.MixturePoissonLinearCombine(encoded, FALSE)

		## save models
		save(model, file=paste(data_file, ".model", sep=""))
	}

	########################################################################
	## Predict
	if (action == "Predict")
	{
		## setting default values for blind regions
		numRow <- nrow(data_table)
		probv <- rep(0.5, numRow)
		probs <- rep(0.5, numRow)
		predv <- rep(-1, numRow)
		preds <- rep(-1, numRow)
		load(file=paste(data_file, ".model", sep=""))
		pack_ids <- round(data_table$index, -2) # package size: default -2
		r2v_dev = 0; r2v_dev0 = 0;
		r2s_dev = 0; r2s_dev0 = 0;

		for(i in unique(pack_ids)) {
			cat("predict package", i, "...\n")
			choose <- pack_ids == i & data_table$tag >= 0
			encoded <- expDataCombine(data_table, choose, left_span, right_span)

			if (combine == "MixPoiSep")
				outcome <- train.MixturePoissonSeperate(encoded, rep_seed=1)
			if (combine == "MixPoiLin") 
				outcome <- predict.MixturePoissonLinearSeperate(model, encoded)
			if (combine == "MPLComSam")
				outcome <- predict.MixturePoissonLinearCombine(model, encoded, TRUE)
			if (combine == "MPLComOpp")
				outcome <- predict.MixturePoissonLinearCombine(model, encoded, FALSE)

			probv[choose] <- outcome$ProbV1
			probs[choose] <- outcome$ProbS1
			predv[choose] <- outcome$PredV1
			preds[choose] <- outcome$PredS1
			r2v_dev <- r2v_dev + outcome$DevV1
			r2s_dev <- r2s_dev + outcome$DevS1
			r2v_dev0 <- r2v_dev0 + outcome$Dev0V1
			r2s_dev0 <- r2s_dev0 + outcome$Dev0S1
		}
		## Save R-squared for V1 and S1
		r2v = 1 - r2v_dev/r2v_dev0
		r2s = 1 - r2s_dev/r2s_dev0

		## Save fitted values into original data table
		data_table$pdv <- predv
		data_table$pds <- preds
		data_table$pbv <- probv
		data_table$pbs <- probs

		write.csv(data_table[,!(names(data_table) %in% c("pdv","pds"))], file=data_file, row.names=FALSE, quote=FALSE)
		write(paste("fit_pars_all.R", action, combine, window, r2v, r2s, r2v_dev, r2s_dev, sep="\t"), file=log_file)
	}

	########################################################################
	## Cross Validation
	if (action == "CV") 
	{
		encoded <- expDataCombine(data_table, data_table$tag >= 0, left_span, right_span)
		r2v = 0; r2s = 0;

		if (combine == "MixPoiSep") {
			r2v = cv.MixturePoisson(encoded$v1)
			r2s = cv.MixturePoisson(encoded$s1)
		}
		if (combine == "MixPoiLin") {
			r2v = cv.MixturePoissonLinear(encoded$v1)
			r2s = cv.MixturePoissonLinear(encoded$s1)
		}
		if (combine == "MPLComSam") {
			tmp = cv.MixturePoissonLinearCombine(encoded, same_dir=TRUE)
			r2v = tmp$r2v
			r2s = tmp$r2s
		}
		if (combine == "MPLComOpp") {
			tmp = cv.MixturePoissonLinearCombine(encoded, same_dir=FALSE)
			r2v = tmp$r2v
			r2s = tmp$r2s
		}

		write(paste("fit_pars_all.R", action, combine, window, r2v, r2s, sep="\t"), file=log_file)
	}
	########################################################################
	## End of all actions
}

## Main Function
main <- function() {
	action <- "Train"
	combine <- "MixPoiLin"
	window <- 2
	data_file <- "../work/PARS_V1_S1_win2.csv"
	args <- commandArgs(TRUE)
	if(length(args) >= 1) {
		action <- args[1]	
	if(length(args) >= 2) {
		combine <- args[2]
	if(length(args) >= 3) {
		window <- as.numeric(args[3])
	if(length(args) >= 4) {
		data_file <- args[4]
	} }}}
	fit_pars(action, combine, window, data_file)
}

main()
