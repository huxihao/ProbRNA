## This script is a part of the supporting source code for paper:             ##
## Hu et. al. Computational identification of protein binding sites on        ##
##            RNA using high-throughput RNA structure-probing data            ##
## It is free for non-commercial academic use.                                ##

# This File is based on the 'mseq' package. Full description is at following:
# Package: mseq
# Type: Package
# Title: Modeling non-uniformity in short-read rates in RNA-Seq data
# Version: 1.1
# Date: 2010-01-27
# Depends: R(>= 2.5),gbm
# Author: Jun Li
# Maintainer: Jun Li <junli07@stanford.edu>
# Description: This package implements all the methods in the paper
#         "Modeling non-uniformity in short-read rates in RNA-Seq data".
#         Especially, it implements both the iterative glm procedure for
#         the Poisson linear model and the training procedure of the MART
#         model. The cross-validation for both of the methods is also
#         implemented. Version 1.1 (current version) also implements the
#         Poisson linear model with dinucleotide composition.
# License: GPL (>= 2)
# LazyData: yes
# Packaged: Wed Jan 27 23:25:43 2010; junli07
# Repository: CRAN
# Date/Publication: 2010-01-28 08:04:22

`CVTrainIndexGene` <-
function(data, fold, seed)
{
	uniq_cat <- unique(data$index)
	train_index <- matrix(TRUE, fold, dim(data)[1])
	
	set.seed(seed)
	cv_cat <- sample(uniq_cat)
	for (j in 1 : length(uniq_cat))
	{
		row_num <- ceiling(j / floor(length(uniq_cat) / fold))
		if (row_num > fold)
		{
			row_num <- fold
		}
		train_index[row_num, data$index == cv_cat[j]] <- FALSE
	}
	
	return(train_index)
}

`expData` <-
function(oriData, llen, rlen)
{
	cat("\nExpanding the surrounding sequences to get the data frame for formulas to use...\n")
	
	seq <- as.character(oriData$seq)

	run_i <- 0
	num_i <- sum(oriData$tag == 0)
	index <- numeric(num_i)
	sseq <- character(num_i * (llen + rlen))
	count <- numeric(num_i)
	for (i in 1 : length(seq))
	{
		if (oriData$tag[i] == 0)
		{
			run_i <- run_i + 1
			index[run_i] <- oriData$index[i]
			sseq[((run_i - 1) * (llen + rlen) + 1) : (run_i * (llen + rlen))] <- seq[(i - llen) : (i + rlen - 1)]
			count[run_i] <- oriData$count[i]
		}
	}

	sseq <- factor(sseq, levels = c('T', 'A', 'C', 'G'))
	sseq <- matrix(sseq, ncol = llen + rlen, byrow = TRUE)
	data <- data.frame(index = index, count = count, sseq)
	cname <- character(2 + llen + rlen)
	cname[1] <- "index"
	cname[2] <- "count"
	for (i in 3 : (2 + llen + rlen))
	{
		j <- i - 3 - llen
		if (j < 0)
		{
			cname[i] <- paste("pM", -j, sep = '')
		}
		else
		{
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

`expData2nt` <-
function(oriData, llen, rlen)
{
	cat("normal expanding...\n")
	data <- expData(oriData, llen, rlen)
	
	c <- llen + rlen
	### the numeric data
	cat("Getting the numberic data...\n")
	data.num <- matrix(0, dim(data)[1], c)
	for (i in 1 : c)
	{
		data.num[, i] <- as.numeric(data[, 2 + i])
	}
	
	### the single nucleotide numeric data
	cat("Getting single nucleotide numeric data...\n")
	data.1d <- matrix(0, dim(data)[1], 4 * c)
	for (i in 1 : c)
	{
		for (j in 1 : 4)
		{
			data.1d[data.num[, i] == j, 4 * (i - 1) + j] <- 1
		}
	}
	
	### the bi-nucleotide numeric data
	cat("Getting bi-nucleotide numeric data...\n")
	data.2d <- matrix(0, dim(data)[1], 16 * (c - 1))
	for (i in 1 : (c - 1))
	{
		for (j in 1 : 4)
		{
			for (k in 1 : 4)
			{
				data.2d[(data.1d[, (i - 1) * 4 + j] == 1) & (data.1d[, i * 4 + k] == 1), (i - 1) * 16 + 4 * (j - 1) + k] <- 1
			}
		}
	}
	
	### delete the unnecessary columns
	cat("Deleting unnecessary columns...\n")
	data.1d <- data.1d[, c(1 : (4 * c)) %% 4 != 1]
	data.2d <- data.2d[, (c(1 : (16 * (c - 1))) %% 16 != 1) & (c(1 : (16 * (c - 1))) %% 16 != 8) & (c(1 : (16 * (c - 1))) %% 16 != 0) & (c(1 : (16 * (c - 1))) %% 16 < 12)]
	
	cat("Number of columns in data.1d = ", dim(data.1d)[2], "\n", sep = '')
	cat("Number of columns in data.2d = ", dim(data.2d)[2], "\n", sep = '')
	
	### combine them
	cat("Combining the single and bi nucleotide numeric data...\n")
	data.f <- data.frame(index = data$index, count = data$count, data.1d, data.2d)
	
	### set the column names
	cname <- character(2 + c * 3 + (c - 1) * 9)
	cname[1] <- "index"
	cname[2] <- "count"
	
	char1 <- c('A', 'C', 'G')
	for (k in 1 : 3)
	{
		for (i in 1 : c)
		{
			j <- i - llen - 1
			if (j < 0)
			{
				cname[2 + (i - 1) * 3 + k] <- paste("pM", -j, '.', char1[k], sep = '')
			}
			else
			{
				cname[2 + (i - 1) * 3 + k] <- paste("p", j, '.', char1[k], sep = '')
			}
		}
	}
	
	char2 <- c('TA', 'TC', 'TG', 'AT', 'AA', 'AC', 'CT', 'CA', 'CC')
	for (k in 1 : 9)
	{
		for (i in 1 : (c - 1))
		{
			j <- i - llen - 1
			if (j < -1)
			{
				cname[2 + 3 * c + (i - 1) * 9 + k] <- paste("pM", -j, '.', "pM", -(j + 1), '.', char2[k], sep = '')
			}
			if (j > -1)
			{
				cname[2 + 3 * c + (i - 1) * 9 + k] <- paste("p", j, '.', "p", j + 1, '.', char2[k], sep = '')
			}
			if (j == -1)
			{
				cname[2 + 3 * c + (i - 1) * 9 + k] <- paste('pM1.p0.', char2[k], sep = '')
			}
		}
	}
	
	colnames(data.f) <- cname

	cat("Number of columns in data.f = ", dim(data.f)[2], "\n", sep = '')
	
	return(data.f)
}

`getDev` <-
function(pred_count, real_count)
{
	dev <- sum(real_count[real_count != 0] * log(real_count[real_count != 0] / pred_count[real_count != 0]))
	dev <- dev - sum(real_count) + sum(pred_count)
	dev <- 2 * dev
	
	return(dev)
}

`getNullCount` <-
function(data)
{
	uniq_cat <- unique(data$index)
	null_count <- rep(0, length(data$count))
	for (k in uniq_cat)
	{
		null_count[data$index == k] <- sum(data$count[data$index == k]) / sum(data$index == k)
	}
	
	return(null_count)
}

`getOriLogPref` <-
function(data, small_count = 0.5)
{
	cat("Calculating the original log preferences...\n")
	
	modi_count <- data$count
	modi_count[data$count == 0] <- small_count
	
	uniq_cat <- unique(data$index)
	log_pref <- rep(0, length(data$index))
	
	for (k in uniq_cat)
	{
		mean_k <- sum(data$count[data$index == k]) / sum(data$index == k)
		log_pref[data$index == k] <- log(modi_count[data$index == k] / mean_k)
	}
	
	return(log_pref)
}

`getPredCount` <-
function(data, pred_pref)
{
	uniq_cat <- unique(data$index)
	pred_count <- rep(0, length(data$count))
	for (k in uniq_cat)
	{
		pred_count[data$index == k] <- sum(data$count[data$index == k]) / sum(pred_pref[data$index == k]) * pred_pref[data$index == k]
	}
	return(pred_count)
}

`getWeight` <-
function(data)
{
	cat("Calculating the weight...\n")
	uniq_cat <- unique(data$index)
	weight <- rep(0, length(data$index))
	for (k in uniq_cat)
	{
		weight[data$index == k] <- sum(data$count[data$index == k]) / sum(data$index == k)
	}
	return(weight)
}

`glmPred` <-
function(train.glm, newdata)
{
	#### get the coefficients, 4 for a position ####
	oriCoef <- train.glm$coefficients[-1]
	p <- length(oriCoef) / 3
	coef <- rep(0, p * 4)
	for (i in 1 : p)
	{
		coef[(i - 1) * 4 + 2 : 4] <- oriCoef[(i - 1) * 3 + 1 : 3]
	}

	#### get the predicted log preferences ####
	log_pref <- rep(0, dim(newdata)[1])
	for (i in 1 : p)
	{
		log_pref <- log_pref + coef[(i - 1) * 4 + as.numeric(newdata[, i + 2])]
	}

	return(log_pref)
}

`glmPred2nt` <-
function(train.glm, newdata)
{
	log_pref <- as.numeric(as.matrix(newdata[, -c(1, 2)]) %*% as.numeric(train.glm$coefficients[-1]))
	
	return(log_pref)
}

`iterGlm` <-
function(data, thrd = 0.01, max_iter = 10)
{
	cat("Begin iterative glm...\n")
	
	fordev <- 0
	log_expr <- setOriOffset(data)
	
	for (k in 1 : max_iter)
	{
		cat("\nIteration", k, "...\n")
		data.glm <- glm(count ~ . - index, data = data, family = poisson(link = "log"), offset = log_expr)
		log_expr <- updateOffset(data, data.glm$fitted.values, log_expr)
			
		cat("coefficients = ", data.glm$coefficients, "\n")
		nowdev <- getDev(data.glm$fitted.values, data$count)
		cat("deviance = ", nowdev, "\n")
		if (k > 1 && abs((nowdev - fordev) / fordev) < thrd)
		{
			break
		}
		fordev <- nowdev
	}

	dev <- getDev(data.glm$fitted.values, data$count)
	null_dev <- getDev(getNullCount(data), data$count)
	cat("deviance =", dev, "\n")
	cat("null deviance =", null_dev, "\n")
	cat("R_squared =", 1 - sum(dev) / sum(null_dev), "\n")
	
	return(data.glm)
}

`iterGlm2nt` <-
function(data, thrd = 0.01, max_iter = 10)
{
	return(iterGlm(data, thrd, max_iter))
}
`iterGlmCV` <-
function(data, fold = 5, seed = 281142, thrd = 0.01, max_iter = 10)
{
	cat("Cross validation using iterative glm...\n")
	train_index <- CVTrainIndexGene(data, fold, seed)
	
	dev <- rep(0, fold)
	null_dev <- rep(0, fold)
	
	for (j in 1 : fold)
	{
		cat("\nFold", j, "...\n")

		log_expr <- setOriOffset(data[train_index[j, ], ])
		fordev <- 0
		for (k in 1 : max_iter)
		{
			cat("\nIteration", k, "...\n")
			train.glm <- glm(count ~ . - index, data = data[train_index[j, ], ], family = poisson(link = "log"), offset = log_expr)
			log_expr <- updateOffset(data[train_index[j, ], ], train.glm$fitted.values, log_expr)
			
			cat("coefficients = ", train.glm$coefficients, "\n")
			nowdev <- getDev(getPredCount(data[train_index[j, ], ], exp(glmPred(train.glm, data[train_index[j, ], ]))), data[train_index[j, ], ]$count)
			cat("deviance = ", nowdev, "\n")
			if (k > 1 && abs((nowdev - fordev) / fordev) < thrd)
			{
				break
			}
			fordev <- nowdev
		}

		dev[j] <- getDev(getPredCount(data[!train_index[j, ], ], exp(glmPred(train.glm, data[!train_index[j, ], ]))), data[!train_index[j, ], ]$count)
		null_dev[j] <- getDev(getNullCount(data[!train_index[j, ], ]), data[!train_index[j, ], ]$count)
		cat("\nFold ", j, ": null deviance = ", null_dev[j], "; deviance = ", dev[j], "; R squared = ", 1 - dev[j] / null_dev[j], "\n", sep = '')
	}
	
	cat("\n")
	cat("deviance =", dev, "\n")
	cat("null deviance =", null_dev, "\n")
	
	R_squared <- 1 - sum(dev) / sum(null_dev)
	cat("R_squared =", R_squared, "\n")
	
	return(R_squared)
}

`iterGlmCV2nt` <-
function(data, fold = 5, seed = 281142, thrd = 0.01, max_iter = 10)
{
	cat("Cross validation using iterative glm...\n")
	train_index <- CVTrainIndexGene(data, fold, seed)
	
	dev <- rep(0, fold)
	null_dev <- rep(0, fold)
	
	for (j in 1 : fold)
	{
		cat("\nFold", j, "...\n")

		log_expr <- setOriOffset(data[train_index[j, ], ])
		fordev <- 0
		for (k in 1 : max_iter)
		{
			cat("\nIteration", k, "...\n")
			train.glm <- glm(count ~ . - index, data = data[train_index[j, ], ], family = poisson(link = "log"), offset = log_expr)
			log_expr <- updateOffset(data[train_index[j, ], ], train.glm$fitted.values, log_expr)
			
			cat("coefficients = ", train.glm$coefficients, "\n")
			nowdev <- getDev(getPredCount(data[train_index[j, ], ], exp(glmPred2nt(train.glm, data[train_index[j, ], ]))), data[train_index[j, ], ]$count)
			cat("deviance = ", nowdev, "\n")
			if (k > 1 && abs((nowdev - fordev) / fordev) < thrd)
			{
				break
			}
			fordev <- nowdev
		}

		dev[j] <- getDev(getPredCount(data[!train_index[j, ], ], exp(glmPred2nt(train.glm, data[!train_index[j, ], ]))), data[!train_index[j, ], ]$count)
		null_dev[j] <- getDev(getNullCount(data[!train_index[j, ], ]), data[!train_index[j, ], ]$count)
		cat("\nFold ", j, ": null deviance = ", null_dev[j], "; deviance = ", dev[j], "; R squared = ", 1 - dev[j] / null_dev[j], "\n", sep = '')
	}
	
	cat("\n")
	cat("deviance =", dev, "\n")
	cat("null deviance =", null_dev, "\n")
	
	R_squared <- 1 - sum(dev) / sum(null_dev)
	cat("R_squared =", R_squared, "\n")
	
	return(R_squared)
}

`martCV` <-
function(data, fold = 5, seed = 281142, shrinkage = 0.06, interaction.depth = 10, n.trees = 2000, small_count = 0.5)
{
	cat("Doing cross-validation for MART...\n")
	train_index <- CVTrainIndexGene(data, fold, seed)
	weight <- getWeight(data)
	
	dev <- rep(0, fold)
	null_dev <- rep(0, fold)
	
	for (j in 1 : fold)
	{
		cat("\nFold", j, "...\n")

		log_pref <- getOriLogPref(data[train_index[j, ], ], small_count)
		
		train.gbm <- gbm(log_pref ~ . - count - index, distribution = "gaussian", data = data[train_index[j, ], ], weights = weight[train_index[j, ]], n.trees = n.trees, shrinkage = shrinkage, interaction.depth = interaction.depth)

		dev[j] <- getDev(getPredCount(data[!train_index[j, ], ], exp(martPred(train.gbm, data[!train_index[j, ], ], n.trees))), data[!train_index[j, ], ]$count)
		null_dev[j] <- getDev(getNullCount(data[!train_index[j, ], ]), data[!train_index[j, ], ]$count)
		
		cat("R_squared for fold", j, "=", 1 - dev[j] / null_dev[j], "\n")
	}
	
	cat("deviance =", dev, "\n")
	cat("null deviance =", null_dev, "\n")
	
	R_squared <- 1 - sum(dev) / sum(null_dev)
	cat("R_squared =", R_squared, "\n")
	
	return(R_squared)
}

`martPred` <-
function(train.gbm, newdata, n.trees = 2000)
{
	# original prediction
	pred <- predict.gbm(train.gbm, newdata, n.trees)
	
	# construct all T data and get its prediction
	T_data <- newdata[1, ]
	T_data[1, 3 : dim(newdata)[2]] <- 'T'
	for (j in 3 : dim(newdata)[2])
	{
		T_data[, j] <- factor(T_data[, j], levels = c('T', 'A', 'C', 'G'))
	}
	pred_T <- predict.gbm(train.gbm, T_data, n.trees)
	
	# modify the prediction
	pred <- pred - pred_T
	
	return(pred)
}

`martTrain` <-
function(data, shrinkage = 0.06, interaction.depth = 10, n.trees = 2000, small_count = 0.5)
{
	cat("Training the data using Gbm...\n")
	weight <- getWeight(data)
	
	log_pref <- getOriLogPref(data, small_count)
	
	train.gbm <- gbm(log_pref ~ . - count - index, distribution = "gaussian", data = data, weights = weight, n.trees = n.trees, shrinkage = shrinkage, interaction.depth = interaction.depth)
	
	null_dev <- getDev(getNullCount(data), data$count)
	dev <- getDev(getPredCount(data, exp(martPred(train.gbm, data, n.trees))), data$count)
	cat("deviance =", dev, "\n")
	cat("null deviance =", null_dev, "\n")
	cat("R_squared =", 1 - sum(dev) / sum(null_dev), "\n")
	
	return(train.gbm)
}

plotCoef <- function(data.glm, llen, rlen)
{
	# get the coefficients
	coef <- data.glm$coefficients[-1]
	p <- length(coef) / 3
	coef_T <- rep(0, p)
	coef_A <- coef[c(1 : p) * 3 - 2]
	coef_C <- coef[c(1 : p) * 3 - 1]
	coef_G <- coef[c(1 : p) * 3]

	# get the positions
	pos <- c((-llen) : (rlen - 1))

	plot(NA, NA, xlim = c(min(pos), max(pos)), ylim = 1.0 * c(min(coef), max(coef)), xlab = "position", ylab = "coefficients")
	title("coefficients, red-T, green-A, blue-C, black-G")
	
	lines(pos, coef_T, type = 'o', col = 'red', cex = 0.5)
	lines(pos, coef_A, type = 'o', col = 'green', cex = 0.5)
	lines(pos, coef_C, type = 'o', col = 'blue', cex = 0.5)
	lines(pos, coef_G, type = 'o', col = 'black', cex = 0.5)
}
plotCounts <- function(counts)
{
	plot(c(1, length(counts)), c(0, 0), type = 'l', ylim = c(0, max(counts)), col = 'white', xlab = "position", ylab = "counts")
	
	for (k in 1 : length(counts))
	{
		if (counts[k] != 0)
		{
			lines(c(k, k), c(0, counts[k]), col = 'black')
		}
	}
}

`setOriOffset` <-
function(data)
{
	uniq_cat <- unique(data$index)
	log_expr <- rep(0, length(data$index))
	
	for (k in uniq_cat)
	{
		log_expr[data$index == k] <- log(sum(data$count[data$index == k]) / sum(data$index == k))
	}

	return(log_expr)
}

`updateOffset` <-
function(data, pred_count, log_expr)
{
	uniq_cat <- unique(data$index)
	
	log_pref <- log(pred_count) - log_expr
	
	for (k in uniq_cat)
	{
		log_expr[data$index == k] <- log(sum(data$count[data$index == k]) / sum(exp(log_pref[data$index == k])))
	}

	return(log_expr)
}

