
library("kernlab")
library("ROCR")

###############################################################################
## Kernel Part

dyn.load("../src/kernel_la.so")
dyn.load("../src/kernel_bpla.so")
dyn.load("../src/kernel_dpla.so")

la_kernel <- function(beta= 0.5, open= -10, extend= -5) 
{
	## subsitution matrix, which can also be a parameter
	## here is just an example with only two characters
	sub_mat = matrix(c(
		## RIBOSUM85-60 matrix
		## Reference:
		##    Robert J Klein and Sean R Eddy. 
		##    RSEARCH: Finding homologs of single structured RNA sequences
		## http://www.biomedcentral.com/1471-2105/4/44/figure/F3?highres=y
		##  A,     C,     G,     U
		 2.22, -1.86, -1.46, -1.39, ## A
		-1.86,  1.16, -2.48, -1.05, ## C
		-1.46, -2.48,  1.03, -1.74, ## G
		-1.39, -1.05, -1.74,  1.65  ## U
		), 4)

	## make it to be a semi-definent
	sub_mat <- sub_mat - diag(min(eigen(sub_mat)$value)-1e-9, nrow(sub_mat))

	## input parameters and initialize
	.C("la_set_para", as.integer(nrow(sub_mat)), as.numeric(sub_mat),
		as.double(beta), as.double(open), as.double(extend))

	## calculate the dot product for given inputs
	kval <- function(x, y = NULL) {
	## format: seq, Pr_l, Pr_r, Pr_v, Pr_s
	len <- length(x)/6
	.C("la_get_norm", 
	   as.integer(len),
	   as.integer(x[1:len]), 
	   as.integer(y[1:len]), 
	   score=as.double(0))$score
	}
	return(new("kernel", .Data=kval, 
	       kpar=list(beta=beta, open=open, extend=extend)))
}

bpla_kernel <- function(alpha= 7, beta= 0.5, open= -10, extend= -5) 
{
	## subsitution matrix, which can also be a parameter
	## here is just an example with only two characters
	sub_mat = matrix(c(
		## RIBOSUM85-60 matrix
		## Reference:
		##    Robert J Klein and Sean R Eddy. 
		##    RSEARCH: Finding homologs of single structured RNA sequences
		## http://www.biomedcentral.com/1471-2105/4/44/figure/F3?highres=y
		##  A,     C,     G,     U
		 2.22, -1.86, -1.46, -1.39, ## A
		-1.86,  1.16, -2.48, -1.05, ## C
		-1.46, -2.48,  1.03, -1.74, ## G
		-1.39, -1.05, -1.74,  1.65  ## U
		), 4)

	## make it to be a semi-definent
	sub_mat <- sub_mat - diag(min(eigen(sub_mat)$value)-1e-9, nrow(sub_mat))

	## input parameters and initialize
	.C("bpla_set_para", as.integer(nrow(sub_mat)), as.numeric(sub_mat),
	   as.double(alpha), as.double(beta), as.double(open), as.double(extend))

	## calculate the dot product for given inputs
	kval <- function(x, y = NULL) {
	## format: Seq, Pr_Left, Pr_Right, Pr_Unpair
	len <- length(x)/6
	.C("bpla_get_norm", 
	   as.integer(len),
	   as.integer(x[1:len]),
	   as.integer(y[1:len]), 
	   as.numeric(x[(len+1):(4*len)]),
	   as.numeric(y[(len+1):(4*len)]),
	   score=as.double(0))$score
	}
	return(new("kernel", .Data=kval, 
	       kpar=list(alpha=alpha, beta=beta, open=open, extend=extend)))
}

dpla_kernel <- function(alpha= 7, beta= 0.5, open= -10, extend= -5) 
{
	## subsitution matrix, which can also be a parameter
	## here is just an example with only two characters
	sub_mat = matrix(c(
		## RIBOSUM85-60 matrix
		## Reference:
		##    Robert J Klein and Sean R Eddy. 
		##    RSEARCH: Finding homologs of single structured RNA sequences
		## http://www.biomedcentral.com/1471-2105/4/44/figure/F3?highres=y
		##  A,     C,     G,     U
		 2.22, -1.86, -1.46, -1.39, ## A
		-1.86,  1.16, -2.48, -1.05, ## C
		-1.46, -2.48,  1.03, -1.74, ## G
		-1.39, -1.05, -1.74,  1.65  ## U
		), 4)

	## make it to be a semi-definent
	sub_mat <- sub_mat - diag(min(eigen(sub_mat)$value)-1e-9, nrow(sub_mat))

	## input parameters and initialize
	.C("dpla_set_para", as.integer(nrow(sub_mat)), as.numeric(sub_mat),
	   as.double(alpha), as.double(beta), as.double(open), as.double(extend))

	## calculate the dot product for given inputs
	kval <- function(x, y = NULL) {
	## format: Seq, Pr_Left, Pr_Right, Pr_Unpair
	len <- length(x)/6
	.C("dpla_get_norm", 
	   as.integer(len),
	   as.integer(x[1:len]),
	   as.integer(y[1:len]), 
	   as.numeric(x[(len+1):(4*len)]),
	   as.numeric(y[(len+1):(4*len)]),
	   as.numeric(x[(4*len+1):(6*len)]),
	   as.numeric(y[(4*len+1):(6*len)]),
	   score=as.double(0))$score
	}
	return(new("kernel", .Data=kval, 
	       kpar=list(alpha=alpha, beta=beta, open=open, extend=extend)))
}


#kernel_for_all <<- matrix(0)
#precomputed_kernel <- function()
#{
#	cat("Size of precomputed kernel matrix is", dim(kernel_for_all), "\n")
#	kval <- function(x, y = NULL) {
#		## we assume that x and y are just indices to be evaluated
#		kernel_for_all[x[1], y[1]]
#	}
#	return(new("kernel", .Data=kval, kpar=list()))
#}

###############################################################################
## Action Part

learn_domain <- function(action, data_file, para) {
data_all <- read.csv(file=data_file)
data_all$gene <- factor(data_all$gene)
cat(names(data_all), "\n")
marks <- c("gene", "pos")
pred <- rep(0, nrow(data_all))

model_name <- "kernlab.model"
set.seed(2011)

###############################################################################
## Use Precomputed Kernel Matrix to speed up
#if(nrow(kernel_for_all)==0) { ## first running
#kernel_function <- la_kernel(0.5, -10, -5)
#kernel_for_all <<- as.matrix(kernelMatrix(kernel_function, 
#	as.matrix(data_all[, !names(data_all) %in% c("label", marks)])))

## only save label and indexes in the kernel matrix
#data_all <- data_all[, c("label", marks)]
#data_all$index_in_kernel_matrix <- 1:nrow(data_all)
#data_all$add_something_just_to_fix_bug <- 0
#}
###############################################################################


#n ACTION ONE
if(action == "CV") {
cat("Cross validatoin on", data_file, para, "\n")
avg_auc <- 0.0
for(name in unique(data_all$gene)) {
	data_train <- data_all[data_all$gene != name, !names(data_all) %in% marks]
	data_test <- data_all[data_all$gene == name, !names(data_all) %in% marks]

	## deal with unbalanced data
	pos_index <- seq(1:nrow(data_train))[data_train$label > 0]
	neg_index <- seq(1:nrow(data_train))[data_train$label <= 0]
	predicted <- rep(0, nrow(data_test))
	model <- as.list(rep(0,bagging))
	for(bag in 1:bagging) {
		random_neg_index <- sample(neg_index, length(pos_index))
		model[[bag]] <- ksvm(
			label ~ ., data=data_train[c(pos_index, random_neg_index),],
			scaled=FALSE, type="eps-svr", 
			#kernel=dpla_kernel, 
			#kernel=bpla_kernel, 
			#kpar=list(alpha=7, beta= 0.07, open= -14, extend= -0.07),
			kernel=la_kernel,
			kpar=list(beta= 0.07, open= -14, extend= -0.07),
			C=para, epsilon=0.1)
		if(bag == 1){
			predicted <- predict(model[[bag]], data_test[,-1])
		}else{
			predicted <- predicted + predict(model[[bag]], data_test[,-1])
		}
	}
	predicted <- predicted / bagging

	pred[data_all$gene == name] <- predicted
	pair <- prediction(predicted, (data_test[,1]>0))
	pref <- performance(pair, measure = "tpr", x.measure = "fpr")
	auc <- performance(pair, "auc")@y.values[[1]]
	avg_auc <- avg_auc + auc
	plot(pref, col=1, type='l', main=paste(name, round(auc,2)))
	cat("Train", nrow(data_train), "->", "Test", nrow(data_test), name, auc, "\n")
	write(paste(name, auc, sep="\t"), file=paste(data_file, ".log", sep=""), append=TRUE)
}
avg_auc <- avg_auc / length(unique(data_all$gene))
write(paste(para, avg_auc, sep="\t"), file=paste(data_file, ".log", sep=""), append=TRUE)
mse <- mean((data_all$label-pred)**2)**0.5
## Plot
plot(data_all$label,col=1,type='l',main=paste("Avg AUC =", avg_auc))
points(pred,col=2,pch=20)
write.csv(data.frame(gene=data_all$gene, pos=data_all$pos, label=data_all$label, pred=pred), file=paste("cv_comp.csv"))
return(avg_auc)

## ACTION TWO
}else if(action == "Train") {
cat("Training on", data_file, para, "\n")
	data_train <- data_all[, !names(data_all) %in% marks]
	pos_index <- seq(1:nrow(data_train))[data_train$label > 0]
	neg_index <- seq(1:nrow(data_train))[data_train$label <= 0]
	model <- as.list(rep(0,bagging))
	for(bag in 1:bagging) {
		random_neg_index <- sample(neg_index, length(pos_index))
		model[[bag]] <- ksvm(label ~ ., 
			data=data_train[c(pos_index, random_neg_index),],
			scaled=FALSE, type="eps-svr",
			kernel=dpla_kernel, 
			kpar=list(alpha=7, beta= 0.07, open= -14, extend= -0.07),
			C=para, epsilon=0.1)
	}
	save(model, file=model_name)
return(0)

## ACTION THREE
} else if(action == "Test") {
cat("Testing on", data_file, "\n")
	data_test <- data_all[, !names(data_all) %in% marks]
	load(file=model_name)

	predicted <- rep(0, nrow(data_test))
	for(name in unique(data_all$gene)){
	choose <- data_all$gene == name
	cat("Predict", name, "\n")
	for(bag in 1:bagging) {
		if(bag == 1){
			predicted[choose] <- predict(model[[bag]], data_test[choose,-1])
		}else{
			predicted[choose] <- predicted[choose] + predict(model[[bag]], data_test[choose,-1])
		}
	}
	}
	predicted <- predicted / bagging

	plot(data_test[,1],col=1,pch=20)
	points(predicted,col=2)
	if(sum(data_test[,1]) > 0) { ## has know positive samples so do comparison
		pair <- prediction(predicted, (data_test[,1]>0))
		pref <- performance(pair, measure = "tpr", x.measure = "fpr")
		auc <- performance(pair, "auc")@y.values[[1]]
		plot(pref, col=1, lty=2, main=paste("AUC =", round(auc,2)))
		return(auc)
	}else{ ## save data
		data_all$label <- predicted
		write.csv(data_all, row.names=FALSE, file=data_file)
	}
return(0)
}
## End of Actions
}

###############################################################################
## Main Function
action <- "CV"
data_file <- "../work/domain_encode.csv"
args <- commandArgs(TRUE)
bagging <- 1 ## global
para <- 0/0

cat("Command:", args, "\n")
if(length(args) >= 2) {
	action <- args[1]	
	data_file <- args[2]
	if(length(args) >= 3) { ## parameter
		bagging <- as.numeric(args[3])
		if(length(args) >= 4) { ## bagging time
			para <- as.numeric(args[4])
		}
	}
}
pdf("../work/learn_domain_kernlab.pdf")

if(action == "CV") { ## search best parameter
write("learn_domain_kernlab_la.R", file=paste(data_file, ".log", sep=""))
best_par <- 1
best_acc <- 0
last_pos <- 0 ## left: -1; begin: 0; right: 1

if(is.nan(para)){ ## search best parameter
para <- 1
repeat{
for(new_par in c(0.01,0.1,1,10,100)*para) {
	acc <- learn_domain(action, data_file, new_par)
	if(acc > best_acc){
		best_acc <- acc
		best_par <- new_par
	}
}
if(best_par < 0.05*para){ ## shift to left
	if(last_pos == -1) {break}
	para <- para * 1e-5
	last_pos <- 1 ## last best from right
}else{ 
if(best_par > 50*para){ ## shift to right
	if(last_pos == 1) {break}
	para <- para * 1e5
	last_pos <- -1 ## last best from left
}else{ break }}
cat("Current best paramter is", best_par, "with", best_acc, "\n")
}
}else{
best_par <- para
}
## recall the best parameter
best_acc <- learn_domain(action, data_file, best_par)
cat("Best paramter is", best_par, "with", best_acc, "\n")

}else{ ## other actions
cat("Parameter", para, "get", learn_domain(action, data_file, para), "\n")
}

dev.off()
