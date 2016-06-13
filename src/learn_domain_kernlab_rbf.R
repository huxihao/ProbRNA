
library("kernlab")
library("ROCR")


learn_domain <- function(action, data_file, para, old=TRUE) {
write(paste("Parameter", para, sep="\t"), file=paste(data_file, ".log", sep=""), append=old)
data_all <- read.csv(file=data_file)
data_all$gene <- factor(data_all$gene)
cat(names(data_all), "\n")
marks <- c("gene", "pos")
pred <- rep(0, nrow(data_all))

model_name <- "kernlab.model"
set.seed(2011)

## ACTION ONE
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
		model[[bag]] <- ksvm(label ~ ., data=data_train[c(pos_index, random_neg_index),],
		         scaled=FALSE, type="eps-svr", kernel="rbfdot", C=para, epsilon=0.1)
		if(bag == 1){
			predicted <- predict(model[[bag]], data_test[,-1])
		}else{
			#predicted <- apply(cbind(predicted, predict(model, data_test[,-1])), 1, min)
			predicted <- predicted + predict(model[[bag]], data_test[,-1])
		}
	}
	predicted <- predicted / as.numeric(bagging)

	pred[data_all$gene == name] <- predicted
	pair <- prediction(predicted, (data_test[,1]>0))
	pref <- performance(pair, measure = "tpr", x.measure = "fpr")
	auc <- performance(pair, "auc")@y.values[[1]]
	avg_auc <- avg_auc + auc
	
	## Plots
	if(FALSE) {
	plot(pref, col=1, type='l', main=paste(name, round(auc,2)))
	gene_pos <- data_all$pos[data_all$gene==name]
	sort_label <- rep(0, max(gene_pos))
	sort_pred <- rep(0, max(gene_pos))
	sort_label[gene_pos] <- data_test[,1] #data_all$label[data_all$gene==name]
	sort_pred[gene_pos] <- predicted
	plot(sort_label,col=1,type='l',main=paste(name, length(pos_index), length(neg_index)))
	points(sort_pred,col=2,pch=20)
	}

	cat("Train", nrow(data_train), "->", "Test", nrow(data_test), name, auc, "\n")
	write(paste(name, auc, sep="\t"), file=paste(data_file, ".log", sep=""), append=TRUE)
}
avg_auc <- avg_auc / length(unique(data_all$gene))
write(paste("Average", avg_auc, sep="\t"), file=paste(data_file, ".log", sep=""), append=TRUE)
mse <- mean((data_all$label-pred)**2)**0.5
#write.csv(data.frame(gene=data_all$gene, pos=data_all$pos, label=data_all$label, pred=pred), file=paste("cv_comp.csv"))
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
		model[[bag]] <- ksvm(label ~ ., data=data_train[c(pos_index, random_neg_index),],
		         scaled=FALSE, type="eps-svr", kernel="rbfdot", C=para, epsilon=0.1)
	}
	save(model, file=model_name)
return(0)

## ACTION THREE
} else if(action == "Predict") {
cat("Predicting on", data_file, "\n")
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
	predicted <- predicted / as.numeric(bagging)

#	plot(data_test[,1],col=1,pch=20)
#	points(predicted,col=2)
	if(sum(data_test[,1]) > 0) { ## has know positive samples so do comparison
		pair <- prediction(predicted, (data_test[,1]>0))
		pref <- performance(pair, measure = "tpr", x.measure = "fpr")
		auc <- performance(pair, "auc")@y.values[[1]]
#		plot(pref, col=1, lty=2, main=paste("AUC =", round(auc,2)))
		return(auc)
	}else{ ## save data
		data_all <- data_all[,names(data_all) %in% marks]
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

if(length(args) >= 2) {
	action <- args[1]	
	data_file <- args[2]
	if(length(args) >= 3) { ## bagging time
		bagging <- as.numeric(args[3])
		if(length(args) >= 4) { ## input default parameter
			para <- as.numeric(args[4])
		}
	}
}
cat("Commands:", action, data_file, bagging, para, "\n")

#pdf("../work/learn_domain_kernlab.pdf")

if(action == "CV") { ## search best parameter
write("learn_domain_kernlab_rbf.R", file=paste(data_file, ".log", sep=""))
best_par <- 1
best_acc <- 0
last_pos <- 0 ## left: -1; begin: 0; right: 1

if(is.nan(para)){ ## search best parameter
para <- 1
#repeat{
for(new_par in c(0.01,0.1,1,10,100)*para) {
	acc <- learn_domain(action, data_file, new_par)
	if(acc > best_acc){
		best_acc <- acc
		best_par <- new_par
	}
}
#if(best_par < 0.05*para){ ## shift to left
#	if(last_pos == -1) {break}
#	para <- para * 1e-5
#	last_pos <- 1 ## last best from right
#}else{ 
#if(best_par > 50*para){ ## shift to right
#	if(last_pos == 1) {break}
#	para <- para * 1e5
#	last_pos <- -1 ## last best from left
#}else{ break }}
#cat("Current best paramter is", best_par, "with", best_acc, "\n")
#} ## end of repeat
#}else{
#best_par <- para
}

## recall the best parameter
best_acc <- learn_domain(action, data_file, best_par, FALSE)
cat("Best parameter is", best_par, "with", best_acc, "\n")

}else{ ## other actions
cat("Parameter", para, "get", learn_domain(action, data_file, para), "\n")
}

#dev.off()
