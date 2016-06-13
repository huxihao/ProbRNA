
library("kernlab")
library("ROCR")


learn_domain <- function(action, data_file, cost) {
data_all <- read.csv(data_file)
data_all$gene <- factor(data_all$gene)
cat(names(data_all), "\n")
marks <- c("gene", "pos")
pred <- rep(0, nrow(data_all))


## ACTION ONE
if(action == "CV") {
cat("Cross validatoin on", data_file, cost, "\n")
avg_auc <- 0.0
for(name in unique(data_all$gene)) {
	data_train <- data_all[data_all$gene != name, !names(data_all) %in% marks]
	data_test <- data_all[data_all$gene == name, !names(data_all) %in% marks]
	model <- ksvm(label ~ ., data=data_train, scaled=FALSE, type="eps-svr", kernel="rbfdot", C=cost, epsilon=0.1)
	predicted  <- predict(model, data_test[,-1])
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
write(paste(cost, avg_auc, sep="\t"), file=paste(data_file, ".log", sep=""), append=TRUE)
mse <- mean((data_all$label-pred)**2)**0.5
## Plot
plot(data_all$label,col=1,type='l',main=paste("Avg AUC =", avg_auc))
points(pred,col=2,pch=20)
write.csv(data.frame(gene=data_all$gene, pos=data_all$pos, label=data_all$label, pred=pred), file=paste("cv_comp.csv"))
return(avg_auc)

## ACTION TWO
}else if(action == "Train") {
cat("Training on", data_file, "\n")
	data_train <- data_all[, !names(data_all) %in% marks]
	model <- ksvm(label ~ ., data=data_train, type="eps-svr", kernel="rbfdot", C=cost, epsilon=0.1)
	save(model, file=paste(data_file, ".model", sep=""))
return(0)

## ACTION THREE
} else if(action == "Test") {
cat("Testing on", data_file, "\n")
	data_test <- data_all[, !names(data_all) %in% marks]
	load(file=paste(data_file, ".model", sep=""))
	predicted  <- predict(model, data_test[,-1])
	plot(data_test[,1],col=1,pch=20)
	points(predicted,col=2)
	if(sum(data_test[,1]>0)) { ## has know positive samples
		pair <- prediction(predicted, (data_test[,1]>0))
		pref <- performance(pair, measure = "tpr", x.measure = "fpr")
		auc <- performance(pair, "auc")@y.values[[1]]
		plot(pref, col=1, lty=2, main=paste("AUC =", round(auc,2)))
	}
return(0)
}
## End of actions
}

###############################################################################
## Main Function
action <- "CV"
data_file <- "../work/domain_encode.csv"
args <- commandArgs(TRUE)
para <- 0/0

cat("Command:", args, "\n")
if(length(args) >= 2) {
	action <- args[1]	
	data_file <- args[2]
	if(length(args) >= 3) { ## parameter
		para <- as.numeric(args[3])
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
cat("Parameter", par, "get", learn_domain(action, data_file, par), "\n")
}

dev.off()
