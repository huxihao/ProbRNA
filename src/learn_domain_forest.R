## This script is a part of the supporting source code for paper:             ##
## Hu et. al. Computational identification of protein binding sites on        ##
##            RNA using high-throughput RNA structure-probing data            ##
## It is free for non-commercial academic use.                                ##


library("randomForest")
library("ROCR")


learn_domain <- function(action, data_file, mtryCoef) 
{
## FUNCTION BEGIN
data_table <- read.csv(data_file)
data_table$gene <- factor(data_table$gene)
cat(names(data_table), "\n")
marks <- c("gene", "pos")
pred <- rep(0, nrow(data_table))


## ACTION ONE
if(action == "CV") {
cat("Cross validation on", data_file, "\n")
avg_auc <- 0.0
for(name in unique(data_table$gene)) {
	data_train <- data_table[data_table$gene != name, !names(data_table) %in% marks]
	data_test <- data_table[data_table$gene == name, !names(data_table) %in% marks]
	model <- randomForest(x=data_train[,-1], y=factor(data_train[,1]), sampsize=min(5000,nrow(data_train)))
	predicted  <- predict(model, data_test[,-1], type="prob")[,2]
	pred[data_table$gene == name] <- predicted
	pair <- prediction(predicted, (data_test[,1]>0))
	pref <- performance(pair, measure = "tpr", x.measure = "fpr")
	auc <- performance(pair, "auc")@y.values[[1]]
	avg_auc <- avg_auc + auc
	#plot(pref, col=1, type='l', main=paste(name, round(auc,2)))
	cat("Train", nrow(data_train), "->", "Test", nrow(data_test), name, auc, "\n")
	write(paste(name, auc, nrow(data_train), nrow(data_test), sep="\t"), file=paste(data_file, ".log", sep=""), append=TRUE)
}
avg_auc <- avg_auc / length(unique(data_table$gene))
pair <- prediction(pred, data_table$label)
all_auc <- performance(pair, "auc")@y.values[[1]]
## calculate bit gain
#pred = data_table[,1] ## for control
cutoff = sort(unique(pred), decreasing=FALSE)
info_bit = rep(0, length(cutoff))
for(i in 1:length(cutoff)) {
	pos_p <- sum(data_table[pred>=cutoff[i],1]>0)/as.double(sum(pred>=cutoff[i]))
	pos_q <- 1 - pos_p
	if(pos_p > 0 & pos_q >0) pos_bit <- - pos_p * log2(pos_p) - pos_q * log2(pos_q)
	else pos_bit = 0
	if(i == 1) {
		## for the random predictor
		info_bit[1] = pos_bit
	}else{
		neg_p <- sum(data_table[pred<cutoff[i],1]<0)/as.double(sum(pred<cutoff[i]))
		neg_q <- 1 - neg_p
		if(neg_p > 0 & neg_q >0) neg_bit <- - neg_p * log2(neg_p) - neg_q * log2(neg_q)
		else neg_bit = 0
		info_bit[i] = as.double(sum(pred>=cutoff[i]))*pos_bit + as.double(sum(pred<cutoff[i]))*neg_bit
		info_bit[i] = info_bit[i] / as.double(length(pred))
	}
}
write(paste("Total", all_auc, nrow(data_table), info_bit[1], min(info_bit), sep="\t"), file=paste(data_file, ".log", sep=""), append=TRUE)
write('Predicted Labels:', file=paste(data_file, ".log", sep=""), append=TRUE)
write.table(pred, file=paste(data_file, ".log", sep=""), sep='\t', append=TRUE, row.names=FALSE, col.names=FALSE)
mse <- mean((data_table$label - pred)**2)**0.5
## Plot
#plot(data_table$label,col=1,type='l',main=paste("Avg AUC =", avg_auc))
#points(pred,col=2,pch=20)
cat("AUC =", avg_auc, "MSE =", mse, "\n")

## ACTION TWO
}else if(action == "Train") {
	cat("Train on", data_file, "\n")
	data_train <- data_table[, !names(data_table) %in% marks]
	model <- randomForest(x=data_train[,-1], y=factor(data_train[,1]), sampsize=min(5000,nrow(data_train)))
	save(model, file=paste(data_file, ".model", sep=""))
	## Plot
	#varImpPlot(model)
	#plot(model$importance)
	#plot(data_train[,1],col=1,type='l')
	#points(model$predicted, col=2)
	#pair <- prediction(model$predicted, (data_train[,1]>0))
	#pref <- performance(pair, measure = "tpr", x.measure = "fpr")
	#auc <- performance(pair, "auc")@y.values[[1]]
	#plot(pref, col=1, type='l', main=paste("AUC =", round(auc,2)))

## ACTION THREE
} else if(action == "Predict") {
	cat("Predict on", data_file, "\n")
	data_test <- data_table[, !names(data_table) %in% marks]
	load(file=paste(data_file, ".model", sep=""))
	predicted  <- predict(model, data_test[,-1], type="prob")[,2]
	if(sum(data_test[,1]) > 0) { ## has know positive samples so do comparison
		pair <- prediction(predicted, (data_test[,1]>0))
		pref <- performance(pair, measure = "tpr", x.measure = "fpr")
		auc <- performance(pair, "auc")@y.values[[1]]
#		plot(pref, col=1, lty=2, main=paste("AUC =", round(auc,2)))
		return(auc)
	}else{ ## save data
		data_table$label <- predicted
		write.csv(data_table, row.names=FALSE, file=data_file)
	}
	return(0)
}
## FUNCTION END
}

## Main Function
action <- "CV"
data_file <- "../work/domain_encode.csv"
args <- commandArgs(TRUE)

if(length(args) >= 2) {
	action <- args[1]	
	data_file <- args[2]
}

begin <- 1
end <- 1
if(action == "CV") {
	write("Model\tlearn_domain_forest.R", file=paste(data_file, ".log", sep=""), append=FALSE)
	if(length(args) >= 4) {
		begin <- as.numeric(args[3])
		end <- as.numeric(args[4])
	}
}
cat("Commands:", action, data_file, begin, end, "\n")

#pdf("learn_domain_forest.pdf")
#for(mtryCoef in begin:end)
#learn_domain(action, data_file, mtryCoef)
learn_domain(action, data_file, 1)
#dev.off()

