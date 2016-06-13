## This script is a part of the supporting source code for paper:             ##
## Hu et. al. Computational identification of protein binding sites on        ##
##            RNA using high-throughput RNA structure-probing data            ##
## It is free for non-commercial academic use.                                ##


data <- read.table("validate_table.txt", sep=",", header=TRUE)

library("ROCR")

AUCs <- rep(0, length(unique(data$index)))
i <- 0
pdf(paste("validate_all_",length(AUCs),".pdf",sep=""))
cat("We have", length(AUCs), "genes\n")
for (index in unique(data$index)) {
	one_cell <- data[data$index == index,]
	pair <- prediction(one_cell$score, one_cell$label>0)
	pref <- performance(pair, measure = "tpr", x.measure = "fpr")
	auc <- performance(pair, "auc")@y.values[[1]]
	i <- i + 1
	AUCs[i] <- auc
}
write(paste("GeneNumber", length(AUCs), sep="\t"), file="log.txt", append=TRUE)
write(paste("Average", mean(AUCs), sep="\t"), file="log.txt", append=TRUE)
h <- hist(AUCs, 40)
write(paste("Breaks", "Counts", "Density", sep="\t"), file="log.txt", append=TRUE)
for (i in 1:length(h$breaks)) {
write(paste(h$breaks[i], h$counts[i], h$density[i], sep="\t"), file="log.txt", append=TRUE)
}
dev.off()
