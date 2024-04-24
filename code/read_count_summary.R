# get the list of genes and associated taxa
# both.filt generated in code/setup.R
#### OUTPUT sum.table.counts
#### OUTPUT mapping_prop.pdf

# generate simple plots of the filtered and unfiltered data and taxonomy barplots

# assign a taxonomy to all rows
tax.vec <- tax.table[rownames(both.filt),2]
names(tax.vec) <- rownames(both.filt)

# lets get the total number of genes in each species
n.genes <- vector()
for(i in levels(factor(tax.vec))){
    n.genes[i] <- length(which(tax.vec == i))
}

# 75 genes is a pretty small number to observe per species
# returns 23 taxa
tax.order <- names(which(n.genes > 75))

count.by.tax <- matrix(data=NA, nrow=length(tax.order), ncol=ncol(both.filt))

for(i in 1:length(tax.order)){
	count.by.tax[i, ] <- colSums(both.filt[which(tax.vec == tax.order[i]),])
}

sum.table.counts <- data.frame(matrix(data=NA, nrow=ncol(both.filt), ncol=6), stringsAsFactors=F)

rownames(sum.table.counts) <- colnames(both.filt)
colnames(sum.table.counts) <- c("total", "filtered", "most.frequent", "prop.unfilt", "prop.filt", "prop.abund")

sum.table.counts$total <- colSums(both.data)
sum.table.counts$filtered <- colSums(both.filt)
sum.table.counts$most.frequent <- colSums(count.by.tax)

sum.table.counts$prop.unfilt <- round((sum.table.counts$total / sum.table.counts$total),3)
sum.table.counts$prop.filt <- round((sum.table.counts$filtered / sum.table.counts$total),3)

sum.table.counts$prop.abund <- round((sum.table.counts$most.frequent / sum.table.counts$total),3)

sum.table.counts <- rbind(sum.table.counts, rep(1, 5))
sum.table.counts <- rbind(sum.table.counts, rep(0, 5))

co.cols <- c(rep("blue", 15), "red", rep("blue", 28), "black", "black")

pdf("mapping_prop.pdf")
parcoord(sum.table.counts[,4:6], var.label=T, col=co.cols, lwd=2, pch=19, lty=c(rep(2, 46), 1,1))
abline(h=0.709, lty=2, col="grey")
text(2.7,0.167, labels="0013B  16.7%", col="red")
text(2.5,0.39, labels="006A  38.9%", col="blue")

dev.off()
