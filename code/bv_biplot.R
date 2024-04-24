source("code/setup.R")
# compare BV only from two datasets
# pull from the merged datasets
library(sva)

bv.data <- data.frame(both.data[,grep('v.', colnames(both.data))])
bv.filt <-  codaSeq.filter(bv.data, min.occurrence=0.30,
    min.prop=0.00005, samples.by.row=F)

btch.bv <- c(rep(1,14), rep(2,14))
bv.combat <- sva::ComBat_seq(as.matrix(bv.filt), batch=btch.bv, group=NULL, full_mod=T)

bv.clr <- apply(bv.combat + 0.5, 2, function(x) log(x) - mean(log(x)))

# assign a taxonomy to all rows
tax.vec.bv <- as.vector(tax.table[rownames(bv.combat),2])
names(tax.vec.bv) <- rownames(bv.combat)

# lets get the total number of genes in each species
n.genes.bv <- vector()
for(i in levels(factor(tax.vec.bv))){
n.genes.bv[i] <- length(which(tax.vec.bv == i))
}

tax.order.bv <- names(which(n.genes.bv > 50))
tax.colours.bv <- c(tax.colors[tax.order.bv,], rgb(0,0,0,0.05))

loadings.group.bv <- list()

for(i in 1:length(tax.order.bv)){
    loadings.group.bv[[tax.order.bv[i]]] <- which(tax.vec.bv == tax.order.bv[i])
}

loadings.group.bv[["other"]] <- setdiff(c(1:nrow(bv.clr)),unlist(loadings.group.bv))
tax.symbol <- c(rep(20, 8), 19, 18, 17, 15, 5, rep(20, 6))

bv.pcx <- prcomp(t(bv.clr))
biplot(bv.pcx, cex=c(0.7,0.1), var.axes=F, col=c('black', rgb(1,0,0,0.1)))
#######works to here

#pdf("BV_only.pdf", height=9, width=9)
codaSeq.PCAplot(bv.pcx, plot.groups=F, plot.circles=F,
    plot.loadings=TRUE, loadings.grp=loadings.group.bv,
    loadings.col=tax.colours.bv, loadings.sym=19, main="", PC=c(1,3))

legend(-320,300, legend=names(loadings.group.bv),
    col=tax.colours.bv, pch=19, cex=0.5)
# can change the size of the legend with cex - small bc it was covering - cc
#dev.off()


####### Count of data that are in common organisms
####### use for heatmap

# adding in many other taxa don't change the proportion mapping very
# much, except for 001A which is wonky anyway
# 013A is mapping only 17-18% of reads
# remove from analysis.

# count.by.tax <- matrix(data=NA, nrow=length(tax.order.bv), ncol=ncol(bv.filt))

#for(i in 1:length(tax.order.bv)){
#	count.by.tax[i, ] <- colSums(bv.filt[which(tax.vec.bv == tax.order.bv[i]),])
#}
#
# apply(count.by.tax, 2, sum)/ colSums(bv.filt)
