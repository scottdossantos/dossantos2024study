# healthy only biplog
# depends only on setup.r

source("code/setup.R")

 egg.both.n <- egg.both.all[,grep('^h', colnames(egg.both.all))]
 
 # don't think i need...
 # 015B appears to be wierd
 #both.n$h.015B <- NULL
 
 egg.both.n.filt <- codaSeq.filter(egg.both.n,min.prop=0.00005, min.occurrence=0.3, samples.by.row=FALSE)
 
 egg.both.n.clr <- apply(egg.both.n.filt +0.5, 2, function(x) log(x) - mean(log(x)))
 
 egg.both.n.pcx <- prcomp(t(egg.both.n.clr))
 
biplot(egg.both.n.pcx, cex=c(0.7,0.2), var.axes=FALSE, col=c("black", rgb(1,0,0,0.3) ))

 
 # assign a taxonomy to all rows
egg.tax.vec.h <- as.vector(tax.table[rownames(egg.both.n.filt),2])
names(egg.tax.vec.h) <- rownames(egg.both.n.filt)


# lets get the total number of genes in each species
egg.n.genes.h <- vector()
for(i in levels(factor(egg.tax.vec.h))){
egg.n.genes.h[i] <- length(which(egg.tax.vec.h == i))
}

egg.tax.order.h <- names(which(egg.n.genes.h > 50))
egg.tax.colours.h <- c(tax.colors[egg.tax.order.h,], rgb(0,0,0,0.075))

egg.loadings.group.h <- list()

for(i in 1:length(egg.tax.order.h)){
    egg.loadings.group.h[egg.tax.order.h[i]] <- which(egg.tax.vec.h == egg.tax.order.h[i])
}

egg.loadings.group.h[["other"]] <- setdiff(c(1:nrow(egg.both.n.clr)),unlist(egg.loadings.group.h))
egg.tax.symbol <- c(20, 19, 18, 17, 15, 20)

#codaSeq.PCAplot(both.n.pcx, plot.groups=F, plot.circles=F,
    #plot.loadings=TRUE, loadings.grp=loadings.group.h,
    #loadings.col=tax.colours.h, loadings.sym=tax.symbol, main="", PC=c(1,2))

#legend(-320,325, legend=names(loadings.group.h),
    #col=tax.colours.h, pch=tax.symbol, cex=0.7)


biplot(egg.both.n.pcx, cex=c(0.7,0.2), var.axes=FALSE, col=c("black", rgb(1,0,0,0.1) ))
# changing transparency (0.1)
