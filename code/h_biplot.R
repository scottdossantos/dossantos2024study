# healthy only biplog
# depends only on setup.r

source("code/setup.R")

 both.n <- both.data[,grep('^h', colnames(both.data))]
 
 # 015B appears to be wierd
 both.n$h.015B <- NULL
 
 both.n.filt <- codaSeq.filter(both.n,min.prop=0.00005, min.occurrence=0.3, samples.by.row=FALSE)
 
 both.n.clr <- apply(both.n.filt +0.5, 2, function(x) log(x) - mean(log(x)))
 
 both.n.pcx <- prcomp(t(both.n.clr))
 
# biplot(both.n.pcx, var.axes=F, cex=c(1,0.5), col=c('black', rgb(1,0,0,0.05)))
 
 
 # assign a taxonomy to all rows
tax.vec.h <- as.vector(tax.table[rownames(both.n.filt),2])
names(tax.vec.h) <- rownames(both.n.filt)

# lets get the total number of genes in each species
n.genes.h <- vector()
for(i in levels(factor(tax.vec.h))){
n.genes.h[i] <- length(which(tax.vec.h == i))
}

tax.order.h <- names(which(n.genes.h > 50))
tax.colours.h <- c(tax.colors[tax.order.h,], rgb(0,0,0,0.075))

loadings.group.h <- list()

for(i in 1:length(tax.order.h)){
    loadings.group.h[[tax.order.h[i]]] <- which(tax.vec.h == tax.order.h[i])
}

loadings.group.h[["other"]] <- setdiff(c(1:nrow(both.n.clr)),unlist(loadings.group.h))
tax.symbol <- c(20, 19, 18, 17, 15, 14)

# loadings.sym parameter not working - fix in codaSeq package
codaSeq.PCAplot(both.n.pcx, plot.groups=F, plot.circles=F,
    plot.loadings=TRUE, loadings.grp=loadings.group.h,
    loadings.col=tax.colours.h, loadings.sym=19, main="", PC=c(1,2))

legend(-320,325, legend=names(loadings.group.h),
    col=tax.colours.h, pch=19, cex=0.7)


biplot(both.n.pcx, cex=c(0.7,0.2), var.axes=FALSE, main="Healthy biplot", col=c("black", rgb(1,0,0,0.1) ))
# changing transparency (0.1)
