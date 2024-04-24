# healthy only biplog
# depends only on setup.r

source("code/setup.R")

ec.both.n <- ec.both.all[,grep('^h', colnames(ec.both.all))]
 
 # don't think i need...
 # 015B appears to be wierd
 #both.n$h.015B <- NULL
 
ec.both.n.filt <- codaSeq.filter(ec.both.n,min.prop=0.00005, min.occurrence=0.3, samples.by.row=FALSE)
 
ec.both.n.clr <- apply(ec.both.n.filt +0.5, 2, function(x) log(x) - mean(log(x)))
 
ec.both.n.pcx <- prcomp(t(ec.both.n.clr))

biplot(ec.both.n.pcx, cex=c(0.7,0.2), var.axes=FALSE, main="EC Healthy biplot", col=c("black", rgb(1,0,0,0.4) ))
# changing transparency (0.1)
