source("code/setup.R")
# compare BV only from two datasets
# pull from the merged datasets

ec.bv.data <- data.frame(ec.both.all[,grep('v.', colnames(ec.both.all))])

ec.bv.filt <-  codaSeq.filter(ec.bv.data, min.occurrence=0.30,
    min.prop=0.00005, samples.by.row=F)
ec.bv.clr <- apply(ec.bv.filt + 0.5, 2, function(x) log(x) - mean(log(x)))


ec.bv.pcx <- prcomp(t(ec.bv.clr))
biplot(ec.bv.pcx, cex=c(0.6,0.15), var.axes=F, main="EC BV biplot", col=c('black', rgb(1,0,0,0.7)))
