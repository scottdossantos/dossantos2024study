# combined plot for both full datasets
# source("code/setup.R")

ec.clr.input <- cmultRepl(t(both.filt), method='CZM', label=0)

ec.both.clr <- apply(t(ec.clr.input) , 2, function(x) log(x) - mean(log(x)))

ec.both.pcx <- prcomp(t(ec.both.clr))
biplot(ec.both.pcx, cex=c(0.5,0.1), var.axes=F, main="EC ALL biplot", col=c('black',rgb(0,0,0,0.1)))


