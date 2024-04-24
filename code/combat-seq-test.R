BiocManager::install("sva")
library(sva)

devtools::load_all('~/Documents/0_git/ALDEx.bioc')
# or library(ALDEx2)

# load test dataset
load('Rdata/egg.both.Rda')

conds <- c(rep('h', 8), rep('b', 28), rep('h',8))
btch <- c(rep('1',22), rep('2', 22))

egg.btch <- ComBat_seq(as.matrix(egg.both), batch=btch, group=conds, full_mod=F)

btch.clr <- apply(egg.btch+0.5, 2, function(x) log(x) - mean(log(x)))

# the biplots look quite different
# would be interesting to see how the species plot on this
biplot(prcomp(t(both.clr)), cex=c(0.5,0.1), var.axes=F)
biplot(prcomp(t(btch.clr)), cex=c(0.5,0.1), var.axes=F)



x.both <- aldex.clr(egg.both, conds=conds)
x.btch <- aldex.clr(egg.btch, conds=conds)

# these plots look very similar
aldex.plot(x.e.both, test='effect')
aldex.plot(x.e.btch, test='effect')


# in fact, the effect sizes are barely changed
# so ALDEx2 must be relatively unaffected by batch effects
plot(x.e.both$effect, x.e.btch$effect)

# example batch correcting
# from biplot_all.R
both.filt.btch <- ComBat_seq(as.matrix(both.filt), batch=btch, group=conds, full_mod=F)
