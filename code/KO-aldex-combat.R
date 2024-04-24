# this shows the CAMP pathway is significantly up in H vs BV
# replicates in EC (only 1 enzyme though)
# replicates for 3 of 4 eggnog orthologs

library(sva)

devtools::load_all('~/Documents/0_git/ALDEx.bioc')
# or library(ALDEx2)

# load test dataset
source('code/setup.R')

conds <- c(rep('h', 8), rep('b', 28), rep('h',8))
btch <- c(rep('1',22), rep('2', 22))

ko.btch <- ComBat_seq(as.matrix(ko.both.all), batch=btch, group=conds, full_mod=F)

gamma=1
mu.vec <- gsub('h', log2(mu1), conds)
mu.vec <- as.numeric(gsub('b', log2(mu2), mu.vec))
mu.mod <- sapply(mu.vec, FUN = function(mu) rlnorm(128, mu, gamma))


x.all <- aldex(ko.btch, conditions=conds, gamma=t(mu.mod))
# FDR 0.05
# eg:  rownames(x.all)[x.all$we.eBH < 0.05 & x.all$diff.btw < 1]

camp.lacto <- c('K03367', 'K03739', 'K03740', 'K14188')
camp.bv <- c("K00677","K01364","K01448","K03585","K03760","K03767","K07264","K12340","K14205")
par(mfrow=c(1:3))
aldex.plot(x.all, type='MA')
 points(x.all[camp.lacto,'rab.all'], (x.all[camp.lacto,'diff.btw']))
 points(x.all[camp.bv,'rab.all'], (x.all[camp.bv,'diff.btw']))
 
 aldex.plot(x.all)
  points(x.all[camp.lacto,'diff.win'], (x.all[camp.lacto,'diff.btw']))
  points(x.all[camp.bv,'diff.win'], (x.all[camp.bv,'diff.btw']))

 aldex.plot(x.all, type='volcano')
 points(x.all[camp.lacto,'diff.btw'], -1*log10(x.all[camp.lacto,'we.eBH']))
 points(x.all[camp.bv,'diff.btw'], -1*log10(x.all[camp.bv,'we.eBH']))
 
# K03739 = COG1696
# K03367 = COG1020
# K03740 = COG3966
# K14188 = COG0236

egg.btch <- ComBat_seq(as.matrix(egg.both.all), batch=btch, group=conds, full_mod=F)

camp.e <- c('COG1696', 'COG1020', 'COG3966', 'COG0236')
 
x.all <- aldex(egg.btch, conditions=conds, gamma=t(mu.mod))
par(mfrow=c(1:3))
aldex.plot(x.all, type='MA')
 points(x.all[camp.e,'rab.all'], (x.all[camp.e,'diff.btw']))
 
 aldex.plot(x.all)
  points(x.all[camp.e,'diff.win'], (x.all[camp.e,'diff.btw']))

 aldex.plot(x.all, type='volcano')
 points(x.all[camp.e,'diff.btw'], -1*log10(x.all[camp.e,'we.eBH']))



ec.btch <- ComBat_seq(as.matrix(ec.both.all), batch=btch, group=conds, full_mod=F)

camp.ec <- c('ec:6.1.1.13')
 
x.all <- aldex(ec.btch, conditions=conds, gamma=t(mu.mod))

par(mfrow=c(1:3))
aldex.plot(x.all, type='MA')
 points(x.all[camp.ec,'rab.all'], (x.all[camp.ec,'diff.btw']))
 
 aldex.plot(x.all)
  points(x.all[camp.ec,'diff.win'], (x.all[camp.ec,'diff.btw']))

 aldex.plot(x.all, type='volcano')
 points(x.all[camp.ec,'diff.btw'], -1*log10(x.all[camp.ec,'we.eBH']))
 
