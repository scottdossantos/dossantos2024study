library(codacore)

site <- c(rep('L', 14), rep('E',14))
bv <- c(rep('H', 8), rep('B', 14), rep('B',14), rep('H', 8))

response <- data.frame(site)

# so codacore is somewhat stable

model.am <- codacore(t(egg.bv + 0.5), bv.conds, logRatioType='SLR')

model.bal <- codacore(t(egg.bv + 0.5), bv.conds)