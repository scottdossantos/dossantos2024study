# identifying signature pathways for health vs. BV and BV subgroups

#################################### setup ####################################

# install the coda4microbiome package
# install.packages("coda4microbiome")
# BiocManager::install("ComplexHeatmap") # required but not installed above

library(coda4microbiome)
library(ALDEx2)

# path to github
user <- "sds/cc"
path.to.github <- switch(user,
                         "sds/cc" = "~/Documents/GitHub/dossantos2024study/",
                         "gg" = "~/Documents/0_git/projects/dossantos2024study/")

# load virginia KO dataset, group info for H/BV and BV1/BV2, lookup tables for
# KO and taxonomy
load(paste0(path.to.github,"Rdata/ko.virginia.filt.Rda"))
load(paste0(path.to.github,"Rdata/virginia.bv.groups.Rda"))
load(paste0(path.to.github,"Rdata/virginia.groups.Rda"))
load(paste0(path.to.github,"Rdata/tax.vec.virginia.Rda"))
load(paste0(path.to.github,"Rdata/virginia.ko.path.filt.Rda"))

# subset dataset for only BV subgroups, transposes and check for zero sum rows
ko.virginia.filt.bv <- ko.virginia.filt[,rownames(virginia.bv.groups)]
any(rowSums(ko.virginia.filt.bv)==0) # returns FALSE

# define group conds
conds.h.bv <- factor(virginia.groups$groups.2)
conds.bv1.bv2 <- factor(virginia.bv.groups$group)

################################ signature KOs ################################

# TAKES AGES DO NOT RUN ------------------------------------------

# # look for microbial signatures by KO number
# 
# # health vs. BV
# sigs.ko.health.bv <- coda_glmnet(x = t(ko.virginia.filt), y = conds.h.bv,
#                               showPlots = TRUE)
# 
# # BV1 vs. BV2
# sigs.ko.bv.subs <- coda_glmnet(x = t(ko.virginia.filt.bv), y = conds.bv1.bv2,
#                             showPlots = TRUE)
# TAKES AGES DO NOT RUN ------------------------------------------

######################## signature functions: H vs. BV ########################

# aggregate tables by functions present in pathway lookup table and set row
# names
func.h.bv <- aggregate(ko.virginia.filt,
                       by = list(virginia.ko.path.filt$pathway), FUN = sum)
rownames(func.h.bv) <- func.h.bv$Group.1
func.h.bv$Group.1 <- NULL

# CLR transform data using low variance, high abundance (housekeeping) functions
# as the denominator when calculating the geometric mean
set.seed(1234)
clr.h.bv <- aldex.clr(reads = func.h.bv, conds = as.character(conds.h.bv),
                      mc.samples = 128, denom = "lvha")

# get median CLR values for all pathways & samples
clr.h.bv.effect <- aldex.effect(clr = clr.h.bv, include.sample.summary = TRUE)

# get t-test results and combine with effect data
clr.h.bv.ttest <- aldex.ttest(clr = clr.h.bv, verbose = TRUE)
clr.h.bv.combined <- cbind(clr.h.bv.effect,clr.h.bv.ttest)

# exponentiate CLR values to undo log2 transform and retrieve pseudocounts
prop.h.bv <- 2^clr.h.bv.effect[,grep("rab.sample.", colnames(clr.h.bv.effect))]

# extract only pathways with low variance and high effect scores, then remove
# 'Unknown' pathways
prop.h.bv.lvhe <- prop.h.bv[clr.h.bv.effect$diff.win <3 & abs(clr.h.bv.effect$effect) >0.5,]
prop.h.bv.lvhe <- prop.h.bv[c(1:40, 42),]

# # check if the data is centred around housekeeping functions: clr transform
# # pseudocounts, then calculate within-group dispersion (standard deviation) and
# # between-group mean CLR difference
# prop.h.bv.clr <- apply(prop.h.bv, 2, function(x){log2(x) - mean(log2(x))})
# prop.h <- conds.h.bv == "Healthy"
# prop.bv <- conds.h.bv == "BV"
# prop.h.bv.dif <- apply(prop.h.bv.clr, 1, function(x){mean(x[prop.h]) - mean(x[prop.bv])})
# prop.h.bv.sd <- apply(prop.h.bv.clr, 1, function(x){sqrt(var(x[prop.h]) + var(x[prop.bv]))})

# plot dispersion and difference: HK functions offset as no scale included!
# plot(x = prop.h.bv.sd, prop.h.bv.dif, pch = 1, cex = 1, xlim = c(0.285,7.5),
#      xlab = expression("Median Log"[2]~"Dispersion"),
#      ylab = expression("Median Log"[2]~"Difference"),
#      main = "CLR-MW plot: centering of low variance/high abundance functions")
# points(prop.h.bv.sd[c("Ribosome","Glycolysis / Gluconeogenesis","Aminoacyl-tRNA biosynthesis", "RNA polymerase")],
#        prop.h.bv.dif[c("Ribosome","Glycolysis / Gluconeogenesis","Aminoacyl-tRNA biosynthesis","RNA polymerase")],
#        col = "red")
# abline(a = 0, b = 1, col = "grey", lty = 2)
# abline(a = 0, b = -1, col = "grey", lty = 2)
# abline(h = 0, col = "blue", lty = 2)

# calculate signature functions for health vs. BV (samples expected in rows and
# features expected in columns). NOTE: although there is no set.seed() call,
# the features that are defined as signatures for health and BV are robust and
# the same functions will be identified when running this command repeatedly
sigs.health.bv <- coda_glmnet(x = t(prop.h.bv.lvhe), y = conds.h.bv,
                              showPlots = T, alpha = 4, coef_threshold = 0.3,
                              lambda = 'lambda.1se', nfolds = 100)

# pull list of functions where the log-contrast coefficients are large
large <- abs(sigs.health.bv$`log-contrast coefficients`) > 0.05

# create MW plot of all functions and overlay functions which are maximally 
# distinguishable
# png("~/virginia_health_bv_DA_maxDist.png", units = "in",
#     width = 6, height = 6, res = 400)
aldex.plot(clr.h.bv.combined, test = "effect", type = "MW", xlim = c(0.315, 8.2))
title("ALDEx2: low var / high abund functions - health vs. BV",
      adj = 0, line = 0.8)
points(clr.h.bv.effect[sigs.health.bv$taxa.name[large],"diff.win"],
       clr.h.bv.effect[sigs.health.bv$taxa.name[large],"diff.btw"],
       col = "blue", cex = 1.2)
legend("bottomleft", legend = c("Rare","Non-sig.","Sig.", "Max dist."),
       col = c("black", rgb(0,0,0,0.2), "red", "blue"), pch = c(19,19,19,1),
       pt.cex = c(0.4,0.4,0.9,0.9), cex = 0.85, ncol = 1)
# dev.off()

# save signature and prediction plots for health vs. BV
# png("~/virginia_health_bv_signaturePlot.png", units = "in",
#     width = 11, height = 6, res = 300)
# sigs.health.bv$`signature plot`
# dev.off()
# 
# png("~/virginia_health_bv_predictionPlot.png", units = "in",
#     width = 11, height = 6, res = 300)
# sigs.health.bv$`predictions plot`
# dev.off()

###################### signature functions: BV subgroups ######################

# aggregate subset of the data (only BV subgroup 1 and 2 samples) by function
func.bv.subs <- aggregate(ko.virginia.filt.bv,
                          by = list(virginia.ko.path.filt$pathway), FUN = sum)
rownames(func.bv.subs) <- func.bv.subs$Group.1
func.bv.subs$Group.1 <- NULL

# set the seed for RNG then CLR transform using low-variance high-abundance
# functions as the denominator for geometric mean
set.seed(1234)
clr.bv.subs <- aldex.clr(reads = func.bv.subs, conds = as.character(conds.bv1.bv2),
                         mc.samples = 128, denom = "lvha")

# get median CLR values for functions
clr.bv.subs.effect <- aldex.effect(clr = clr.bv.subs, verbose = TRUE,
                                   include.sample.summary = TRUE)

# calculate t-test statistics and combine with effect data
clr.bv.subs.ttest <- aldex.ttest(clr = clr.bv.subs, verbose = TRUE)
clr.bv.subs.combined <- cbind(clr.bv.subs.effect, clr.bv.subs.ttest)

# exponentiate CLR values to undo log2 transform and retrieve pseudocounts, then
# extract only functions with low variance and high effect scores
prop.bv.subs <- 2^clr.bv.subs.effect[,grep("rab.sample.", colnames(clr.bv.subs.effect))]
prop.bv.subs.lvhe <- prop.bv.subs[clr.bv.subs.effect$diff.win <4 & abs(clr.bv.subs.effect$effect) >0.5,]

# check if the data is centred around housekeeping functions: clr transform
# pseudocounts, then calculate within-group dispersion (standard deviation) and
# between-group mean CLR difference
prop.bv.subs.clr <- apply(prop.bv.subs, 2, function(x){log2(x) - mean(log2(x))})
prop.bv1 <- conds.bv1.bv2 == "BV Subgroup 1"
prop.bv2 <- conds.bv1.bv2 == "BV Subgroup 2"
prop.bv.subs.dif <- apply(prop.bv.subs.clr, 1, function(x){mean(x[prop.bv1]) - mean(x[prop.bv2])})
prop.bv.subs.sd <- apply(prop.bv.subs.clr, 1, function(x){sqrt(var(x[prop.bv1]) + var(x[prop.bv2]))})

# plot dispersion and difference: HK functions offset as no scale included!
plot(x = prop.bv.subs.sd, prop.bv.subs.dif, pch = 1, cex = 0.75,
     xlab = expression("Median Log"[2]~"Dispersion"),
     ylab = expression("Median Log"[2]~"Difference"),
     main = "CLR-MW plot: centering of low variance/high abundance functions")
points(prop.bv.subs.sd[c("Ribosome","Glycolysis / Gluconeogenesis","Aminoacyl-tRNA biosynthesis", "RNA polymerase")],
       prop.bv.subs.dif[c("Ribosome","Glycolysis / Gluconeogenesis","Aminoacyl-tRNA biosynthesis","RNA polymerase")],
       col = "red")
abline(a = 0, b = 1, col = "grey", lty = 2)
abline(a = 0, b = -1, col = "grey", lty = 2)
abline(h = 0, col = "blue", lty = 2)

# calculate signature functions for BV subgroups
sigs.bv.subs <- coda_glmnet(x = t(prop.bv.subs.lvhe), y = conds.bv1.bv2,
                            showPlots = TRUE, alpha = 4, coef_threshold = 0.3,
                            lambda = 'lambda.1se', nfolds = 100)

# create MW plot of all functions and overlay functions which are maximally 
# distinguishable
# png("~/virginia_bvSubgroups_DA_maxDist.png", units = "in",
#     width = 6, height = 6, res = 400)
aldex.plot(clr.bv.subs.combined, test = "effect", type = "MW", xlim = c(0.3075,7.5))
title("ALDEx2: low var / high abund functions - BV1/BV2",
      adj = 0, line = 0.8)
points(clr.bv.subs.effect[sigs.bv.subs$taxa.name,"diff.win"],
       clr.bv.subs.effect[sigs.bv.subs$taxa.name,"diff.btw"],
       col = "blue", cex = 1.2)
legend("bottomleft", legend = c("Rare","Non-sig.","Sig.", "Max dist."),
       col = c("black",rgb(0,0,0,0.2), "red", "blue"), pch = c(19,19,19,1),
       pt.cex = c(0.4,0.4,0.9,0.9), cex = 0.85, ncol = 1)
# dev.off()


# save signature and prediction plots for BV subgroups
# png("~/virginia_bvSubgroups_signaturePlot.png", units = "in",
#     width = 11, height = 6, res = 300)
# sigs.bv.subs$`signature plot`
# dev.off()
# 
# png("~/virginia_bvSubgroups_predictionPlot.png", units = "in",
#     width = 11, height = 6, res = 300)
# sigs.bv.subs$`predictions plot`
# dev.off()

