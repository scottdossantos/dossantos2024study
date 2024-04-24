# comparing ALDEx2 with and without scale simulation: london/europe dataset

#################################### setup ####################################

# load ALDEx2
library(ALDEx2) # v1.35 (bioconductor v3.19)

# NOTE - to install the latest version of ALDEx2 that incorporates scale 
#        simulation, use devtools to install ALDEx2 from the GitHub repository:
#        devtools::install_github("ggloor/ALDEx_bioc")

# set path to github repository (change to match your machine)
user <- "gg"
path.to.github <- switch(user,
                         "sds/cc" = "~/Documents/GitHub/metatranscriptome/",
                         "gg" = "~/Documents/0_git/projects/metatranscriptome/")

# load in KO-aggregated lon/eur feature table and KO to pathway table
load(paste0(path.to.github, "Rdata/ko.both.Rda"))
load(paste0(path.to.github, "Rdata/ko.both.path.Rda"))

# remove two samples from the filtered london/europe feature table aggregated by
# K0 number (classed as BV but almost no BV organisms):
#   -  v.001A: close to 100% L. gasseri with practically no BV organisms
#   -  v.019A: around 80 % iners, ~ 20% crispatus and a tiny bit of Gardnerella
ko.both <- ko.both[,-c(9,22)]

# make conditions vector (corresponds to lon/eur samples with 2 removed)
lon.eur.conds <- c(rep('H',8), rep('B',12), rep('B',14), rep('H', 8)) 

# remove non-bacterial KOs from the K number-aggregated, filtered feature table
# (these were discovered during curation of 'Unknown' pathways)
#    - K03364: Eukaryotic cell division cycle 20-like protein 1
#    - K13963: Serpin B (eukaryotic serine protease inhibitor)
#    - K01173: Mitochondrial endonuclease G
#    - K12373: Human lysosomal hexosaminidase
#    - K14327: Eukaryotic regulator of nonsense transcripts 2
#    - K00863: Human triose/dihydroxyacetone kinase
#    - K00599: Eukaryotic tRNA N(3)-methylcytidine methyltransferase
#    - K13993: Human HSP20
#    - K00811: Chloroplastic aspartate aminotransferase
#    - K03260: Eukaryotic translation initiation factor 4G
#    - K00985: Enterovirus RNA-directed RNA polymerase

ko.both <- ko.both[-which(grepl(paste("K03364","K13963","K01173","K12373",
                                      "K14327","K00863","K00599","K13993",
                                      "K00811","K03260","K00985", sep = "|"),
                                rownames(ko.both))),]

# pull housekeeping functions
path.ribo <- which(ko.both.path$pathway == "Ribosome")
path.trna <- which(ko.both.path$pathway == "Aminoacyl-tRNA biosynthesis")
path.glyc <- which(ko.both.path$pathway == "Glycolysis / Gluconeogenesis")
path.hk <- c(path.ribo, path.trna, path.glyc)

############################### aldex2: no scale ###############################

# set seed for RNG
set.seed(2023)

# run aldex.clr
# gg changed scale to 1e-3 so we can get the background scale
aldex.le.ns.clr <- aldex.clr(reads = ko.both, conds = lon.eur.conds,
                             mc.samples = 128, denom = "all", gamma = 1e-3,
                             verbose = TRUE)

aldex.le.ns.clr.e <- aldex.effect(clr = aldex.le.ns.clr, verbose = TRUE,
                                  include.sample.summary = TRUE)

aldex.le.ns.clr.t <- aldex.ttest(clr = aldex.le.ns.clr,
                                 verbose = TRUE)

aldex.le.ns.clr.all <- cbind(aldex.le.ns.clr.e,
                             aldex.le.ns.clr.t)

## scale is in aldex.le.ns.clr@scaleSamps
bg.scale <- aldex.le.ns.clr@scaleSamps
# mean(rowMeans(bg.scale)[lon.eur.conds=="B"]) 
# = 14.25, 4.53 from the final scaled data
# mean(rowMeans(bg.scale)[lon.eur.conds=="H"])
# = 17.20, 5.2 from the final scaled data
############################ aldex2: default scale ############################

# set seed for RNG
set.seed(2023)

# run aldex.clr with gamma set to 0.5
aldex.le.g.clr <- aldex.clr(reads = ko.both, conds = lon.eur.conds,
                            mc.samples = 128, denom = "all", gamma = 0.5,
                            verbose = TRUE)

aldex.le.g.clr.e <- aldex.effect(clr = aldex.le.g.clr, verbose = TRUE,
                                  include.sample.summary = TRUE)

aldex.le.g.clr.t <- aldex.ttest(clr = aldex.le.g.clr,
                                 verbose = TRUE)

aldex.le.g.clr.all <- cbind(aldex.le.g.clr.e,
                            aldex.le.g.clr.t)

############################## aldex2: full scale ##############################

# full scale model of lon/eur data with scale difference of 15% and gamma = 0.5 
# was generated in 'lon_eur_health_vs_BV.R' as part of the analysis
load(paste0(path.to.github, "Rdata/ko.both.clr.all.Rda"))

# rename for consistency with non-scaled data
aldex.le.s.clr.all <- ko.both.clr.all
rm(ko.both.clr.all)

########################## visualise effect of scale ##########################

# MW plots (effect plots) ------------------------------------------------------

# make three MW plots side-by-side showing the results of: 
#     A) unscaled ALDEx2
#     B) scaled ALDEx2 with gamma of 0.5 (default scale model)
#     C) scaled ALDEx2 with gamma of 0.5 and scale difference between groups of
#        15 % (full scale model)

# NOTE: data points are coloured as below:
#       grey = non-significant feature (BH-corrected P-value >0.01)
#       black = rare feature (CLR abundance < 0)
#       orange = housekeeping gene (ribosome, aatRNA biosynthesis, glycolysis)
#       red = significant feature (BH-corrected P-value <0.01)
#       blue outline = significant feature with absolute effect size >1

# png(paste0(path.to.github,"figs_for_paper/R_FigsLondonEurope/suppl_lon_eur_scaleMWs.png"),
#     units = "in", height = 6, width = 12, res = 600)

par(mfrow= c(1,3))

# panel A: unscaled aldex2
aldex.plot(aldex.le.ns.clr.all, type = "MW", cutoff.pval = 0.01, xlim = c(0.325,9))
title("A)  ALDEx2 unscaled - (gamma = 0, ~3% difference)", adj = 0, line = 0.8)
points(aldex.le.ns.clr.all$diff.win[path.hk],
       aldex.le.ns.clr.all$diff.btw[path.hk],
       pch = 19, col = rgb(1,0.7,0), cex = 0.5)
points(aldex.le.ns.clr.all$diff.win[abs(aldex.le.ns.clr.all$effect) > 1 & aldex.le.ns.clr.all$we.eBH <0.01],
       aldex.le.ns.clr.all$diff.btw[abs(aldex.le.ns.clr.all$effect) > 1 & aldex.le.ns.clr.all$we.eBH <0.01],
       pch = 21, col = rgb(0,0,1,0.5), cex = 0.5)
abline(h = 0, col = rgb(0,0,1,0.5), lty = 2)


# panel B: aldex2 w/ default scale model (gamma = 0.5)
aldex.plot(aldex.le.g.clr.all, type = "MW", cutoff.pval = 0.01, xlim = c(0.325,9),
           ylab = "")
title("B)  ALDEx2 - scaled (gamma = 0.5, ~3% difference)", adj = 0, line = 0.8)
points(aldex.le.g.clr.all$diff.win[path.hk],
       aldex.le.g.clr.all$diff.btw[path.hk],
       pch = 19, col = rgb(1,0.7,0), cex = 0.5)
points(aldex.le.g.clr.all$diff.win[abs(aldex.le.g.clr.all$effect) > 1 & aldex.le.g.clr.all$we.eBH <0.01],
       aldex.le.g.clr.all$diff.btw[abs(aldex.le.g.clr.all$effect) > 1 & aldex.le.g.clr.all$we.eBH <0.01],
       pch = 21, col = rgb(0,0,1,0.5), cex = 0.5)
abline(h = 0, col = rgb(0,0,1,0.5), lty = 2)


# panel C: aldex2 w/ full scale model (gamma = 0.5, scale difference of 15%)
aldex.plot(aldex.le.s.clr.all, type = "MW", cutoff.pval = 0.01, xlim = c(0.325,9),
           ylab = "")
title("C)  ALDEx2 - scaled (gamma = 0.5, 15% difference)", adj = 0, line = 0.8)
points(aldex.le.s.clr.all$diff.win[path.hk],
       aldex.le.s.clr.all$diff.btw[path.hk],
       pch = 19, col = rgb(1,0.7,0), cex = 0.5)
points(aldex.le.s.clr.all$diff.win[abs(aldex.le.s.clr.all$effect) > 1 & aldex.le.s.clr.all$we.eBH <0.01],
       aldex.le.s.clr.all$diff.btw[abs(aldex.le.s.clr.all$effect) > 1 & aldex.le.s.clr.all$we.eBH <0.01],
       pch = 21, col = rgb(0,0,1,0.5), cex = 0.5)
abline(h = 0, col = rgb(0,0,1,0.5), lty = 2)

# dev.off()

# volcano plots ----------------------------------------------------------------

# also explore volcano plots of the same sets of data:
#     A) unscaled ALDEx2
#     B) scaled ALDEx2 with gamma of 0.5 (default scale model)
#     C) scaled ALDEx2 with gamma of 0.5 and scale difference between groups of
#        15 % (full scale model)

# NOTE: horizontal grey lines = -1 * log10(p.value cutoff)
#       vertical grey lines = 1.5 (represents 2^1.5)
#       vertical blue line = 0 (represents log2 difference between groups of 0)

par(mfrow=c(1,3))

# panel A: unscaled aldex2
p.add <- min(aldex.le.ns.clr.all$we.eBH[aldex.le.ns.clr.all$we.eBH > 0])/10
all.p <- aldex.le.ns.clr.all$we.eBH + p.add

aldex.plot(aldex.le.ns.clr.all, type = "volcano", cutoff.pval = 0.01)
title("A)  ALDEx2 unscaled", adj = 0, line = 0.8)
points(aldex.le.ns.clr.all$diff.btw[path.hk],
       -1*log10(all.p)[path.hk],
       pch = 19, col = rgb(1,0.7,0), cex= 0.5)
points(aldex.le.ns.clr.all$diff.btw[abs(aldex.le.ns.clr.all$effect) > 1 & aldex.le.ns.clr.all$we.eBH <0.01],
       -1*log10(all.p)[abs(aldex.le.ns.clr.all$effect) > 1 & aldex.le.ns.clr.all$we.eBH <0.01],
       pch = 1, col = rgb(0,0,1,0.15), cex = 0.85)
abline(v = 0, col = rgb(0,0,1,0.5), lty = 1)


# panel B: aldex2 w/ default scale model (gamma = 0.5)
p.add <- min(aldex.le.g.clr.all$we.eBH[aldex.le.g.clr.all$we.eBH > 0])/10
all.p <- aldex.le.g.clr.all$we.eBH + p.add

aldex.plot(aldex.le.g.clr.all, type = "volcano", cutoff.pval = 0.01, ylab = "")
title("B)  ALDEx2 - scaled (gamma = 0.5)", adj = 0, line = 0.8)
points(aldex.le.g.clr.all$diff.btw[path.hk],
       -1*log10(all.p)[path.hk],
       pch = 19, col = rgb(1,0.7,0), cex= 0.5)
points(aldex.le.g.clr.all$diff.btw[abs(aldex.le.g.clr.all$effect) > 1 & aldex.le.g.clr.all$we.eBH <0.01],
       -1*log10(all.p)[abs(aldex.le.g.clr.all$effect) > 1 & aldex.le.g.clr.all$we.eBH <0.01],
       pch = 1, col = rgb(0,0,1,0.15), cex = 0.85)
abline(v = 0, col = rgb(0,0,1,0.5), lty = 1)


# panel C: aldex2 w/ full scale model (gamma = 0.5, scale difference of 15%)
p.add <- min(aldex.le.s.clr.all$we.eBH[aldex.le.s.clr.all$we.eBH > 0])/10
all.p <- aldex.le.s.clr.all$we.eBH + p.add

aldex.plot(aldex.le.s.clr.all, type = "volcano", cutoff.pval = 0.01, ylab = "")
title("C)  ALDEx2 - scaled (gamma = 0.5, 15% difference)", adj = 0, line = 0.8)
points(aldex.le.s.clr.all$diff.btw[path.hk],
       -1*log10(all.p)[path.hk],
       pch = 19, col = rgb(1,0.7,0), cex= 0.5)
points(aldex.le.s.clr.all$diff.btw[abs(aldex.le.s.clr.all$effect) > 1 & aldex.le.s.clr.all$we.eBH <0.01],
       -1*log10(all.p)[abs(aldex.le.s.clr.all$effect) > 1 & aldex.le.s.clr.all$we.eBH <0.01],
       pch = 1, col = rgb(0,0,1,0.15), cex = 0.85)
abline(v = 0, col = rgb(0,0,1,0.5), lty = 1)


