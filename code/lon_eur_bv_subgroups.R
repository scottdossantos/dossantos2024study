# differential abundance analysis for two BV subgroups

#################################### setup ####################################

library(dplyr) # for data manipulation; MAKE SURE 'plyr' ISN'T ALSO LOADED!!!
library(forcats) # for reordering factors
library(ggplot2) # for plotting
library(pheatmap) # for making heatmaps
library(viridisLite) # for making colour blind-friendly palettes
library(vegan) # for procrustes rotation

# set path to the github repository
user <- "sds/cc"
path.to.github <- switch(user,
                         "sds/cc" = "~/Documents/GitHub/metatranscriptome/",
                         "gg" = "~/Documents/0_git/projects/metatranscriptome/")

# run setup code to get batch-corrected feature tables and taxa vector
source(paste(path.to.github, "code/setup.R", sep=""))

# load in manually curated KO-> pathway lookup table for filtered london & 
# europe dataset
load(paste(path.to.github, "Rdata/ko.both.path.Rda", sep = ""))

# remove two samples from the filtered london/europe feature table aggregated by
# KO number (classed as BV but almost no BV organisms):
#   -  v.001A: close to 100% L. gasseri with practically no BV organisms
#   -  v.019A: around 80 % iners, ~ 20% crispatus and a tiny bit of Gardnerella
ko.both<-ko.both[,-c(9,22)]

# remove non-bacterial KOs from the K number-aggregated, filtered feature table
# (discovered during curation of 'Unknown' pathways)
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

# load in cutree output from clustering of scale sim'd CLR-transformed z-score
# data, aggregated by function (from VIRGO pathway table)
load(paste(path.to.github, "Rdata/lon.eur.bv.groups.Rda", sep = ""))

# filter for BV subgroups 1 & 2, and a healthy sample classed as BV ('h.015B')
# with no symptoms but dominated by Prevotella bivia
lon.eur.bv.groups <- lon.eur.bv.groups %>% 
  filter(group %in% c("BV1","BV2") & !rownames(lon.eur.bv.groups) %in% c("h.015B","v.008A",
                                                                         "v.a8_3","v.a6_4",
                                                                         "v.014B","v.a4_1"))

# load in previously generated cst and dataset info for london/europe samples,
# subset for BV samples used in this sub-analysis and bind with subgroup data
load(paste(path.to.github, "Rdata/hm.metadata.Rda", sep = ""))
lon.eur.bv.meta <- hm.metadata[rownames(lon.eur.bv.groups),]
lon.eur.bv.meta <- data.frame(Dataset = lon.eur.bv.meta$Dataset,
                              Subgroup = lon.eur.bv.groups$group,
                              row.names = rownames(lon.eur.bv.groups))

##################### BV subgroups: differential abundance #####################

# if ALDEx2 .clr object has previously been generated, load it in. Otherwise,
# subset the KO-aggregated feature table for the BV subgroups, make a matrix of
# scale values from a log-normal distribution (stdev 0.5, 15% difference), run
# ALDEx2 with scale-simulation, calculate effect sizes & t-test statistics, then
# combine the outputs into a summary data frame

if(file.exists(paste(path.to.github,"Rdata/ko.both.2bv.clr.all.Rda",sep = ""))){
  load(paste(path.to.github,"Rdata/ko.both.2bv.clr.all.Rda",sep = ""))
} else{
  
# subset ko.both and ko.both.all by BV1 & BV2 samples
  ko.both.2bv <- ko.both[,(rownames(lon.eur.bv.groups))]
  any(rowSums(ko.both.2bv)==0) # returns FALSE
  
# set seed for RNG
  set.seed(2023)
  
# generate matrix of scale values with mean difference of 10%, std dev of
# log-normal distribution = 0.5
  mu.mat <- aldex.makeScaleMatrix(gamma = 0.5, mu = c(1, 1.10),
                                  lon.eur.bv.groups$group)
  
# transform
  ko.both.2bv.clr <- aldex.clr(ko.both.2bv, lon.eur.bv.groups$group,
                               gamma = mu.mat, verbose = T)
  
# calculate effect sizes, plus median CLR values (all MC instances) of all KOs
# for each sample
  ko.both.2bv.clr.e <- aldex.effect(ko.both.2bv.clr, include.sample.summary = T)
  
# calculate t-test statistics
  ko.both.2bv.clr.t <- aldex.ttest(ko.both.2bv.clr, verbose = T)
  
# combine effect size and t-test data into summary data frame
  ko.both.2bv.clr.all <- cbind(ko.both.2bv.clr.e, ko.both.2bv.clr.t) 
  
# save summary ALDEx2 output as .Rda object
  save(ko.both.2bv.clr.all, 
       file = paste(path.to.github,"Rdata/ko.both.2bv.clr.all.Rda",sep = ""))
  
}

# extract indices of K0s assigned to 'housekeeping' pathways (corresponds to K0
# feature table rows)
path.ribo <- which(ko.both.path$pathway == "Ribosome")
path.trna <- which(ko.both.path$pathway == "Aminoacyl-tRNA biosynthesis")
path.glyc <- which(ko.both.path$pathway == "Glycolysis / Gluconeogenesis")

# plot dispersion of K0s vs. mean CLR difference between groups
# png(paste(path.to.github,
#           "figs_for_paper/lon_eur_BVSubgroups_diffAbund_MW.png", sep = ""),
#     units = "in", height = 7, width = 14, res = 400)

aldex.plot(ko.both.2bv.clr.all, xlim=c(0.346,9), cutoff.pval=0.01)
title('ALDEx2 w/ ScaleSim (K0-aggregated London & Europe - BV subgroups): mu = 1 : 1.10, gamma = 0.5, p <0.01',
      adj=0, line= 0.8)

points(ko.both.2bv.clr.all[path.ribo ,"diff.win"],
       ko.both.2bv.clr.all[path.ribo ,"diff.btw"],
       pch= 19, col= "blue3", cex= 0.5)

points(ko.both.2bv.clr.all[path.trna,"diff.win"],
       ko.both.2bv.clr.all[path.trna,"diff.btw"],
       pch= 19, col= "royalblue1", cex= 0.5)

points(ko.both.2bv.clr.all[path.glyc,"diff.win"],
       ko.both.2bv.clr.all[path.glyc,"diff.btw"],
       pch= 19, col= "lightskyblue1", cex= 0.5)

points(ko.both.2bv.clr.all[which(abs(ko.both.2bv.clr.all$effect) >1),"diff.win"],
       ko.both.2bv.clr.all[which(abs(ko.both.2bv.clr.all$effect) >1),"diff.btw"],
       pch= 19, col= "purple3", cex= 0.5)

legend("bottomleft", pch = 19, cex = 0.75,
       col = c( "black", "purple3", "blue3", "royalblue1", "lightskyblue1"),
       legend = c("Non-significant", "Absolute effect size >1",
                  "Ribosome", "Aminoacyl-tRNA biosynthesis", "Glycolysis / Gluconeogenesis"))

# dev.off()

# pull all K0s with absolute effect size > 1 and the corresponding effect sizes
bv2.sig.ko <- ko.both.path[which(abs(ko.both.2bv.clr.all$effect) >1), ,
                           drop = FALSE]
bv2.sig.ko$effect <- ko.both.2bv.clr.all[rownames(bv2.sig.ko), 'effect']

# filter out K0s with unknown pathway
bv2.sig.ko <- bv2.sig.ko %>% 
  filter(pathway != "Unknown")

# append " (BV1)" or " (BV2)" to pathways if there is >1 K0 for each pathway
# AND there are positive and negative effect sizes for that pathway
# note: this code is disgusting to read and I am well aware of it. Sorry.
for(i in levels(factor(bv2.sig.ko$pathway))){
  if(length(which(bv2.sig.ko$pathway == i)) >1 & length(unique(sign(bv2.sig.ko$effect[which(bv2.sig.ko$pathway == i)]))) != 1){
    bv2.sig.ko$pathway <- case_when(bv2.sig.ko$pathway == i & sign(bv2.sig.ko$effect) == 1 ~ paste(i, "(BV2)", sep = " "),
                                            bv2.sig.ko$pathway == i & sign(bv2.sig.ko$effect) == -1 ~ paste(i, "(BV1)", sep = " "),
                                            .default = bv2.sig.ko$pathway)
  }
}

# get vector indicating whether effect size is positive or negative
bv2.sig.ko$group <- sign(bv2.sig.ko$effect)
bv2.sig.ko$group[which(bv2.sig.ko$group == 1)] <- "BV2"
bv2.sig.ko$group[which(bv2.sig.ko$group == -1)] <- "BV1"

# collapse data frame by function, calculate median effect size and count number
# of significantly different KO terms per function, and filter out 'Unassigned'
bv2.sig.ko.plot<- bv2.sig.ko %>% 
  group_by(pathway) %>% 
  summarise(median.effect = median(effect),counts = n()) %>% 
  mutate(bar.cols = case_when(sign(median.effect) == 1 ~ "goldenrod4",
                              sign(median.effect) == -1 ~ "goldenrod2"))

# make title column
bv2.sig.ko.plot$title <- "Differentially abundant KEGG functions (vs. BV2)"

# plot median effect sizes for all functions other than 'Unassigned', with
# functions re-ordered by decreasing effect size. Numbers in middle of bars
# represent number of significantly different KO terms per function

# png(paste(path.to.github, "figs_for_paper/lon_eur_BVSubgroups_diffAbund_effectSize.png",
#           sep = ""),
#     units = "in", width = 10, height = 8, res = 400)

bv2.sig.ko.plot%>%
  ggplot(aes(x=median.effect, y=fct_reorder(pathway,median.effect)))+
  geom_col(col="black", linewidth= 0.25, fill= bv2.sig.ko.plot$bar.cols)+
  geom_hline(yintercept = 18.5, linewidth=0.5, color= "grey80", linetype= 2)+
  geom_text(aes(-1.60,19.5),label= "Up in BV2", color= "goldenrod4")+
  geom_text(aes(-1.60,17.5),label= "Up in BV1", color= "goldenrod2")+
  geom_text(aes(0.5*median.effect, fct_reorder(pathway,median.effect)),
            label=bv2.sig.ko.plot$counts,
            col="white", fontface="bold")+
  geom_text(aes(-0.19*sign(median.effect),
                fct_reorder(pathway,median.effect)),
            label= format(signif(bv2.sig.ko.plot$median.effect, digits = 3),
                          nsmall = 3),
            size=3.5, color= bv2.sig.ko.plot$bar.cols)+
  facet_grid(.~title)+
  xlab("ALDEx2 effect size")+ ylab("")+theme_bw()+
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(colour = "black",size = 11),
        axis.text.x = element_text(colour = "black",size = 12),
        axis.title.x = element_text(size = 12),
        strip.text = element_text(size = 12))

# dev.off()

################## BV subgroups: Z-scoring & pathway heatmap ##################

# subset summary ALDEx2 dataframe for significant K0s and remove 'rab.sample.'
# from sample names
bv2.clr.ko <- ko.both.2bv.clr.all[rownames(bv2.sig.ko),
                                  grep("rab.sample.", colnames(ko.both.2bv.clr.all))]

colnames(bv2.clr.ko) <- gsub("rab.sample.", "", colnames(bv2.clr.ko))

# make matrix for aggregating differential K0s by pathways
bv2.clr.path <- matrix(data = NA, 
                       nrow = nrow(bv2.sig.ko.plot),ncol = ncol(bv2.clr.ko),
                       dimnames = list(row = bv2.sig.ko.plot$pathway,
                                       col = colnames(bv2.clr.ko)))

# fill in matrix with pathway means
for(i in rownames(bv2.clr.path)){
  bv2.clr.path[i, 1:ncol(bv2.clr.path)] <- colMeans(bv2.clr.ko[which(bv2.sig.ko$pathway == i),])
}

# create function for calculating z-scores
zscore <- function(x){
  (x - mean(x)) / sd(x)
}

# convert to z-score
bv2.clr.path.z <- apply(bv2.clr.path, 2, zscore)

# make list for column colour bars
lon.eur.bv.column.cols <- list(Dataset = c(London = "snow",
                                           Europe = "grey10"),
                               Subgroup = c(BV1 = "goldenrod2",
                                            BV2 = "goldenrod4"))

# get min and max z scores to centre scale around 0
max(bv2.clr.path.z) # 2.1
min(bv2.clr.path.z) # -3.03

# plot heatmap
# png(paste(path.to.github,
#           "figs_for_paper/lon_eur_BVSubgroups_diffAbund_K0.png",sep = ""),
#     units = "in", height = 6.25, width = 8.5, res = 400)

pheatmap(bv2.clr.path.z, cutree_cols = 2, cutree_rows = 4,
         treeheight_row = 0, treeheight_col = 25, breaks = seq(-3,3,0.06),
         show_colnames = FALSE, annotation_col = lon.eur.bv.meta,
         annotation_colors = lon.eur.bv.column.cols, fontsize_row = 9.5,
         cellheight = 9, cellwidth = 20, border_color = rgb(0,0,0,0.2))

# dev.off()

# ========== FOR INVESTIGATING PATHWAYS OF INTEREST (START) ==========

# filter master VIRGO vNum -> KO lookup table for signicant KOs between BV1 and
# BV2
bv2.vnum.ko <- KO %>% 
  filter(V2 %in% rownames(bv2.sig.ko))

# subset bv2 vNum -> KO table for presence in lon/eur dataset
bv2.vnum.ko <- bv2.vnum.ko[rownames(bv2.vnum.ko) %in% rownames(new.both.data),]

# subset vNum -> KO table for several genes of interest
bv2.vnum.expolsynth <- bv2.vnum.ko %>% 
  filter(V2 == "K11936")
bv2.vnum.chemotax <- bv2.vnum.ko %>% 
  filter(V2 %in% c("K00575","K03406","K03407","K03408","K03410","K03411","K03412","K03415"))
bv2.vnum.cofactor <- bv2.vnum.ko %>% 
  filter(V2 %in% c("K01719","K00798","K01698","K02230"))
bv2.vnum.flagella <- bv2.vnum.ko %>% 
  filter(V2 %in% c("K02390", "K02392", "K02395", "K02396", "K02397", "K02398",
                   "K02400", "K02406", "K02407", "K02409", "K02410", "K02411",
                   "K02412", "K02414", "K02416", "K02417", "K02422", "K02556",
                   "K06603", "K13626"))

# bind dataframes of interest together and add pathway & taxonomy info
bv2.vnum.interest <- rbind(bv2.vnum.expolsynth, bv2.vnum.chemotax,
                           bv2.vnum.cofactor, bv2.vnum.flagella)

bv2.vnum.interest$pathway <- rep(c("Exopolysaccharide biosynthesis",
                                   "Bacterial chemotaxis",
                                   "Porphyrin metabolism",
                                   "Flagellar asssembly"),
                                 c(nrow(bv2.vnum.expolsynth),
                                   nrow(bv2.vnum.chemotax),
                                   nrow(bv2.vnum.cofactor),
                                   nrow(bv2.vnum.flagella)))

bv2.vnum.interest$taxonomy <- tax.table[rownames(bv2.vnum.interest),2]

# subset entire lon/eur dataset by vNumbers of interest and BV1/2 samples
bv2.vnum.interest.df <- new.both.data[rownames(bv2.vnum.interest),
                                      colnames(bv2.clr.ko)]

# pull row sums for each vNumber
bv2.vnum.interest$sum.reads <- rowSums(bv2.vnum.interest.df)

# relocate column for KO to the right of pathway
bv2.vnum.interest <- bv2.vnum.interest %>% 
  relocate(V2, .after = pathway)

# check total reads for each vNumber, grouped by pathway and taxonomy
bv2.vnum.interest.sum <- bv2.vnum.interest %>% 
  group_by(pathway, taxonomy) %>% 
  summarise(total.reads = sum(sum.reads))

# considering reads with known taxonomy, most of the reads assigned to
# exopolysaccharide biosynthesis, flagellar assembly and chemotaxis come from
# BVAB1, while porphyrin metabolism reads mostly contributed by prevotella bivia
# (for the KO up in BV1) and timonensis (for the KOs up in BV2).

# considering all reads, genes whose taxonomy is unknown contribute the largest
# proportion of reads to all categories of interest.

# ========== FOR INVESTIGATING PATHWAYS OF INTEREST (END) ==========

#################### BV subgroups: biplot - significant K0s ####################

# perform PCA on subset of CLR-transformed feature table containing only 
# significantly different KO terms
pca.ko.sig <- prcomp(t(bv2.clr.ko))

# make list of significantly different functions (n = 23) containing indices
# of the significantly different KO terms (n = 66)
ind.sig.ko <- list()
for(func in levels(factor(bv2.sig.ko$pathway))){
  ind.sig.ko[[func]] <- which(bv2.sig.ko$pathway == func)
}

# make vector of function names sorted by their mean effect size (as in
# 'bv.dif.eff.path.plot'), in ascending order
sort.path <- bv2.sig.ko.plot %>% 
  arrange(median.effect)
sort.path <- sort.path$pathway

# sort list of function indices by sorted vector of function names
ind.sig.ko <- ind.sig.ko[sort.path]

# create colour vector for loadings (use two different viridis palettes)
pca.load.cols <- turbo(39)

# make list of grouping indices
ind.grp <- list()
for(group in c("BV1", "BV2")){
  ind.grp[[group]] <- which(lon.eur.bv.groups$group == group)
}

# make the biplot with codaSeq.PCAplot()
# NOTE: trying to add density plots will cause ggplot2's geom_density function
#       to throw a wobbler as it won't plot densities of groups with <2 data
#       points and will instead drop those groups from the plot (here, that
#       means curves for 'Flagellar assembly' and 'Bacterial chemotaxis' are
#       drawn)

# png(paste(path.to.github,
#           "figs_for_paper/lon_eur_BVSubgroups_biplot_sigK0.png", sep = ""),
#     units = "in", height = 9, width = 12, res = 400)

codaSeq.PCAplot(pcx = pca.ko.sig , plot.groups = TRUE, plot.loadings = TRUE,
                plot.ellipses = NULL, plot.density = NULL, 
                grp = ind.grp, grp.col = c("goldenrod2", "goldenrod4"),
                grp.sym = "text", grp.cex = 1,
                load.grp = ind.sig.ko, load.col =  pca.load.cols, 
                load.sym = c(rep(c(15,17), c(18,19)),21,20) , load.cex = 1,
                PC = c(1,2), title = "PCA - BV subgroups: significantly different KO terms",
                plot.legend = "loadings", leg.position = "bottom",
                leg.cex = 0.7282, leg.columns = 5)

# dev.off()

################## biplots: differentially abundant vNumbers ##################

# identify all vNumbers in 'new.both.filt' which correspond to the significantly
# different KO terms
vnum.sig.ko <- KO %>%                       # pull ALL vnumbers corresponding
  filter(V2 %in% rownames(bv2.sig.ko))      # to significantly different KO terms

vnum.sig.ko.filt <- vnum.sig.ko %>%                           # filter to just
  filter(rownames(vnum.sig.ko) %in% rownames(new.both.filt))  # vnums in dataset

# subset 'new.both.filt' by the vNumbers corresponding to the significantly
# different KO terms (n = 285)
new.both.filt.sig <- new.both.filt[rownames(vnum.sig.ko.filt),] # rows
new.both.filt.sig <- new.both.filt.sig[,colnames(bv2.clr.ko)] # cols

# check for and correct columns/rows with zero sums
any(colSums(new.both.filt.sig) == 0)  # returns FALSE
any(rowSums(new.both.filt.sig) == 0)  # returns TRUE
new.both.filt.sig <- new.both.filt.sig[-which(rowSums(new.both.filt.sig) == 0),]

# set seed for RNG
set.seed(2023)

# generate matrix of scale values with mu = 1 and 1.15, std dev of log normal
# distribution = 0.5
mu.matrix <- aldex.makeScaleMatrix(gamma = 0.5, mu = c(1, 1.10),
                                lon.eur.bv.groups$group)

# CLR transform data and calculate median CLR values for all features (vNum)

vnum.sig.ko.clr <- aldex.clr(new.both.filt.sig,
                             conds = lon.eur.bv.groups$group,
                             gamma = mu.matrix, verbose = TRUE)

vnum.sig.ko.clr.e <- aldex.effect(vnum.sig.ko.clr, verbose = TRUE,
                                  include.sample.summary = TRUE)

# extract log-ratio table and perform PCA
vnum.sig.ko.clr.df <- vnum.sig.ko.clr.e[, grep('rab.sample.',
                                               colnames(vnum.sig.ko.clr.e))]
colnames(vnum.sig.ko.clr.df) <- gsub("rab.sample.", "", colnames(vnum.sig.ko.clr.df))
pca.vnum.sig <- prcomp(t(vnum.sig.ko.clr.df))

# get named vector of taxa present in the vNumber subset table
vnum.taxa <- tax.table[rownames(vnum.sig.ko.clr.df),2]
names(vnum.taxa) <- rownames(vnum.sig.ko.clr.df)

# make list of vNumber indices for each species
ind.sig.vnum <- list()
for(i in levels(factor(vnum.taxa))){
  ind.sig.vnum[[i]] <- which(vnum.taxa == i)
}

# get colours corresponding to taxa in the list of vNumbers corresponding to 
# significantly different KO terms, plus transparent grey for 'Other'
pca.load.cols.vnum <- c(tax.colors[names(ind.sig.vnum),], rgb(0,0,0,0.1))

# add colours for taxa that weren't included in 'tax.colors'
pca.load.cols.vnum[8] <- "blue3"     # Lactobacillus johnsonii
pca.load.cols.vnum[11] <- "seagreen4"    # Mycoplasma hominis
pca.load.cols.vnum[12] <- "chocolate4"    # Peptoniphilus harei

#add vNumbers with unknown taxa to 'Other'
ind.sig.vnum[["Other"]] <- setdiff(1:nrow(new.both.filt.sig),
                                   unlist(ind.sig.vnum))

# make species biplot of vNumbers corresponding to significantly different KO 
# terms

# png(paste(path.to.github,
#           "figs_for_paper/lon_eur_BVSubgroups_biplot_sigvNum.png", sep = ""),
#     units = "in", height = 9, width = 12, res = 400)

codaSeq.PCAplot(pca.vnum.sig, plot.groups = FALSE, plot.loadings = TRUE,
                plot.ellipses = NULL, plot.density = NULL,
                load.grp = ind.sig.vnum, load.col = pca.load.cols.vnum,
                load.sym = rep(19,22), load.cex = 1, PC = c(1,2),
                plot.legend = "loadings", leg.position = "topleft", 
                leg.cex = 0.85, leg.columns = 2,
                title = "PCA - BV subgroups: significantly different vNumbers")

# dev.off()

######################### biplots: procrustes rotation #########################

# A procrustes rotation will take  two matrices or ordination objects, overlay 
# them (i.e. superimpose the two plots) and rotate the second such that the 
# similarity between them is maximised. Mathematically, the sum of the squared
# deviations (denoted m2) is minimised between samples. The deviation refers to 
# the difference between the location of a sample on the first PCA plot and the
# location of the same sample on the rotated PCA plot, and is also called the 
# vector residual. A small vector residual indicates a close agreement between
# the location of a sample on the two plots.

# perform procrustes rotation using the KO term PCA as the 'target', rotating 
# the vNumber PCA (first two axes only)
procrust.vnum.ko <- procrustes(X = pca.ko.sig, Y = pca.vnum.sig,
                               scores = "sites", choices= c(1,2),
                               scale = FALSE,  symmetric = FALSE)

# extract PCA co-ordinates from KO term prcomp object and add BV groupings &
# colours
pca.ko.sig.vals <- as.data.frame(pca.ko.sig$x[,c(1,2)])
pca.ko.sig.vals$group <- lon.eur.bv.groups$group
pca.ko.sig.vals$colour <- case_when(pca.ko.sig.vals$group == "BV1" ~ "goldenrod2",
                                    pca.ko.sig.vals$group == "BV2" ~ "goldenrod4")

# plot the procrustes rotation (target = K0 PCA, rotated = vNumber PCA) and
# save as .png
# png(paste(path.to.github,
#           "figs_for_paper/lon_eur_BVSubgroups_biplot_procrustes.png", sep = ""),
#     units = "in", height = 9, width = 12, res = 400)

plot(procrust.vnum.ko, kind=1, type= "text", length=0,
     xlab= "PC1", ylab= "PC2", cex = 0.8, ar.col= "dodgerblue3",
     main= "Procrustes rotation: target = K0 PCA, rotated = vNumber PCA")

text(pca.ko.sig.vals$PC1, pca.ko.sig.vals$PC2, cex = 0.8,
     labels = rownames(pca.ko.sig.vals), col = pca.ko.sig.vals$colour)

legend("bottomleft", legend = c("BV1","BV2"), pch = "vv",
       col = c("goldenrod2","goldenrod4"), pt.cex = 1.5)

# dev.off()
