# DESeq2 analysis of Lon/Eur data at KO level without SRI

# load deseq2
library(DESeq2)
library(dplyr)

# gh path(adjust if needed)
path.to.github <- "~/Documents/GitHub/dossantos2024study/"


#################### DESeq2 analysis of London/Europe data ####################
# load in KO-aggregated counts table and KO -> pathway lookup table
load(paste0(path.to.github, "Rdata/ko.both.Rda"))
load(paste0(path.to.github, "Rdata/ko.both.path.Rda"))

# remove two samples from the filtered london/europe feature table aggregated by
# K0 number (classed as BV but almost no BV organisms):
#   -  v.001A: close to 100% L. gasseri with practically no BV organisms
#   -  v.019A: around 80 % iners, ~ 20% crispatus and a tiny bit of Gardnerella
ko.both <- ko.both[,-c(9,22)]

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

# load aldex2 results for lon/eur health vs. bv and pull significant KOs
load(paste0(path.to.github,"Rdata/ko.both.clr.all.Rda"))
diff.sig.ko <- ko.both.path[which(abs(ko.both.clr.all$effect) >1), , drop=FALSE]
diff.sig.ko$effect <- ko.both.clr.all$effect[which(abs(ko.both.clr.all$effect) >1)]
diff.sig.ko <- diff.sig.ko %>% 
  filter(pathway != 'Unknown')

# make colData (42 rows corresponding to samples, columns = H vs. BV and batch)
conds <- data.frame(status = factor(c(rep('H',8), rep('B',12), rep('B',14), rep('H', 8))),
                    batch = factor(c(rep("London", 20), rep("Europe", 22))),
                    row.names = colnames(ko.both))

# make deseqdataset
dds.le.hbv <- DESeqDataSetFromMatrix(countData = ko.both,
                                     colData = conds,
                                     design = ~ batch + status)

# run deseq2
dds.le.hbv <- DESeq(dds.le.hbv)

# get results names and results for H vs. BV
resultsNames(dds.le.hbv)
res.le.hbv <- results(dds.le.hbv, name="status_H_vs_B", alpha = 0.05)

# plot MA
plotMA(res.le.hbv)

# get results data and filter for P adj <0.01
res.le.hbv.data <- data.frame(res.le.hbv)
res.le.hbv.data <- res.le.hbv.data[which(res.le.hbv.data$padj <0.01),]

# check which aldex2 significant KOs in deseq2 significant KOs
rownames(diff.sig.ko) %in% rownames(res.le.hbv.data)

#################### DESeq2 analysis of Virginia data ####################

# load in KO-aggregated counts table and KO -> pathway lookup table
load(paste0(path.to.github, "Rdata/ko.virginia.filt.Rda"))
load(paste0(path.to.github, "Rdata/virginia.ko.path.filt.Rda"))
load(paste0(path.to.github, "Rdata/virginia.groups.Rda"))

# aldex results
load(paste(path.to.github, "Rdata/virginia.filt.ko.clr.c.Rda", sep = ""))

# subset the virginia K0 -> pathway table based on 567 K0 terms with absolute
# effect size >1 (note: 'drop = FALSE' when subsetting retains row names)
diff.virginia.ko.path <- virginia.ko.path.filt[which(abs(virginia.filt.ko.clr.c$effect)>1), , drop=F]

# add effect sizes to this subset
diff.virginia.ko.path$effect <- virginia.filt.ko.clr.c[which(abs(virginia.filt.ko.clr.c$effect)>1),"effect"]

# remove K0 terms with "Unknown" pathway
diff.virginia.ko.path <- diff.virginia.ko.path %>% 
  filter(pathway != "Unknown")

# conds for virginia ko
conds.virg <- data.frame(status=virginia.groups$groups.2)
rownames(conds.virg) <- colnames(ko.virginia.filt)

# make deseqdataset
dds.virg.hbv <- DESeqDataSetFromMatrix(countData = ko.virginia.filt,
                                       colData = conds.virg,
                                       design = ~ status)

# run deseq2
dds.virg.hbv <- DESeq(dds.virg.hbv)

# get results names and results for H vs. BV
resultsNames(dds.virg.hbv)
res.virg.hbv <- results(dds.virg.hbv, name="status_Healthy_vs_BV", alpha = 0.05)

# plot MA
plotMA(res.virg.hbv)

# get results data and filter for P adj <0.01
res.virg.hbv.data <- data.frame(res.virg.hbv)
res.virg.hbv.data <- res.virg.hbv.data[which(res.virg.hbv.data$padj <0.01),]

# check which aldex2 significant KOs in deseq2 significant KOs
rownames(diff.virginia.ko.path) %in% rownames(res.virg.hbv.data)

################### DESeq2 analysis of lon/eur BV subgroups ###################

# count data
load(paste0(path.to.github, "Rdata/bvSubgroup_names.Rda"))
loneur.bvsubg <- ko.both[,rownames(lon.eur.bv.groups)]

# aldex results
load(paste0(path.to.github, "Rdata/ko.both.2bv.clr.all.Rda"))

# pull all K0s with absolute effect size > 1 and the corresponding effect sizes
bv2.sig.ko <- ko.both.path[which(abs(ko.both.2bv.clr.all$effect) >1), ,
                           drop = FALSE]
bv2.sig.ko$effect <- ko.both.2bv.clr.all[rownames(bv2.sig.ko), 'effect']

# filter out K0s with unknown pathway
bv2.sig.ko <- bv2.sig.ko %>% 
  filter(pathway != "Unknown")

# make deseqdataset
dds.le.bvsub <- DESeqDataSetFromMatrix(countData = loneur.bvsubg,
                                       colData = lon.eur.bv.groups,
                                       design = ~ group)

# run deseq2
dds.le.bvsub <- DESeq(dds.le.bvsub)

# get results names and results for H vs. BV
resultsNames(dds.le.bvsub)
res.le.bvsub <- results(dds.le.bvsub, name="group_BV2_vs_BV1", alpha = 0.05)

# plot MA
plotMA(res.le.bvsub)

# get results data and filter for P adj <0.01
res.le.bvsub.data <- data.frame(res.le.bvsub)
res.le.bvsub.data <- res.le.bvsub.data[which(res.le.bvsub.data$padj <0.01),]

rownames(bv2.sig.ko) %in% rownames(res.le.bvsub.data)

################### DESeq2 analysis of Virginia BV subgroups ###################

# count data
load(paste0(path.to.github, "Rdata/ko.virginia.filt.Rda"))
load(paste0(path.to.github, "Rdata/virginia.bv.groups.Rda"))
virg.bvsubg <- ko.virginia.filt[,rownames(virginia.bv.groups)] 

# aldex results
load(paste0(path.to.github, "Rdata/sub.virginia.ko.clr.all.Rda"))

# get pathways associated with K0s showing absolute effect size >1
sub.virginia.sig.ko <- virginia.ko.path.filt[which(abs(sub.virginia.ko.clr.all$effect) >1), , drop = FALSE]

# ABC transporters, bacterial chemotaxis, flagellar assembly, porphyrin
# biosynthesis and two-component systems are all on there!

# pull effect sizes for these differential K0s (should all be negative)
sub.virginia.sig.ko$effect <- sub.virginia.ko.clr.all[rownames(sub.virginia.sig.ko),"effect"]

# remove K0 terms with 'Unknown' pathway
sub.virginia.sig.ko <- sub.virginia.sig.ko %>% 
  filter(pathway != "Unknown")


# make deseqdataset
dds.virg.bvsub <- DESeqDataSetFromMatrix(countData = virg.bvsubg,
                                         colData = virginia.bv.groups,
                                         design = ~ group)

# run deseq2
dds.virg.bvsub <- DESeq(dds.virg.bvsub)

# get results names and results for H vs. BV
resultsNames(dds.virg.bvsub)
res.virg.bvsub <- results(dds.virg.bvsub, name="group_BV.Subgroup.2_vs_BV.Subgroup.1", alpha = 0.05)

# plot MA
plotMA(res.virg.bvsub)

# get results data and filter for P adj <0.01
res.virg.bvsub.data <- data.frame(res.virg.bvsub)
res.virg.bvsub.data <- res.virg.bvsub.data[which(res.virg.bvsub.data$padj <0.01),]

rownames(sub.virginia.sig.ko) %in% rownames(res.virg.bvsub.data)

#################### DESeq2 analysis figures for reviewers ####################

# png("~/deseq2_diffAbund_healthBV.png", units = "in",
#     height = 4, width = 8, res = 400)

par(mfrow = c(1,2))

plotMA(res.le.hbv)
title(main = "DESeq2 - Health vs. BV (Lon/Eur)")

plotMA(res.virg.hbv)
title(main = "DESeq2 - Health vs. BV (Virginia)")

# dev.off()
