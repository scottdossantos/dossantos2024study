# assessing effect of removing BVAB1 genes on differential expression results

################################################################################
#################################### setup ####################################
################################################################################

library(ALDEx2) # differential abundance analysis
library(dplyr) # data manipulation
library(viridisLite) # colourblind-friendly palettes
library(pheatmap) # nicer heatmaps
library(CoDaSeq) # biplots

# set working directory
setwd("~/Documents/GitHub/dossantos2024study/")
wd <- paste0(getwd(), "/")

# load in following .Rda objects:
#    - Filtered virginia feature table (gene-level)
#    - Filtered virginia taxa vector
#    - Filtered virginia KO -> pathway lookup table (edited)
#    - Grouping vector for virginia samples (health vs. BV)
#    - Grouping vector for virginia BV samples (BV1 vs. BV2)
#    - CST vector for virginia samples
#    - VIRGO vNumber -> KO table
#    - Taxa colour vector
#    - Renamed virginia taxa vector (G vaginalis edits)
load(paste0(wd,"Rdata/virginia.filt.Rda"))
load(paste0(wd,"Rdata/virginia.tax.vec.Rda"))
load(paste0(wd,"Rdata/virginia.ko.path.filt.Rda"))
load(paste0(wd,"Rdata/virginia.groups.Rda"))
load(paste0(wd,"Rdata/virginia.bv.groups.Rda"))
load(paste0(wd,"Rdata/virginia.cst.Rda"))
vnum.to.ko <- read.table(paste0(wd,"1_VIRGO/8.A.kegg.ortholog.txt"), header = F,
                         sep = "\t", quote = '', row.names = 1, 
                         col.names = c("vnum", "ko", "code", "function"))
taxa.cols <- read.table(paste0(wd, "Rdata/species_colors_updated.txt"), header = T,
                        sep = "\t", quote = '', row.names = 1)

# make vector of KO numbers corresponding to vNumbers in filtered virginia data
# (or load final object if it's already been made)
ko.vec.virginia <- vnum.to.ko[rownames(virginia.filt),1]
names(ko.vec.virginia) <- rownames(virginia.filt) 

# remove all vNumbers corresponding to BVAB1 from feature table and ko vector
noBVAB1.filt <- virginia.filt[-c(which(tax.vec.virginia == "BVAB1")),]
noBVAB1.ko.vec <- ko.vec.virginia[-which(tax.vec.virginia == "BVAB1")]

# aggregate filtered virginia feature table by KO term, remove suspected
# eukaryotic KO terms, and filter ko -> pathway object by no BVAB1 feature table
noBVAB1.filt.ko <- aggregate(noBVAB1.filt, by = list(noBVAB1.ko.vec), FUN = sum)
rownames(noBVAB1.filt.ko) <- noBVAB1.filt.ko$Group.1
noBVAB1.filt.ko$Group.1 <- NULL
noBVAB1.filt.ko <- noBVAB1.filt.ko[-which(grepl(paste("K03260","K06237","K12373",
                                                      "K00863","K13993","K03648",
                                                      "K01408","K00599",sep = "|"),
                                                  rownames(noBVAB1.filt.ko))),]
noBVAB1.ko.path <- virginia.ko.path.filt[rownames(noBVAB1.filt.ko), , drop = F]

# filter feature table for only BV samples (defined by clustering on KO terms)
noBVAB1.filt.ko <- noBVAB1.filt.ko[,virginia.groups$groups.2=="BV"]

# sanity check for 0 sum rows/cols
any(rowSums(noBVAB1.filt.ko)==0)  # returns FALSE
any(colSums(noBVAB1.filt.ko)==0)  # returns FALSE

################################################################################
######################### virginia: no BVAB1 - ALDEx2 #########################
################################################################################

# filter feature table for only samples in original BV1 vs. BV2 groups
bv12.ko <- noBVAB1.filt.ko[,rownames(virginia.bv.groups)]
any(rowSums(bv12.ko) == 0)  # returns FALSE

# set seed for RNG, make scale matrix, then CLR transform with ALDEx2 inc. scale
# and get effect sizes/sample summary and t-test results (or load final object
# if it's already been made)
set.seed(2024)
  
bv12.mu <- aldex.makeScaleMatrix(gamma = 0.5, mu = c(1,1.15), mc.samples = 128,
                                 conditions = virginia.bv.groups$group)
  
bv12.clr <- aldex.clr(reads = bv12.ko, conds = virginia.bv.groups$group,
                      mc.samples = 128, denom = "all", gamma = bv12.mu,
                      verbose = TRUE)
  
bv12.clr.e <- aldex.effect(bv12.clr, verbose = TRUE,
                           include.sample.summary = TRUE)
  
bv12.clr.t <- aldex.ttest(bv12.clr, verbose = TRUE)
  
bv12.clr.all <- cbind(bv12.clr.e, bv12.clr.t)


# pull pathways with absolute effect size >1 and P <0.01
bv12.sig.ko <- noBVAB1.ko.path[which(abs(bv12.clr.all$effect) >1 & bv12.clr.all$we.eBH <0.01), , drop = F]
bv12.sig.ko$effect <- bv12.clr.all$effect[which(abs(bv12.clr.all$effect) >1 & bv12.clr.all$we.eBH <0.01)]

# edit pathway for K03309 from 'Unknown' to 'Dicarboxylate transport' (see also
# https://pubmed.ncbi.nlm.nih.gov/19581365/)
bv12.sig.ko["K03309", "pathway"] <- "Dicarboxylate transport"

# sum reads in BV1 and BV2 samples (these differences are very real)
bv12.sig.ko$reads.bv1 <- NA
bv12.sig.ko$reads.bv2 <- NA
for(i in rownames(bv12.sig.ko)){
  bv12.sig.ko[i,"reads.bv1"] <- rowSums(bv12.ko[i,which(virginia.bv.groups$group=="BV Subgroup 1")])
  bv12.sig.ko[i,"reads.bv2"] <- rowSums(bv12.ko[i,which(virginia.bv.groups$group=="BV Subgroup 2")])
}

# pull median CLR values from aldex results
bv12.clr.df <- bv12.clr.all[,grep('rab.sample.', colnames(bv12.clr.all))]
colnames(bv12.clr.df) <- gsub('rab.sample.', '', colnames(bv12.clr.df))

# make matrix for significantly different pathways and add row for biofilm form
bv12.path.df <- as.data.frame(matrix(data = NA, ncol = ncol(bv12.clr.df),
                                     nrow = length(levels(factor(bv12.sig.ko$pathway)))+1,
                                     dimnames = list(rows = c(levels(factor(bv12.sig.ko$pathway)), "Biofilm formation"),
                                                     cols = colnames(bv12.clr.df))))


# fill in data frame with mean CLR values per pathway
for(i in levels(factor(bv12.sig.ko$pathway))){
  bv12.path.df[i,] <- colMeans(bv12.clr.df[rownames(bv12.sig.ko)[which(bv12.sig.ko$pathway==i)],])
}
bv12.path.df["Biofilm formation",] <- colMeans(bv12.clr.df["K11936",])

# make z-score function and calculate zscores
zscore <- function(x){
  (x -mean(x)) / sd(x)
}

bv12.path.z <- t(as.data.frame(apply(bv12.path.df, 1, zscore)))

# get metadata objects and make list object for heatmap
bv12.meta <- data.frame(virginia.cst[colnames(bv12.path.z), , drop = F],
                        virginia.bv.groups$group)
colnames(bv12.meta) <- c("CST", "Subgroup")
bv12.meta$Subgroup <- gsub(" Subgroup ", "", bv12.meta$Subgroup)

bv12.hm.cols <- list(CST= c(III = viridis(7)[3],
                            `III / IV` = viridis(7)[4],
                            IV = viridis(7)[6]),
                     Subgroup = c(BV1 = "goldenrod2",
                                  BV2 = "goldenrod4"))

# save(bv12.meta, file = paste0(wd,"Rdata/bv12.meta.Rda"))
# save(bv12.hm.cols, file = paste0(wd,"Rdata/bv12.hm.cols.Rda"))

# make heatmap 
# png("virginia_noBVAB1_bvSubgroups_difAbund.png",
#     units = "in", height = 2.5, width = 11, res = 400)

pheatmap(mat = bv12.path.z, clustering_method = "ward.D2", cutree_cols = 2,
         annotation_colors = bv12.hm.cols, annotation_col = bv12.meta,
         show_colnames = F, fontsize_row = 10.5, border_color = rgb(0,0,0,0.1),
         cellwidth = 5, cellheight = 10.5, treeheight_row = 0, breaks=seq(-2.5,2.5,0.05))

# dev.off()

################################################################################
########################### virginia: no BVAB1 - PCA ###########################
################################################################################

# subset taxa vector for no bvab1 
tax.vec.nobvab1 <- tax.vec.virginia[rownames(noBVAB1.filt)]

# count number of genes for each species and pull only those w/ >75
n.genes <- vector()
for(i in levels(factor(tax.vec.nobvab1))){
  n.genes[i] <- length(which(tax.vec.nobvab1 == i))
}

n.genes.75 <- n.genes[which(n.genes >75)]

# add taxa colours for f. magna, l. coleohominis and p. harei
taxa.cols[25:27,] <- c("yellow3", "blue4", "chocolate4")
rownames(taxa.cols)[25:27] <- c("Finegoldia_magna",
                                "Lactobacillus_coleohominis",
                                "Peptoniphilus_harei")

# make index list and colour vector for species and add 'Other'
ind.spp <- list()
for(i in names(n.genes.75)){
  ind.spp[[i]] <- which(tax.vec.nobvab1 == i)
}
ind.spp[["Other"]] <- setdiff(1:length(tax.vec.nobvab1),
                              unlist(ind.spp))
cols.spp <- c(taxa.cols[names(n.genes.75),], rgb(0,0,0,0.1))

# make index list and colour vector for groups
ind.bv <- list(BV1 = which(virginia.bv.groups=="BV Subgroup 1"),
               BV2 = which(virginia.bv.groups=="BV Subgroup 2"))

cols.grp <- c("goldenrod1", "black")

# filter vnumber feature table by only BV species, CLR transform and pull
# median CLR values
noBVAB1.filt.bv <- noBVAB1.filt[, rownames(virginia.bv.groups)]
any(rowSums(noBVAB1.filt.bv) == 0) # returns FALSE

set.seed(2024)
noBVAB1.filt.bv.clr <- aldex.clr(reads = noBVAB1.filt.bv, mc.samples = 128,
                                 conds = virginia.bv.groups$group, denom = "all",
                                 gamma = bv12.mu, verbose = TRUE)
  
noBVAB1.filt.bv.clr.e <- aldex.effect(clr = noBVAB1.filt.bv.clr, verbose = TRUE,
                                      include.sample.summary = TRUE)

# extract median CLR values and edit sample names to remove 'rab.sample.'
noBVAB1.filt.bv.clr.df <- noBVAB1.filt.bv.clr.e[,grep('rab.sample.',
                                                      colnames(noBVAB1.filt.bv.clr.e))]

colnames(noBVAB1.filt.bv.clr.df) <- gsub("rab.sample.", "",
                                         colnames(noBVAB1.filt.bv.clr.df))

# pca on vnumber CLR values
pca.noBVAB1.vnum <- prcomp(t(noBVAB1.filt.bv.clr.df))

# plot biplot - density = loadings
# png("virginia_noBVAB1_vnum_load.png",
#     units = "in", height = 6, width = 12, res = 400)

codaSeq.PCAplot(pcx = pca.noBVAB1.vnum, plot.groups = T, plot.loadings = T,
                plot.density = "loadings", grp = ind.bv, grp.col = cols.grp,
                grp.cex = 0.4, load.grp = ind.spp, load.col = cols.spp,
                load.sym = rep(19,17), load.cex = 0.2, PC = c(1,2),
                leg.position = "bottomleft", plot.legend = "groups",
                leg.cex = 0.75, title = "PCA: BV subgroups (no BVAB1 genes)")
legend(x = 0.805, y = 1.065, legend = gsub("_", " ", names(ind.spp)),
       col = cols.spp, pch = 19, pt.cex = 0.4, cex = 0.6, ncol = 2, box.lwd = 0)

# dev.off()

# plot biplot - density = groups
# png("virginia_noBVAB1_vnum_group.png",
#     units = "in", height = 6, width = 12, res = 400)

codaSeq.PCAplot(pcx = pca.noBVAB1.vnum, plot.groups = T, plot.loadings = T,
                plot.density = "groups", grp = ind.bv, grp.col = cols.grp,
                grp.cex = 0.4, load.grp = ind.spp, load.col = cols.spp,
                load.sym = rep(19,17), load.cex = 0.2, PC = c(1,2),
                leg.position = "bottomleft", plot.legend = "groups",
                leg.cex = 0.75, title = "PCA: BV subgroups (no BVAB1 genes)")
legend(x = 0.805, y = 1.065, legend = gsub("_", " ", names(ind.spp)),
       col = cols.spp, pch = 19, pt.cex = 0.4, cex = 0.6, ncol = 2, box.lwd = 0)

# dev.off()
