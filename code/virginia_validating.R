# validating previous observations in Virginia dataset

#################################### setup ####################################

library(dplyr) # data manipulation
library(viridisLite) # colourblind-friendly palettes
library(pheatmap) # nicer heatmaps

user <- "sds/cc"
path.to.github <- switch(user,
                         "sds/cc" = "~/Documents/GitHub/dossantos2024study/",
                         "gg" = "~/Documents/0_git/projects/dossantos2024study/")

# run 'setup.R'
source(paste(path.to.github,"code/setup.R", sep = ""))

# load in following .Rda objects:
#    - filtered and unfiltered virginia feature tables (297 samples)
#    - taxa vector corresponding to the filtered virginia feature table
#    - clinical and demographic metadata for 297 virginia samples (omitted from GH)
#    - hclust output for 297 virginia samples clustered by K0 numbers
#    - cutree output for 297 virginia samples clustered by K0 numbers
#    - grouping variables (healthy subgroups 1 & 2 vs. molecular BV )
#    - filtered virginia feature table aggregated by K0 number
#    - K0 -> pathway table for filtered virginia dataset (manually edited)
#    - summary data frame from running ALDEx2 w/ scale sim on virginia K0 table 
load(paste(path.to.github, "Rdata/virginia.data.Rda", sep = ""))
load(paste(path.to.github, "Rdata/virginia.filt.Rda", sep = ""))
load(paste(path.to.github, "Rdata/tax.vec.virginia.Rda", sep = ""))
load("[LOCAL-PATH-TO]/virginia.meta.Rda") # this object contains dbGaP metadata and is omitted from the github repo
load(paste(path.to.github, "Rdata/virginia.hclust.Rda", sep = ""))
load(paste(path.to.github, "Rdata/virginia.cst.Rda", sep = ""))
load(paste(path.to.github, "Rdata/virginia.groups.Rda", sep = ""))
load(paste(path.to.github, "Rdata/ko.virginia.filt.Rda", sep = ""))
load(paste(path.to.github, "Rdata/virginia.ko.path.filt.Rda", sep = ""))
load(paste(path.to.github, "Rdata/virginia.filt.ko.clr.c.Rda", sep = ""))

############### health vs. molecular BV: differential abundance ################

# subset the virginia K0 -> pathway table based on 567 K0 terms with absolute
# effect size >1 (note: 'drop = FALSE' when subsetting retains row names)
diff.virginia.ko.path <- virginia.ko.path.filt[which(abs(virginia.filt.ko.clr.c$effect)>1), , drop=F]

# add effect sizes to this subset
diff.virginia.ko.path$effect <- virginia.filt.ko.clr.c[which(abs(virginia.filt.ko.clr.c$effect)>1),"effect"]

# remove K0 terms with "Unknown" pathway (down to 325 K0 terms)
diff.virginia.ko.path <- diff.virginia.ko.path %>% 
  filter(pathway != "Unknown")

# many instances where K0 terms in same pathway have effect sizes with opposite
# (i.e. -ve and +ve) effect sizes- need to edit pathway names to distinguish

# append " (H)" or " (BV)" to all pathways only if there is more than one K0 for
# each pathway AND there are positive and negative effect sizes for that pathway
for(i in levels(factor(diff.virginia.ko.path$pathway))){
  if(length(which(diff.virginia.ko.path$pathway == i)) >1 & length(unique(sign(diff.virginia.ko.path$effect[which(diff.virginia.ko.path$pathway == i)]))) != 1){
    diff.virginia.ko.path$pathway <- case_when(diff.virginia.ko.path$pathway == i & sign(diff.virginia.ko.path$effect) == 1 ~ paste(i, "(H)", sep = " "),
                                               diff.virginia.ko.path$pathway == i & sign(diff.virginia.ko.path$effect) == -1 ~ paste(i, "(BV)", sep = " "),
                                               .default = diff.virginia.ko.path$pathway)
  }
}

# check how many pathways there are among the 325 K0s (93 unique pathway names)
diff.virginia.ko.path.grp <- diff.virginia.ko.path %>% 
  group_by(pathway) %>% 
  count()

# =========== INVESTIGATING PREVIOUSLY ID'D KOs OF INTEREST (START) ===========

# load the KEGGREST package to allow access to KEGG database
library(KEGGREST)

# load unfiltered virginia KO-aggregated dataset
load(paste0(path.to.github, "RData/ko.virginia.Rda"))

# pull all KOs associated with the CAMP pathway (map01503) using 'keggLink()'
# and remove "ko:" from all entries
camp.ko <- keggLink("ko","map01503") # 54 KOs
camp.ko <- gsub("ko:", "", camp.ko)

# subset for only the CAMP KOs in the virginia dataset
camp.ko.virginia <- rownames(ko.virginia)[which(rownames(ko.virginia) %in% camp.ko)]

# copy this into a text editor, then copy/paste into KEGG KO to check functions
cat(camp.ko.virginia)

# load taxa vector for unfiltered virginia dataset if it exists, or make it if
# it doesn't
if(file.exists(paste0(path.to.github, "Rdata/tax.vec.virginia.all.Rda"))){
  load(paste0(path.to.github, "Rdata/tax.vec.virginia.all.Rda"))
} else{
  tax.vec.virginia.all <- tax.table[rownames(virginia.data),2] # takes >10 min
  
  save(tax.vec.virginia.all, file = paste0(path.to.github,
                                           "Rdata/tax.vec.virginia.all.Rda")) 
}

# filter vnum -> KO table for KOs in virginia dataset, then vNumbers in 
# virginia dataset
camp.vnum.ko.virg <- KO %>% 
  filter(V2 %in% rownames(ko.virginia))

camp.vnum.ko.virg <- camp.vnum.ko.virg %>% 
  filter(rownames(camp.vnum.ko.virg) %in% rownames(virginia.data))

# check if K01354 (African trypanosomiasis/PtrB - oligopeptidase B) is present
# in virginia dataset
"K01354" %in% rownames(ko.virginia) # returns TRUE

# add PtrB KO to the CAMP KO vector and reorder
camp.ko.virginia <- append(camp.ko.virginia, "K01354", after = length(camp.ko.virginia))
camp.ko.virginia <- sort(camp.ko.virginia)

# for each camp KO in virginia dataset, pull vNumbers, taxonomy info and 
# corresponding read counts across virignia dataset
camp.res.genes.tax.virg <- list()
for(i in camp.ko.virginia){
  df <- camp.vnum.ko.virg %>% 
    filter(V2 == i)
  taxa <- tax.table[rownames(df), 2]
  names(taxa) <- rownames(df)
  camp.ko.sums <- rowSums(virginia.data[rownames(df),])
  camp.res.genes.tax.virg[[i]] <- data.frame(vnum= rownames(df),
                                             ko = rep(i, nrow(df)),
                                             tax.id = taxa,
                                             total.reads =camp.ko.sums,
                                             row.names = NULL)
}

# save as .Rda object
# save(camp.res.genes.tax.virg, file = paste0(path.to.github,
#                                             "Rdata/camp.res.genes.tax.virg.Rda"))

# collapse list of CAMP resistance genes in virginia dataset and their counts
# to dataframe and set rownames to numbers
camp.res.genes.tax.virg.df <- do.call(rbind, camp.res.genes.tax.virg)
rownames(camp.res.genes.tax.virg.df) <- NULL

# group by KO and then sum reads
camp.res.genes.tax.virg.sum <- camp.res.genes.tax.virg.df %>% 
  group_by(ko, tax.id) %>% 
  summarise(sum=sum(total.reads))

# many of the new KOs not seen in Lon/Eur dataset are represented by <250
# reads in total. Threshold for a 'believeable' result looks to be about 1,000
# reads (stuff that makes sense and isn't just 1 vNumber corresponding to
# E. coli)

# now do the same thing but for specific iron acquisition KOs, collagen deg KO
# and urease KOs (2012,2217,7243 = iron, 8303 = collagenase, 3191 = ureI)

# make vector of KO terms to pull (googled the KO terms for UreA / UreB)
ir.co.ur.ko <- c("K02012","K02217","K07243",
                 "K08303",
                 "K03191", "K01429", "K01430")

#check if all KOs are in virginia dataset 
ir.co.ur.ko %in% rownames(ko.virginia) # all TRUE

# filter vnum -> KO table for these KOs, then for vNums in database
ir.co.ur.ko.virg <- KO %>% 
  filter(V2 %in% ir.co.ur.ko)

ir.co.ur.ko.virg <- ir.co.ur.ko.virg %>% 
  filter(rownames(ir.co.ur.ko.virg) %in% rownames(virginia.data))

# make summary list object of vNumbers, their corresponding KOs and counts in
# the virginia dataset
ir.co.ur.tax.virg <- list()
for(i in ir.co.ur.ko){
  df <- ir.co.ur.ko.virg %>% 
    filter(V2 == i)
  taxa <- tax.table[rownames(df), 2]
  names(taxa) <- rownames(df)
  ir.co.ur.ko.sums <- rowSums(virginia.data[rownames(df),])
  ir.co.ur.tax.virg[[i]] <- data.frame(vnum= rownames(df),
                                             ko = rep(i, nrow(df)),
                                             tax.id = taxa,
                                             total.reads =ir.co.ur.ko.sums,
                                             row.names = NULL)
}

# save as .Rda object
# save(ir.co.ur.tax.virg, file = paste0(path.to.github,
#                                             "Rdata/ir.co.ur.tax.virg.Rda"))

# collapse list of genes in virginia dataset and their counts to dataframe and
# set rownames to numbers
ir.co.ur.tax.virg.df <- do.call(rbind, ir.co.ur.tax.virg)
rownames(ir.co.ur.tax.virg.df) <- NULL

# group by KO and then sum reads
ir.co.ur.tax.virg.sum <- ir.co.ur.tax.virg.df %>% 
  group_by(ko,tax.id) %>% 
  summarise(sum=sum(total.reads))

# ============ INVESTIGATING PREVIOUSLY ID'D KOs OF INTEREST (END) ============

# pull out columns corresponding to median CLR values (all samples) and remove
# 'rab.sample.' from column names
diff.virginia.ko.clr <- virginia.filt.ko.clr.c[, grep("rab.sample.",
                                                      colnames(virginia.filt.ko.clr.c))]

colnames(diff.virginia.ko.clr) <- gsub("rab.sample.", "",
                                       colnames(diff.virginia.ko.clr))

# subset the scaled, CLR feature table of K0 terms table to retain only the 
# differential, non-Unknown K0 terms
diff.virginia.ko.clr <- diff.virginia.ko.clr[rownames(diff.virginia.ko.path),]

# create a matrix of NAs of dimensions corresponding to the number of samples
# and renamed pathways
diff.virginia.path.clr <- matrix(data = NA, 
                                  ncol = ncol(diff.virginia.ko.clr),
                                  nrow = nrow(diff.virginia.ko.path.grp),
                                  dimnames = list(diff.virginia.ko.path.grp$pathway,
                                                  colnames(diff.virginia.ko.clr)))

# make named vector of pathways present in above feature table
diff.path.vec <- diff.virginia.ko.path$pathway
names(diff.path.vec) <- rownames(diff.virginia.ko.path)

# calculate mean CLR values per-sample for each pathway, using the scaled, CLR
# feature table of K0 terms
for(i in levels(factor(diff.path.vec))){
  diff.virginia.path.clr[i,1:297] <- colMeans(diff.virginia.ko.clr[which(diff.path.vec == i),])  
}

# create a function for calculating the z-score of an object
zscore <- function(x){
  (x -mean(x)) / sd(x)
}

# calculate z-scores for each pathway, normalising per sample (by column)
diff.virginia.path.clr.z <- apply(diff.virginia.path.clr, 2, zscore)

# create new metadata for differential abundance heatmap and edit group names
diff.meta <- as.data.frame(cbind(CST= virginia.cst$cst,
                                 Group= virginia.groups$groups.3), 
                           row.names = colnames(diff.virginia.path.clr.z))

diff.meta$Group <- case_when(diff.meta$Group == "Healthy-G1" ~ "H-1",
                             diff.meta$Group == "Healthy-G2" ~ "H-2",
                             .default = diff.meta$Group)

#rename diff.meta and save as .Rda for making heatmaps of individual K0s
# virginia.diff.meta <- diff.meta
# save(virginia.diff.meta,
#      file = paste(path.to.github, "Rdata/virginia.diff.meta.Rda",sep = ""))

# create colour bar list for CST & Group
diff.column.cols <- list(CST = c(I = viridis(7)[1], II = viridis(7)[2],
                                 III = viridis(7)[3], `III / IV` = viridis(7)[4],
                                 `III / V` = viridis(7)[5], IV = viridis(7)[6],
                                 V = viridis(7)[7]),
                         Group = c(`H-1` = "steelblue1",
                                   `H-2` = "royalblue3",
                                   `BV` = "goldenrod2"))

# make heatmap of z-scores for all differentially abundant abs(effect size >1)
# pathways (healthy vs. bv), using previous hclust object: CLR values derived
# from counts of K0 terms (health vs. bv; Ward's method, euclidean distance)
# png(paste(path.to.github,
#           "figs_for_paper/virginia_healthBV_diffAbund_K0.png", sep = ""),
#     units = "in", height = 11.8, width = 16.6, res = 400)

pheatmap(mat = diff.virginia.path.clr.z, cluster_cols = virginia.hclust,
         cutree_cols = 2, cutree_rows = 3, show_colnames = FALSE,
         annotation_colors = diff.column.cols, annotation_col = diff.meta,
         fontsize = 9, fontsize_row = 9, border_color = rgb(0,0,0,0.1),
         treeheight_row = 0, cellwidth = 3.2, cellheight = 7.75)

# dev.off()

# make the same heatmap, but set the breaks argument so that the heatmap colour
# gradient is centred around 0
# png(paste(path.to.github,
#           "figs_for_paper/virginia_healthBV_diffAbund_K0_scale.png", sep = ""),
#     units = "in", height = 11.8, width = 16.6, res = 400)

pheatmap(mat = diff.virginia.path.clr.z, cluster_cols = virginia.hclust,
         cutree_cols = 2, cutree_rows = 3, show_colnames = F, breaks = seq(-4,4,0.08),
         annotation_colors = diff.column.cols, annotation_col = diff.meta,
         fontsize = 9, fontsize_row = 9, border_color = rgb(0,0,0,0.1),
         treeheight_row = 0, cellwidth = 3.2, cellheight = 7.75)

# dev.off()

# subset scale sim'd, CLR-transformed health only data for differential K0s and
# make tibble of K0 and diff.btw (column_to_rownames() compatibility)
diff.path.20 <- virginia.filt.ko.clr.c[rownames(diff.virginia.ko.path),]
diff.path.20 <- tibble(`diff.btw` = as.numeric(diff.path.20$diff.btw),
                         path = diff.virginia.ko.path$pathway)

# get top 10 pathways most differential pathways in each group (20 largest 
# and 20 smallest diff.btw values) to plot on smaller version of heatmap
diff.path.20 <- diff.path.20 %>% 
  group_by(path) %>%
  summarise(mean.diff.btw= mean(diff.btw)) %>% 
  arrange(mean.diff.btw) %>% 
  filter(row_number() %in% c(1:10,(n()-9):n())) %>% 
  tibble::column_to_rownames("path")

# subset health only Z-scores for the most differential pathways in each group
diff.virginia.path.clr.z.sub <- diff.virginia.path.clr.z[rownames(diff.path.20),]
diff.virginia.path.clr.z.sub <- diff.virginia.path.clr.z.sub[-which(rownames(diff.virginia.path.clr.z.sub) %in% c("DNA replication","RNA polymerase (BV)","Pentose phosphate pathway (H)")),]
diff.virginia.path.clr.z.sub <- rbind(diff.virginia.path.clr.z.sub,
                                      diff.virginia.path.clr.z[c("Butanoate metabolism","Iron acquisition","CAMP resistance (H)"),1:297])
rownames(diff.virginia.path.clr.z.sub)[18:20] <- c("Butanoate metabolism","Iron acquisition","CAMP resistance (H)")

# plot z-score heatmap for most differential pathways clustering samples based
# on CLR values derived from counts of K0 terms (health vs. bv; Ward's method, 
# euclidean distance)
# png(paste(path.to.github,
#           "figs_for_paper/virginia_healthBV_diffAbund_K0_sub.png", sep = ""),
#     units = "in", height = 5.1, width = 16.75, res = 400)

pheatmap(mat = diff.virginia.path.clr.z.sub, cluster_cols = virginia.hclust,
         cutree_cols = 2, cutree_rows = 3, show_colnames = FALSE,
         annotation_colors = diff.column.cols, annotation_col = diff.meta,
         border_color = rgb(0,0,0,0.1), treeheight_row = 0, fontsize_row = 10.5,
         cellwidth = 3.2, cellheight = 13.5, breaks = seq(-4,4,0.08))

# dev.off()

############### molecular BV subgroups: differential abundance ################

# subset the K0-aggregated virginia feature table to only include molecular BV
# samples and check for rows with 0 sums
ko.virginia.bv <- ko.virginia.filt[,which(virginia.groups$groups.3 == "BV")]
any(rowSums(ko.virginia.bv)==0) # returns FALSE

# set seed for RNG
set.seed(2023)

# zero replacement for just bv samples in virginia dataset (samples expected
# as rows, so input dataframe transposed)
bv.virginia.ko.zr <- cmultRepl(t(ko.virginia.bv),label = 0, 
                               method = "CZM", z.warning = 0.99)

# CLR transform of BV only, zero-replaced, virgina K0 dataset (apply CLR across
# columns, so input dataframe transposed)
bv.virginia.ko.clr <- apply(t(bv.virginia.ko.zr), 2, function(x){log(x) - mean(log(x))})

# cluster euclidean distances using Ward's method
bv.virginia.hclust <- hclust(dist(t(bv.virginia.ko.clr),method = "euclidean"),
                    method = "ward.D2")

# plot dendrogram of BV samples 
plot(bv.virginia.hclust, hang= -1, cex = 0.7)

# there are three clear subgroups, though one only has around 20 samples. Split
# samples into these three groups and count number of samples in each
bv.virginia.cutree <- as.data.frame(cutree(bv.virginia.hclust, k = 3), nm = "group")
bv.virginia.cutree %>% 
  group_by(group) %>% 
  count()

# marker samples for each group:
# 1 = SRR6744999 (57 samples)
# 2 = SRR6743927 (56 samples)
# 3 = SRR6744444 (21 samples)

# for now, limit analysis to groups 1 and 2 (could run 3-way analysis later);
# filter out group 3 and rename groups 1 and 2
# NOTE: I realised that groups are arbitrarily swapped relative to lon/eur BV 
#       subgroup analysis, such that cutree group 1 here is BV2 (non-motile) and
#       cutree group 2 is BV1 (motile). Swap these in the case_when call to get
#       consistent group labels (SDS, 2024-04-24).
bv.virginia.cutree <- bv.virginia.cutree %>% 
  filter(group == 1 | group ==  2)

bv.virginia.cutree$group <- case_when(bv.virginia.cutree$group == 2 ~ "BV Subgroup 1",
                             bv.virginia.cutree$group == 1 ~ "BV Subgroup 2")

# save bv subgroup vector for stacked bar chart
# virginia.bv.groups <- bv.virginia.cutree
# save(virginia.bv.groups,
#      file = paste0(path.to.github,"Rdata/virginia.bv.groups.Rda"))

# subset K0-aggregated virginia feature table for only these samples and check
# for zero row sums
sub.virginia.ko <- ko.virginia.bv[,rownames(bv.virginia.cutree)]
any(rowSums(sub.virginia.ko) == 0) # returns FALSE

# make scale matrix using sd of 0.5 and difference of 15%
mu.matrix <- aldex.makeScaleMatrix(gamma = 0.5, mu = c(1, 1.15),
                                   conditions = bv.virginia.cutree$group)

# CLR transform K0-aggregated virgina bv-only feature table using scale sim
sub.virginia.ko.clr <- aldex.clr(reads = sub.virginia.ko, mc.samples = 128,
                                 conds = bv.virginia.cutree$group,
                                 gamma = mu.matrix, verbose = TRUE)

# rename aldex.clr object as .Rda for network analysis later
graph.clr <- sub.virginia.ko.clr

# calculate effect sizes, run t-tests and combine data frames
sub.virginia.ko.clr.e <- aldex.effect(sub.virginia.ko.clr, verbose = TRUE,
                                      include.sample.summary = TRUE)
sub.virginia.ko.clr.t <- aldex.ttest(sub.virginia.ko.clr,
                                     verbose = TRUE)
sub.virginia.ko.clr.all <- as.data.frame(cbind(sub.virginia.ko.clr.e,
                                               sub.virginia.ko.clr.t))

# get names of K0 terms aligning to 'Aminoacyl-tRNA biosynthesis', Ribosome'
# and 'Glycolysis / Gluconeogenesis' for checking centering on the MW plot
cent.trna <- rownames(virginia.ko.path.filt)[which(virginia.ko.path.filt$pathway == "Aminoacyl-tRNA biosynthesis")]
cent.ribo <- rownames(virginia.ko.path.filt)[which(virginia.ko.path.filt$pathway == "Ribosome")]
cent.glyc <- rownames(virginia.ko.path.filt)[which(virginia.ko.path.filt$pathway == "Glycolysis / Gluconeogenesis")]

# make MW plot and overlay points corresponding to K0 numbers for absolute 
# effect size  >1, 'Ribosome' and 'Aminoacyl-tRNA biosynthesis'
# png(paste(path.to.github,
#           "figs_for_paper/virginia_BVSubgroups_MW.png", sep = ""),
#     units = "in", height =7,  width = 14, res = 400)

aldex.plot(sub.virginia.ko.clr.all, xlim=c(0.4225,11), ylim = c(-12, 12), cutoff.pval=0.01)
title("ALDEx2 w/ ScaleSim (Virginia BV-only filtered K0): mu = 1 : 1.15, gamma = 0.5, p <0.01",
      adj=0, line= 0.8)

points(sub.virginia.ko.clr.all[which(abs(sub.virginia.ko.clr.all$effect) >1), "diff.win"],
       sub.virginia.ko.clr.all[which(abs(sub.virginia.ko.clr.all$effect) >1), "diff.btw"],
       pch= 19, col= "purple3", cex= 0.5)

points(sub.virginia.ko.clr.all[cent.ribo, "diff.win"],
       sub.virginia.ko.clr.all[cent.ribo, "diff.btw"],
       pch= 19, col= "blue3", cex= 0.5)

points(sub.virginia.ko.clr.all[cent.trna, "diff.win"],
       sub.virginia.ko.clr.all[cent.trna, "diff.btw"],
       pch= 19, col= "royalblue1", cex= 0.5)

points(sub.virginia.ko.clr.all[cent.glyc, "diff.win"],
       sub.virginia.ko.clr.all[cent.glyc, "diff.btw"],
       pch= 19, col= "lightskyblue1", cex= 0.5)

legend("bottomleft", pch = 19, cex = 0.75,
       col = c( "black", "red", "purple3", "blue3", "royalblue1", "lightskyblue1"),
       legend = c("Non-significant", "P <0.01", "Absolute effect size >1",
                  "Ribosome", "Aminoacyl-tRNA biosynthesis", "Glycolysis / Gluconeogenesis"))

# dev.off()

# looks like there aren't any K0s up in subgroup 2, but several in subgroup 1

# get pathways associated with K0s showing absolute effect size >1
sub.virginia.sig.ko <- virginia.ko.path.filt[which(abs(sub.virginia.ko.clr.all$effect) >1), , drop = FALSE]

# ABC transporters, bacterial chemotaxis, flagellar assembly, porphyrin
# biosynthesis and two-component systems are all on there!

# pull effect sizes for these differential K0s (should all be negative)
sub.virginia.sig.ko$effect <- sub.virginia.ko.clr.all[rownames(sub.virginia.sig.ko),"effect"]

# remove K0 terms with 'Unknown' pathway
sub.virginia.sig.ko <- sub.virginia.sig.ko %>% 
  filter(pathway != "Unknown")

# ========== FOR INVESTIGATING PATHWAYS OF INTEREST (START) ==========

# pull KO numbers for differential pathways of interest (BV subgroups)
rownames(sub.virginia.sig.ko %>% 
           filter(pathway == "Bacterial chemotaxis"))

# K00575 K03406 K03407 K03408 K03410 K03411 K03412 K03413 K03415 K06595

rownames(sub.virginia.sig.ko %>% 
           filter(pathway == "Porphyrin metabolism"))

# K01698

rownames(sub.virginia.sig.ko %>% 
           filter(pathway == "Flagellar assembly"))
# K02389 K02390 K02392 K02396 K02397 K02398 K02400 K02405 K02406 K02407
# K02408 K02409 K02410 K02411 K02412 K02414 K02416 K02417 K02419 K02422
# K02556 K06603 K13626

rownames(sub.virginia.sig.ko %>% 
           filter(pathway == "Biofilm formation"))
# K11936; note that this label is changed in the final heatmap figure to
#         'exopolysaccharide biosynthesis' to match the KO term label in
#         KEGG and the corresponding heatmap from the lon/eur dataset. Did this
#         in inkscape rather than R because I am lazy.

# almost all found in lon/eur dataset; exceptions are chemotaxis response
# regulator CheY (K03413), heme-based oxygen sensor hemAT (K06595), flagellar
# basal body rod modification protein flgD (K02389), flagellar sigma factor fliA
# (K02405) and flagellar hook-basal body complex subunit fliE

# filter master VIRGO vNum -> KO lookup table for signicant KOs between BV1 and
# BV2
sub.interest.vnum.ko <- KO %>% 
  filter(V2 %in% rownames(sub.virginia.sig.ko))

# subset vNum -> KO table for presence in lon/eur dataset
sub.interest.vnum.ko <- sub.interest.vnum.ko[rownames(sub.interest.vnum.ko) %in% rownames(virginia.data),]

# subset vNum -> KO table for several genes of interest
sub.interest.exopoly <- sub.interest.vnum.ko %>% 
  filter(V2 %in% c("K11936")) # NOTE: this is the same KO as 'exopolysaccharide
                              #       biosynthesis' in lon/eur dataset. Edited
                              #       in Inkscape because I am lazy.

sub.interest.porph <- sub.interest.vnum.ko %>% 
  filter(V2 %in% c("K01698"))

sub.interest.chemotax <- sub.interest.vnum.ko %>% 
  filter(V2 %in% c("K00575","K03406","K03407","K03408","K03410","K03411",
                   "K03412","K03413","K03415","K06595"))

sub.interest.flagella <- sub.interest.vnum.ko %>% 
  filter(V2 %in% c("K02389","K02390","K02392","K02396","K02397","K02398",
                   "K02400","K02405","K02406","K02407","K02408","K02409",
                   "K02410","K02411","K02412","K02414","K02416","K02417",
                   "K02419","K02422","K02556","K06603","K13626"))

# bind dataframes of interest together and add pathway & taxonomy info
sub.interest <- rbind(sub.interest.chemotax, sub.interest.exopoly,
                      sub.interest.flagella, sub.interest.porph)

sub.interest$pathway <- rep(c("Bacterial chemotaxis",
                              "Exopolysaccharide biosynthesis",
                              "Flagellar assembly", "Porphyrin metabolism"),
                            c(nrow(sub.interest.chemotax),
                              nrow(sub.interest.exopoly),
                              nrow(sub.interest.flagella),
                              nrow(sub.interest.porph)))

sub.interest$taxonomy <- tax.table[rownames(sub.interest),2]

# subset entire virginia dataset by vNumbers of interest and BV1/2 samples
sub.interest.df <- virginia.data[rownames(sub.interest),
                                 colnames(sub.virginia.ko)]

# pull row sums for each vNumber
sub.interest$sum.reads <- rowSums(sub.interest.df)

# relocate column for KO to the right of pathway
sub.interest <- sub.interest %>% 
  relocate(V2, .after = pathway)

# check total reads for each vNumber, grouped by pathway and taxonomy
sub.interest.sum <- sub.interest %>% 
  group_by(pathway, taxonomy) %>% 
  summarise(total.reads = sum(sum.reads))

# check to see if L. iners flagellar gene expression is specific to BV1.Pull
# bv1 and bv2 samples
sub.bv1 <- rownames(bv.virginia.cutree %>% filter(group=="BV Subgroup 1"))
sub.bv2 <- rownames(bv.virginia.cutree %>% filter(group=="BV Subgroup 2"))

# subset virginia data for vNumbers of interest for BV1 and BV2 samples
sub.bv1.df <- virginia.data[rownames(sub.interest), sub.bv1]
sub.bv2.df <- virginia.data[rownames(sub.interest), sub.bv2]

# pull row sums for each vNumber
sub.interest$total.bv1 <- rowSums(sub.bv1.df)
sub.interest$total.bv2 <- rowSums(sub.bv2.df)

# check total reads for each vNumber, grouped by pathway and taxonomy
sub.interest.sum <- sub.interest %>% 
  group_by(pathway, taxonomy) %>% 
  summarise(total.reads = sum(sum.reads),
            bv1.reads = sum(total.bv1),
            bv2.reads = sum(total.bv2))

# ========== FOR INVESTIGATING PATHWAYS OF INTEREST (END) ==========

####################### molecular BV subgroups: heatmap #######################

# extract matrix of CLR values from combined aldex object and edit column names
# to remove 'rab.sample'
sub.virginia.ko.clr <- sub.virginia.ko.clr.all[, 4:116]
colnames(sub.virginia.ko.clr) <- gsub("rab.sample.", "", colnames(sub.virginia.ko.clr))

# create empty matrix for calculating mean pathway CLR values
sub.virginia.path.clr <- matrix(data = NA, 
                                nrow = length(levels(factor(sub.virginia.sig.ko$pathway))),
                                ncol = ncol(sub.virginia.ko),
                                dimnames = list(levels(factor(sub.virginia.sig.ko$pathway)),
                                                colnames(sub.virginia.ko.clr)))

# generate pathway vector for all K0 terms in virginia dataset
virginia.path.vec <- virginia.ko.path.filt$pathway
names(virginia.path.vec) <- rownames(virginia.ko.path.filt)

# calculate mean CLR values per sample for each pathway
for(i in levels(factor(sub.virginia.sig.ko$pathway))){
  sub.virginia.path.clr[i,1:113] <- colMeans(sub.virginia.ko.clr[which(virginia.path.vec == i),])
}

# convert mean CLR values for each pathway to z-scores
sub.virginia.path.z <- apply(sub.virginia.path.clr, 2, zscore)

# make dataframe for metadata (CST and BV subgroup)  and remove 'subgroup' from
# the Subgroup column
sub.meta <- as.data.frame(cbind(CST= diff.meta[colnames(sub.virginia.path.z),1],
                                Subgroup=bv.virginia.cutree$group),
                          row.names = colnames(sub.virginia.path.z))

sub.meta$Subgroup <- gsub(" Subgroup ", "-", sub.meta$Subgroup)

# check which CSTs are present in these subgroups and remind myself which
# viridis colours are associated with these CSTs
levels(factor(sub.meta$CST))

# III = viridis(7)[3]
# III / IV = viridis(7)[4]
# IV = viridis(7)[6]

# make list for metadata columns (CST & BV subgroup)
sub.column.cols <- list(CST= c(III = viridis(7)[3],
                               `III / IV` = viridis(7)[4],
                               IV = viridis(7)[6]),
                        Subgroup = c(`BV-1` = "goldenrod2",
                                     `BV-2` = "goldenrod4"))

# plot z-score heatmap using Ward's method for clustering (non-scaled)
# png(paste(path.to.github,
#           "figs_for_paper/virginia_BVSubgroups_diffAbund_K0.png", sep = ""),
#     units = "in", height = 3.15, width = 8.4, res = 400)

pheatmap(mat = sub.virginia.path.z, clustering_method = "ward.D2",
         cutree_rows = 2, cutree_cols = 2, show_colnames = FALSE,
         annotation_colors = sub.column.cols, annotation_col = sub.meta,
         fontsize_row = 10.5, border_color = rgb(0,0,0,0.1),
         cellwidth = 3.2, cellheight = 13.5, treeheight_row = 0, scale = "none")

# dev.off()

# plot z-score heatmap using Ward's method for clustering (scaled)
# png(paste(path.to.github,
#           "figs_for_paper/virginia_BVSubgroups_diffAbund_K0_scaled.png", sep = ""),
#     units = "in", height = 3.15, width = 8.4, res = 400)

pheatmap(mat = sub.virginia.path.z, clustering_method = "ward.D2",
         cutree_rows = 2, cutree_cols = 2, show_colnames = FALSE,
         annotation_colors = sub.column.cols, annotation_col = sub.meta,
         fontsize_row = 10.5, border_color = rgb(0,0,0,0.1),
         cellwidth = 3.2, cellheight = 13.5, treeheight_row = 0, scale = "column")

# dev.off()

# make reproducible hclust object of clustering for reordering stacked bar graph
# virginia.bv.groups.hclust <- hclust(dist(t(sub.virginia.path.z),
#                                          method = "euclidean"),
#                                     method = "ward.D2")
# 
# save(virginia.bv.groups.hclust, file = paste0(path.to.github,
#                                               "Rdata/virginia.bv.groups.hclust.Rda"))

####################### molecular BV subgroups: biplots #######################

# KOs

# perform PCA on scale sim'd, CLR transformed BV data - ALL K0 terms
pca.virginia.bv.sub <- prcomp(t(sub.virginia.ko.clr))

# get list of indices for differential pathways in the pathway vector for all
# filtered dataset K0 numbers (will be highlighted on PCA plots) and add "Other"
sub.pca.path.ind <- list()
for(i in levels(factor(rownames(sub.virginia.path.clr)))){
  sub.pca.path.ind[[i]] <- which(virginia.path.vec == i)
}

sub.pca.path.ind[["Other"]] <- setdiff(1:length(virginia.path.vec),
                                          unlist(sub.pca.path.ind))

# make list of grouping indices
sub.pca.grp.ind <- list()
for(i in levels(factor(bv.virginia.cutree$group))){
  sub.pca.grp.ind[[i]] <- which(bv.virginia.cutree$group == i)
}

# make vector for colouring pathways
sub.pca.path.cols <- c("olivedrab2","tomato3","dodgerblue3","purple3",
                       "brown4","bisque2","yellow4","yellow2","black",
                       "cyan3",rgb(0,0,0,0.035))

# draw pca plots for Virginia dataset: all K0s and differential K0s
# png(paste(path.to.github,
#           "figs_for_paper/virginia_BVSubgroups_biplot_K0all.png", sep = ""),
#     units = "in", height = 6, width = 12, res = 400)

codaSeq.PCAplot(pca.virginia.bv.sub, plot.groups = TRUE, plot.loadings = TRUE,
                plot.ellipses = NULL, plot.density = "groups", PC = c(1,2),
                grp = sub.pca.grp.ind, grp.col = c("goldenrod2", "goldenrod4"), 
                grp.sym = "text", grp.cex = 0.5, load.grp = sub.pca.path.ind,
                load.col = sub.pca.path.cols, load.sym = rep(19, 11), 
                load.cex = 0.5, plot.legend = "loadings", 
                leg.position = "bottomright", leg.cex = 0.65, leg.columns = 3,
                title = "Virginia - K0 terms: BV subgroups (all K0 numbers)")

legend(x = 0.69485, y = 0.6045, legend = c("BV Subgroup 1", "BV Subgroup 2"),
       pch = "SS", col = c("goldenrod2","goldenrod4"), cex = 0.55, pt.cex = 0.65)

# dev.off()

# vNumbers by species

# identify all vNumbers in 'new.both.filt' which correspond to the significantly
# different KO terms
bv.vnum.sig.ko <- KO %>%                        # pull ALL vnumbers corresponding
  filter(V2 %in% rownames(sub.virginia.sig.ko)) # to significant KOs

bv.vnum.sig.ko <- bv.vnum.sig.ko %>%                            # filter to just
  filter(rownames(bv.vnum.sig.ko) %in% rownames(virginia.filt)) # vnums in dataset

# subset 'new.both.filt' by the vNumbers corresponding to the significantly
# different KO terms (n = 114)
bv.virginia.filt <- virginia.filt[rownames(bv.vnum.sig.ko),]
bv.virginia.filt <- bv.virginia.filt[,colnames(sub.virginia.ko)]

# check for row and column sum zeros
any(colSums(bv.virginia.filt)==0) # FALSE
any(colSums(bv.virginia.filt)==0) # FALSE

# set seed for RNG
set.seed(2023)

# generate matrix of scale values with mu = 1 and 1.15, std dev of log normal
# distribution = 0.5
mu.matrix <- aldex.makeScaleMatrix(gamma = 0.5, mu = c(1, 1.15),
                                   conditions = bv.virginia.cutree$group)

# CLR transform data and calculate median CLR values for all features (vNum)
bv.virginia.filt.clr <- aldex.clr(bv.virginia.filt,
                                  conds = bv.virginia.cutree$group,
                                  gamma = mu.matrix, verbose = TRUE)
bv.virginia.filt.clr.e <-aldex.effect(bv.virginia.filt.clr, verbose = TRUE,
                                      include.sample.summary = TRUE)

bv.virginia.filt.clr.df <- bv.virginia.filt.clr.e[, grep("rab.sample.",
                                                         colnames(bv.virginia.filt.clr.e))]

colnames(bv.virginia.filt.clr.df) <- gsub("rab.sample.", "", colnames(bv.virginia.filt.clr.df))

# run PCA
pca.bv.sig.vnum <- prcomp(t(bv.virginia.filt.clr.df))

# get named vector of taxa present in the vNumber subset table
bv.sig.vnum.taxa <- tax.table[rownames(bv.virginia.filt.clr.df),2]
names(bv.sig.vnum.taxa) <- rownames(bv.virginia.filt.clr.df)

# make list of vNumber indices for each species
bv.ind.sig.vnum <- list()
for(i in levels(factor(bv.sig.vnum.taxa))){
  bv.ind.sig.vnum[[i]] <- which(bv.sig.vnum.taxa == i)
}

# get colours corresponding to taxa in the list of vNumbers corresponding to 
# significantly different KO terms, plus transparent grey for 'Other'
pca.load.cols.vnum <- c(tax.colors[names(bv.ind.sig.vnum),], rgb(0,0,0,0.1))

#add vNumbers with unknown taxa to 'Other'
bv.ind.sig.vnum[["Other"]] <- setdiff(1:nrow(bv.virginia.filt.clr.df),
                                   unlist(bv.ind.sig.vnum))

# make PCA plot of significantly different vNumbers by species in BV subgroups
codaSeq.PCAplot(pca.bv.sig.vnum, plot.groups = FALSE, plot.loadings = TRUE,
                plot.ellipses = NULL, plot.density = NULL, grp.cex = 0.3,
                load.grp = bv.ind.sig.vnum, load.col = pca.load.cols.vnum,
                load.sym = rep(19,5), load.cex = 1, PC = c(1,2),
                plot.legend = NULL, leg.position = "topleft", 
                leg.cex = 0.85, leg.columns = 2,
                title = "PCA - Virginia BV subgroups: significantly different vNumbers")

################# molecular BV subgroups: correlation network #################
# NOTE: the 'Get started' page on the igraph website is really useful for
#       understanding the syntax of igraph functions & how to edit graphs; see:
#       https://r.igraph.org/articles/igraph.html

library(igraph)

# set seed for reproducibility
set.seed(2023)

# get expected values of Phi from the dirichlet log-ratio distribution
graph.phi <- codaSeq.phi(graph.clr)

# subset for Phi <0.025
graph.phi.sub <- subset(graph.phi, phi <=0.025)

# pull pathways for all K0s in 'row' and 'col' columns
# NOTE: pulling pathways by subsetting the K0->pathway table returns incorrect 
#       pathways (returns pathway from the row above what is intended) for
#       several K0s, so had to do it via a loop. Possibly due to duplicate K0 
#       numbers in the phi subset dataframe?
for(i in levels(factor(graph.phi.sub$row))){
  index <- grep(i, rownames(virginia.ko.path.filt))
  path <- virginia.ko.path.filt$pathway[index]
  graph.phi.sub$pathr <- case_when(graph.phi.sub$row == i ~ path,
                                   .default = graph.phi.sub$pathr)
}

for(i in levels(factor(graph.phi.sub$col))){
  index <- grep(i, rownames(virginia.ko.path.filt))
  path <- virginia.ko.path.filt$pathway[index]
  graph.phi.sub$pathc <- case_when(graph.phi.sub$col == i ~ path,
                                   .default = graph.phi.sub$pathc)
}

# filter phi subset by pathway to remove any correlations where one or both 
# pathways are unknown
graph.phi.sub <- graph.phi.sub %>% 
  filter(pathc != "Unknown" & pathr != "Unknown")

# create a non-redundant vector of K0 -> pathway in alphabetical order
graph.vert.pathways <-as.character(c(graph.phi.sub$row, graph.phi.sub$col))
names(graph.vert.pathways) <-c(graph.phi.sub$pathr, graph.phi.sub$pathc)
graph.vert.pathways <- graph.vert.pathways[-which(duplicated(graph.vert.pathways))]

# remove 'pathr' and 'pathc' columns from phi dataframe
graph.phi.sub <- graph.phi.sub %>% 
  dplyr::select(c(row,col,phi))

# generate the igraph graphical object and set vertex label cex to 0.8
graph.sub <- graph.data.frame(graph.phi.sub, directed = FALSE)
V(graph.sub)$label.cex <- 0.8

# confirm graph vertex names and pathway names are equal
length(which(V(graph.sub)$name == graph.vert.pathways)) # 105x TRUE

# set vertex names to pathways represented by these K0s
V(graph.sub)$pathway <- names(graph.vert.pathways)

# set layout for graphs using Fruchterman-Reingold (force-directed) algorithm
graph.layout <- layout_with_fr(graph.sub)

# plot graphs with colour automatically set
plot(graph.sub, layout = graph.layout,
     vertex.color = as.factor(V(graph.sub)$pathway),
     vertex.size = 5)

# decompose graph as there is only one with  a large number of (non-ribosomal)
# vertices
graph.sub.d <- decompose.graph(graph.sub) # located in 3rd list element

# vertex names for decomposed graph #3 and corresponding pathways
graph.vert.pathways.3 <- V(graph.sub.d[[3]])$name
names(graph.vert.pathways.3) <- names(graph.vert.pathways[which(graph.vert.pathways %in% graph.vert.pathways.3)])

# reassign pathways for decomposed graph
V(graph.sub.d[[3]])$pathway <- names(graph.vert.pathways.3)

# set seed for layout
set.seed(100)

# set layout for decomposed graph using FR algorithm
graph.layout.3<- layout_with_gem(graph.sub.d[[3]])

# make vector of colours for decomposed graph vertices
graph.vert.cols.3 <- as.character(factor(names(graph.vert.pathways.3),
                                         labels = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdYlBu")))(16)))

# make character vector of pathways for legend and add two-letter abbreviations
graph.path.abr <- levels(factor(V(graph.sub.d[[3]])$pathway))

graph.path.abr <- paste(graph.path.abr, sep = " ",c("(At)","(Ap)","(Cr)","(Dd)",
                                                    "(Fa)","(Ly)","(Nm)","(Os)",
                                                    "(Op)","(Pp)","(Pg)","(Pu)",
                                                    "(Py)","(Sd)","(Tp)","(Tc)"))

# plot correlation network for the largest, non-ribosomal graph
# png(paste(path.to.github,
#           "figs_for_paper/virginia_BVSubgroups_corrNetwork.png", sep = ""),
#     units = "in", height = 10, width = 10, res = 400)

plot(graph.sub.d[[3]], layout = graph.layout.3,
     vertex.color = graph.vert.cols.3, vertex.label.cex = 0.95,
     vertex.size = 17, vertex.label.color= "black", vertex.label.family = "sans")

legend(x = -1.15, y = -0.135, legend = graph.path.abr, pch = 21, cex = 1,
       pt.cex = 2, box.lwd = 0, bg = rgb(0,0,0,0), ncol = 1,
       pt.bg = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdYlBu")))(16))

# dev.off()

################## healthy subgroups: differential abundance ###################

# subset virginia grouping dataframe for only healthy samples
health.groups <- virginia.groups %>% 
  filter(groups.3 != "BV") %>% 
  dplyr::select(-groups.2)

# subset K0-aggregagted, filtered virginia data for only healthy samples and
# check for zero-sum rows
health.ko.filt <- ko.virginia.filt[,rownames(health.groups)]
any(rowSums(health.ko.filt)==0) # FALSE

# set seed again to ensure consistency of results
set.seed(2023)

# make matrix of scale distributions for healthy samples: sd of 0.5 and
# difference of 1: 1.15
mu.matrix.h <- aldex.makeScaleMatrix(gamma = 0.5, mu = c(1, 1.15), 
                                     conditions = health.groups$groups.3)

# CLR-transform healthy-only virginia K0 dataset with scale sim
health.ko.filt.clr <- aldex.clr(health.ko.filt, conds = health.groups$groups.3,
                                gamma = mu.matrix.h, verbose = TRUE)

# calculate effect size and t-test statistics, then merge dataframes
health.ko.filt.clr.e <- aldex.effect(health.ko.filt.clr, verbose = TRUE,
                                     include.sample.summary = TRUE)
health.ko.filt.clr.t <- aldex.ttest(health.ko.filt.clr,
                                    verbose = TRUE)
health.ko.filt.clr.c <- cbind(health.ko.filt.clr.e, health.ko.filt.clr.t)

# plot intra-group CLR dispersion vs. median CLR difference between groups (MW)
# png(paste(path.to.github,
#           "figs_for_paper/virginia_HSubgroups_MW.png", sep = ""),
#     units = "in", height = 7, width = 14, res = 400)

aldex.plot(health.ko.filt.clr.c, xlim = c(0.346, 8.75), cutoff.pval = 0.01)
title("ALDEx2 w/ ScaleSim (Virginia Healthy-only filtered K0): mu = 1 : 1.15, gamma = 0.5, p <0.01",
      adj=0, line= 0.8)

points(health.ko.filt.clr.c[cent.ribo, "diff.win"],
       health.ko.filt.clr.c[cent.ribo, "diff.btw"],
       pch = 19, col = "blue3", cex = 0.5)

points(health.ko.filt.clr.c[cent.trna, "diff.win"],
       health.ko.filt.clr.c[cent.trna, "diff.btw"],
       pch = 19, col = "royalblue1", cex = 0.5)

points(health.ko.filt.clr.c[cent.glyc, "diff.win"],
       health.ko.filt.clr.c[cent.glyc, "diff.btw"],
       pch = 19, col = "lightskyblue1", cex = 0.5)

points(health.ko.filt.clr.c[which(abs(health.ko.filt.clr.c$effect) >1), "diff.win"],
       health.ko.filt.clr.c[which(abs(health.ko.filt.clr.c$effect) >1), "diff.btw"],
       pch = 19, col = "purple3", cex = 0.5)

legend("bottomleft", pch = 19, cex = 0.75,
       col = c( "black", "red", "purple3", "blue3", "royalblue1", "lightskyblue1"),
       legend = c("Non-significant", "P <0.01", "Absolute effect size >1",
                  "Ribosome", "Aminoacyl-tRNA biosynthesis", "Glycolysis / Gluconeogenesis"))

# dev.off()

# stuff upregulated in both groups

# subset K0 -> pathway table for terms with abs(effect) >1, pull efffect sizes
# and remove 'Unknown' K0s
health.sig.ko.path <- virginia.ko.path.filt[which(abs(health.ko.filt.clr.c$effect)>1), , drop = F]
health.sig.ko.path$effect <- health.ko.filt.clr.c[rownames(health.sig.ko.path),'effect']
health.sig.ko.path <- health.sig.ko.path %>% 
  filter(pathway != "Unknown")

# append " (H1)" or " (H2)" to pathways ONLY if there is >1 K0 for each pathway
# AND there are positive and negative effect sizes for that pathway
for(i in levels(factor(health.sig.ko.path$pathway))){
  if(length(which(health.sig.ko.path$pathway == i)) >1 & length(unique(sign(health.sig.ko.path$effect[which(health.sig.ko.path$pathway == i)]))) != 1){
    health.sig.ko.path$pathway <- case_when(health.sig.ko.path$pathway == i & sign(health.sig.ko.path$effect) == 1 ~ paste(i, "(H2)", sep = " "),
                                            health.sig.ko.path$pathway == i & sign(health.sig.ko.path$effect) == -1 ~ paste(i, "(H1)", sep = " "),
                                            .default = health.sig.ko.path$pathway)
  }
}

# count number of K0s per pathway (total 93 pathways)
health.sig.ko.path.grp <- health.sig.ko.path %>% 
  group_by(pathway) %>% 
  count()

# subset scale-sim'd, CLR-transformed data for only differential K0s and remove
# 'rab.sample.' from column headers
health.ko.filt.clr.df <- health.ko.filt.clr.c[rownames(health.sig.ko.path),4:166]
colnames(health.ko.filt.clr.df) <- gsub("rab.sample.", "", colnames(health.ko.filt.clr.df))

# =========== INVESTIGATING PREVIOUSLY ID'D KOs OF INTEREST (START) ===========

# filter KO lookup tables for significance and presence in virginia healthy
# subgroup subset
health.interest.ko.vnum <- KO %>% 
  filter(V2 %in% rownames(health.sig.ko.path))

health.interest.ko.vnum <- health.interest.ko.vnum %>% 
  filter(rownames(health.interest.ko.vnum) %in% rownames(virginia.data[,colnames(health.ko.filt.clr.df)]))

# filter for putrescine/spermidine transporter and dioxin degradation KOs
# (K02055 and K01821), make vnumber column and sort by it
health.interest.filt <- health.interest.ko.vnum %>% 
  filter(V2 %in% c("K02055","K01821"))

health.interest.filt <- health.interest.filt %>% 
  mutate(vnum=rownames(health.interest.filt)) %>% 
  relocate(vnum, .after = V2) %>% 
  arrange(V2,vnum)

# make dataframe of vNumbers, KO term, taxonomy and read totals in virginia data
health.interest.filt.df <- data.frame(vnum=health.interest.filt$vnum,
                                      ko=health.interest.filt$V2,
                                      taxid=tax.table[health.interest.filt$vnum,2],
                                      read.totals= rowSums(virginia.data[health.interest.filt$vnum,
                                                                         colnames(health.ko.filt.clr.df)]))

# summarise read totals by KO and species
health.interest.filt.sum <- health.interest.filt.df %>% 
  group_by(taxid, ko) %>% 
  summarise(sum=sum(read.totals))

# lactobacilli (specifically crispatus, but also jensenii and gasseri to a far
# lesser extent) account for almost all reads assigned to K01821 and K02055

# ============ INVESTIGATING PREVIOUSLY ID'D KOs OF INTEREST (END) ============

# make pathway vector for differential pathways and add K0s as names
health.sig.path.vec <- health.sig.ko.path$pathway
names(health.sig.path.vec) <- rownames(health.sig.ko.path)

# make matrix for calculating mean CLR values per pathway
health.sig.path <- matrix(data = NA,
                          ncol = ncol(health.ko.filt.clr.df),
                          nrow = nrow(health.sig.ko.path.grp),
                          dimnames = list(row=levels(factor(health.sig.ko.path.grp$pathway)),
                                          col=colnames(health.ko.filt.clr.df)))
 
# calculate mean CLR values for all samples, per pathway
for(i in levels(factor(health.sig.path.vec))){
  health.sig.path[i,1:163] <- colMeans(health.ko.filt.clr.df[which(health.sig.path.vec == i),])
}

# apply Z-score function over rows
health.sig.path.z <- apply(health.sig.path, 2, zscore)

# subset virginia CST vector for only healthy samples and check which CSTs are
# present (6/7 CSTs are present, but only 3x samples with CST IV; remember these
# are based on molecular function)
virginia.cst.h <- virginia.cst[colnames(health.sig.path), , drop = FALSE]
virginia.cst.h %>% 
  group_by(cst) %>% 
  count()

# make dataframe of metadata for healthy samples and remove 'G' from group names
health.meta <- as.data.frame(cbind(CST = virginia.cst.h$cst,
                                   Group = health.groups$groups.3),
                             row.names = colnames(health.ko.filt.clr.df))

health.meta$Group <- gsub("G", "", health.meta$Group)

# make list for column colour bars
health.column.cols <- list(CST= c(I = viridis(7)[1],
                                  II = viridis(7)[2],
                                  III = viridis(7)[3],
                                  `III / V` = viridis(7)[5],
                                  IV = viridis(7)[6],
                                  V = viridis(7)[7]),
                           Group = c(`Healthy-1` = "steelblue1",
                                     `Healthy-2` = "royalblue3"))

# make heatmap of Z-scores per pathway for healthy subgroups (non-scaled)
# png(paste(path.to.github,
#           "figs_for_paper/virginia_HSubgroups_diffAbund_K0.png", sep = ""),
#     units = "in", height = 10.6, width = 12.8, res = 400)

pheatmap(mat = health.sig.path.z, cutree_cols = 3, cutree_rows = 3,
         annotation_col = health.meta, annotation_colors = health.column.cols,
         show_colnames = FALSE, treeheight_row = 0, treeheight_col = 15,
         fontsize_row = 8.75, cellwidth = 4, cellheight = 7.65, scale = "none",
         border_color = rgb(0,0,0,0.2), clustering_method = "ward.D2")

# dev.off()

# make heatmap of Z-scores per pathway for healthy subgroups (scaled)
# png(paste(path.to.github,
#           "figs_for_paper/virginia_HSubgroups_diffAbund_K0_scale.png", sep = ""),
#     units = "in", height = 10.6, width = 12.8, res = 400)

pheatmap(mat = health.sig.path.z, cutree_cols = 3, cutree_rows = 3,
         annotation_col = health.meta, annotation_colors = health.column.cols,
         show_colnames = FALSE, treeheight_row = 0, treeheight_col = 15,
         fontsize_row = 8.75, cellwidth = 4, cellheight = 7.65, scale = "column",
         border_color = rgb(0,0,0,0.2), clustering_method = "ward.D2")

# dev.off()

# cluster healthy subgroup samples on Z-scores of mean CLR per pathway as in
# above heatmaps (for clustering of a subset of pathways)
health.path.20.hclust <- hclust(dist(t(health.sig.path.z), method = "euclidean"),
                                method = "ward.D2")

# subset scale sim'd, CLR-transformed health only data for differential K0s and
# make tibble of K0 and diff.btw (data.frame() doesn't play well with numerics)
health.path.20 <- health.ko.filt.clr.c[rownames(health.sig.ko.path),]
health.path.20 <- tibble(`diff.btw` = as.numeric(health.path.20$diff.btw),
                         path = health.sig.ko.path$pathway,)

# get top 20 pathways most differential pathways in each group (20 largest 
# and 20 smallest diff.btw values) to plot on smaller version of heatmap
health.path.20 <- health.path.20 %>% 
  group_by(path) %>%
  summarise(mean.diff.btw= mean(diff.btw)) %>% 
  arrange(mean.diff.btw) %>% 
  filter(row_number() %in% c(1:10,(n()-9):n())) %>% 
  tibble::column_to_rownames("path")

# subset health only Z-scores for the most differential pathways in each group
health.sig.path.z.sub <- health.sig.path.z[rownames(health.path.20),]
health.sig.path.rem <- c("Base excision repair","Protein digestion & absorption","Protein export","Phe-Tyr-Trp metabolism","Purine metabolism (H1)")
health.sig.path.z.sub <- health.sig.path.z.sub[-which(rownames(health.sig.path.z.sub) %in% health.sig.path.rem),]
health.sig.path.add <- c("Pentose & glucuronate interconversions","Putrescine/Spermidine transport","Extracellular DNA uptake","Flagellar assembly","Vitamin B3 metabolism (H1)")
health.sig.path.z.sub <- rbind(health.sig.path.z.sub,
                               health.sig.path.z[health.sig.path.add,1:163])
rownames(health.sig.path.z.sub) <- gsub("Pentose & glucuronate interconversions",
                                        "Pentose-glucuronate interconversion",
                                        rownames(health.sig.path.z.sub))

# make scaled heatmap of Z-scores for the top 20 differential pathways
# png(paste(path.to.github,
#           "figs_for_paper/virginia_HSubgroups_diffAbund_K0_sub.png", sep = ""),
#     units = "in", height = 4.75, width = 11.3, res = 400)

pheatmap(mat = health.sig.path.z.sub, cluster_cols = health.path.20.hclust,
         cutree_cols = 3, cutree_rows = 2, show_colnames = FALSE,
         annotation_col = health.meta, annotation_colors = health.column.cols,
         treeheight_row = 0, treeheight_col = 35, fontsize_row = 10.5,
         cellwidth = 3.2, cellheight = 13.5, breaks = seq(-2.5,2.5,0.05),
         border_color = rgb(0,0,0,0.1), clustering_method = "ward.D2")

# dev.off()
