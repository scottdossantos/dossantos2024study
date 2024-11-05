# no SRI: london/europe differential abundance

# load required packages
library(dplyr)
library(pheatmap)
library(ALDEx2) # MUST be version >1.35 for scale simulation (see note below)

# NOTE: Install the current version of ALDEx2 from GitHub using the following
#       command: devtools::install_github("ggloor/ALDEx_bioc")

# set path to github directory (edit 'path.to.github' to reflect your machine!)
user <- "sds/cc"
path.to.github <- switch(user,
                         "sds/cc" = "~/Documents/GitHub/metatranscriptome/",
                         "gg" = "~/Documents/0_git/projects/metatranscriptome/")

# run 'setup.R' to load feature tables etc.
source(paste(path.to.github, "code/setup.R", sep = ""))

# load KO -> pathway table for london/europe
load(paste(path.to.github, "Rdata/ko.both.path.Rda", sep = ""))

# load in vector containing cst info from london/europe species heatmap
load(paste(path.to.github, 'Rdata/hm.metadata.Rda', sep = ""))

# load in london/europe heatmap column colour bar list
load(paste(path.to.github, "Rdata/hm.column.cols.Rda",sep = ""))

# remove two samples from the filtered london/europe feature table aggregated by
# K0 number (classed as BV but almost no BV organisms):
#   -  v.001A: close to 100% L. gasseri with practically no BV organisms
#   -  v.019A: around 80 % iners, ~ 20% crispatus and a tiny bit of Gardnerella
ko.both<-ko.both[,-c(9,22)]

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

# make grouping vector for ALDEx2 with no SRI
ko.conds <- c(rep('H',8), rep('B',12), rep('B',14), rep('H', 8)) 


# set seed
set.seed(2023)

# run ALDEx2 but DO NOT USE SRI
ko.both.clr <- aldex.clr(reads = ko.both, conds = ko.conds,
                         gamma = 0.001, verbose = TRUE)

# calculate effect size
ko.both.clr.e <- aldex.effect(ko.both.clr, verbose = TRUE, 
                              include.sample.summary = TRUE)
# run t-tests
ko.both.clr.t <- aldex.ttest(ko.both.clr, verbose = TRUE)

# combine effect and t-test dataframes
ko.both.clr.all <- cbind(ko.both.clr.e, ko.both.clr.t)

# get indices for genes assigned to 'Ribosome', 'Aminoacyl-tRNA biosynthesis' 
# and 'Glycolysis / Gluconeogenesis'
path.ribo <- which(ko.both.path$pathway == "Ribosome")
path.trna <- which(ko.both.path$pathway == "Aminoacyl-tRNA biosynthesis")
path.glyc <- which(ko.both.path$pathway == "Glycolysis / Gluconeogenesis")

# plot difference vs. dispersion of median CLR values for health vs. BV and
# add coloured points for K0s with effect size >1, K0s assigned to 'Ribosome',
# 'Aminoacyl-tRNA biosynthesis' and 'Glycolysis / Gluconeogenesis' (housekeep.)

# png(filename = paste(path.to.github,
#                      "figs_for_paper/lon_eur_healthBV_MW.png", sep = ""),
#     units = "in", height = 7, width = 14, res = 400)

aldex.plot(ko.both.clr.all, xlim=c(0.346,9), cutoff.pval=0.01)
title('ALDEx2 w/ ScaleSim (London/Europe filtered K0): mu = 1: 1.15, gamma = 0.75, p <0.01',
      adj=0, line= 0.8)

points(ko.both.clr.all[which(abs(ko.both.clr.all$effect) >1),"diff.win"],
       ko.both.clr.all[which(abs(ko.both.clr.all$effect) >1),"diff.btw"],
       pch= 19, col= "purple3", cex= 0.5)

points(ko.both.clr.all[path.ribo ,"diff.win"],
       ko.both.clr.all[path.ribo ,"diff.btw"],
       pch= 19, col= "blue3", cex= 0.5)

points(ko.both.clr.all[path.trna,"diff.win"],
       ko.both.clr.all[path.trna,"diff.btw"],
       pch= 19, col= "royalblue1", cex= 0.5)

points(ko.both.clr.all[path.glyc,"diff.win"],
       ko.both.clr.all[path.glyc,"diff.btw"],
       pch= 19, col= "lightskyblue1", cex= 0.5)

legend("bottomleft", pch = 19, cex = 0.75,
       col = c( "black", "purple3", "blue3", "royalblue1", "lightskyblue1"),
       legend = c("Non-significant", "Absolute effect size >1",
                  "Ribosome", "Aminoacyl-tRNA biosynthesis", "Glycolysis / Gluconeogenesis"))

# dev.off()

# data are not centered as can be seen by location of housekeeping genes!

# pull KOs with absolute effect size >1 and their corresponding effect sizes,
# and filter out K0s whose pathways are unknown
diff.sig.ko <- ko.both.path[which(abs(ko.both.clr.all$effect) >1), , drop=FALSE]
diff.sig.ko$effect <- ko.both.clr.all$effect[which(abs(ko.both.clr.all$effect) >1)]
diff.sig.ko <- diff.sig.ko %>% 
  filter(pathway != 'Unknown')

# for pathways where there is >1 KO AND there are both +ve and -ve effect sizes,
# add a suffix of ' (H)' or ' (BV)' as applicable
for(i in levels(factor(diff.sig.ko$pathway))){
  if(length(which(diff.sig.ko$pathway == i)) >1 & length(unique(sign(diff.sig.ko$effect[which(diff.sig.ko$pathway == i)]))) != 1){
    diff.sig.ko$pathway <- case_when(diff.sig.ko$pathway == i & sign(diff.sig.ko$effect) == 1 ~ paste(i, "(H)", sep = " "),
                                     diff.sig.ko$pathway == i & sign(diff.sig.ko$effect) == -1 ~ paste(i, "(BV)", sep = " "),
                                     .default = diff.sig.ko$pathway)
  }
}

# collapse list of differential KOs by pathway
diff.sig.ko.path <- diff.sig.ko %>% 
  group_by(pathway) %>% 
  count()

# extract median CLR values (across all 128 MC instances) for all samples from
# summary ALDEx2 output, remove 'rab.sample.' from column headers and subset
# for only K0s with absolute effect size >1
ko.both.clr.df <- ko.both.clr.all[,4:45]
colnames(ko.both.clr.df) <- gsub('rab.sample.', '', colnames(ko.both.clr.df))
ko.both.clr.df <- ko.both.clr.df[rownames(diff.sig.ko),]

# make matrix for lon/eur samples, one row per differential pathway
ko.both.clr.path <- matrix(data = NA,
                           nrow = nrow(diff.sig.ko.path),
                           ncol = ncol(ko.both.clr.df),
                           dimnames = list(row = diff.sig.ko.path$pathway,
                                           col = colnames(ko.both.clr.df)))

# fill in matrix rows by calculating mean CLR values per pathway for all samples
for(i in rownames(ko.both.clr.path)){
  ko.both.clr.path[i,1:42] <- colMeans(ko.both.clr.df[which(diff.sig.ko$pathway == i),])
}

# create function for calculating z-scores
zscore <- function(x){
  (x - mean(x)) / sd(x)
}

# calculate z-scores for each pathway, normalising per sample (by column)
ko.both.clr.path.z <- apply(ko.both.clr.path, 2, zscore)

# subset the ALDEx2 summary output for only differential KOs and pull pathways
# and diff.btw into tibble (column_to_rownames() compatibility)
diff.sig.path.20 <- ko.both.clr.all[rownames(diff.sig.ko),]

diff.sig.path.20 <- tibble(diff.btw = diff.sig.path.20$diff.btw,
                           path = diff.sig.ko$pathway)

# group K0s by pathway, calculate mean diff.btw and take pathways with the 10
# highest and lowest mean diff.btw values (only 9 pathways with +ve diff.btw)
diff.sig.path.20 <- diff.sig.path.20 %>% 
  group_by(path) %>% 
  summarise(mean.diff.btw = mean(diff.btw)) %>% 
  arrange(mean.diff.btw) %>% 
  filter(row_number() %in% c(1:10,(n()-8):n())) %>% 
  tibble::column_to_rownames("path")

# subset mean Z-score/pathway dataframe by top 10 pathways
ko.both.clr.path.z.sub <- ko.both.clr.path.z[rownames(diff.sig.path.20),]

# plot PCA
# png("~/Documents/GitHub/dossantos2024study/figs_for_paper/R_FigsLondonEurope/suppl_lon_eur_diffAbund_noSRI.png",
#     units = "in", res = 400, height = 4.25, width = 10)

pheatmap(ko.both.clr.path.z.sub, cutree_cols = 2, cutree_rows = 2,
         annotation_col = hm.metadata[,-3], annotation_colors = hm.column.cols,
         show_colnames = FALSE, treeheight_row = 0, treeheight_col = 50,
         clustering_method = "ward.D2", border_color = rgb(0,0,0,0.2), 
         breaks = seq(-4.5,4.5,0.09), fontsize = 10, cellwidth = 10, cellheight = 10)

# dev.off()