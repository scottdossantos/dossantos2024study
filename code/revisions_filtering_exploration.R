# explore less stringent filtering for taxonomic biplots

library(ALDEx2)
library(CoDaSeq)
library(dplyr)

#################################### setup ####################################

# gh path
setwd("~/Documents/GitHub/dossantos2024study/")
github <- "~/Documents/GitHub/dossantos2024study/"

# load in complete feature tables for lon/eur and virginia datasets
load(paste0(github, "Rdata/new.both.data.Rda")) # batch-corrected
load(paste0(github, "Rdata/virginia.data.Rda"))

# species lookup table for VIRGO
lookup.spp <- read.table(paste0(github, "1_VIRGO/1.taxon.tbl.txt"),
                         header = F, sep = "\t", quote = "", row.names = 2)

# # load tax vec for complete virginia dataset and make one for lon/eur
# load(paste0(github, "Rdata/virginia.tax.vec.all.Rda"))
# tax.vec.loneur <- lookup.spp[rownames(new.both.data),2]
# names(tax.vec.loneur) <- rownames(new.both.data)

# filter with less stringent criteria or load if already exist
if(all(file.exists(paste0(github, "Rdata/filt.ls.loneur.Rda")),
       file.exists(paste0(github, "Rdata/filt.ls.virg.Rda")))){
  load(paste0(github,"Rdata/filt.ls.loneur.Rda"))
  load(paste0(github,"Rdata/filt.ls.virg.Rda"))
} else{
  filt.ls.loneur <- codaSeq.filter(new.both.data, samples.by.row=F,
                                   min.occurrence=0.10, min.prop=0.00005)
  
  filt.ls.virg <- codaSeq.filter(virginia.data, samples.by.row=F,
                                 min.occurrence=0.10, min.prop=0.00005)
  
  save(filt.ls.loneur, file = paste0(github, "Rdata/filt.ls.loneur.Rda"))
  save(filt.ls.virg, file = paste0(github, "Rdata/filt.ls.virg.Rda"))
}

# make new taxa vectors for these filtered datasets (load in if already exist)
if(all(file.exists(paste0(github, "Rdata/tax.vec.loneur.ls.Rda")),
       file.exists(paste0(github, "Rdata/tax.vec.virg.ls.Rda")))){
  load(paste0(github,"Rdata/tax.vec.loneur.ls.Rda"))
  load(paste0(github,"Rdata/tax.vec.virg.ls.Rda"))
} else{
  tax.vec.loneur.ls <- lookup.spp[rownames(filt.ls.loneur),2]
  names(tax.vec.loneur.ls) <- rownames(filt.ls.loneur)
  save(tax.vec.loneur.ls, file = paste0(github, "Rdata/tax.vec.loneur.ls.Rda"))
  
  tax.vec.virg.ls <- lookup.spp[rownames(filt.ls.virg),2]
  names(tax.vec.virg.ls) <- rownames(filt.ls.virg)
  save(tax.vec.virg.ls, file = paste0(github, "Rdata/tax.vec.virg.ls.Rda"))
}

# load filtered feature tables and taxa vectors
load(paste0(github, "Rdata/virginia.filt.Rda"))
load(paste0(github, "Rdata/virginia.tax.vec.Rda"))
load(paste0(github, "Rdata/new.both.filt.Rda"))
load(paste0(github, "Rdata/new.tax.vec.Rda"))

# remove two janky samples from london datasets and check for zero sum rows
#   - v.001A - close to 100% L. gasseri with practically no BV organisms
#   - v.019A - around 80 % iners, ~ 20% crispatus and a tiny bit of Gardnerella
filt.ls.loneur <- filt.ls.loneur[,-c(9,22)]
any(rowSums(filt.ls.loneur) == 0) # FALSE

# get grouping vectors for london/europe and virginia
loneur.conds <- c(rep('H',8), rep('B',12), rep('B',14), rep('H', 8)) 
load(paste0(github, "Rdata/virginia.groups.Rda"))

# load functional lookup table (if it doesn't exist, go to final section to
# make from scratch!)
load(paste0(github, "Rdata/lookup.func.Rda"))

# load species colour table
lookup.cols <- read.table(paste0(github,"Rdata/species_colors_updated.txt"),
                          header = F, sep = "\t", quote = "", row.names = 1)

################################### lon/eur ###################################

# zero replacement and CLR transformation of count data, then PCA
filt.ls.loneur.zr <- cmultRepl(t(filt.ls.loneur), label = 0,
                               method = "CZM", z.warning = 0.95)

filt.ls.loneur.clr <- apply(t(filt.ls.loneur.zr), 2, function(x){log(x) - mean(log(x))})

pca.ls.loneur <- prcomp(t(filt.ls.loneur.clr))

# calculate number of genes per spp
genes.le <- vector()
for(i in levels(factor(tax.vec.loneur.ls))){
  genes.le[i] <- length(which(tax.vec.loneur.ls == i))
}

# @@@@@@@@@@@@@@@@@@@@@@ for reviewer 1, re: comment # 5 @@@@@@@@@@@@@@@@@@@@

# calculate total reads for species represented by <75 genes
genes.le.ind <- list()
for(i in levels(factor(tax.vec.loneur.ls))){
  genes.le.ind[[i]] <- which(tax.vec.loneur.ls == i)
}

genes.le.sum <- vector()
for(i in names(which(genes.le <75))){
  genes.le.sum[i] <- sum(rowSums(filt.ls.loneur[genes.le.ind[[i]],]))
}

# express these read counts as a proportion of all reads (most abundant spp.
# is  ~0.15% of all reads in this filtered (10%) dataset; total is ~0.75% of the
# filtered dataset)
(genes.le.sum/sum(rowSums(filt.ls.loneur)))*100

# check sum of all genes across samples added by reducing the threshold from
# 30% to 10%
newgenes.le <- rowSums(filt.ls.loneur[setdiff(rownames(filt.ls.loneur),
                                              rownames(new.both.filt)),])
# get quantiles of these row sums
quantile(newgenes.le, probs = seq(0,1,0.05))

# check top 20 most abundant genes
tail(sort(newgenes.le), n = 20)

# top gene averages out to 2012 reads per sample
84524/42

# @@@@@@@@@@@@@@@@@@@@@@ for reviewer 1, re: comment # 5 @@@@@@@@@@@@@@@@@@@@

# count species represented by >75 genes
genes.le.75 <- names(which(genes.le >=75))

# get indices of genes corresponding to these spp. and add "other"
ind.le.75 <- list()
for(i in genes.le.75){
  ind.le.75[[i]] <- which(tax.vec.loneur.ls == i)
}

ind.le.75[["Other"]] <- setdiff(1:nrow(filt.ls.loneur),
                                unlist(ind.le.75))

cols.le <- c(lookup.cols[genes.le.75,1], rgb(0,0,0,0.05))

# png(paste0("~/filter_lon_eur_summarySpecies_biplot.png"),
#     units = "in", height=9, width=12,res = 400)

codaSeq.PCAplot(pca.ls.loneur, plot.groups = FALSE, plot.loadings = TRUE, 
                plot.ellipses = NULL, plot.density = "loadings",
                load.grp = ind.le.75, load.col = cols.le, load.sym = 19,
                load.cex = 0.5, PC = c(1,2), plot.legend = "loadings",
                leg.position = "bottomright", leg.columns = 3, leg.cex = 0.65,
                title = "London & Europe: filter = 10% prev. / 0.005% abund.")

# dev.off()

################################### virginia ###################################

# zero replacement and CLR transformation of count data, then PCA
filt.ls.virg.zr <- cmultRepl(t(filt.ls.virg), label = 0,
                               method = "CZM", z.warning = 0.95)

filt.ls.virg.clr <- apply(t(filt.ls.virg.zr), 2, function(x){log(x) - mean(log(x))})

pca.ls.virg <- prcomp(t(filt.ls.virg.clr))

# calculate number of genes per spp
genes.virg <- vector()
for(i in levels(factor(tax.vec.virg.ls))){
  genes.virg[i] <- length(which(tax.vec.virg.ls == i))
}

# @@@@@@@@@@@@@@@@@@@@@@ for reviewer 1, re: comment # 5 @@@@@@@@@@@@@@@@@@@@

# calculate total reads for species represented by <75 genes
genes.virg.ind <- list()
for(i in levels(factor(tax.vec.virg.ls))){
  genes.virg.ind[[i]] <- which(tax.vec.virg.ls == i)
}

genes.virg.sum <- vector()
for(i in names(which(genes.virg <75))){
  genes.virg.sum[i] <- sum(rowSums(filt.ls.virg[genes.virg.ind[[i]],]))
}

# express these read counts as a proportion of all reads (most abundant spp.
# is  ~0.15% of all reads in this filtered (10%) dataset; total is ~0.75% of the
# filtered dataset)
(genes.virg.sum/sum(rowSums(filt.ls.virg)))*100

# check sum of all genes across samples added by reducing the threshold from
# 30% to 10%
newgenes.virg <- rowSums(filt.ls.virg[setdiff(rownames(filt.ls.virg),
                                              rownames(virginia.filt)),])
# get quantiles of these row sums
quantile(newgenes.virg, probs = seq(0,1,0.05))

# check top 20 most abundant genes
tail(sort(newgenes.virg), n = 20)

# top gene averages out to 3686 reads per sample
1094605/297

# @@@@@@@@@@@@@@@@@@@@@@ for reviewer 1, re: comment # 5 @@@@@@@@@@@@@@@@@@@@

# count species represented by >75 genes
genes.virg.75 <- names(which(genes.virg >=75))

# get indices of genes corresponding to these spp. and add "other"
ind.virg.75 <- list()
for(i in genes.virg.75){
  ind.virg.75[[i]] <- which(tax.vec.virg.ls == i)
}

ind.virg.75[["Other"]] <- setdiff(1:nrow(filt.ls.virg),
                                unlist(ind.virg.75))

cols.virg <- c(lookup.cols[genes.virg.75,1], rgb(0,0,0,0.05))

# 33 spp., of which, 13 don't have colour. Add these in now
cols.virg[c(1,4,6:8,10,13,19,20,25,28,32,33)] <- c("deeppink4","gold4","cyan3",
                                                   "cyan","forestgreen","blue",
                                                   "dodgerblue","palevioletred1",
                                                   "chocolate4","purple2",
                                                   "deepskyblue3","sandybrown",
                                                   "gold")

# png(paste0("~/filter_virginia_summarySpecies_biplot.png"),
#     units = "in", height=9, width=12,res = 400)

codaSeq.PCAplot(pca.ls.virg, plot.groups = FALSE, plot.loadings = TRUE, 
                plot.ellipses = NULL, plot.density = "loadings", grp.cex = 0.2,
                load.grp = ind.virg.75, load.col = cols.virg, load.sym = 19,
                load.cex = 0.2, PC = c(1,2), plot.legend = "loadings",
                leg.position = "bottom", leg.columns = 6, leg.cex = 0.5,
                title = "Virginia: filter = 10% prev. / 0.005% abund.")

# dev.off()

# also want a KO-aggregated biplot as this was the corresponding figure
lookup.ko <- read.table(paste0(github, "1_VIRGO/8.A.kegg.ortholog.txt"),
                        header = F, sep = "\t", quote = "", row.names = 1)

virg.ko.vec <- lookup.ko[rownames(filt.ls.virg),1]
names(virg.ko.vec) <- rownames(filt.ls.virg)
ko.ls.virg <- aggregate(filt.ls.virg, by = list(virg.ko.vec), FUN = sum)
rownames(ko.ls.virg) <- ko.ls.virg$Group.1
ko.ls.virg$Group.1 <- NULL

# clr transform and check that HK genes are centred
set.seed(2023)

mu.matrix <- aldex.makeScaleMatrix(gamma = 0.5, mu = c(1,1.15),
                                   conditions = virginia.groups$groups.2)

ko.ls.virg.clr <- aldex.clr(reads = ko.ls.virg, conds = virginia.groups$groups.2,
                            gamma = mu.matrix, verbose = TRUE, mc.samples = 128)

ko.ls.virg.clr.e <- aldex.effect(ko.ls.virg.clr, verbose = TRUE,
                                 include.sample.summary = TRUE)

ko.ls.virg.clr.t <- aldex.ttest(ko.ls.virg.clr, verbose = TRUE)

ko.ls.virg.clr.all <- cbind(ko.ls.virg.clr.e, ko.ls.virg.clr.t)

# subset ko-> pathway lookup for KOs in virginia dataset (remove eukaryotic
# first!)
ko.ls.virg <- ko.ls.virg[-which(grepl(paste("K03364","K13963","K01173",
                                                "K12373","K14327","K00863",
                                                "K00599","K13993","K00811",
                                                "K03260","K00985", sep = "|"),
                                          rownames(ko.ls.virg))),]
lookup.func.virg <- lookup.func[rownames(ko.ls.virg), , drop = F]

# sanity check that order is equal, then pull housekeeping genes
length(which(rownames(lookup.func.virg)==rownames(ko.ls.virg)))
path.ribo <- rownames(ko.ls.virg)[which(lookup.func.virg$pathway == "Ribosome")]
path.trna <- rownames(ko.ls.virg)[which(lookup.func.virg$pathway == "Aminoacyl-tRNA biosynthesis")]
path.glyc <- rownames(ko.ls.virg)[which(lookup.func.virg$pathway == "Glycolysis / Gluconeogenesis")]

aldex.plot(ko.ls.virg.clr.all, xlim=c(0.346,9), cutoff.pval = 0.01)
title('ALDEx2 w/ ScaleSim (Virginia filtered K0): mu = 1 : 1.15, gamma = 0.5, p <0.01',
      adj=0, line= 0.8)

points(ko.ls.virg.clr.all[path.ribo ,"diff.win"],
       ko.ls.virg.clr.all[path.ribo ,"diff.btw"],
       pch= 19, col= "blue", cex= 0.5)

points(ko.ls.virg.clr.all[path.trna,"diff.win"],
       ko.ls.virg.clr.all[path.trna,"diff.btw"],
       pch= 19, col= "dodgerblue", cex= 0.5)

points(ko.ls.virg.clr.all[path.glyc,"diff.win"],
       ko.ls.virg.clr.all[path.glyc,"diff.btw"],
       pch= 19, col= "lightskyblue1", cex= 0.5)

points(ko.ls.virg.clr.all[which(abs(ko.ls.virg.clr.all$effect) >1),"diff.win"],
       ko.ls.virg.clr.all[which(abs(ko.ls.virg.clr.all$effect) >1),"diff.btw"],
       pch= 19, col= "purple3", cex= 0.5)

abline(h = 0, col = "blue", lty = 2)

legend("bottomleft", pch = 19, cex = 0.75,
       col = c( "black", "red", "purple3", "blue3", "royalblue1", "lightskyblue1"),
       legend = c("Non-significant", "P <0.01", "Absolute effect size >1",
                  "Ribosome", "Aminoacyl-tRNA biosynthesis", "Glycolysis / Gluconeogenesis"))

# HK genes look nice and centred- great. Extract CLR values and do PCA
ko.ls.virg.clr.df <- ko.ls.virg.clr.all[,grep("rab.sample.", colnames(ko.ls.virg.clr.all))]
colnames(ko.ls.virg.clr.df) <- gsub("rab.sample.", "", colnames(ko.ls.virg.clr.df))

pca.ls.virg.ko <- prcomp(t(ko.ls.virg.clr.df))

# make list of pathways to highlight
highlight.path <- c("Ribosome", "Aminoacyl-tRNA biosynthesis",
                    "CAMP resistance","Flagellar assembly",
                    "Bacterial chemotaxis", "Porphyrin metabolism",
                    "Butanoate metabolism", "Two-component system",
                    "Starch & sucrose metabolism", "ABC transporters")

ind.path <- list()
for(i in levels(factor(highlight.path))){
  ind.path[[i]] <- which(lookup.func.virg$pathway == i)
}

# add all other pathways to list under "Other"
ind.path[["Other"]] <- setdiff(1:nrow(ko.ls.virg.clr.df),
                                      unlist(ind.path))

# make vector of pathway colours ("Other" as transparent black)
cols.path <- c("olivedrab2", "lightskyblue1", "tomato3", "gold",
                        "chocolate4", "purple3", "bisque2", "black",
                        "deeppink3", "cyan3", rgb(0,0,0,0.05))

# make list of group indices
ind.grp <- list()
for(i in levels(factor(virginia.groups$groups.3))){
  ind.grp[[i]] <- which(virginia.groups$groups.3 == i)
}

# make vector of group colours
cols.virg.grp= c("goldenrod2", "steelblue1", "royalblue3")

# plot PCA for virginia metatranscriptome functions
# png(paste0("~/filter_virginia_summaryKO_biplot.png"),
#     units = "in", height=6, width=12,res = 400)

codaSeq.PCAplot(pca.ls.virg.ko, plot.groups = TRUE, plot.loadings = TRUE,
                plot.ellipses = NULL, plot.density = NULL, PC = c(1,2),
                grp = ind.grp, grp.col = cols.virg.grp, grp.sym = "text",
                grp.cex = 0.3, load.grp = ind.path, load.col = cols.path, 
                load.sym = rep(19,11), load.cex = 0.5, plot.legend = "loadings",
                leg.position = "top", leg.cex = 0.6, leg.columns = 4,
                title = "Virginia - K0 terms: Health vs. Molecular BV")

legend("bottomleft", legend = c("Healthy (Subgroup 1)", "Healthy (Subgroup 2)", "Molecular BV"),
       col = c("steelblue1", "royalblue3", "goldenrod2"), pch = "SS", cex = 0.6)

# dev.off()

############################## functional lookups ##############################

# # load in vnum -> KO and vnum -> pathway for l/e, merge and retain only
# # Vnums, KOs and pathways
# lookup.ko <- read.table(paste0(github, "1_VIRGO/8.A.kegg.ortholog.txt"),
#                         header = F, sep = "\t", quote = "")
# 
# lookup.path <- read.table(paste0(github, "1_VIRGO/8.C.kegg.pathway.txt"),
#                           header = F, sep = "\t", quote = "")
# 
# lookup.func <- left_join(x = lookup.ko, y = lookup.path, by = "V1")
# lookup.func <- lookup.func[,c(2,6)]
# colnames(lookup.func) <- c("ko", "pathway")
# 
# # replace NAs in pathway with 'Unknown'
# lookup.func$pathway <- case_when(is.na(lookup.func$pathway) ~ "Unknown",
#                                  .default = lookup.func$pathway)
# 
# # arrange by KO and remove duplicates
# lookup.func <- lookup.func[which(!duplicated(lookup.func$ko)),]
# rownames(lookup.func) <- lookup.func$ko
# lookup.func <- lookup.func %>%
#   arrange(ko) %>%
#   dplyr::select(pathway)
# 
# # remove eukaryotic KOs (discovered during curation of 'Unknown' pathways and
# # reflects changes to KO aggregated dataframe)
# lookup.func <- lookup.func[-which(grepl(paste("K03364","K13963","K01173",
#                                                 "K12373","K14327","K00863",
#                                                 "K00599","K13993","K00811",
#                                                 "K03260","K00985", sep = "|"),
#                                           rownames(lookup.func))), ,
#                              drop = FALSE]
# 
# # use info from lon/eur and virginia filtering to replace unknowns where
# # applicable
# lookup.func$pathway <- case_when(rownames(lookup.func) %in% c("K01580","K01915","K00812","K11358","K01775",
#                                                                 "K01776") ~ "Alanine aspartate and glutamate metabolism",
#                                   rownames(lookup.func) == "K01710" ~ "Amino sugar metabolism",
#                                   rownames(lookup.func) %in% c("K01358","K03544") ~ "ATP-dependent proteolysis",
#                                   rownames(lookup.func) == "K03648" ~ "Base excision repair",
#                                   rownames(lookup.func) %in% c("K01715","K14534") ~ "Butanoate metabolism",
#                                   rownames(lookup.func) %in% c("K03367","K03739","K03740","K14188","K14205") ~ "CAMP resistance",
#                                   rownames(lookup.func) %in% c("K03531","K03588","K03589","K03590") ~ "Cell division",
#                                   rownames(lookup.func) %in% c("K00031","K01679","K00029","K00024","K00174",
#                                                                 "K00175","K00176","K00177","K00239","K00240",
#                                                                 "K00241","K00625","K00925","K01676","K01681",
#                                                                 "K01847","K01848","K01895","K01902","K01903",
#                                                                 "K01966","K03737","K05606","K13788","K13581") ~ "Citrate cycle (TCA cycle)",
#                                   rownames(lookup.func) == "K12234" ~ "Coenzyme F biosynthesis",
#                                   rownames(lookup.func) == "K11031" ~ "Cytolysis",
#                                   rownames(lookup.func) %in% c("K02313","K02314") ~ "DNA replication",
#                                   rownames(lookup.func) == "K01897" ~ "Fatty acid metabolism",
#                                   rownames(lookup.func) %in% c("K02406","K03092") ~ "Flagellar assembly",
#                                   rownames(lookup.func) %in% c("K01632","K07407") ~ "Galactose metabolism",
#                                   rownames(lookup.func) == "K00864" ~ "Glycerolipid metabolism",
#                                   rownames(lookup.func) %in% c("K00058","K00600","K01079") ~ "Glycine serine and threonine metabolism",
#                                   rownames(lookup.func) %in% c("K00134", "K13953","K01835","K00873","K00927",
#                                                                 "K01006","K01610","K01624","K01803","K02446",
#                                                                 "K04041","K01007","K00850","K01834","") ~ "Glycolysis / Gluconeogenesis",
#                                   rownames(lookup.func) == "K00018" ~ "Glyoxylate and dicarboxylate metabolism",
#                                   rownames(lookup.func) == "K00817" ~ "Histidine metabolism",
#                                   rownames(lookup.func) == "K01338" ~ "Lon protease",
#                                   rownames(lookup.func) == "K01582" ~ "Lysine metabolism",
#                                   rownames(lookup.func) == "K00123" ~ "Methanoate metabolism",
#                                   rownames(lookup.func) %in% c("K04079","K04043","K04077") ~ "Molecular chaperones",
#                                   rownames(lookup.func) == "K01354" ~ "Oligopeptidase activity",
#                                   rownames(lookup.func) %in% c("K00297","K01491","K01938") ~ "One carbon pool by folate",
#                                   rownames(lookup.func) %in% c("K02108","K02109","K02110","K02111","K02112",
#                                                                 "K02113","K02114","K02115","K02117","K02118",
#                                                                 "K02119","K02120","K02122","K02123","K02124") ~ "Oxidative phosphorylation",
#                                   rownames(lookup.func) == "K01195" ~ "Pentose and glucuronate interconversions",
#                                   rownames(lookup.func) %in% c("K01621","K01783","K01807","K01808","K08094",
#                                                                 "K00615") ~ "Pentose phosphate pathway",
#                                   rownames(lookup.func) == "K02563" ~ "Peptidoglycan biosynthesis",
#                                   rownames(lookup.func) %in% c("K00832","K00483","K01451") ~ "Phenylalanine tyrosine and tryptophan biosynthesis",
#                                   rownames(lookup.func) == "K00889" ~ "Phosphatidylinositol signaling system",
#                                   rownames(lookup.func) %in% c("K00088","K00857","K01951","K01488") ~ "Purine metabolism",
#                                   rownames(lookup.func) %in% c("K00757","K00760","K00876") ~ "Pyrimidine metabolism",
#                                   rownames(lookup.func) %in% c("K01759", "K01759","K01595") ~ "Pyruvate metabolism",
#                                   rownames(lookup.func) == "K11749" ~ "Regulated intramembrane proteolysis",
#                                   rownames(lookup.func) %in% c("K13288","K03685") ~ "RNA degradation",
#                                   rownames(lookup.func) == "K04564" ~ "ROS degradation",
#                                   rownames(lookup.func) == "K01186" ~ "Sialidase activity",
#                                   rownames(lookup.func) == "K00688" ~ "Starch and sucrose metabolism",
#                                   rownames(lookup.func) %in% c("K03313","K03315") ~ "Sodium-proton antiport",
#                                   rownames(lookup.func) == "K00869" ~ "Terpenoid backbone biosynthesis",
#                                   rownames(lookup.func) == "K01078" ~ "Thiamine metabolism",
#                                   rownames(lookup.func) == "K02358" ~ "Translation",
#                                   rownames(lookup.func) == "K02040" ~ "Two-component system",
#                                   rownames(lookup.func) == "K08303" ~ "Type I collagen degradation",
#                                   rownames(lookup.func) %in% c("K00568","K01661","K02548","K03183","K03688") ~ "Ubiquinone biosynthesis",
#                                   rownames(lookup.func) == "K03191" ~ "Urea transport",
#                                   rownames(lookup.func) == "K01428" ~ "Urease activity",
#                                   rownames(lookup.func) %in% c("K01066","K00532") ~ "Unknown",
#                                   rownames(lookup.func) == "K01703" ~ "Valine leucine and isoleucine biosynthesis",
#                                   rownames(lookup.func) == "K02361" ~ "Vitamin K2 biosynthesis",
#                                   .default = lookup.func$pathway)
# 
# # abbreviate long pathways to shorter versions
# lookup.func$pathway <- case_when(lookup.func$pathway == "Alanine aspartate and glutamate metabolism" ~ "Ala-Asp-Glu metabolism",
#                                   lookup.func$pathway == "Amino sugar and nucleotide sugar metabolism" ~ "Amino sugar metabolism",
#                                   lookup.func$pathway == "Arginine and proline metabolism" ~ "Arg-Pro metabolism",
#                                   lookup.func$pathway == "Biosynthesis of unsaturated fatty acids" ~ "Fatty acid biosynthesis",
#                                   lookup.func$pathway == "Cysteine and methionine metabolism" ~ "Cys-Met metabolism",
#                                   lookup.func$pathway == "Glycine serine and threonine metabolism" ~ "Gly-Ser-Thr metabolism",
#                                   lookup.func$pathway == "Histidine metabolism" ~ "His metabolism",
#                                   lookup.func$pathway == "Lipopolysaccharide biosynthesis" ~ "LPS biosynthesis",
#                                   lookup.func$pathway == "Lysine biosynthesis" ~ "Lys metabolism",
#                                   lookup.func$pathway == "Lysine degradation" ~ "Lys metabolism",
#                                   lookup.func$pathway == "Lysine metabolism" ~ "Lys metabolism",
#                                   lookup.func$pathway == "Nicotinate and nicotinamide metabolism" ~ "Vitamin B3 metabolism",
#                                   lookup.func$pathway == "Phenylalanine tyrosine and tryptophan biosynthesis" ~ "Phe-Tyr-Trp metabolism",
#                                   lookup.func$pathway == "Phosphotransferase system (PTS)" ~ "Phosphotransferase system",
#                                   lookup.func$pathway == "Porphyrin and chlorophyll metabolism" ~ "Porphyrin metabolism",
#                                   lookup.func$pathway == "Terpenoid backbone biosynthesis" ~ "Terpenoid biosynthesis",
#                                   lookup.func$pathway == "Valine leucine and isoleucine biosynthesis" ~ "Val-Leu-Ile metabolism",
#                                   .default = lookup.func$pathway)
# 
# # assign formerly 'Unknown' KOs a function based on current KEGG database and
# # literature (see end of script for references)
# lookup.func$pathway <- case_when(rownames(lookup.func) == "K03281" ~ "Acid survival",
#                                   rownames(lookup.func) == "K03320" ~ "Ammonium transport",
#                                   rownames(lookup.func) == "K00561" ~ "Antimicrobial resistance",
#                                   rownames(lookup.func) == "K03402" ~ "Arg-Pro metabolism",
#                                   rownames(lookup.func) == "K03758" ~ "Arginine uptake",
#                                   rownames(lookup.func) %in% c("K01419","K03695","K03696","K03697") ~ "ATP-dependent proteolysis",
#                                   rownames(lookup.func) %in% c("K02477","K02478") ~ "Bacteriocin production",
#                                   rownames(lookup.func) == "K01448" ~ "CAMP resistance",
#                                   rownames(lookup.func) == "K03798" ~ "Cell division",
#                                   rownames(lookup.func) %in% c("K07012","K09952") ~ "CRISPR/Cas system",
#                                   rownames(lookup.func) %in% c("K01361","K01364") ~ "Cytokine degradation",
#                                   rownames(lookup.func) %in% c("K06189","K11068") ~ "Cytolysis",
#                                   rownames(lookup.func) %in% c("K03631","K04485","K01356") ~ "DNA damage repair",
#                                   rownames(lookup.func) %in% c("K02346","K02469","K02621","K02622","K03168",
#                                                                 "K03169","K03346") ~ "DNA replication",
#                                   rownames(lookup.func) == "K11936" ~ "Exopolysaccharide biosynthesis",
#                                   rownames(lookup.func) %in% c("K02237","K02238","K02240","K02242","K02243",
#                                                                 "K02244","K02248") ~ "Extracellular DNA uptake",
#                                   rownames(lookup.func) %in% c("K00777", "K02395","K02481","K06603","K13626") ~ "Flagellar assembly",
#                                   rownames(lookup.func) %in% c("K01795","K02444","K03436") ~ "Fructose and mannose metabolism",
#                                   rownames(lookup.func) %in% c("K02429","K02431","K07046") ~ "Fucose metabolism",
#                                   rownames(lookup.func) %in% c("K01854","K02082","K02529") ~ "Galactose metabolism",
#                                   rownames(lookup.func) == "K02445" ~ "Glycerol uptake",
#                                   rownames(lookup.func) == "K07029" ~ "Glycerolipid metabolism",
#                                   rownames(lookup.func) == "K02437" ~ "Gly-Ser-Thr metabolism",
#                                   rownames(lookup.func) %in% c("K02012","K02217","K07243") ~ "Iron acquisition",
#                                   rownames(lookup.func) == "K03303" ~ "Lactate metabolism",
#                                   rownames(lookup.func) %in% c("K04744","K06041") ~ "LPS biosynthesis",
#                                   rownames(lookup.func) %in% c("K03284","K06213","K07507") ~ "Magnesium transport",
#                                   rownames(lookup.func) %in% c("K03322","K12950") ~ "Manganese transport",
#                                   rownames(lookup.func) %in% c("K03686","K04078") ~ "Molecular chaperones",
#                                   rownames(lookup.func) == "K03297" ~ "Multidrug resistance - efflux pump",
#                                   rownames(lookup.func) %in% c("K03282","K03442") ~ "Osmotic stress adaptation",
#                                   rownames(lookup.func) %in% c("K02279","K02282","K02283","K02652") ~ "Pilus assembly",
#                                   rownames(lookup.func) == "K01991" ~ "Polysaccharide transport",
#                                   rownames(lookup.func) == "K01529" ~ "Purine metabolism",
#                                   rownames(lookup.func) == "K03445" ~ "Purine efflux",
#                                   rownames(lookup.func) == "K03756" ~ "Putrescine transport",
#                                   rownames(lookup.func) == "K02824" ~ "Pyrimidine metabolism",
#                                   rownames(lookup.func) %in% c("K01153","K01154","K01155","K01156","K03427") ~ "Restriction endonuclease",
#                                   rownames(lookup.func) %in% c("K01200", "K01215","K02438","K05992","K06896","K00701") ~ "Starch and sucrose metabolism",
#                                   rownames(lookup.func) == "K05946" ~ "Teichoic acid biosynthesis",
#                                   rownames(lookup.func) %in% c("K02483","K02484") ~ "Two-component system",
#                                   rownames(lookup.func) %in% c("K03187","K03188") ~ "Urease activity",
#                                   .default = lookup.func$pathway)
# 
# lookup.func$pathway <- case_when(lookup.func$pathway == "Adipocytokine signaling pathway" ~ "Fatty acid biosynthesis",
#                                            lookup.func$pathway == "African trypanosomiasis" ~ "Prolyl oligopeptidase",
#                                            lookup.func$pathway == "Alanine aspartate and glutamate metabolism" ~ "Ala-Asp-Glu metabolism",
#                                            lookup.func$pathway == "Amino sugar and nucleotide sugar metabolism" ~ "Amino sugar metabolism",
#                                            lookup.func$pathway == "Arginine and proline metabolism" ~ "Arg-Pro metabolism",
#                                            lookup.func$pathway == "Bacterial invasion of epithelial cells" ~ "Intracellular invasion",
#                                            lookup.func$pathway == "Biosynthesis of vancomycin group antibiotics" ~ "Amino sugar metabolism",
#                                            lookup.func$pathway == "Biosynthesis of unsaturated fatty acids" ~ "Fatty acid biosynthesis",
#                                            lookup.func$pathway == "Cationic antimicrobial peptide (CAMP) resistance" ~ "CAMP resistance",
#                                            lookup.func$pathway == "Cysteine and methionine metabolism" ~ "Cys-Met metabolism",
#                                            lookup.func$pathway == "Drug metabolism - cytochrome P450" ~ "Glycolysis / Gluconeogenesis",
#                                            lookup.func$pathway == "Drug metabolism - other enzymes" ~ "Drug metabolism",
#                                            lookup.func$pathway == "Fructose and mannose metabolism" ~ "Glycolysis / Gluconeogenesis",
#                                            lookup.func$pathway == "Glycerolipid metabolism" ~ "Glycerophospholipid metabolism",
#                                            lookup.func$pathway == "Glycine serine and threonine metabolism" ~ "Gly-Ser-Thr metabolism",
#                                            lookup.func$pathway == "Glycosphingolipid biosynthesis - globo series" ~ "Glycosphingolipid biosynthesis",
#                                            lookup.func$pathway == "Histidine metabolism" ~ "His metabolism",
#                                            lookup.func$pathway == "Lysine biosynthesis" ~ "Lys metabolism",
#                                            lookup.func$pathway == "Lysine degradation" ~ "Lys metabolism",
#                                            lookup.func$pathway == "MAPK signaling pathway - yeast" ~ "Pyruvate metabolism",
#                                            lookup.func$pathway == "Nicotinate and nicotinamide metabolism" ~ "Vitamin B3 metabolism",
#                                            lookup.func$pathway == "Phenylalanine metabolism" ~ "Phe-Tyr-Trp metabolism",
#                                            lookup.func$pathway == "Phenylalanine tyrosine and tryptophan biosynthesis" ~ "Phe-Tyr-Trp metabolism",
#                                            lookup.func$pathway == "Phosphotransferase system (PTS)" ~ "Phosphotransferase system",
#                                            lookup.func$pathway == "Photosynthesis" ~ "Oxidative phosphorylation",
#                                            lookup.func$pathway == "Porphyrin and chlorophyll metabolism" ~ "Porphyrin metabolism",
#                                            lookup.func$pathway == "Ribosome biogenesis in eukaryotes" ~ "Oligoribonuclease activity",
#                                            lookup.func$pathway == "Terpenoid backbone biosynthesis" ~ "Terpenoid biosynthesis",
#                                            lookup.func$pathway == "Valine leucine and isoleucine biosynthesis" ~ "Val-Leu-Ile metabolism",
#                                            lookup.func$pathway == "Valine leucine and isoleucine degradation" ~ "Val-Leu-Ile metabolism",
#                                            .default = lookup.func$pathway)
# 
# # further editing of specific KO terms to correct pathway (searched on KEGG)
# lookup.func$pathway <-case_when(rownames(lookup.func) == "K01580" ~ "ABC transporters",
#                                           rownames(lookup.func) == "K01580" ~ "Acid survival",
#                                           rownames(lookup.func) %in% c("K00812","K11358","K01915",
#                                                                                  "K01775","K01776","K00811") ~ "Ala-Asp-Glu metabolism",
#                                           rownames(lookup.func) == "K14195" ~ "Bacterial adherence",
#                                           rownames(lookup.func) == "K06595" ~ "Bacterial chemotaxis",
#                                           rownames(lookup.func) == "K11936" ~ "Biofilm formation",
#                                           rownames(lookup.func) == "K14534" ~ "Butanoate metabolism",
#                                           rownames(lookup.func) %in% c("K03367","K03739","K03740","K14188","K14205") ~ "CAMP resistance",
#                                           rownames(lookup.func) %in% c("K03531","K03588","K03589","K03590","K13581") ~ "Cell division",
#                                           rownames(lookup.func) %in% c("K00024","K00174","K00175","K00239","K00240",
#                                                                                  "K00241","K01676","K01681","K01847","K01848",
#                                                                                  "K01849","K01902","K01903","K01966","K03737",
#                                                                                  "K05606","K13788","K00031","K00177","K00625",
#                                                                                  "K00925","K01595","K01895","K05884","K01679") ~ "Citrate cycle (TCA cycle)",
#                                           rownames(lookup.func) %in% c("K01358","K03544") ~ "Clp-dependent proteolysis",
#                                           rownames(lookup.func) == "K11031" ~ "Cytolysis",
#                                           rownames(lookup.func) %in% c("K02313","K02314") ~ "DNA replication",
#                                           rownames(lookup.func) == "K00864" ~ "Glycerolipid metabolism",
#                                           rownames(lookup.func) == "K00981" ~ "Glycerophospholipid metabolism",
#                                           rownames(lookup.func) %in% c("K02446","K04041","K00134","K01803","K00873",
#                                                                                  "K00927","K01610","K01624","K00850","K01834","") ~ "Glycolysis / Gluconeogenesis",
#                                           rownames(lookup.func) == "K00018" ~ "Glyoxylate and dicarboxylate metabolism",
#                                           rownames(lookup.func) %in% c("K00058","K00600","K01079") ~ "Gly-Ser-Thr metabolism",
#                                           rownames(lookup.func) == "K01715" ~ "Fatty acid biosynthesis",
#                                           rownames(lookup.func) %in% c("K02405", "K02406", "K03092","K06603","K13626") ~ "Flagellar assembly",
#                                           rownames(lookup.func) == "K00817" ~ "His metabolism",
#                                           rownames(lookup.func) %in% c("K11749","K02012","K02217","K07243") ~ "Iron acquisition",
#                                           rownames(lookup.func) == "K01338" ~ "Lon protease",
#                                           rownames(lookup.func) == "K01582" ~ "Lys metabolism",
#                                           rownames(lookup.func) %in% c("K05299","K00123","K00124","K00127") ~ "Methanoate metabolism",
#                                           rownames(lookup.func) %in% c("K04079","K04043","K04077") ~ "Molecular chaperones",
#                                           rownames(lookup.func) %in% c("K00297","K01491","K01938") ~ "One carbon pool by folate",
#                                           rownames(lookup.func) %in% c("K02117","K02118","K02119","K02120","K02121",
#                                                                                  "K02122","K02123","K02124") ~ "Oxidative phosphorylation",
#                                           rownames(lookup.func) == "K00832" ~ "Phe-Tyr-Trp metabolism",
#                                           rownames(lookup.func) == "K01195" ~ "Pentose and glucuronate interconversions",
#                                           rownames(lookup.func) %in% c("K01621","K01632","K01783","K01807",
#                                                                                  "K01808","K00615") ~ "Pentose phosphate pathway",
#                                           rownames(lookup.func) == "K02563" ~ "Peptidoglycan biosynthesis",
#                                           rownames(lookup.func) == "K01488" ~ "Purine metabolism",
#                                           rownames(lookup.func) == "K00757" ~ "Pyrimidine metabolism",
#                                           rownames(lookup.func) %in% c("K00029","K01006","K01007") ~ "Pyruvate metabolism",
#                                           rownames(lookup.func) == "K01661" ~ "Quinone biosynthesis",
#                                           rownames(lookup.func) %in% c("K14051","K07173") ~ "Quorum sensing",
#                                           rownames(lookup.func) == "K01186" ~ "Sialidase activity",
#                                           rownames(lookup.func) %in% c("K03313","K03315") ~ "Sodium-proton antiport",
#                                           rownames(lookup.func) == "K00688" ~ "Starch and sucrose metabolism",
#                                           rownames(lookup.func) == "K04564" ~ "ROS degradation",
#                                           rownames(lookup.func) == "K00869" ~ "Terpenoid biosynthesis",
#                                           rownames(lookup.func) == "K02358" ~ "Translation",
#                                           rownames(lookup.func) == "K02040" ~ "Two-component system",
#                                           rownames(lookup.func) == "K08303" ~ "Type I collagen degradation",
#                                           rownames(lookup.func) == "K00532" ~ "Unknown",
#                                           rownames(lookup.func) == "K01428" ~ "Urease activity",
#                                           rownames(lookup.func) %in% c("K03183","K03688") ~ "Ubiquinone biosynthesis",
#                                           rownames(lookup.func) %in% c("K02548","K02548") ~ "Vitamin K2 biosynthesis",
#                                           .default = lookup.func$pathway)
# 
# # even further editing of specific KO terms to correct pathway (discovered
# # after analysing data downstream- these were all 'Unknown'; see end of file)
# lookup.func$pathway <- case_when(rownames(lookup.func) %in% c("K00435","K01772") ~ "Porphyrin metabolism",
#                                            rownames(lookup.func) == "K07777" ~ "Two-component system",
#                                            rownames(lookup.func) == "K01005" ~ "Teichoic acid biosynthesis",
#                                            rownames(lookup.func) == "K00903" ~ "Exopolysaccharide biosynthesis",
#                                            rownames(lookup.func) %in% c("K01356","K03631","K03700") ~ "DNA damage repair",
#                                            rownames(lookup.func) == "K02055" ~ "Putrescine/Spermidine transport",
#                                            rownames(lookup.func) %in% c("K02068","K02069") ~ "Iron homeostasis",
#                                            rownames(lookup.func) %in% c("K02237","K02244","K02248") ~ "Extracellular DNA uptake",
#                                            rownames(lookup.func) %in% c("K02444","K03436") ~ "Fructose metabolism",
#                                            rownames(lookup.func) %in% c("K03282","K03442") ~ "Osmotic stress adaptation",
#                                            rownames(lookup.func) == "K02530" ~ "Lactose metabolism",
#                                            rownames(lookup.func) == "K02859" ~ "Riboflavin metabolism",
#                                            rownames(lookup.func) == "K03284" ~ "Magnesium transport",
#                                            rownames(lookup.func) %in% c("K03297","K08161") ~ "Multidrug resistance - efflux pump",
#                                            rownames(lookup.func) == "K03303" ~ "Lactate metabolism",
#                                            rownames(lookup.func) == "K03322" ~ "Manganese acquisition",
#                                            rownames(lookup.func) %in% c("K03484", "K05992") ~ "Starch and sucrose metabolism",
#                                            rownames(lookup.func) == "K04075" ~ "Aminoacyl-tRNA biosynthesis",
#                                            rownames(lookup.func) == "K07043" ~ "Pyrimidine metabolism",
#                                            rownames(lookup.func) == "K01595" ~ "Pyruvate metabolism",
#                                            rownames(lookup.func) == "K01066" ~ "Unknown",
#                                            .default = lookup.func$pathway)
# 
# lookup.func$pathway <- gsub(pattern = "and", replacement = "&",
#                                       lookup.func$pathway)
# 
# # identify pathways that can be abbreviated in line with the above modifications
# levels(factor(lookup.func$pathway))
# 
# lookup.func$pathway <- case_when(lookup.func$pathway == "Alanine, aspartate & glutamate metabolism" ~ "Ala-Asp-Glu metabolism",
#                                  lookup.func$pathway == "Biosynthesis of 12-, 14- & 16-membered macrolides" ~ "Macrolide biosynthesis",
#                                  lookup.func$pathway == "Fatty acid biosynthesis" ~ "Fatty acid metabolism",
#                                  lookup.func$pathway == "Glycine, serine & threonine metabolism" ~ "Gly-Ser-Thr metabolism",
#                                  lookup.func$pathway == "Glycosphingolipid biosynthesis - lacto & neolacto series" ~ "Glycosphingolipid biosynthesis",
#                                  lookup.func$pathway == "Metabolism of xenobiotics by cytochrome P450" ~ "Xenobiotic metabolism",
#                                  lookup.func$pathway == "Phenylalanine, tyrosine & tryptophan biosynthesis" ~ "Phe-Tyr-Trp metabolism",
#                                  lookup.func$pathway == "Putrescine transport" ~ "Putrescine/Spermidine transport",
#                                  lookup.func$pathway == "Ubiquinone & other terpenoid-quinone biosynthesis" ~ "Ubiquinone biosynthesis",
#                                  lookup.func$pathway == "Valine, leucine & isoleucine biosynthesis" ~ "Val-Leu-Ile metabolism",
#                                  lookup.func$pathway == "Valine, leucine & isoleucine degradation" ~ "Val-Leu-Ile metabolism",
#                                  .default = lookup.func$pathway)
# 
# # save final KO -> pathway lookup
# save(lookup.func, file = paste0(github, "Rdata/lookup.func.Rda"))
