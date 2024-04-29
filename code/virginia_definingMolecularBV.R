# processing Virginia dataset & molecular definition of BV

#################################### setup ####################################

library(dplyr) # for data manipulation
library(pheatmap) # for heatmap plotting
library(viridisLite) # for colourblind-friendly palettes
library(dendextend) # for dendrogram plotting

user <- "sds/cc"
path.to.github <- switch(user,
                         "sds/cc" = "~/Documents/GitHub/dossantos2024study/",
                         "gg" = "~/Documents/0_git/projects/dossantos2024study/")

# run 'setup.R'
source(paste(path.to.github,"code/setup.R", sep = ""))

# load the .Rda object containing the unfiltered virginia feature table if it
# exists in 'Rdata/', OR read it in from a .zip archive (519418 genes in 297 
# samples), remove the trailing '_1 from every sample name' and save as .Rda.
# NOTE: unz()  requires the complete file path  of the .zip archive for
# 'description', but only the name for 'filename')
if(file.exists(paste(path.to.github, "Rdata/virginia.data.Rda", sep = ""))){
  load(paste(path.to.github,"Rdata/virginia.data.Rda", sep = ""))
} else{
  virginia.data <-read.table(unz(description = paste(path.to.github,
                                                     "Rdata/t2t_summary.virginia.NR.abundance.txt.zip",
                                                     sep = ""),
                                 filename = "t2t_summary.virginia.NR.abundance.txt"),
                             header = T, row.names = 1, quote = '', sep = "\t")
  
  colnames(virginia.data) <- gsub("_1", "", colnames(virginia.data))
  
  save(virginia.data,
       file = paste(path.to.github, "Rdata/virginia.data.Rda", sep = ""))
}

# load the.Rda object containing the filtered virginia feature table if it 
# exists in 'Rdata/', OR filter the virginia dataset using parameters used for
# london & europe datasets: (29059 genes in 297 samples) and save as .Rda
if(file.exists(paste(path.to.github, "Rdata/virginia.filt.Rda", sep = ""))){
  load(paste(path.to.github, "Rdata/virginia.filt.Rda", sep = ""))
} else{
  virginia.filt <- codaSeq.filter(virginia.data, samples.by.row = FALSE,
                                  min.occurrence=0.30, min.prop=0.00005)
  
  save(virginia.filt,
       file = paste(path.to.github, "Rdata/virginia.filt.Rda", sep = ""))
}

# load in KEGG pathway table (provided as part of VIRGO)
pathway.table <- read.table(paste(path.to.github,
                                  '1_VIRGO/8.C.kegg.pathway.copy.txt', sep=""),
                            sep="\t", header=T, row.names=1, fill=TRUE)

# load in virginia metadata (output of 'Rdata/processing_virginia_metadata.R')
# NOTE: This is one of the files omitted from the repo as it contains clinical 
#       metadata hosted on dbGaP.
load("[LOCAL-PATH-TO]/virginia.meta.Rda")

# sanity check on feature table colnames vs. metadata 'sra.id' for equality
length(which(colnames(virginia.data) == virginia.meta$sra_id))  # all TRUE

######################### heatmap: species composition #########################

# load in the .Rda object containing the taxa vector for the filtered virginia
# dataset if it exists in 'Rdata/', OR make it from scratch andsave as .Rda
if(file.exists(paste(path.to.github, "Rdata/tax.vec.virginia.Rda", sep = ""))){
  load(paste(path.to.github, "Rdata/tax.vec.virginia.Rda", sep = ""))
} else{
  tax.vec.virginia <- tax.table[rownames(virginia.filt),2]
  names(tax.vec.virginia) <- rownames(virginia.filt)
  save(tax.vec.virginia, file = paste(path.to.github,
                                      "Rdata/tax.vec.virginia.Rda", sep = "")) 
}

# get total number of reads per species
virginia.reads.sp <- list()
for(i in levels(factor(tax.vec.virginia))){
  virginia.reads.sp[[i]] <-colSums(virginia.filt[which(tax.vec.virginia==i),])
}

# count number of reads for gene with no taxonomy information in VIRGO
virginia.reads.sp[["No_tax_info"]] <- colSums(virginia.filt[is.na(tax.vec.virginia),])

# convert summary list of species read counts to data frame
virginia.reads.sp.df <- as.data.frame(do.call(rbind, virginia.reads.sp))

# count number of genes representing each species
virginia.n.genes <- list()
for(i in levels(factor(tax.vec.virginia))){
  virginia.n.genes[[i]] <- length(which(tax.vec.virginia == i))
}

# count number of genes with no taxonomy information in VIRGO
virginia.n.genes[["No_tax_info"]] <- length(which(is.na(tax.vec.virginia)))

# calculate total, mean, median, min and max reads for each species (from non-
# zero values), then summarise in data frame
reads.total <- rowSums(virginia.reads.sp.df)
reads.mean <- apply (virginia.reads.sp.df, 1, function(x) mean(x[x!=0]))
reads.median <- apply (virginia.reads.sp.df, 1, function(x) median(x[x!=0]))
reads.min <- apply (virginia.reads.sp.df, 1, function(x) min(x[x!=0]))
reads.max <-apply (virginia.reads.sp.df, 1, max)

virginia.read.summary <- as.data.frame(cbind(no.genes= unlist(virginia.n.genes),
                                             reads.total= reads.total,
                                             reads.mean= reads.mean,
                                             reads.median= reads.median,
                                             reads.min= reads.min,
                                             reads.max= reads.max))

# get vector of taxa represented by >75 genes
virginia.n.genes.75 <- names(which(virginia.n.genes >75))

# subset data frame for only taxa with >75 genes and sanity check rows/cols
virginia.reads.sp.df.75 <- virginia.reads.sp.df[virginia.n.genes.75,]
any(rowSums(virginia.reads.sp.df.75) == 0) # returns FALSE
any(colSums(virginia.reads.sp.df.75) == 0) # returns FALSE

# generate matrix for checking what % of reads are covered by species
# represented by various read cut-off values and >75 genes
virginia.top.sp <- matrix(ncol = ncol(virginia.filt), nrow = 7, data = NA)
colnames(virginia.top.sp) <- colnames(virginia.filt)
rownames(virginia.top.sp) <- c("total",paste0(rep("10^",5),c(3:7)),">75 genes")

# count total reads, reads present in dataset using various count thresholds,
# and reads from by species represented by >75 genes
virginia.top.sp[1,] <- colSums(virginia.reads.sp.df[])
virginia.top.sp[2,] <- colSums(virginia.reads.sp.df[which(rowSums(virginia.reads.sp.df) >10^3),])
virginia.top.sp[3,] <- colSums(virginia.reads.sp.df[which(rowSums(virginia.reads.sp.df) >10^4),])
virginia.top.sp[4,] <- colSums(virginia.reads.sp.df[which(rowSums(virginia.reads.sp.df) >10^5),])
virginia.top.sp[5,] <- colSums(virginia.reads.sp.df[which(rowSums(virginia.reads.sp.df) >10^6),])
virginia.top.sp[6,] <- colSums(virginia.reads.sp.df[which(rowSums(virginia.reads.sp.df) >10^7),])
virginia.top.sp[7,] <- colSums(virginia.reads.sp.df.75)

# convert summary matrix of read counts at various thresholds to proportions
virginia.top.sp.cut <- t(t(virginia.top.sp)/virginia.top.sp[1,])

# calculate number of samples retaining >95% of their reads when using >75 gene
# threshold
sort(virginia.top.sp.cut[7,],decreasing = TRUE) # only a handful with low prop
length(which(virginia.top.sp.cut[7,] >0.945)) # 284 samples >95 %

# ordinarily, I would remove these 13 samples with <95 % reads retained after
# filtering for taxa represented by >75 genes, but I need to know what CST they
# are (in order to include that info when aggregating by KO term for clustering)

# convert matrix of counts into proportions (excluding genes with no taxonomy 
# information) and save as .Rda 
virginia.filt.prop <- t(t(virginia.reads.sp.df.75[-18,])/colSums(virginia.reads.sp.df.75[-18,]))

# make vectors of metadata for vaginal pH, and vaginal discharge, odor, itching
# and whether patient had bv since last visit, replacing NA with "Unknown"
meta.ph <- virginia.meta$vaginal_ph

meta.discharge <- case_when(is.na(virginia.meta$current_vag_abnorm_discharge) ~ "Unknown",
                            .default = as.character(virginia.meta$current_vag_abnorm_discharge))

meta.odor <- case_when(is.na(virginia.meta$current_vag_bad_odor) ~ "Unknown",
                       .default = as.character(virginia.meta$current_vag_bad_odor))

meta.itching <- case_when(is.na(virginia.meta$current_vag_itching) ~ "Unknown",
                          .default = as.character(virginia.meta$current_vag_itching))

# cbind metadata rows for vaginal pH, and vaginal discharge, odor and itching
hm.meta <- as.data.frame(cbind(ph = meta.ph,
                               discharge = meta.discharge,
                               odor = meta.odor,
                               itching = meta.itching))

# set row names to sample IDs and edit column names for uppercase letters at
# start of word
rownames(hm.meta) <- colnames(virginia.filt)
colnames(hm.meta) <- c("Vaginal pH", "Discharge", "Odor", "Itching")

# replace underscores with spaces and decapitalise 'S' for 'Not_Sure'
hm.meta$Discharge <- gsub("_", " ", hm.meta$Discharge)
hm.meta$Discharge <- gsub("Sure", "sure", hm.meta$Discharge)

hm.meta$Odor <- gsub("_", " ", hm.meta$Odor)
hm.meta$Odor <- gsub("Sure", "sure", hm.meta$Odor)

hm.meta$Itching <- gsub("_", " ", hm.meta$Itching)
hm.meta$Itching <- gsub("Sure", "sure", hm.meta$Itching)

# get levels of factor objects from heatmap metadata
levels(factor(hm.meta$`Vaginal pH`))
levels(factor(hm.meta$Discharge))
levels(factor(hm.meta$Odor))
levels(factor(hm.meta$Itching))

# make list of colours for heatmap column bars using above factor levels
hm.column.cols <- list(`Vaginal pH` = c("4" = "red4", "5" = "red3", 
                                        "5.5" = "indianred",
                                        "6" = "lightskyblue1", 
                                        "6.5" = "steelblue3", 
                                        "7" = "royalblue3", "7.5" = "blue4",
                                        "Unknown" = "grey"),
                       Discharge = c("No" = "blue", `Not sure` = "magenta4", 
                                     "Unknown" = "grey", "Yes" = "red"),
                       Odor = c("No" = "blue", `Not sure` = "magenta4",
                                "Unknown" = "grey",  "Yes" = "red"),
                       Itching = c("No" = "blue", `Not sure` = "magenta4",
                                   "Unknown" = "grey", "Yes" = "red"))

# unsupervised clustering of species-level metatranscriptomes using Euclidean 
# distances and Ward's method
hm.dendrogram <- hclust(dist(t(virginia.filt.prop), method = "euclidean"),
                        method = "ward.D2")

# make vector of colourblind-friendly heatmap colours
hm.colours<-rev(magma(100))

# # make list of row names to be italicised
hm.row.italics<-lapply(gsub("_", " ", rownames(virginia.filt.prop)),
                       function(x){
                         bquote(italic(.(x)))
                       })

# edit list items for which italics should not be used
hm.row.italics[2]<-"BVAB1"
hm.row.italics[10]<-expression(italic("Megasphaera")~"genomosp.")

# plotting the heatmap at this point showed that all five 'canonical CSTs are 
# present: I (crispatus), II (gasseri) III, (iners), IV (BV organisms, mostly 
# G. vaginalis), and V (jensenii). Additionally, there are two further groups 
# consisting of a combintion of iners / BV organisms, and iners/jensenii (i.e.
# III / IV and III / V). To get this information from cutree, one needs to split
# into 12 groups, then edit groups as required

# plot the dendrogram and see what the clusters look like when altering 'k' in
# rect.hclust() (NOTE: as noted above, all five CSTs, plus two 'combination
# state types, are present. Need k= 12 to split these easily)
plot(hm.dendrogram, hang= -1, cex = 0.4)
rect.hclust(hm.dendrogram, k = 12)

# cut tree into 12 groups and rename column
hm.cutree <- as.data.frame(cutree(hm.dendrogram, k = 12))
colnames(hm.cutree) <- "cst"

# look at heatmap (the version produced the first time I explored this data, and
# NOT the one produced below) and dendrogram from rect.hclust(k = 12) output to
# identify marker samples for each group (to be used to cross-check cluster
# numbers)

#   1 = SRR6744107 = CST III (~90-100 % iners)
#   2 = SRR6744164 = CST IV (~70-100 % G. vag, plus other BVs or iners)
#   3 = SRR6744423 = CST IV (>60-75 % G. vag, plus other BVs or iners)
#   4 = SRR6744778 = CST III / V (relatively even split iners & jensenii)
#   5 = SRR6744485 = CST I (~40-70 % crisp, plus other lacto; crisp highest spp)
#   6 = SRR6744777 = CST I (~100% L. crispatus w/ small amounts of other lacto)
#   7 = SRR6743971 = CST IV (2x mostly Prevotella spp, 1x mostly A. vaginae)
#   8 = SRR6744809 = CST III (~75-90 % iners, plus mostly other lactobacilli)
#   9 = SRR6744443 = CST III (>55-80 % L. iners, plus G. vag & other lactos)
#  10 = SRR6744029 = CST III /IV (relatively even split iners & G.vag/other BV)
#  11 = SRR6744099 = CST V (~50-100 % jensenii, plus other lactos/bit of G. vag)
#  12 = SRR6744630 = CST II (~50-85 % gasseri, plus other lactos/bit of G. vag)

# make edits to cluster membership using the above information
hm.cutree$cst <- case_when(hm.cutree$cst == 1 ~ "III",
                           hm.cutree$cst == 2 ~ "IV",
                           hm.cutree$cst == 3 ~ "IV",
                           hm.cutree$cst == 4 ~ "III / V",
                           hm.cutree$cst == 5 ~ "I",
                           hm.cutree$cst == 6 ~ "I",
                           hm.cutree$cst == 7 ~ "IV",
                           hm.cutree$cst == 8 ~ "III",
                           hm.cutree$cst == 9 ~ "III",
                           hm.cutree$cst == 10 ~ "III / IV",
                           hm.cutree$cst == 11 ~ "V",
                           hm.cutree$cst == 12 ~ "II")

# change SRR6744877 to CST IV (very diverse; around 35% lacto, but most dominant
# taxon is E. faecalis; other CST IV spp. present, too)

hm.cutree["SRR6744877",1] <- "IV"

# summarise CST information
hm.cutree %>% 
  group_by(cst) %>% 
  count()

# rename cutree object and save as .Rda
# virginia.cst <- hm.cutree
# save(virginia.cst, file = paste(path.to.github,
#                                 "Rdata/virginia.cst.Rda", sep = ""))

#     cst             n
#     <chr>       <int>
# 1   I              66
# 2   II              4
# 3   III            88
# 4   III / IV       20
# 5   III / V        17
# 6   IV             94
# 7   V               8
# .....................
#     TOTAL         297
# .....................

# add CST info from cutree to heatmap metadata and move the column so it appears
# below 'Vaginal pH' on the heatmap
hm.meta$CST <- hm.cutree$cst
hm.meta <- hm.meta %>% 
  relocate(CST, .before = `Vaginal pH`)

# define colours for each CST, add these to the column colours list and reorder
# list elements so CST is closest to the heatmap
hm.column.cols[["CST"]] <- c(I = viridis(7)[1],
                             II = viridis(7)[2],
                             III = viridis(7)[3],
                             `III / IV` = viridis(7)[4],
                             `III / V` = viridis(7)[5],
                             IV = viridis(7)[6],
                             V = viridis(7)[7])

hm.column.cols <- hm.column.cols[c("CST", "Vaginal pH",
                                   "Discharge", "Odor", "Itching")]

# make heatmap for 297 virginia samples based on proportional abundances of 
# species within vaginal metatranscriptomes, with colour bars for vaginal pH,
# current vaginal symptoms and whether the patient had BV since their last visit
# png(filename = paste(path.to.github,
#                      "figs_for_paper/virginia_summarySpecies_heatmap.png", sep = ""),
#     units = "in", height = 8.75, width = 16, res = 400)

pheatmap(mat = virginia.filt.prop, color = hm.colours, show_colnames = FALSE,
         cluster_cols = hm.dendrogram, cluster_rows = FALSE,
         annotation_colors = hm.column.cols, annotation_col = hm.meta,
         legend_breaks = c(0,0.2,0.4,0.6,0.8,0.9939),
         legend_labels = c("0","0.2","0.4","0.6","0.8","1.0"),
         labels_row = as.expression(hm.row.italics),
         cellwidth = 3, cellheight = 20, border_color = rgb(0, 0, 0, 0.1))

# dev.off()

##################### molecular BV definition: dendrogram ######################

# load in the .Rda object containing a named vector of KO terms in the filtered
# virginia dataset if it exists in 'Rdata/', OR make it and save as .Rda
if(file.exists(paste(path.to.github, "Rdata/ko.virginia.filt.vec.Rda", sep = ""))){
  load(paste(path.to.github, "Rdata/ko.virginia.filt.vec.Rda", sep = ""))
} else{
  ko.virginia.filt.vec <- KO[rownames(virginia.filt),1] # subset vNum->KO lookup
  names(ko.virginia.filt.vec) <- rownames(virginia.filt) # assign vNums as names

  save(ko.virginia.filt.vec, file = paste(path.to.github,
                                     "Rdata/ko.virginia.filt.vec.Rda", sep = "")) 
}

# aggregate vNumber feature tables by KO term
ko.virginia.filt <- aggregate(virginia.filt, by=list(ko.virginia.filt.vec), FUN=sum)

# set row names to the KO number and remove the now-redundant column from the
# dataframe
rownames(ko.virginia.filt) <- ko.virginia.filt$Group.1
ko.virginia.filt$Group.1 <- NULL

# remove the following eukaryotic KOs (found these during curation of 'Unknown'
# pathways)
# K03260 (eukaryotic translation initiation factor 4G)
# K06237 (type IV collagen)
# K12373 (human lysosomal hexosaminidase)
# K00863 (human triose/dihydroxyacetone kinase)
# K13993 (human HSP20)
# K03648 (human uracil-DNA glycosylase)
# K01408 (eukaryotic insulin-degrading enzyme)
# K00599 (eukaryotic tRNA N(3)-methylcytidine methyltransferase)
ko.virginia.filt <- ko.virginia.filt[-which(grepl(paste("K03260","K06237","K12373","K00863",
                                                        "K13993","K03648","K01408","K00599",sep = "|"),
                                                  rownames(ko.virginia.filt))),]

# save filtered KO feature table for 297 virginia samples as a .Rda object
# save(ko.virginia.filt,
#      file = paste(path.to.github, "Rdata/ko.virginia.filt.Rda", sep = ""))

# do the same for the unfiltered dataset (this takes ~20 minutes): load in the
# .Rda objects if they exist in 'Rdata/', OR, make them from scratch and save
# as .Rda
if(file.exists(paste(path.to.github, "Rdata/ko.virginia.vec.Rda", sep = ""))){
  load(paste(path.to.github, "Rdata/ko.virginia.vec.Rda", sep = ""))
} else{
  ko.virginia.vec <- KO[rownames(virginia.data),1]
  names(ko.virginia.vec) <- rownames(virginia.data)
  
  save(ko.virginia.vec, 
       file = paste(path.to.github, "Rdata/ko.virginia.vec.Rda", sep = ""))
}

if(file.exists(paste(path.to.github, "Rdata/ko.virginia.Rda", sep = ""))){
  load(paste(path.to.github, "Rdata/ko.virginia.Rda", sep = ""))
} else{
  ko.virginia <- aggregate(virginia.data, by=list(ko.virginia.vec), FUN=sum)
  rownames(ko.virginia) <- ko.virginia$Group.1
  ko.virginia$Group.1 <- NULL
  ko.virginia <- ko.virginia[-which(grepl(paste("K03260","K06237","K12373","K00863",
                                                "K13993","K03648","K01408","K00599",sep = "|"),
                                          rownames(ko.virginia))),]
  save(ko.virginia,
       file = paste(path.to.github, "Rdata/ko.virginia.Rda", sep = ""))
}

# set seed for RNG
set.seed(2023)

# zero replacement with cMultRepl (transposed as samples expected to be in rows)
ko.virginia.filt.zr <- cmultRepl(t(ko.virginia.filt), method = "CZM", 
                                 label = "0", z.warning = 0.99)

# CLR transform virginia KO dataset (transposed to account for above transpose)
ko.virginia.filt.clr <- apply(t(ko.virginia.filt.zr), 2, function(x) log(x) - mean(log(x)))

# unsupervised clustering of CLR KO data using Euclidean distances & Ward's
# method
virginia.hclust <- dist(t(ko.virginia.filt.clr), method = "euclidean") %>% 
  hclust(method = "ward.D2")
  
# save hclust object for later clustering in 'virginia_validating.R'
# save(virginia.hclust,
#      file = paste(path.to.github,"Rdata/virginia.hclust.Rda",sep = ""))

# convert hclust output to dendrogram object for dendextend
virginia.hclust <- as.dendrogram(virginia.hclust)

# make vectors for additional, possibly confounding, metadata variables (to be
# added to dendrogram colour bars)
meta.abx <- case_when(is.na(virginia.meta$antibiotics_current) ~ "Unknown",
                       .default = as.character(virginia.meta$antibiotics_current))

meta.abx.since.last.visit <- case_when(is.na(virginia.meta$antibiotics_since_last_visit) ~ "Unknown",
                      .default = as.character(virginia.meta$antibiotics_since_last_visit))

meta.prior.rec.bv <- case_when(is.na(virginia.meta$prior_diagnosis_recurrent_bv) ~ "Unknown",
                      .default = as.character(virginia.meta$prior_diagnosis_recurrent_bv))

meta.bv.since.last.visit <- case_when(is.na(virginia.meta$bv_since_last_visit) ~ "Unknown",
                      .default = as.character(virginia.meta$bv_since_last_visit))

meta.prior.rec.yeast <- case_when(is.na(virginia.meta$prior_diagnosis_recurrent_yeast) ~ "Unknown",
                      .default = as.character(virginia.meta$prior_diagnosis_recurrent_yeast))

meta.yeast.since.last.visit <- case_when(is.na(virginia.meta$yeast_infection_since_last_visit) ~ "Unknown",
                      .default = as.character(virginia.meta$yeast_infection_since_last_visit))

meta.uti.since.last.visit <- case_when(is.na(virginia.meta$uti_since_last_visit) ~ "Unknown",
                      .default = as.character(virginia.meta$uti_since_last_visit))

meta.diabetes <- case_when(is.na(virginia.meta$diabetes) ~ "Unknown",
                      .default = as.character(virginia.meta$diabetes))

meta.trimester <- case_when(is.na(virginia.meta$visit_trimester) ~ "Unknown",
                            virginia.meta$visit_trimester == "1" ~ "First",
                            virginia.meta$visit_trimester == "2" ~ "Second",
                            virginia.meta$visit_trimester == "3" ~ "Third",)

meta.ethnicity <- case_when(is.na(virginia.meta$ethnicity) ~ "Unknown",
                      .default = as.character(virginia.meta$ethnicity))

levels(factor(meta.abx))
levels(factor(meta.abx.since.last.visit))
levels(factor(meta.prior.rec.bv))
levels(factor(meta.bv.since.last.visit))
levels(factor(meta.prior.rec.yeast))
levels(factor(meta.yeast.since.last.visit))
levels(factor(meta.uti.since.last.visit))
levels(factor(meta.diabetes))
levels(factor(meta.trimester))
levels(factor(meta.ethnicity))

# combine metadata vectors for virginia samples from previously generated 
# vectors (to be used for plotting coloured bars underneath dendrogram). Note 
# that dendextend plots the last column as the first bar, so order is reversed
# for cbind()
hm.meta <- as.data.frame(cbind(ethnicity = meta.ethnicity,
                               trimester = meta.trimester,
                               diabetes = meta.diabetes,
                               yeast.prior.rec = meta.prior.rec.yeast,
                               bv.prior.rec = meta.prior.rec.bv,
                               uti.last.vis = meta.uti.since.last.visit,
                               yeast.last.vis = meta.yeast.since.last.visit,
                               bv.last.vis = meta.bv.since.last.visit,
                               abx.last.vis = meta.abx.since.last.visit,
                               abx = meta.abx,
                               itching = meta.itching,
                               odor = meta.odor,
                               discharge = meta.discharge,
                               cst = hm.cutree$cst, 
                               ph = meta.ph))


# change factors to R colours (needed for dendextend's colored_bars() function)
hm.meta$ph <- case_when(hm.meta$ph == "4" ~ "red4",
                        hm.meta$ph == "5" ~ "red3",
                        hm.meta$ph == "5.5" ~ "indianred",
                        hm.meta$ph == "6" ~ "lightskyblue1",
                        hm.meta$ph == "6.5" ~ "steelblue3",
                        hm.meta$ph == "7" ~ "royalblue3",
                        hm.meta$ph == "7.5" ~ "blue4",
                        hm.meta$ph == "Unknown" ~ "grey")

hm.meta$cst <- case_when(hm.meta$cst == "I" ~ viridis(7)[1],
                         hm.meta$cst == "II" ~ viridis(7)[2],
                         hm.meta$cst == "III" ~ viridis(7)[3],
                         hm.meta$cst == "III / IV" ~ viridis(7)[4],
                         hm.meta$cst == "III / V" ~ viridis(7)[5],
                         hm.meta$cst == "IV" ~ viridis(7)[6],
                         hm.meta$cst == "V" ~ viridis(7)[7])

hm.meta$discharge <- case_when(hm.meta$discharge == "Yes" ~ "red",
                               hm.meta$discharge == "No" ~ "blue",
                               hm.meta$discharge == "Not_Sure" ~ "magenta4",
                               hm.meta$discharge == "Unknown" ~ "grey")

hm.meta$odor <- case_when(hm.meta$odor == "Yes" ~ "red",
                               hm.meta$odor == "No" ~ "blue",
                               hm.meta$odor == "Not_Sure" ~ "magenta4",
                               hm.meta$odor == "Unknown" ~ "grey")

hm.meta$itching <- case_when(hm.meta$itching == "Yes" ~ "red",
                               hm.meta$itching == "No" ~ "blue",
                               hm.meta$itching == "Not_Sure" ~ "magenta4",
                               hm.meta$itching == "Unknown" ~ "grey")

hm.meta$abx <- case_when(hm.meta$abx == "Yes" ~ "red",
                          hm.meta$abx == "No" ~ "blue",
                          hm.meta$abx == "Unknown" ~ "grey")

hm.meta$abx.last.vis <- case_when(hm.meta$abx.last.vis == "Yes" ~ "red",
                                  hm.meta$abx.last.vis == "No" ~ "blue",
                                  hm.meta$abx.last.vis == "Unknown" ~ "grey")

hm.meta$bv.prior.rec <- case_when(hm.meta$bv.prior.rec == "Yes" ~ "red",
                                  hm.meta$bv.prior.rec == "No" ~ "blue",
                                  hm.meta$bv.prior.rec == "Not_Sure" ~ "magenta4",
                                  hm.meta$bv.prior.rec == "Unknown" ~ "grey")

hm.meta$bv.last.vis <- case_when(hm.meta$bv.last.vis == "Yes" ~ "red",
                                 hm.meta$bv.last.vis == "No" ~ "blue",
                                 hm.meta$bv.last.vis == "Unknown" ~ "grey")

hm.meta$yeast.prior.rec <- case_when(hm.meta$yeast.prior.rec == "Yes" ~ "red",
                                     hm.meta$yeast.prior.rec == "No" ~ "blue",
                                     hm.meta$yeast.prior.rec == "Not_Sure" ~ "magenta4",
                                     hm.meta$yeast.prior.rec == "Unknown" ~ "grey")

hm.meta$yeast.last.vis <- case_when(hm.meta$yeast.last.vis == "Yes" ~ "red",
                                    hm.meta$yeast.last.vis == "No" ~ "blue",
                                    hm.meta$yeast.last.vis == "Unknown" ~ "grey")

hm.meta$uti.last.vis <- case_when(hm.meta$uti.last.vis == "Yes" ~ "red",
                                  hm.meta$uti.last.vis == "No" ~ "blue",
                                  hm.meta$uti.last.vis == "Unknown" ~ "grey")

hm.meta$diabetes <- case_when(hm.meta$diabetes == "Yes" ~ "red",
                              hm.meta$diabetes == "No" ~ "blue",
                              hm.meta$diabetes == "Not_Sure" ~ "magenta4",
                              hm.meta$diabetes == "Unknown" ~ "grey")

hm.meta$trimester <- case_when(hm.meta$trimester == "First" ~ "snow2",
                               hm.meta$trimester == "Second" ~ "grey40",
                               hm.meta$trimester == "Third" ~ "black")

hm.meta$ethnicity <- case_when(hm.meta$ethnicity == "African-American" ~ "yellow2",
                               hm.meta$ethnicity == "African-American / Caucasian" ~ "yellow3",
                               hm.meta$ethnicity == "African-American / Hispanic or Latino" ~ "yellow4",
                               hm.meta$ethnicity == "American Indian or Alaska Native" ~ "cyan2",
                               hm.meta$ethnicity == "American Indian or Alaska Native / Caucasian" ~ "cyan4",
                               hm.meta$ethnicity == "Asian" ~ "darkseagreen",
                               hm.meta$ethnicity == "Asian / Caucasian" ~ "darkseagreen1",
                               hm.meta$ethnicity == "Caucasian" ~ "deepskyblue",
                               hm.meta$ethnicity == "Caucasian / Hispanic or Latino" ~ "deepskyblue3",
                               hm.meta$ethnicity == "Hispanic or Latino" ~ "firebrick",
                               hm.meta$ethnicity == "Unknown" ~ "grey")

# set row & col names, then re-order so pH is first colour bar & itching is last
# (last column is displayed as the first colour bar with the dendextend package)
rownames(hm.meta) <- colnames(virginia.filt)
colnames(hm.meta) <- c("Ethnicity", "Trimester sampled", "Diabetes",
                       "Prior recurrent yeast infection", "Prior recurrent BV", 
                       "UTI since last visit",
                       "Yeast infection since last visit",
                       "BV since last visit",
                       "Antibiotics since last visit", "Antibiotics currently",
                       "Itching", "Odor", "Discharge", "CST", "Vaginal pH")

# plot dendrogram & colour bars corresponding to vaginal pH, CST and vaginal
# symptoms, hiding sample names
# png(filename = paste(path.to.github,
#                      "figs_for_paper/virginia_summaryK0_dendrogram.png", sep = ""),
#     units = "in", width = 10, height = 5, res = 400)

virginia.hclust %>%
  set("labels_colors", "white") %>% 
  plot(cex.axis = 0.6)

colored_bars(colors = hm.meta, dend = virginia.hclust, 
             y_shift = -20, y_scale = 515, cex.rowLabels = 0.5,
             rowLabels = c("Ethnicity", "Trimester sampled", "Diabetes",
                           "Prior recurrent yeast infection", "Prior recurrent BV", 
                           "UTI since last visit",
                           "Yeast infection since last visit",
                           "BV since last visit",
                           "Antibiotics since last visit", "Antibiotics currently",
                           "Itching", "Odor", "Discharge", "CST", "Vaginal pH"))

legend("topright", legend = c("4.0", "5.0", "5.5",
                              "6.0","6.5","7.0","7.5"),
       fill = c("red4", "red3", "indianred", "lightskyblue1",
                "steelblue3", "royalblue3", "blue4"),
       border = rgb(0,0,0,0), box.lwd = 0.5, cex = 0.475, title = "Vaginal pH", inset = c(0.005, 0.035))

legend("topright", legend = c("Yes", "No", "Not sure", "Unknown"),
       fill = c("red", "blue", "magenta4", "grey"),
       border = rgb(0,0,0,0), box.lwd = 0.5, cex = 0.475, title = "Clinical factors", inset = c(0.125, 0.085))

legend("topright", legend = c("I", "II", "III", "III / IV", "III / V", "IV", "V"),
       fill = viridis(7),
       border = rgb(0,0,0,0), box.lwd = 0.5, cex = 0.475, title = "CST", inset = c(0.06375, 0.035))

legend("topright", legend = c("First", "Second", "Third"),
       fill = c("snow2", "grey40", "black"),
       border = rgb(0,0,0,0), box.lwd = 0.5, cex = 0.475, title = "Trimester", inset = c(0.198, 0.1))

legend("topright", legend = c("African-American", "African-American / Caucasian",
                              "African-American / Hispanic or Latino",
                              "American Indian or Alaska Native", 
                              "American Indian or Alaska Native / Caucasian",
                              "Asian", "Asian / Caucasian", "Caucasian",
                              "Caucasian / Hispanic or Latino",
                              "Hispanic or Latino", "Unknown"),
       fill = c("yellow2", "yellow3", "yellow4", "cyan2", "cyan4", 
                "darkseagreen", "darkseagreen1", "deepskyblue", "deepskyblue3",
                "firebrick", "grey"),
       border = rgb(0,0,0,0), box.lwd = 0.5, cex = 0.45, title = "Ethnicity", inset = c(0.005, 0.34), ncol = 1)

# dev.off()

# quite a clear split into two main branches- great. Main points:
#   - CST IV associated with higher pH and CST I and V associated with low pH.
#   - Within 'healthy', crispatus is very easy to identify and, as expected,
#     appears to be associated with Caucasian ethnicity
#   - Only three CST IV samples on the 'healthy' side: SRR6744826 (60% G. vag,
#     ~20% crispatus, ~15% iners), SRR6744176 (~50% G. vag, ~20% iners, ~25% 
#     jensenii), SRR6744338 (~75% G. vag, 15 % iners and ~5% jensenii).
#   - CST III is split relatively evenly between 'healthy' and BV sides, as
#     expected (slightly more on low pH side?)
#   - Samples do NOT cluster by vaginal symptoms, other clinical factors or
#     trimester


################## molecular BV definition: editing pathways ###################

# load in .Rda object containing KO -> pathway table for unfiltered virginia
# dataset if it exists in 'Rdata/', OR, make it from scratch by: subsetting
# pathway and KO tables for unfiltered virginia vNumbers, removing vNumbers
# without an assigned KO term, merging the two dataframes, keeping KOs with
# no information (replace with 'Unknown'), retaining only 'ko' / 'pathway'
# columns, removing duplicate KOs and saving as a .Rda object

# NOTE 1: aggregate() removes rows which have no corresponding entry in the 'by='
# input, so any vNumbers which have not been assigned a KO number are removed

# NOTE 2: if making from scratch, the first two commands take a total of ~45
#         minutes!

if(file.exists(paste(path.to.github, "Rdata/virginia.ko.path.Rda", sep = ""))){
  load(paste(path.to.github, "Rdata/virginia.ko.path.Rda", sep = ""))
} else{
  virginia.vnum.ko <- KO[rownames(virginia.data),] # vnum -> KO
  virginia.vnum.path <-pathway.table[rownames(virginia.data),] # vnum -> pathway
  length(which(is.na(virginia.vnum.ko$V2))) # 256,581 without a KO number
  
  virginia.vnum.ko <-na.omit(virginia.vnum.ko) # remove vNumbers without a KO
  virginia.ko.path <- merge(virginia.vnum.ko, virginia.vnum.path,
                            by = 0, all.x = TRUE) # merge lookup tables
  
  colnames(virginia.ko.path)[2] <- "ko"
  rownames(virginia.ko.path) <- virginia.ko.path$Row.names
  virginia.ko.path$pathway <- case_when(is.na(virginia.ko.path$pathway) ~ "Unknown",
                                        .default = virginia.ko.path$pathway)
  
  virginia.ko.path <- virginia.ko.path %>% # retain only 'ko' & 'pathway'
    dplyr::select(ko, pathway)
  
  # de-duplicate KOs
  virginia.ko.path <- virginia.ko.path[which(!duplicated(virginia.ko.path$ko)),]
  rownames(virginia.ko.path) <- virginia.ko.path$ko 
  virginia.ko.path <- virginia.ko.path %>% 
    arrange(ko) %>% 
    dplyr::select(-ko)
  
  # remove eukaryotic KOs
  virginia.ko.path <- virginia.ko.path[rownames(ko.virginia), , drop = FALSE]
  
  save(virginia.ko.path,
       file = paste(path.to.github, "Rdata/virginia.ko.path.Rda", sep = "")) 
}

# load in .Rda object containing KO -> pathway table for filtered virginia
# dataset if it exists in 'Rdata/', OR, make it from scratch by: subsetting the
# unfiltered K0 -> pathway table to retain only the K0s found in the filtered
# virginia dataset, manually editing many K0s and saving as a .Rda object
if(file.exists(paste(path.to.github, "Rdata/virginia.ko.path.filt.Rda", sep = ""))){
  load(paste(path.to.github, "Rdata/virginia.ko.path.filt.Rda", sep = ""))
} else{
  virginia.ko.path.filt <- virginia.ko.path[rownames(ko.virginia.filt), ,drop = FALSE]
  
# manually edit pathway names (many are misassigned or too long)
  virginia.ko.path.filt$pathway <- case_when(virginia.ko.path.filt$pathway == "Adipocytokine signaling pathway" ~ "Fatty acid biosynthesis",
                                             virginia.ko.path.filt$pathway == "African trypanosomiasis" ~ "Prolyl oligopeptidase", 
                                             virginia.ko.path.filt$pathway == "Alanine aspartate and glutamate metabolism" ~ "Ala-Asp-Glu metabolism",
                                             virginia.ko.path.filt$pathway == "Amino sugar and nucleotide sugar metabolism" ~ "Amino sugar metabolism",
                                             virginia.ko.path.filt$pathway == "Arginine and proline metabolism" ~ "Arg-Pro metabolism",
                                             virginia.ko.path.filt$pathway == "Bacterial invasion of epithelial cells" ~ "Intracellular invasion",
                                             virginia.ko.path.filt$pathway == "Biosynthesis of vancomycin group antibiotics" ~ "Amino sugar metabolism",
                                             virginia.ko.path.filt$pathway == "Biosynthesis of unsaturated fatty acids" ~ "Fatty acid biosynthesis",
                                             virginia.ko.path.filt$pathway == "Cationic antimicrobial peptide (CAMP) resistance" ~ "CAMP resistance",
                                             virginia.ko.path.filt$pathway == "Cysteine and methionine metabolism" ~ "Cys-Met metabolism",
                                             virginia.ko.path.filt$pathway == "Drug metabolism - cytochrome P450" ~ "Glycolysis / Gluconeogenesis",
                                             virginia.ko.path.filt$pathway == "Drug metabolism - other enzymes" ~ "Drug metabolism",
                                             virginia.ko.path.filt$pathway == "Fructose and mannose metabolism" ~ "Glycolysis / Gluconeogenesis",
                                             virginia.ko.path.filt$pathway == "Glycerolipid metabolism" ~ "Glycerophospholipid metabolism",
                                             virginia.ko.path.filt$pathway == "Glycine serine and threonine metabolism" ~ "Gly-Ser-Thr metabolism",
                                             virginia.ko.path.filt$pathway == "Glycosphingolipid biosynthesis - globo series" ~ "Glycosphingolipid biosynthesis",
                                             virginia.ko.path.filt$pathway == "Histidine metabolism" ~ "His metabolism",
                                             virginia.ko.path.filt$pathway == "Lysine biosynthesis" ~ "Lys metabolism",
                                             virginia.ko.path.filt$pathway == "Lysine degradation" ~ "Lys metabolism",
                                             virginia.ko.path.filt$pathway == "MAPK signaling pathway - yeast" ~ "Pyruvate metabolism",
                                             virginia.ko.path.filt$pathway == "Nicotinate and nicotinamide metabolism" ~ "Vitamin B3 metabolism",
                                             virginia.ko.path.filt$pathway == "Phenylalanine metabolism" ~ "Phe-Tyr-Trp metabolism",
                                             virginia.ko.path.filt$pathway == "Phenylalanine tyrosine and tryptophan biosynthesis" ~ "Phe-Tyr-Trp metabolism",
                                             virginia.ko.path.filt$pathway == "Phosphotransferase system (PTS)" ~ "Phosphotransferase system",
                                             virginia.ko.path.filt$pathway == "Photosynthesis" ~ "Oxidative phosphorylation",
                                             virginia.ko.path.filt$pathway == "Porphyrin and chlorophyll metabolism" ~ "Porphyrin metabolism",
                                             virginia.ko.path.filt$pathway == "Ribosome biogenesis in eukaryotes" ~ "Oligoribonuclease activity",
                                             virginia.ko.path.filt$pathway == "Terpenoid backbone biosynthesis" ~ "Terpenoid biosynthesis",
                                             virginia.ko.path.filt$pathway == "Valine leucine and isoleucine biosynthesis" ~ "Val-Leu-Ile metabolism",
                                             virginia.ko.path.filt$pathway == "Valine leucine and isoleucine degradation" ~ "Val-Leu-Ile metabolism",
                                             .default = virginia.ko.path.filt$pathway)

# further editing of specific KO terms to correct pathway (searched on KEGG)
  virginia.ko.path.filt$pathway <-case_when(rownames(virginia.ko.path.filt) == "K01580" ~ "ABC transporters",
                                            rownames(virginia.ko.path.filt) == "K01580" ~ "Acid survival",
                                            rownames(virginia.ko.path.filt) %in% c("K00812","K11358","K01915",
                                                                                   "K01775","K01776","K00811") ~ "Ala-Asp-Glu metabolism",
                                            rownames(virginia.ko.path.filt) == "K14195" ~ "Bacterial adherence",
                                            rownames(virginia.ko.path.filt) == "K06595" ~ "Bacterial chemotaxis",
                                            rownames(virginia.ko.path.filt) == "K11936" ~ "Biofilm formation",
                                            rownames(virginia.ko.path.filt) == "K14534" ~ "Butanoate metabolism",
                                            rownames(virginia.ko.path.filt) %in% c("K03367","K03739","K03740","K14188","K14205") ~ "CAMP resistance",
                                            rownames(virginia.ko.path.filt) %in% c("K03531","K03588","K03589","K03590","K13581") ~ "Cell division",
                                            rownames(virginia.ko.path.filt) %in% c("K00024","K00174","K00175","K00239","K00240",
                                                                                   "K00241","K01676","K01681","K01847","K01848",
                                                                                   "K01849","K01902","K01903","K01966","K03737",
                                                                                   "K05606","K13788","K00031","K00177","K00625",
                                                                                   "K00925","K01595","K01895","K05884","K01679") ~ "Citrate cycle (TCA cycle)",
                                            rownames(virginia.ko.path.filt) %in% c("K01358","K03544") ~ "Clp-dependent proteolysis",
                                            rownames(virginia.ko.path.filt) == "K11031" ~ "Cytolysis",
                                            rownames(virginia.ko.path.filt) %in% c("K02313","K02314") ~ "DNA replication",
                                            rownames(virginia.ko.path.filt) == "K00864" ~ "Glycerolipid metabolism",
                                            rownames(virginia.ko.path.filt) == "K00981" ~ "Glycerophospholipid metabolism",
                                            rownames(virginia.ko.path.filt) %in% c("K02446","K04041","K00134","K01803","K00873",
                                                                                   "K00927","K01610","K01624","K00850","K01834","") ~ "Glycolysis / Gluconeogenesis",
                                            rownames(virginia.ko.path.filt) == "K00018" ~ "Glyoxylate and dicarboxylate metabolism",
                                            rownames(virginia.ko.path.filt) %in% c("K00058","K00600","K01079") ~ "Gly-Ser-Thr metabolism",
                                            rownames(virginia.ko.path.filt) == "K01715" ~ "Fatty acid biosynthesis",
                                            rownames(virginia.ko.path.filt) %in% c("K02405", "K02406", "K03092","K06603","K13626") ~ "Flagellar assembly",
                                            rownames(virginia.ko.path.filt) == "K00817" ~ "His metabolism",
                                            rownames(virginia.ko.path.filt) %in% c("K11749","K02012","K02217","K07243") ~ "Iron acquisition",
                                            rownames(virginia.ko.path.filt) == "K01338" ~ "Lon protease",
                                            rownames(virginia.ko.path.filt) == "K01582" ~ "Lys metabolism",
                                            rownames(virginia.ko.path.filt) %in% c("K05299","K00123","K00124","K00127") ~ "Methanoate metabolism",
                                            rownames(virginia.ko.path.filt) %in% c("K04079","K04043","K04077") ~ "Molecular chaperones",
                                            rownames(virginia.ko.path.filt) %in% c("K00297","K01491","K01938") ~ "One carbon pool by folate",
                                            rownames(virginia.ko.path.filt) %in% c("K02117","K02118","K02119","K02120","K02121",
                                                                                   "K02122","K02123","K02124") ~ "Oxidative phosphorylation",
                                            rownames(virginia.ko.path.filt) == "K00832" ~ "Phe-Tyr-Trp metabolism",
                                            rownames(virginia.ko.path.filt) == "K01195" ~ "Pentose and glucuronate interconversions",
                                            rownames(virginia.ko.path.filt) %in% c("K01621","K01632","K01783","K01807",
                                                                                   "K01808","K00615") ~ "Pentose phosphate pathway",
                                            rownames(virginia.ko.path.filt) == "K02563" ~ "Peptidoglycan biosynthesis",
                                            rownames(virginia.ko.path.filt) == "K01488" ~ "Purine metabolism",
                                            rownames(virginia.ko.path.filt) == "K00757" ~ "Pyrimidine metabolism",
                                            rownames(virginia.ko.path.filt) %in% c("K00029","K01006","K01007") ~ "Pyruvate metabolism",
                                            rownames(virginia.ko.path.filt) == "K01661" ~ "Quinone biosynthesis",
                                            rownames(virginia.ko.path.filt) %in% c("K14051","K07173") ~ "Quorum sensing",
                                            rownames(virginia.ko.path.filt) == "K01186" ~ "Sialidase activity",
                                            rownames(virginia.ko.path.filt) %in% c("K03313","K03315") ~ "Sodium-proton antiport",
                                            rownames(virginia.ko.path.filt) == "K00688" ~ "Starch and sucrose metabolism",
                                            rownames(virginia.ko.path.filt) == "K04564" ~ "ROS degradation",
                                            rownames(virginia.ko.path.filt) == "K00869" ~ "Terpenoid biosynthesis",
                                            rownames(virginia.ko.path.filt) == "K02358" ~ "Translation",
                                            rownames(virginia.ko.path.filt) == "K02040" ~ "Two-component system",
                                            rownames(virginia.ko.path.filt) == "K08303" ~ "Type I collagen degradation",
                                            rownames(virginia.ko.path.filt) == "K00532" ~ "Unknown",
                                            rownames(virginia.ko.path.filt) == "K01428" ~ "Urease activity",
                                            rownames(virginia.ko.path.filt) %in% c("K03183","K03688") ~ "Ubiquinone biosynthesis",
                                            rownames(virginia.ko.path.filt) %in% c("K02548","K02548") ~ "Vitamin K2 biosynthesis",
                                            .default = virginia.ko.path.filt$pathway)

  # even further editing of specific KO terms to correct pathway (discovered
  # after analysing data downstream- these were all 'Unknown'; see end of file)
  virginia.ko.path.filt$pathway <- case_when(rownames(virginia.ko.path.filt) %in% c("K00435","K01772") ~ "Porphyrin metabolism",
                                          rownames(virginia.ko.path.filt) == "K07777" ~ "Two-component system",
                                          rownames(virginia.ko.path.filt) == "K01005" ~ "Teichoic acid biosynthesis",
                                          rownames(virginia.ko.path.filt) == "K00903" ~ "Exopolysaccharide biosynthesis",
                                          rownames(virginia.ko.path.filt) %in% c("K01356","K03631","K03700") ~ "DNA damage repair",
                                          rownames(virginia.ko.path.filt) == "K02055" ~ "Putrescine/Spermidine transport",
                                          rownames(virginia.ko.path.filt) %in% c("K02068","K02069") ~ "Iron homeostasis",
                                          rownames(virginia.ko.path.filt) %in% c("K02237","K02244","K02248") ~ "Extracellular DNA uptake",
                                          rownames(virginia.ko.path.filt) %in% c("K02444","K03436") ~ "Fructose metabolism",
                                          rownames(virginia.ko.path.filt) %in% c("K03282","K03442") ~ "Osmotic stress adaptation",
                                          rownames(virginia.ko.path.filt) == "K02530" ~ "Lactose metabolism",
                                          rownames(virginia.ko.path.filt) == "K02859" ~ "Riboflavin metabolism",
                                          rownames(virginia.ko.path.filt) == "K03284" ~ "Magnesium transport",
                                          rownames(virginia.ko.path.filt) %in% c("K03297","K08161") ~ "Multidrug resistance - efflux pump",
                                          rownames(virginia.ko.path.filt) == "K03303" ~ "Lactate metabolism",
                                          rownames(virginia.ko.path.filt) == "K03322" ~ "Manganese acquisition",
                                          rownames(virginia.ko.path.filt) %in% c("K03484", "K05992") ~ "Starch and sucrose metabolism",
                                          rownames(virginia.ko.path.filt) == "K04075" ~ "Aminoacyl-tRNA biosynthesis",
                                          rownames(virginia.ko.path.filt) == "K07043" ~ "Pyrimidine metabolism",
                                          rownames(virginia.ko.path.filt) == "K01595" ~ "Pyruvate metabolism",
                                          rownames(virginia.ko.path.filt) == "K01066" ~ "Unknown",
                                          .default = virginia.ko.path.filt$pathway)
  
  virginia.ko.path.filt$pathway <- gsub(pattern = "and", replacement = "&",
                                      virginia.ko.path.filt$pathway)

  save(virginia.ko.path.filt, file = paste(path.to.github,
                                    "Rdata/virginia.ko.path.filt.Rda", sep = ""))
}

# group by pathway and count number of KO numbers per pathway: 126 x pathways
pathways.ko <- virginia.ko.path.filt %>% 
  group_by(pathway) %>% 
  count() %>% 
  arrange(desc(n))

####################### molecular BV definition: biplot #######################

# perform PCA on zero-replaced, CLR-transformed KO feature table
pca.virginia.ko <- prcomp(t(ko.virginia.filt.clr))

# make vector of pathways to highlight on the biplot
highlight.path <- c("Ribosome", "Aminoacyl-tRNA biosynthesis",
                    "CAMP resistance","Flagellar assembly",
                    "Bacterial chemotaxis", "Porphyrin metabolism",
                    "Butanoate metabolism", "Two-component system",
                    "Starch & sucrose metabolism", "ABC transporters")

# get list of indices for each KO pathway
indices.pathway <- list()
for(i in levels(factor(highlight.path))){
  indices.pathway[[i]] <- which(virginia.ko.path.filt$pathway == i)
}

# add all other pathways to list under "Other"
indices.pathway[["Other"]] <- setdiff(1:nrow(virginia.ko.path.filt),
                                            unlist(indices.pathway))

# make vector of pathway colours ("Other" as transparent black)
highlight.path.col <- c("olivedrab2", "lightskyblue1", "tomato3", "gold",
                        "chocolate4", "purple3", "bisque2", "black",
                        "deeppink3", "cyan3", rgb(0,0,0,0.05))

# split samples into two groups based on KO dendrogram and edit names (1 is 
# healthy, 2 is molecular BV)
virginia.groups <- as.data.frame(cutree(virginia.hclust, k = 2),nm = "group")
virginia.groups$group <- case_when(virginia.groups$group == 1 ~ "Healthy",
                                   virginia.groups$group == 2 ~ "BV")

# make list of group indices
indices.group <- list()
for(i in levels(factor(virginia.groups$group))){
  indices.group[[i]] <- which(virginia.groups$group == i)
}

# make biplot of 297 virginia samples from PCA of KO terms, highlighting 
# specific pathways
# png(filename = paste(path.to.github,
#                      "figs_for_paper/virginia_healthBV2_biplot_K0.png", sep = ""),
#     units = "in", height = 6, width = 12, res = 400)

codaSeq.PCAplot(pca.virginia.ko, plot.groups = TRUE, plot.loadings = TRUE,
                plot.ellipses = NULL, plot.density = NULL, PC = c(1,2),
                grp = indices.group, grp.col = c("goldenrod2", "steelblue1"),
                grp.sym = "text", grp.cex = 0.3, load.grp = indices.pathway,
                load.col = highlight.path.col, load.sym = rep(19,11), load.cex = 0.5,
                plot.legend = "loadings", leg.position = "top", leg.cex = 0.6,
                leg.columns = 4, title = "Virginia - K0 terms: Health vs. Molecular BV")

legend("bottomleft", legend = c("Healthy", "Molecular BV"),
       col = c("steelblue1","goldenrod2"), pch = "SS", cex = 0.6)

# dev.off()

# biplot shows two clearly distinct 'healthy' subgroups, based on functional
# clustering. These correspond to the three clusters identified by cutree. Split
# the KO dendrogram into 3, rename groups, remake index list and re-plot
virginia.groups <- as.data.frame(cbind(groups.2= virginia.groups$group,
                                       groups.3= cutree(virginia.hclust, k = 3)))

virginia.groups$groups.3 <- case_when(virginia.groups$groups.3 == 1 ~ "Healthy-G1",
                                      virginia.groups$groups.3 == 2 ~ "BV",
                                      virginia.groups$groups.3 == 3 ~ "Healthy-G2")

# save healthy vs. molecular BV designations
# save(virginia.groups, file = paste(path.to.github,
#                                    "Rdata/virginia.groups.Rda", sep = ""))

indices.group.3 <- list()
for(i in levels(factor(virginia.groups$groups.3))){
  indices.group.3[[i]] <- which(virginia.groups$groups.3 == i)
}

# png(filename = paste(path.to.github,
#                      "figs_for_paper/virginia_healthBV3_biplot_K0.png", sep = ""),
#     units = "in", height = 6, width = 12, res = 400)

codaSeq.PCAplot(pca.virginia.ko, plot.groups = TRUE, plot.loadings = TRUE,
                plot.ellipses = NULL, plot.density = NULL, PC = c(1,2),
                grp = indices.group.3, grp.col = c("goldenrod2", "steelblue1", "royalblue3"),
                grp.sym = "text", grp.cex = 0.3, load.grp = indices.pathway,
                load.col = highlight.path.col, load.sym = rep(19,11), load.cex = 0.5,
                plot.legend = "loadings", leg.position = "top", leg.cex = 0.6,
                leg.columns = 4, title = "Virginia - K0 terms: Health vs. Molecular BV")

legend("bottomleft", legend = c("Healthy (Subgroup 1)", "Healthy (Subgroup 2)", "Molecular BV"),
       col = c("steelblue1", "royalblue3", "goldenrod2"), pch = "SS", cex = 0.6)

# dev.off()

################## molecular BV definition: biplot (centred) ###################

# load in .Rda object containing combined ALDEx2 output for filtered virginia 
# K0s if it exists in 'Rdata/', OR run aldex2 with scale simulation and save 
# the combined output as a. Rda object
if(file.exists(paste(path.to.github, "Rdata/virginia.filt.ko.clr.c.Rda", sep = ""))){
  load(paste(path.to.github, "Rdata/virginia.filt.ko.clr.c.Rda", sep = ""))
} else{
  
  # making this from scratch takes a minute or so
  # set seed for RNG
  set.seed(2023)
  
  # make a scale matrix for running aldex w/ scale simulation: 15 % difference
  # sd for distribution of 0.5
  mu.matrix <- aldex.makeScaleMatrix(gamma=0.5, mu=c(1.15, 1),
                                     virginia.groups$groups.2)
  
  # CLR transformation with scale-simulation
  virginia.filt.ko.clr <- aldex.clr(reads = ko.virginia.filt,
                                   conds = virginia.groups$groups.2,
                                   mc.samples = 128, denom = "all",
                                   gamma = mu.matrix, verbose = TRUE)

  # calculate effect sizes and t-test statistics, then combine data frames
  virginia.filt.ko.clr.e <- aldex.effect(virginia.filt.ko.clr,
                                         include.sample.summary = TRUE)
  virginia.filt.ko.clr.t <- aldex.ttest(virginia.filt.ko.clr,
                                        verbose = TRUE)
  virginia.filt.ko.clr.c <-cbind(virginia.filt.ko.clr.e,
                               virginia.filt.ko.clr.t)

  # save the combined CLR output & summary as a .Rda object
  save(virginia.filt.ko.clr.c,
       file = paste(path.to.github, "Rdata/virginia.filt.ko.clr.c.Rda", sep = ""))
}

# pull out columns corresponding to median CLR values (all samples)
virginia.filt.ko.clr.df <- virginia.filt.ko.clr.c[,grep("rab.sample.", colnames(virginia.filt.ko.clr.c))]

# sanity check that K0 order in clr object is equal to virginia ko -> pathway
# table (all 1690 are TRUE)
length(which(rownames(virginia.ko.path.filt) == rownames(virginia.filt.ko.clr.df)))

# make vectors of K0 numbers which correspond to "Ribosome", "Aminoacyl-tRNA
# biosynthesis", "Glycolysis / Gluconeogenesis"
path.ribo <- rownames(virginia.ko.path.filt)[which(virginia.ko.path.filt$pathway == "Ribosome")]
path.trna <- rownames(virginia.ko.path.filt)[which(virginia.ko.path.filt$pathway == "Aminoacyl-tRNA biosynthesis")]
path.glyco <- rownames(virginia.ko.path.filt)[which(virginia.ko.path.filt$pathway == "Glycolysis / Gluconeogenesis")]

# plot difference vs. dispersion, add title and plot above pathways in blue (all
# are nice and centred after scale-sim)
# png(filename = paste(path.to.github,
#                      "figs_for_paper/virginia_healthBV_MW.png", sep = ""),
#     units = "in", height = 7, width = 14, res = 400)

aldex.plot(virginia.filt.ko.clr.c, xlim=c(0.346,9), cutoff.pval=0.01)
title('ALDEx2 w/ ScaleSim (Virginia filtered K0): mu = 1 : 1.15, gamma = 0.5, p <0.01',
      adj=0, line= 0.8)

points(virginia.filt.ko.clr.c[path.ribo ,"diff.win"],
       virginia.filt.ko.clr.c[path.ribo ,"diff.btw"],
       pch= 19, col= "blue3", cex= 0.5)

points(virginia.filt.ko.clr.c[path.trna,"diff.win"],
       virginia.filt.ko.clr.c[path.trna,"diff.btw"],
       pch= 19, col= "royalblue1", cex= 0.5)

points(virginia.filt.ko.clr.c[path.glyco,"diff.win"],
       virginia.filt.ko.clr.c[path.glyco,"diff.btw"],
       pch= 19, col= "lightskyblue1", cex= 0.5)

points(virginia.filt.ko.clr.c[which(abs(virginia.filt.ko.clr.c$effect) >1),"diff.win"],
       virginia.filt.ko.clr.c[which(abs(virginia.filt.ko.clr.c$effect) >1),"diff.btw"],
       pch= 19, col= "purple3", cex= 0.5)

legend("bottomleft", pch = 19, cex = 0.75,
       col = c( "black", "red", "purple3", "blue3", "royalblue1", "lightskyblue1"),
       legend = c("Non-significant", "P <0.01", "Absolute effect size >1",
                  "Ribosome", "Aminoacyl-tRNA biosynthesis", "Glycolysis / Gluconeogenesis"))

# dev.off()

# looks like almost all K0 terms are 'significantly' different between health and 
# BV... tried multiple mu & gamma value combinations, but picture is the same

# check how many K0 numbers have an absolute effect size >1 (purple points on
# MW plot)
length(which(abs(virginia.filt.ko.clr.c$effect) >1 )) # 567 K0 numbers

# remove 'rab.sample.' from column names
colnames(virginia.filt.ko.clr.df) <- gsub("rab.sample.", "", 
                                          colnames(virginia.filt.ko.clr.df))

# generate new PCA plot using 'centred' CLR values
pca.clr.virginia.ko <- prcomp(t(virginia.filt.ko.clr.df))

# plot new, centred PCA, splitting healthy into two subgroups
# png(filename = paste(path.to.github,
#                      "figs_for_paper/virginia_healthBV3_biplot_K0_centred.png", sep = ""),
#     units = "in", height = 6, width = 12, res = 400)

codaSeq.PCAplot(pca.clr.virginia.ko, plot.groups = TRUE, plot.loadings = TRUE,
                plot.ellipses = NULL, plot.density = NULL, PC = c(1,2),
                grp = indices.group.3, grp.col = c("goldenrod2", "steelblue1", "royalblue3"),
                grp.sym = "text", grp.cex = 0.3, load.grp = indices.pathway,
                load.col = highlight.path.col, load.sym = rep(19,11), load.cex = 0.5,
                plot.legend = "loadings", leg.position = "top", leg.cex = 0.6,
                leg.columns = 4, title = "Virginia - Centred K0 terms: Health vs. Molecular BV")

legend("bottomleft", legend = c("Healthy (Subgroup 1)", "Healthy (Subgroup 2)", "Molecular BV"),
       col = c("steelblue1", "royalblue3", "goldenrod2"), pch = "SSS", cex = 0.6)

# dev.off()

############### assigning Unknown K0 terms: literature evidence ################

# K00435 = hydrogen peroxide-dependent heme synthase: Porphyrin metabolism
# K07777 = two-component system, NarL family, sensor histidine kinase DegS: Two-component system
# K01005 = polyisoprenyl-teichoic acid--peptidoglycan teichoic acid transferase: Teichoic acid biosynthesis
# K00903 = protein-tyrosine kinase: Exopolysaccharide biosynthesis (https://doi.org/10.1101/gad.246397.114.)
# K01356 = repressor LexA: DNA damage repair (https://doi.org/10.1073/pnas.78.7.4204)
# K02055 = putative spermidine/putrescine transport system substrate-binding protein: Putrescine/Spermidine transport
# K02068 = UDP-glucose/iron transport system ATP-binding protein: Iron homeostasis (https://doi.org/10.1128/AEM.02322-13)
# K02069 = UDP-glucose/iron transport system permease proteinL Iron homeostasis (https://doi.org/10.1128/AEM.02322-13)
# K02237 = competence protein ComEA: Extracellular DNA uptake (https://doi.org/10.1128/jb.180.1.41-45.1998)
# K02244 = competence protein ComGB: Extracellular DNA uptake (https://doi.org/10.1128/jb.180.1.41-45.1998)
# K02248 = competence protein ComGF: Extracellular DNA uptake (https://doi.org/10.1128/jb.180.1.41-45.1998)
# K02444 = glycerol-3-phosphate regulon repressor: Fructose metabolism (https://doi.org/10.1128/jb.00244-18)
# K02530 = lactose phosphotransferase system repressor: Lactose metabolism (https://doi.org/10.1128/AEM.01370-14)
# K02859 = riboflavin biosynthesis RibT protein: Riboflavin metabolism (https://doi.org/10.3389/fmicb.2022.856820)
# K03282 = large conductance mechanosensitive channel: Osmotic stress adaptation (https://doi.org/10.1093/emboj/18.7.1730)
# K03284 = Magnesium transporter: Magnesium transport (https://doi.org/10.1080/09687680701441883)
# K03297 = small multidrug resistance pump: Multidrug resistance - efflux pump
# K03303 = lactate permease: Lactate metabolism (https://doi.org/10.1128/AEM.00672-19)
# K03322 = manganese transport protein: Manganese acquisition (https://doi.org/10.1128/AEM.65.11.4746-4752.1999)
# K03436 = fructose operon transcriptional repressor: Fructose metabolism (https://doi.org/10.1111/j.1365-2958.1995.tb02339.x)
# K03442 = small conductance mechanosensitive channel: Osmotic stress adaptation (https://doi.org/10.1093/emboj/18.7.1730)
# K03484 = sucrose operon repressor: Starch & sucrose metabolism (https://doi.org/10.3389/fmicb.2021.754464)
# K03631 = DNA repair protein RecN: DNA damage repair (https://doi.org/10.1038/ncomms15282)
# K03700 = recombination protein U: DNA damage repair (https://doi.org/10.1016/j.str.2005.05.011)
# K04075 = tRNA(Ile)-lysidine synthase: Aminoacyl-tRNA biosynthesis (https://doi.org/10.1074/jbc.M809013200)
# K05992 = maltogenic alpha-amylase: Starch & sucrose metabolism (https://doi.org/10.1111/j.1574-6968.1988.tb03149.x)
# K07043 = UTP pyrophosphatase: Pyrimidine metabolism (https://doi.org/10.1186/s12867-019-0127-x)
# K08161 = MFS transporter, DHA1 family, multidrug resistance protein: Multidrug resistance - efflux pump (https://doi.org/10.1128/jb.184.15.4161-4167.2002)
# K11068 = hemolysin III: Cytolysis (https://doi.org/10.1016/S0005-2736(96)00168-X)
# K01066 = acetyl esterase- Unclassified: metabolism: Unknown

# several additional K0s classed as 'Unknown' have entries on the KEGG website!
# pulled the KEGG web pages for the 13 'Unknowns' (*** = edit manually)
#     - K01303: Acylaminoacyl-peptidase (no other info)
#     - K02647: Carbohydrate diacid regulator (transcription factor)
#     - K05992: Maltogenic alpha-amylase (no other info)
#     - K06412: Stage V sporulation protein (unclassified signalling protein)
#     - K06595: Heme-based aerotactic transducer (bacterial chemotaxis protein) ***
#     - K06603: Flagellar protein FlaG (flagellar assembly) ***
#     - K07814: Cyclic di-GMP phosphodiesterase (no other info)
#     - K09749: Uncharacterised protein (function unknown)
#     - K09766: Uncharacterised protein (function unknown)
#     - K11936: Poly-beta-1,6-N-acetyl-D-glucosamine synthase (biofilm formation) ***
#     - K13244: Cyclic-di-GMP-specific phosphodiesterase (no other info)
#     - K13626: Flagellar assembly factor FliW (flagellar assembly) ***
#     - K14051: Cyclic-di-GMP phosphodiesterase Gmr (quorum sensing) ***