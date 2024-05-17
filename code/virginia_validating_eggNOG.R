# repeating BV subgroup analysis using eggNOG terms - Virginia dataset

#################################### setup ####################################

#load required packages
library(dplyr)
library(CoDaSeq)
library(ALDEx2)
library(pheatmap)

# set and get wd (adjust if needed)
setwd("~/Documents/GitHub/dossantos2024study")
wd <- paste0(getwd(),"/")

# load in virginia datasets and taxa vectors, virginia (sub)group dfs, heatmap
# column colours list/df for BV subgroups, and KO & eggNOG lookup tables
load(paste0(wd, "Rdata/virginia.filt.Rda"))
load(paste0(wd, "Rdata/virginia.groups.Rda"))
load(paste0(wd, "Rdata/virginia.bv.groups.Rda"))
load(paste0(wd, "Rdata/bv12.meta.Rda"))
load(paste0(wd, "Rdata/bv12.hm.cols.Rda"))

eggnog <- read.table(paste0(wd, "1_VIRGO//3.eggnog.NOG.txt"), header = FALSE,
                     sep = "\t",quote = '', row.names = 2)

ko <- read.table(paste0(wd, "1_VIRGO/8.A.kegg.ortholog.txt"), header = FALSE,
                 sep = "\t",quote = '', row.names = 1)

# subset the eggnog lookup table by the vNumbers present in the filtered 
# virginia dataset

if(file.exists(paste0(wd,"Rdata/virginia.eggnog.filt.Rda"))){
  load(paste0(wd, "Rdata/virginia.eggnog.filt.Rda"))
} else{
  virginia.eggnog.filt <- eggnog %>% 
    filter(rownames(eggnog) %in% rownames(virginia.filt))
  
  save(virginia.eggnog.filt, file = paste0(wd, "Rdata/virginia.eggnog.filt.Rda"))
}

######################### editing eggNOG descriptions #########################

# aggregate vNumbers by eggNOG term, or load in final object if already made
# NOTE: 2,868 vNumbers without any eggNOG annotation info
if(file.exists(paste0(wd,"Rdata/egg.virg.filt.Rda"))){
  load(paste0(wd,"Rdata/virginia.eggnog.filt.Rda"))
  load(paste0(wd,"Rdata/egg.virg.filt.Rda"))
} else{
  # subset eggNOG table by vNumbers in fc (filtered/conj) virginia dataset
  virginia.eggnog.filt <- eggnog[rownames(virginia.filt),]
  
  # every vNumber included in 'COG1215', which corresponds to the exopolysacch.
  # biosynthesis KO 'K11936' from the previoius analysis, maps to a different
  # KO. COG1215 is designated as 'glycosyl transferase' which is non-specific.
  # I will create a special term just for the vNumber which corresponds to the
  # KO 'K11936' which was identified in the previous BV subgroup analyses
  # ('V1363555')
  virginia.eggnog.filt["V1363555", "V7"] <- "COGexopol"
  
  # aggregate fc virginia dataset by eggNOG terms
  egg.virg.filt <- aggregate(virginia.filt, by = list(virginia.eggnog.filt$V7), FUN = sum)
  rownames(egg.virg.filt) <- egg.virg.filt$Group.1
  egg.virg.filt$Group.1 <- NULL 
  
  save(virginia.eggnog.filt, file = paste0(wd, "Rdata/virginia.eggnog.filt.Rda"))
  save(egg.virg.filt, file = paste0(wd, "Rdata/egg.virg.filt.Rda"))
}

# make an eggNOG COG -> pathway data frame
# egg.path <- data.frame (matrix(data = NA, nrow = length(rownames(egg.virg.fc)),
#                                ncol = 1, dimnames = list(r=rownames(egg.virg.fc),
#                                                          c = "pathway")))
# egg.path$pathway <- virginia.eggnog.fc$V6[match(rownames(egg.path), virginia.eggnog.fc$V7)]

# these pathway names are horrendously non-specific or are nonsensical (e.g.
# pathway for COG1620 is "L-Lactate", but manual search of eggNOG v6 shows it
# should be "L-Lactate permease"). I will manually edit these to standardise
# and verify every one using the EggNOG v4.5 which was used to build VIRGO 
# NOTE: After trying to do this, I quickly changed my mind. There are ~2500 
#       terms and the assigned EggNOG descriptions don't necessarily match the
#       assignment given by VIRGO, so I have decided to run ALDEx2 and then
#       match up the terms with their assignments for the significant terms
#
# write.csv(egg.path, file = paste0(wd, "Rdata/egg.path.original.csv"))

###################### virginia: ALDEx2 w/ eggNOG - BV1/2 ######################

# subset virginia fc data for only BV1/BV2 samples and check for rows with 0 sum
egg.virg.bv <- egg.virg.filt[, rownames(virginia.bv.groups)]
any(rowSums(egg.virg.bv) == 0) # returns FALSE

# do the CLR transformation and check data is centred, or load in from Rda 
# object if it already exists
if(file.exists(paste0(wd, "Rdata/egg.virg.bv.clr.all.Rda"))){
  load(paste0(wd,"Rdata/egg.virg.bv.clr.all.Rda"))
} else {
  # generate scale matrix with difference of 15% and stdev of 0.5
  set.seed(2024)
  scale.bv12 <- aldex.makeScaleMatrix(gamma = 0.5, mu = c(1, 1.15),
                                      conditions = virginia.bv.groups$group,
                                      mc.samples = 128)
  
  # clr transform eggnog-aggregated virginia fc data
  egg.virg.bv.clr <- aldex.clr(reads = egg.virg.bv, mc.samples = 128, 
                               conds = virginia.bv.groups$group, denom = "all",
                               gamma = scale.bv12, verbose = T)
  
  # calcuate effect size and sample summary
  egg.virg.bv.clr.e <- aldex.effect(clr = egg.virg.bv.clr, verbose = TRUE,
                                    include.sample.summary = TRUE)
  
  # calculate pvals and bind dataframes
  egg.virg.bv.clr.t <- aldex.ttest(clr = egg.virg.bv.clr, verbose = T)
  
  egg.virg.bv.clr.all <- cbind(egg.virg.bv.clr.e, egg.virg.bv.clr.t)
  
  save(egg.virg.bv.clr.all, file = paste0(wd, "Rdata/egg.virg.bv.clr.all.Rda"))
}

# make MW plot
# png(paste0(wd,"figs_for_paper/R_FigsVirginia/eggNOG_virg_bv12_difAbund_MW.png"), 
#     units = "in", height = 6, width =12 , res = 400)

aldex.plot(egg.virg.bv.clr.all, type = "MW", test = "welch", cutoff.pval = 0.01)
title("ALDEx2 - Virginia dataset: eggNOG terms (BV1 vs. BV2, 15% scale diff.)",
      adj=0, line= 0.75)
abline(h = 0, lty = 2, col = "blue")
points(egg.virg.bv.clr.all$diff.win[which(egg.virg.bv.clr.all$we.eBH <0.01 & abs(egg.virg.bv.clr.all$effect) >1)],
       egg.virg.bv.clr.all$diff.btw[which(egg.virg.bv.clr.all$we.eBH <0.01 & abs(egg.virg.bv.clr.all$effect) >1)],
       pch = 19, cex = 0.4, col = "purple")
# points(egg.virg.bv.clr.all["COGexopol","diff.win"],
#        egg.virg.bv.clr.all["COGexopol","diff.btw"],
#        pch = 1, col = "green3", cex = 0.75)
legend("topleft", legend = c("Non-sig.", "P <0.01", "Abs. effect >1"),
       col = c("black", "red", "purple"), pch = 19, cex = 0.6)

# dev.off()

# Looks like data is well-centred and we get same picture as when using KO terms

# extract significant eggNOG terms (P <0.01 and abs. effect >1)
egg.virg.sig <- data.frame(eggnog = rownames(egg.virg.bv.clr.all)[which(egg.virg.bv.clr.all$we.eBH <0.01 & abs(egg.virg.bv.clr.all$effect) >1)])
egg.virg.sig$effect <- egg.virg.bv.clr.all$effect[which(egg.virg.bv.clr.all$we.eBH <0.01 & abs(egg.virg.bv.clr.all$effect) >1)]

# write to csv so I can pull the corrected terms from the .xlsx file
# write.csv(egg.virg.sig, file = paste0(wd,"Rdata/egg.virg.sig.csv"))

# read back in the edited significant eggNOG terms -> pathway table (still in
# the same order, so can pull effect sizes from old data frame)
egg.virg.sig.ed <- read.csv(paste0(wd,"Rdata/egg.virg.sig.edited.csv"), 
                            row.names = 1, header = T)

# check that edited df has same number of eggNOG terms and in the same order
length(which(egg.virg.sig$eggnog == rownames(egg.virg.sig.ed))) # all TRUE

# calculate total reads in each BV subgroup for every significant eggNOG term
egg.virg.sig.ed$reads.bv1 <- rowSums(egg.virg.bv[rownames(egg.virg.sig.ed),
                                                 which(virginia.bv.groups$group == "BV Subgroup 1")])

egg.virg.sig.ed$reads.bv2 <- rowSums(egg.virg.bv[rownames(egg.virg.sig.ed),
                                                 which(virginia.bv.groups$group == "BV Subgroup 2")])

# there is an unknown eggNOG term '0YEGX' with ~6.85 million reads in BV2 and
# 18 thousand in BV1. Explore which vNumbers correspond to this term and pull
# their kegg KO terms and number of reads in BV1/2
egg.0yegx <- virginia.eggnog.filt[which(virginia.eggnog.filt$V7=="0YEGX"),]
egg.0yegx$ko <- ko[rownames(egg.0yegx),1]
egg.0yegx$reads.bv1 <- rowSums(virginia.filt[rownames(egg.0yegx),rownames(virginia.bv.groups)[which(virginia.bv.groups$group == "BV Subgroup 1")]])
egg.0yegx$reads.bv2 <- rowSums(virginia.filt[rownames(egg.0yegx),rownames(virginia.bv.groups)[which(virginia.bv.groups$group == "BV Subgroup 2")]])
  
# looks like most 0YEGX vNumbers have no known KO term, but two are K01227
# (mannosyl-glycoprotein endo-beta-N-acetylglucosaminidase) and one is K01448
# (N-acetylmuramoyl-L-alanine amidase; amiABC). Difference in reads between bv1
# and bv2 spans 1-4 orders of magnitude for these three. Others are larger, but
# lack a KO term.

# extract median clr values summary from aldex object
egg.virg.bv.clr.df <- egg.virg.bv.clr.all[,grep("rab.sample.", colnames(egg.virg.bv.clr.all))]
colnames(egg.virg.bv.clr.df) <- gsub("rab.sample.", "", colnames(egg.virg.bv.clr.df))

# make matrix to hold mean clr values per pathway (and remove unknowns)
egg.virg.bv.clr.path <- as.data.frame(matrix(data = NA, ncol = ncol(egg.virg.bv.clr.df),
                                             nrow = length(levels(factor(egg.virg.sig.ed$pathway))),
                                             dimnames = list(r=levels(factor(egg.virg.sig.ed$pathway)),
                                                             c=colnames(egg.virg.bv.clr.df))))
egg.virg.bv.clr.path <- egg.virg.bv.clr.path[-21,]

# subset median clr value df for only signifcant terms
egg.virg.bv.clr.df.sig <- egg.virg.bv.clr.df[rownames(egg.virg.sig.ed),]

# calculate mean clr values per pathway
for(i in rownames(egg.virg.bv.clr.path)){
  egg.virg.bv.clr.path[i,] <- colMeans(egg.virg.bv.clr.df.sig[which(egg.virg.sig.ed$pathway==i),])
}

# make z-score function and apply to mean pathway df
zscore <- function(x){
  (x -mean(x)) / sd(x)
}

egg.virg.bv.clr.path.z <- t(apply(egg.virg.bv.clr.path, 1, zscore))

# get max and min mean values
min(egg.virg.bv.clr.path.z) # -3.61
max(egg.virg.bv.clr.path.z) #  2.56

# make heatmap for significantly different eggNOG terms
# png(paste0(wd,"figs_for_paper/R_FigsVirginia/eggNOG_virg_bv12_difAbund.png"),
#     units = "in", height = 4, width =12 , res = 400)

pheatmap(egg.virg.bv.clr.path.z, clustering_method = "ward.D2", cutree_cols = 2,
         annotation_colors = bv12.hm.cols, annotation_col = bv12.meta, 
         border_color = rgb(0,0,0,0.1), show_colnames = F, cutree_rows = 2,
         treeheight_row = 0, fontsize_row = 10, breaks = seq(-3,3,0.06),
         cellheight = 10, cellwidth = 5.1)

# dev.off()

###################### virginia: ALDEx2 w/ eggNOG - H/BV ######################

# CLR transform virginia filtered dataset, or load in final object if it already
# exists
if(file.exists(paste0(wd, "Rdata/egg.virg.filt.clr.all.Rda"))){
  load(paste0(wd, "Rdata/egg.virg.filt.clr.all.Rda"))
} else{
  # set seed for RNG
  set.seed(2024)
  
  # make scale matrix for 297 samples
  scale.h.bv <- aldex.makeScaleMatrix(gamma = 0.5, c(1,1.15), mc.samples = 128,
                                      conditions = virginia.groups$groups.2)
  
  # clr transform virginia eggNOG-aggregated data with 15% scale offset
  egg.virg.filt.clr <- aldex.clr(reads = egg.virg.filt, mc.samples = 128,
                                 conds = virginia.groups$groups.2, denom = "all",
                                 gamma = scale.h.bv, verbose = TRUE)
  
  # calculate effect sizes and t-test vals, then bind dfs and save object
  egg.virg.filt.clr.e <- aldex.effect(clr = egg.virg.filt.clr, verbose = TRUE,
                                      include.sample.summary = TRUE)
  
  egg.virg.filt.clr.t <- aldex.ttest(clr = egg.virg.filt.clr,
                                     verbose = TRUE)
  
  egg.virg.filt.clr.all <- cbind(egg.virg.filt.clr.e, egg.virg.filt.clr.t)
  
  save(egg.virg.filt.clr.all, 
       file = paste0(wd, "Rdata/egg.virg.filt.clr.all.Rda"))
}

# make MW plot
# png(paste0(wd,"figs_for_paper/R_FigsLondonEurope/eggNOG_virg_HBV_difAbund_MW.png"),
#     units = "in", height = 6, width =12 , res = 400)

aldex.plot(egg.virg.filt.clr.all, type = "MW", test = "welch", cutoff.pval = 0.01)
title("ALDEx2 - Virginia dataset: eggNOG terms (H vs. BV, 15% scale diff.)",
      adj=0, line= 0.75)
abline(h = 0, lty = 2, col = "blue")
points(egg.virg.filt.clr.all$diff.win[which(egg.virg.filt.clr.all$we.eBH <0.01 & abs(egg.virg.filt.clr.all$effect) >1)],
       egg.virg.filt.clr.all$diff.btw[which(egg.virg.filt.clr.all$we.eBH <0.01 & abs(egg.virg.filt.clr.all$effect) >1)],
       pch = 19, cex = 0.4, col = "purple")
legend("bottomleft", legend = c("Non-sig.", "P <0.01", "Abs. effect >1"),
       col = c("black", "red", "purple"), pch = 19, cex = 0.6)

# dev.off()

# Looks like data is, again, well-centred and we get same picture as when using 
# KO terms with health and bv

# extract significant eggNOG terms (P <0.01 and abs. effect >1), their effect
# sizes and the read counts in healthy and bv samples
egg.virg.filt.sig <- data.frame(eggnog = rownames(egg.virg.filt.clr.all)[which(egg.virg.filt.clr.all$we.eBH <0.01 & abs(egg.virg.filt.clr.all$effect) >1)])
egg.virg.filt.sig$effect <- egg.virg.filt.clr.all$effect[which(egg.virg.filt.clr.all$we.eBH <0.01 & abs(egg.virg.filt.clr.all$effect) >1)]
egg.virg.filt.sig$reads.h <- rowSums(egg.virg.filt[egg.virg.filt.sig$eggnog,
                                                   which(virginia.groups$groups.2 == "Healthy")])
egg.virg.filt.sig$reads.bv <- rowSums(egg.virg.filt[egg.virg.filt.sig$eggnog,
                                                   which(virginia.groups$groups.2 == "BV")])

# pull the description field of their eggNOG term
egg.unique <- eggnog[-which(duplicated(eggnog$V7)),5:6]
egg.virg.filt.sig$description <- NA
for(i in 1:length(egg.virg.filt.sig$eggnog)){
  egg.virg.filt.sig$description[i] <- egg.unique$V6[grep(egg.virg.filt.sig$eggnog[i], egg.unique$V7)]
}

# If a reviewer is reading this, please don't make me repeat the above analysis
# for the health vs. BV comparison. There are several hundred significantly 
# different KO terms and I am only mildly exaggerating when I say that if I had
# to manually check every single one, I might actually die from the tediousness 
# of the task. Instead, here's evidence that the most important functions that 
# we identified in the KEGG KO analysis (e.g. CAMP resistance) are also
# significantly different between health and BV when we use EggNOG terms:

# CAMP resistance genes
# DltD: a D-alanyl lipoteichoic acid transferase up in health
egg.virg.filt.sig[807,] 

# Peptidoglycan biosynthesis
# D-alanyl-D-alanine carboxypeptidase involved in PG crosslinking, up in health
egg.virg.filt.sig[648,]

# amiABC: the N-acetylmuramoyl-L-alanine amidase up in BV
egg.virg.filt.sig[507,]

# arnT: the aminodeoxyarabinose transferase up in BV
egg.virg.filt.sig[665,]

# Butanoate metabolism
# hmgL: one of the butanoate metabolism genes up in BV
egg.virg.filt.sig[319,]

# LPS biosynthesis genes, all up in BV
egg.virg.filt.sig[c(482,484,525,526,633),]

# On the off chance this doesn't convince you, here is a picture of a giraffe
# named Geoffrey. You couldn't possibly disappoint him, could you?

##                         ~~~~~~~~~~~~~~~~
##                        |  Please don't  |
##             O O        |   make Scott   |
##             |_|      / |   repeat the   |
##           <(. .)>   /  |    analysis    |
##            ( u )   /    ~~~~~~~~~~~~~~~~
##                \\
##                 \\
##                  \\               )
##                   \\             /
##                    \\___________/
##                    /|          /|
##                   //||      ( /||
##                  // ||______//_||
##                 //  ||     //  ||
##                //   ||    //   ||
##                \\   ||    \\   ||
##                 \\  ||     \\  ||
##                 //  ||     //  ||
##                /_\ /__\   /_\ /__\