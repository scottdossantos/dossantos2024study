# Making a heatmap
#Heatmap_Europe_London

user <- "sds/cc"
path.to.github <- switch(user,
                         "sds/cc" = "~/Documents/GitHub/metatranscriptome/",
                         "gg" = "~/Documents/0_git/projects/metatranscriptome/")

#since we need to go "back" a folder to metatranscriptome I had to re-assign working directory
# you will need to change this to your computer
setwd("/Users/claracopeland/Documents/GitHub/metatranscriptome")

################################ Setup #######################################
source("code/setup.R")
library("plyr")
library("dendextend")
library("pheatmap")

############################### ADLEX Scaling #################################

# if the ADLEX scaling has been done, and is saved as xt.m.all, than that is loaded
# if not, then the analysis is done - it just takes a little bit of time
if(file.exists(paste(path.to.github,'Rdata/xt.m.all.Rda',sep = ""))){
  load(paste(path.to.github,'Rdata/xt.m.all.Rda',sep = ""))
} else{
  devtools::install_github('ggloor/ALDEx2_bioc')
  library(ALDEx2)
  library(CoDaSeq)
  
  # set seed and load data table as K#'s
  set.seed(2023)
  load('~/Documents/GitHub/metatranscriptome/Rdata/ko.both.Rda')
  
  # v.001A - close to 100% L. gasseri with practically no BV organisms
  # v.019A - around 80 % iners, ~ 20% crispatus and a tiny bit of Gardnerella
  conds.ko <- c(rep('H',8), rep('B',12), rep('B',14), rep('H', 8)) 
  
  ko.both<-ko.both[,-c(9,22)]
  
  mu = c(1,1.15)
  mu.mod <- aldex.makeScaleMatrix(gamma=0.75, mu=mu, conds.ko)
  xt.m <- aldex.clr(ko.both, conds.ko, gamma=mu.mod)
  xt.m.e <- aldex.effect(xt.m, include.sample.summary = TRUE)
  hist(xt.m.e$diff.btw, breaks=99)
  plot(density(xt.m.e$diff.btw))
  abline(v=0, lty=2, lwd=2, col='red')
  xt.m.t <- aldex.ttest(xt.m)
  xt.m.all <- cbind(xt.m.e, xt.m.t)
  
  aldex.plot(xt.m.all, xlim=c(0.3,9), cutoff.pval=0.01)
  title('C: ALDEx2 both scaled', adj=0, line= 0.8)
  
  save(xt.m.all, file='Rdata/xt.m.all.Rda')
}

#subset a list of K#'s by effect size conditions
names.KO.bv <- rownames(xt.m.all)[xt.m.all$effect < -1]
effect.KO.bv <- xt.m.all$effect < -1

names.KO.h <- rownames(xt.m.all)[xt.m.all$effect > 1]
effect.KO.h <- xt.m.all$effect > 1

# Note: Need to run this code to get annotation bar in heat map.
#source("code/Fig2A_heatmap_LondonEurope.R")
# if the metadata file exists, load it
if(file.exists(paste(path.to.github,'Rdata/hm.metadata.Rda',sep = ""))){
  load(paste(path.to.github,'Rdata/hm.metadata.Rda',sep = ""))
}

############################ Reference Tables ################################

if(file.exists(paste(path.to.github,'Rdata/master_table.Rda',sep = ""))){
  load(paste(path.to.github,'Rdata/master_table.Rda',sep = ""))
} else{
  #Master table
  # This code produces a master table saved as file='Rdata/master_table.Rda'
  
  KO.table <- read.table(paste(locn,"1_VIRGO/8.A.kegg.ortholog.txt", sep=""), 
                         header=F, row.names=1, sep="\t", stringsAsFactors=F, quote='')
  KO.table$row_names <- row.names(KO.table)
  
  path.table <- read.table(paste(locn,'1_VIRGO/8.C.kegg.pathway.copy.txt', sep=""), 
                           sep="\t", header=T, row.names=1, fill=TRUE)
  path.table$row_names <- row.names(path.table)
  
  tax.table$row_names <- row.names(tax.table)
  
  #Make table of K# with taxon
  KO.tax <- merge(KO.table, tax.table, by.x = "row_names", by.y = "row_names", all = T)
  
  # ref.table should have 3 variables - ko, pathway, map
  master_table <- merge(KO.table, path.table, by.x = "row_names", by.y = "row_names", all = T)
  master_table <- master_table[!duplicated(master_table$V2),]
  
  #################### Re-assigning pathways ###########################
  master_table$pathway <- replace(master_table$pathway, master_table$V2 %in% c("K03367", "K03739","K03740", "K14188" ,"K14198", "K01401", "K14205"),
                                  "Cationic antimicrobial peptide (CAMP) resistance")
  master_table$pathway <- replace(master_table$pathway,master_table$V2 %in% c("K00688", "K00693"),"Starch and sucrose metabolism")
  master_table$pathway <- replace(master_table$pathway, master_table$V2 %in% c("K03841"),"Glycolysis / Gluconeogenesis")
  master_table$pathway <- replace(master_table$pathway, master_table$V2 %in% c("K11262"),"Fatty acid biosynthesis")
  master_table$pathway <- replace(master_table$pathway, master_table$V2 %in% c("K02991"),"Ribosome")
  master_table$pathway <- replace(master_table$pathway, master_table$V2 %in% c("K01897"),"Fatty acid metabolism")
  master_table$pathway <- replace(master_table$pathway, master_table$V2 %in% c("K07359"),"AMPK signaling pathway")
  master_table$pathway <- replace(master_table$pathway, master_table$V2 %in% c("K01759"),"Pyruvate metabolism")
  master_table$pathway <- replace(master_table$pathway, master_table$V2 %in% c("K13288"),"RNA degradation")
  master_table$pathway <- replace(master_table$pathway, master_table$V2 %in% c("K01580"),"Alanine aspartate and glutamate metabolism")
  master_table$pathway <- replace(master_table$pathway, master_table$V2 %in% c("K03364"),"Cell Cycle")
  master_table$pathway <- replace(master_table$pathway, master_table$V2 %in% c("K01078"),"Thiamine metabolism")
  master_table$pathway <- replace(master_table$pathway, master_table$V2 %in% c("K04079"),"Protein processing in endoplasmic reticulum")
  master_table$pathway <- replace(master_table$pathway, master_table$V2 %in% c("K00863"),"Fructose and mannose metabolism")
  
  # save the table as a table
  # then we can call it easy
  save(master_table, file='Rdata/master_table.Rda')
}

# this line is used to extract the value we want to make the heatmap with
# Choose your adventure:
# your heatmap can either be the K0 numbers, or pathway level
# ref.table <- subset(master_table, select = c(V2))
ref.table <- subset(master_table, select = c(V2, pathway))

rownames(ref.table) <- ref.table$V2
ref.table <- ref.table[ -c(1) ]


######################## MAKING THE HEAT MAP ##########################

# make a column called row_names
ref.table$row_names <- rownames(ref.table)

# Subset path.table by K#'s in lists
table <- ref.table[ref.table$row_names %in% c(names.KO.bv, names.KO.h),]
#table <- ref.table

data <- data.frame(xt.m.all[,grep('rab.sample.', colnames(xt.m.all))])

# fix column names by getting rid of rab.sample. at the start
colClean <- function(x){
  colnames(x) <- gsub("^rab.sample.", "", colnames(x))
} 
data.names <- colClean(data)
colnames(data) <- data.names

# create z-scoring function
zscore <- function(x){
  (x - mean(x)) / sd(x)
}

# for first level z scoring only (KO)
# get z scores by column and turn into new dataset
#data.normed <- apply(data, 2, zscore)

# for second level z scoring only (pathway)
# renaming for multi-functionality purposes
data.normed <- data

#turn rownames to column
rn.data.normed <- as.data.frame(data.normed)
rn.data.normed$row_names <- row.names(rn.data.normed) 

# match by the column row_names
match <- match_df(rn.data.normed, table, on = "row_names")
match.in <- match[ -c(43) ]

# vector matches pathways and V number
vector <- as.vector(table[rownames(match.in),1])
names(vector) <- rownames(match.in)

#Mean gene change per pathway - colSums instead of means (counts instead of clr)
n.genes <- vector()
path.in <- levels(factor(vector))
counts <- matrix(data=NA, nrow=length(path.in), ncol=42)
rownames(counts) <- path.in
for(i in 1:length(path.in)){
  n.genes[path.in[i]] <- length(which(vector == path.in[i]))
  counts[i,1:42] <- as.numeric(colMeans(match.in[which(vector == path.in[i]),1:42]))
}

matrix <- na.omit(as.matrix(counts))
colnames(matrix) <- colnames(data.normed)

# for second level z scoring only (pathway)
matrix <- apply(matrix, 2, zscore)

hclust.normed <- hclust(dist(matrix), method = "complete")

pheatmap(matrix, annotation_col = hm.metadata[,-3], main="Heatmap")

######################## Saving plot as png ##############################
library(ggplot2)
p <- pheatmap(matrix, annotation_col = hm.metadata[,-3], main="Heatmap")
# need to change this to your computer
ggsave(filename="/Users/claracopeland/Desktop/heatmap.png", plot=p, width = 14, height = 12, units = "in")


# modified pathways:

# K00688 Insulin signaling pathway -> Insulin signaling pathway -> Starch and sucrose metabolism (this enzyme converts glycogen to a-D-Glucose-1p (glycogen metabolism))
# K03841 Insulin signaling pathway -> Glycolysis / Gluconeogenesis (this enzyme is fructose-1,6-biphosphatase I and is involved in glycolysis)
# K11262 Insulin signaling pathway -> Fatty acid biosynthesis (this enzyme converts Acetyl-CoA to Malonyl-CoA which is the universal elongation unit in fatty acid biosynthesis)
# K00693 Insulin signaling pathway -> Starch and sucrose metabolism (this enzyme is glycogen synthase which converts UDP-glycose into Amylose (glycogenesis))
# K02991 Insulin signaling pathway -> Ribosome (this protein is a component of the small subunit of the ribosome (S6e))

# K03367 Staphylococcus aureus infection -> Cationic antimicrobial peptide (CAMP) resistance (this protein is involved in modified teichoic acid biosynthesis which is involved in CAMP resistance)
# K03739 Staphylococcus aureus infection -> Cationic antimicrobial peptide (CAMP) resistance (this protein is involved in modified teichoic acid biosynthesis which is involved in CAMP resistance)
# K03740 Staphylococcus aureus infection -> Cationic antimicrobial peptide (CAMP) resistance (this protein is involved in modified teichoic acid biosynthesis which is involved in CAMP resistance)
# K14188 Staphylococcus aureus infection -> Cationic antimicrobial peptide (CAMP) resistance (this protein is involved in modified teichoic acid biosynthesis which is involved in CAMP resistance)
# K14198 Staphylococcus aureus infection -> Cationic antimicrobial peptide (CAMP) resistance (this enzyme is involved in CAMP trapping)
# K01401 Staphylococcus aureus infection -> Cationic antimicrobial peptide (CAMP) resistance (This enzyme is involved in CAMP degredation)
# K14205 Staphylococcus aureus infection -> Cationic antimicrobial peptide (CAMP) resistance (This enzyme is directly invovled in CAMP resistance)

# K01897 Adipocytokine signaling pathway -> Fatty acid metabolism (This enzyme is a long-chain acyl-CoA synthetase which is directly involved in fatty acid metabolism)
# K07359 Adipocytokine signaling pathway -> AMPK signaling pathway (This enzyme phosphorylates AMPK, activating it.)

# K01759 MAPK signaling pathway - yeast -> Pyruvate metabolism (The enzyme is lactoylglutathoine lysase which is directly invovled in pyruvate metabolism)
# K13288 Ribosome biogenesis in eukaryotes -> RNA degradation (oligoribonuclease are more universally involved in RNA degredation. It is being expressed by prokaryotes)
# K01580 Type I diabetes mellitus -> Alanine aspartate and glutamate metabolism (this enzyme is directly involved in alanine, aspartate, and glutamate metabolism)
# K03364 Progesterone-mediated oocyte maturation -> Cell cycle (This enzyme is directly involved in cell cycle)

# K01354 African trypanosomiasis -> Cationic antimicrobial peptide (CAMP) resistance (this enzyme hydrolyzes peptide bonds in Arg, Lys, and Proline-rich Antimicrobial peptides. This means it is also involved in Cationic antimicrobial peptide (CAMP) resistance)
# K01078 Tuberculosis -> Thiamine metabolism (This acid phosphatase directly synthesises thiamine)
# K04079 Prostate cancer -> Protein processing in endoplasmic reticulum (This chaperone protein is directly invovled in Protein processing in endoplasmic reticulum)
# K00863 RIG-I-like receptor signaling pathway -> Fructose and mannose metabolism (This enzyme is directly involved in Frustose and mannose metabolism)


