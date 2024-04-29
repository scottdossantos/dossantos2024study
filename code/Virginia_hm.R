# Making a heatmap
################################ Setup #######################################
user <- "cc"
path.to.github <- switch(user,
                         "cc" = "~/Documents/GitHub/dossantos2024study/",
                         "gg" = "~/Documents/0_git/projects/dossantos2024study/")

#since we need to go "back" a folder to metatranscriptome I had to re-assign working directory
# you will need to change this to your computer
setwd("/Users/claracopeland/Documents/GitHub/metatranscriptome")

source("code/setup.R")
library("plyr")
library("dplyr")
library("dendextend")
library("pheatmap")

load(paste(path.to.github, "Rdata/virginia.ko.path.Rda", sep = ""))
load(paste(path.to.github, "Rdata/virginia.filt.ko.clr.c.Rda", sep = ""))
virginia.filt.ko.clr <- virginia.filt.ko.clr.c

####################### Master table for virginia ############################

if(file.exists(paste(path.to.github,'Rdata/master_table_virginia.Rda',sep = ""))){
  load(paste(path.to.github,'Rdata/master_table_virginia.Rda',sep = ""))
} else{
  
  # make master table with taxon, ko, and pathway for each vnumber
  load(paste(locn,'Rdata/tax.vec.virginia.Rda',sep = ""))
  tax_virginia <- data.frame(tax.vec.virginia)
  tax_virginia$v_number <- row.names(tax_virginia)
  
  load(paste(locn,'Rdata/ko.virginia.vec.Rda',sep = ""))
  ko_virginia <- data.frame(ko.virginia.vec)
  ko_virginia$v_number <- row.names(ko_virginia)
  
  ko_tax_table_virginia <- merge(ko_virginia, tax_virginia, by.x = "v_number", by.y = "v_number", all = T)
  
  load(paste(locn,'Rdata/virginia.ko.path.Rda',sep = ""))
  path_virginia <- data.frame(virginia.ko.path)
  path_virginia$k_number <- row.names(path_virginia)
  
  master_table_virginia <- merge(ko_tax_table_virginia, path_virginia, by.x = "ko.virginia.vec", by.y = "k_number", all = T)
  
  ####### Re-assigning pathways that are mis-assigned###########
  
  master_table_virginia$pathway <- replace(master_table_virginia$pathway, master_table_virginia$ko.virginia.vec %in% c("K03367", "K03739","K03740", "K14188" ,"K14198", "K01401", "K14205"),
                                           "Cationic antimicrobial peptide (CAMP) resistance")
  master_table_virginia$pathway <- replace(master_table_virginia$pathway,master_table_virginia$ko.virginia.vec %in% c("K00688", "K00693"),"Starch and sucrose metabolism")
  master_table_virginia$pathway <- replace(master_table_virginia$pathway, master_table_virginia$ko.virginia.vec %in% c("K03841"),"Glycolysis / Gluconeogenesis")
  master_table_virginia$pathway <- replace(master_table_virginia$pathway, master_table_virginia$ko.virginia.vec %in% c("K11262"),"Fatty acid biosynthesis")
  master_table_virginia$pathway <- replace(master_table_virginia$pathway, master_table_virginia$ko.virginia.vec %in% c("K02991"),"Ribosome")
  master_table_virginia$pathway <- replace(master_table_virginia$pathway, master_table_virginia$ko.virginia.vec %in% c("K01897"),"Fatty acid metabolism")
  master_table_virginia$pathway <- replace(master_table_virginia$pathway, master_table_virginia$ko.virginia.vec %in% c("K07359"),"AMPK signaling pathway")
  master_table_virginia$pathway <- replace(master_table_virginia$pathway, master_table_virginia$ko.virginia.vec %in% c("K01759"),"Pyruvate metabolism")
  master_table_virginia$pathway <- replace(master_table_virginia$pathway, master_table_virginia$ko.virginia.vec %in% c("K13288"),"RNA degradatio")
  master_table_virginia$pathway <- replace(master_table_virginia$pathway, master_table_virginia$ko.virginia.vec %in% c("K01580"),"Alanine aspartate and glutamate metabolism")
  master_table_virginia$pathway <- replace(master_table_virginia$pathway, master_table_virginia$ko.virginia.vec %in% c("K03364"),"Cell Cycle")
  master_table_virginia$pathway <- replace(master_table_virginia$pathway, master_table_virginia$ko.virginia.vec %in% c("K01078"),"Thiamine metabolism")
  master_table_virginia$pathway <- replace(master_table_virginia$pathway, master_table_virginia$ko.virginia.vec %in% c("K04079"),"Protein processing in endoplasmic reticulum")
  master_table_virginia$pathway <- replace(master_table_virginia$pathway, master_table_virginia$ko.virginia.vec %in% c("K00863"),"Fructose and mannose metabolism")
  
  
  # save the table as a table
  # then we can call it easy
  save(master_table_virginia, file='Rdata/master_table_virginia.Rda')
}


######################### Heatmap of selected K0's ############################
# this code is to create a heatmap with specific K0's you want to look at
# virginia.filt.ko.clr is the table of data I want to look at

#this is to clean up the names of my columns (they all started with rab.sample.)
colClean <- function(x){
  colnames(x) <- gsub("^rab.sample.", "", colnames(x))
} 
data.names <- colClean(virginia.filt.ko.clr)
colnames(virginia.filt.ko.clr) <- data.names

########################## Make annotation bar ################################

load(paste(path.to.github, "Rdata/virginia.diff.meta.Rda", sep = ""))

diff.column.cols <- list(CST = c(I = viridis(7)[1], II = viridis(7)[2],
                                 III = viridis(7)[3], `III / IV` = viridis(7)[4],
                                 `III / V` = viridis(7)[5], IV = viridis(7)[6],
                                 V = viridis(7)[7]),
                         Group = c(`H-1` = "steelblue1",
                                   `H-2` = "royalblue3",
                                   `BV` = "red3"))


######## this is where you select what K0's you want to look at############


##################### To look at CAMP resistance #########################
list <- c("K01448", "K03767", "K14205", "K00677", "K14188", "K03739", "K03367", "K03740", "K01354", "K07264", "K03760", "K01364", "K12340", "K03585")
#This is where I define matrix, which is what is going into the heatmap
# this is to get rid of the columns that don't contain data I want to look at
drop <- c("diff.btw","diff.win", "effect", "overlap", "we.ep", "we.eBH", "wi.ep", "wi.eBH", "rab.all", "rab.win.BV", "rab.win.Healthy")
matrix <- virginia.filt.ko.clr[,!(colnames(virginia.filt.ko.clr) %in% drop)]

# this is normalizing all the columns using z score
matrix <- apply(matrix, 2, zscore)

#clustering the columns
hclust.normed <- hclust(dist(matrix), method = "complete")

#sub-setting matrix to the K0's selected above
matrix <- matrix[rownames(matrix) %in% list , ]

# rename K0's to the actual names
rownames(matrix)[rownames(matrix) == "K01448"] <- "amiABC"
rownames(matrix)[rownames(matrix) == "K03767"] <- "ppiA"
rownames(matrix)[rownames(matrix) == "K14205"] <- "MprF"
rownames(matrix)[rownames(matrix) == "K00677"] <- "LpxA"
rownames(matrix)[rownames(matrix) == "K14188"] <- "DltC"
rownames(matrix)[rownames(matrix) == "K03739"] <- "DltB"
rownames(matrix)[rownames(matrix) == "K03367"] <- "DltA"
rownames(matrix)[rownames(matrix) == "K03740"] <- "DltD"
rownames(matrix)[rownames(matrix) == "K01354"] <- "PtrB"
rownames(matrix)[rownames(matrix) == "K07264"] <- "AmT"
rownames(matrix)[rownames(matrix) == "K03760"] <- "EptA"
rownames(matrix)[rownames(matrix) == "K01364"] <- "SpeB"
rownames(matrix)[rownames(matrix) == "K12340"] <- "TolC"
rownames(matrix)[rownames(matrix) == "K03585"] <- "AcrA"

pheatmap(matrix, annotation_col = virginia.diff.meta, annotation_colors = diff.column.cols, main="Virginia Heatmap - CAMP resistance")

#Saving plot as png
# need to change file saving location to your computer
library(ggplot2)
p <- pheatmap(matrix, annotation_col = virginia.diff.meta, annotation_colors = diff.column.cols, main="Virginia Heatmap - CAMP resistance")
ggsave(filename="/Users/claracopeland/Desktop/CAMP_virginia_heatmap.png", plot=p, width = 14, height = 6, units = "in")



############ To look at Iron transport, Collagenase, and Urease ##########
list <- c("K02012", "K02217", "K08303", "K14048", "K03191")
#This is where I define matrix, which is what is going into the heatmap
# this is to get rid of the columns that don't contain data I want to look at
drop <- c("diff.btw","diff.win", "effect", "overlap", "we.ep", "we.eBH", "wi.ep", "wi.eBH", "rab.all", "rab.win.BV", "rab.win.Healthy")
matrix <- virginia.filt.ko.clr[,!(colnames(virginia.filt.ko.clr) %in% drop)]

# this is normalizing all the columns using z score
matrix <- apply(matrix, 2, zscore)

#clustering the columns
hclust.normed <- hclust(dist(matrix), method = "complete")

#sub-setting matrix to the K0's selected above
matrix <- matrix[rownames(matrix) %in% list , ]

# rename K0's to the actual names
rownames(matrix)[rownames(matrix) == "K08303"] <- "PrtC (Collagenase)"
rownames(matrix)[rownames(matrix) == "K02012"] <- "AfuA (Iron transport)"
rownames(matrix)[rownames(matrix) == "K02217"] <- "FtnA (non-heme ferritin)"
rownames(matrix)[rownames(matrix) == "K14048"] <- "UreAB (Urease)"
rownames(matrix)[rownames(matrix) == "K03191"] <- "UreI (Urea channel)"

#Saving plot as png
# need to change file saving location to your computer
library(ggplot2)
p <- pheatmap(matrix, annotation_col = virginia.diff.meta, annotation_colors = diff.column.cols, main="Virginia Heatmap - Iron transport, Collagenase, and Urease")
ggsave(filename="/Users/claracopeland/Desktop/Various_virginia_heatmap.png", plot=p, width = 14, height = 4, units = "in")



##################### To look at flagellar assembly #########################
list <- c("K00575", "K02389", "K02390", "K02392", "K02396", "K02397", "K02398", "K02400", "K02405", "K02406", "K02407", "K02408", "K02409")
#This is where I define matrix, which is what is going into the heatmap
# this is to get rid of the columns that don't contain data I want to look at
drop <- c("diff.btw","diff.win", "effect", "overlap", "we.ep", "we.eBH", "wi.ep", "wi.eBH", "rab.all", "rab.win.BV", "rab.win.Healthy")
matrix <- virginia.filt.ko.clr[,!(colnames(virginia.filt.ko.clr) %in% drop)]

# this is normalizing all the columns using z score
matrix <- apply(matrix, 2, zscore)

#clustering the columns
hclust.normed <- hclust(dist(matrix), method = "complete")

#sub-setting matrix to the K0's selected above
matrix <- matrix[rownames(matrix) %in% list , ]

# rename K0's to the actual names
rownames(matrix)[rownames(matrix) == "K00575"] <- "CheR (chemotaxis)"
rownames(matrix)[rownames(matrix) == "K02389"] <- "FlgD (Rod,P/L ring, hook)"
rownames(matrix)[rownames(matrix) == "K02390"] <- "FlgE (Rod,P/L ring, hook)"
rownames(matrix)[rownames(matrix) == "K02392"] <- "FlgG (Rod,P/L ring, hook)"
rownames(matrix)[rownames(matrix) == "K02396"] <- "FlgK (Hook-filament junction)"
rownames(matrix)[rownames(matrix) == "K02397"] <- "FlgL (Hook-filament junction)"
rownames(matrix)[rownames(matrix) == "K02398"] <- "FlgM (Regulator)"
rownames(matrix)[rownames(matrix) == "K02400"] <- "FlhA (Type III secretion system)"
rownames(matrix)[rownames(matrix) == "K02405"] <- "FliA (Filament Stator RNA polymerase)"
rownames(matrix)[rownames(matrix) == "K02406"] <- "FliC (Filament)"
rownames(matrix)[rownames(matrix) == "K02407"] <- "FliD (Filament cap, chaperone)"
rownames(matrix)[rownames(matrix) == "K02408"] <- "FliE (MS/C ring, Type III secretion system)"
rownames(matrix)[rownames(matrix) == "K02409"] <- "FliF (MS/C ring, Type III secretion system)"

pheatmap(matrix, annotation_column = virginia.diff.meta, annotation_colors = diff.column.cols, main="Virginia Heatmap - Flagellar Assembly")

#Saving plot as png
# need to change file saving location to your computer
library(ggplot2)
p <- pheatmap(matrix, annotation_col = virginia.diff.meta, annotation_colors = diff.column.cols, main="Virginia Heatmap - Flagellar Assembly")
ggsave(filename="/Users/claracopeland/Desktop/Flagella_virginia_heatmap.png", plot=p, width = 14, height = 6, units = "in")














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


