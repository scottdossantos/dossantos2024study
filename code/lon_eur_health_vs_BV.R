# london & europe datasets: functional heatmap (K0 -> pathway)

#################################### setup ####################################

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

# load in vNumber -> KEGG pathway lookup table from VIRGO
path.table <- read.table(paste(locn,'1_VIRGO/8.C.kegg.pathway.copy.txt', sep=""), 
                         sep="\t", header=T, row.names=1, fill=TRUE)

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

################## health vs. BV: edit lookup table pathways ##################

# if K0 -> pathway lookup table already exists as a .Rda object in 'Rdata/', 
# load it in. Otherwise, make it from scratch (takes a little while)
if(file.exists(paste(path.to.github, "Rdata/ko.both.path.Rda", sep = ""))){
  load(paste(path.to.github, "Rdata/ko.both.path.Rda", sep = ""))
}else {
  
# subset vNum -> KO and vNum -> pathway lookup tables for only those present in
# the filtered London/Europe dataset
  vnum.to.K0 <- KO[rownames(new.both.filt),]
  vnum.to.path <- path.table[rownames(new.both.filt),]
  
# check how many vNumbers lack an assigned K number
  length(which(is.na(vnum.to.K0$V2))) # 5,959 without a K number
  
# remove rows with unassigned K number
  vnum.to.K0 <- na.omit(vnum.to.K0)
  
# merge K0 & pathway lookup tables by row names
  ko.both.path <- merge(x = vnum.to.K0, y = vnum.to.path,
                        by = 'row.names', all.x = TRUE)
  
# change K number header to 'ko' and replace all NA values in the 'pathway'
# column with 'Unknown'
  colnames(ko.both.path)[2] <- 'ko'
  ko.both.path$pathway <- case_when(is.na(ko.both.path$pathway) ~ "Unknown",
                                    .default = ko.both.path$path)
  
# remove duplicate KO entries, set row names to K numbers, arrange by row and
# keep only pathway column
  ko.both.path <- ko.both.path[which(!duplicated(ko.both.path$ko)),]
  rownames(ko.both.path) <- ko.both.path$ko
  ko.both.path <- ko.both.path %>% 
    arrange(ko) %>% 
    dplyr::select(pathway)
  
# remove eukaryotic KOs (discovered during curation of 'Unknown' pathways and
# reflects changes to KO aggregated dataframe)
  ko.both.path <- ko.both.path[-which(grepl(paste("K03364","K13963","K01173",
                                                  "K12373","K14327","K00863",
                                                  "K00599","K13993","K00811",
                                                  "K03260","K00985", sep = "|"),
                                            rownames(ko.both.path))), ,
                               drop = FALSE]
  
# edit pathways which look misassigned (e.g. eukaryotic/out of place) based on
# current data in KEGG: literature evidence at end of script
  ko.both.path$pathway <- case_when(rownames(ko.both.path) %in% c("K01580","K01915","K00812","K11358","K01775",
                                                                  "K01776") ~ "Alanine aspartate and glutamate metabolism",
                                    rownames(ko.both.path) == "K01710" ~ "Amino sugar metabolism",
                                    rownames(ko.both.path) %in% c("K01358","K03544") ~ "ATP-dependent proteolysis",
                                    rownames(ko.both.path) == "K03648" ~ "Base excision repair",
                                    rownames(ko.both.path) %in% c("K01715","K14534") ~ "Butanoate metabolism",
                                    rownames(ko.both.path) %in% c("K03367","K03739","K03740","K14188","K14205") ~ "CAMP resistance",
                                    rownames(ko.both.path) %in% c("K03531","K03588","K03589","K03590") ~ "Cell division",
                                    rownames(ko.both.path) %in% c("K00031","K01679","K00029","K00024","K00174",
                                                                  "K00175","K00176","K00177","K00239","K00240",
                                                                  "K00241","K00625","K00925","K01676","K01681",
                                                                  "K01847","K01848","K01895","K01902","K01903",
                                                                  "K01966","K03737","K05606","K13788","K13581") ~ "Citrate cycle (TCA cycle)",
                                    rownames(ko.both.path) == "K12234" ~ "Coenzyme F biosynthesis",
                                    rownames(ko.both.path) == "K11031" ~ "Cytolysis",
                                    rownames(ko.both.path) %in% c("K02313","K02314") ~ "DNA replication",
                                    rownames(ko.both.path) == "K01897" ~ "Fatty acid metabolism",
                                    rownames(ko.both.path) %in% c("K02406","K03092") ~ "Flagellar assembly",
                                    rownames(ko.both.path) %in% c("K01632","K07407") ~ "Galactose metabolism",
                                    rownames(ko.both.path) == "K00864" ~ "Glycerolipid metabolism",
                                    rownames(ko.both.path) %in% c("K00058","K00600","K01079") ~ "Glycine serine and threonine metabolism",
                                    rownames(ko.both.path) %in% c("K00134", "K13953","K01835","K00873","K00927",
                                                                  "K01006","K01610","K01624","K01803","K02446",
                                                                  "K04041","K01007","K00850","K01834","") ~ "Glycolysis / Gluconeogenesis",
                                    rownames(ko.both.path) == "K00018" ~ "Glyoxylate and dicarboxylate metabolism",
                                    rownames(ko.both.path) == "K00817" ~ "Histidine metabolism",
                                    rownames(ko.both.path) == "K01338" ~ "Lon protease",
                                    rownames(ko.both.path) == "K01582" ~ "Lysine metabolism",
                                    rownames(ko.both.path) == "K00123" ~ "Methanoate metabolism",
                                    rownames(ko.both.path) %in% c("K04079","K04043","K04077") ~ "Molecular chaperones",
                                    rownames(ko.both.path) == "K01354" ~ "Oligopeptidase activity",
                                    rownames(ko.both.path) %in% c("K00297","K01491","K01938") ~ "One carbon pool by folate",
                                    rownames(ko.both.path) %in% c("K02108","K02109","K02110","K02111","K02112",
                                                                  "K02113","K02114","K02115","K02117","K02118",
                                                                  "K02119","K02120","K02122","K02123","K02124") ~ "Oxidative phosphorylation",
                                    rownames(ko.both.path) == "K01195" ~ "Pentose and glucuronate interconversions",
                                    rownames(ko.both.path) %in% c("K01621","K01783","K01807","K01808","K08094",
                                                                  "K00615") ~ "Pentose phosphate pathway",
                                    rownames(ko.both.path) == "K02563" ~ "Peptidoglycan biosynthesis",
                                    rownames(ko.both.path) %in% c("K00832","K00483","K01451") ~ "Phenylalanine tyrosine and tryptophan biosynthesis",
                                    rownames(ko.both.path) == "K00889" ~ "Phosphatidylinositol signaling system",
                                    rownames(ko.both.path) %in% c("K00088","K00857","K01951","K01488") ~ "Purine metabolism",
                                    rownames(ko.both.path) %in% c("K00757","K00760","K00876") ~ "Pyrimidine metabolism",
                                    rownames(ko.both.path) %in% c("K01759", "K01759","K01595") ~ "Pyruvate metabolism",
                                    rownames(ko.both.path) == "K11749" ~ "Regulated intramembrane proteolysis",
                                    rownames(ko.both.path) %in% c("K13288","K03685") ~ "RNA degradation",
                                    rownames(ko.both.path) == "K04564" ~ "ROS degradation",
                                    rownames(ko.both.path) == "K01186" ~ "Sialidase activity",
                                    rownames(ko.both.path) == "K00688" ~ "Starch and sucrose metabolism",
                                    rownames(ko.both.path) %in% c("K03313","K03315") ~ "Sodium-proton antiport",
                                    rownames(ko.both.path) == "K00869" ~ "Terpenoid backbone biosynthesis",
                                    rownames(ko.both.path) == "K01078" ~ "Thiamine metabolism",
                                    rownames(ko.both.path) == "K02358" ~ "Translation",
                                    rownames(ko.both.path) == "K02040" ~ "Two-component system",
                                    rownames(ko.both.path) == "K08303" ~ "Type I collagen degradation",
                                    rownames(ko.both.path) %in% c("K00568","K01661","K02548","K03183","K03688") ~ "Ubiquinone biosynthesis",
                                    rownames(ko.both.path) == "K03191" ~ "Urea transport",
                                    rownames(ko.both.path) == "K01428" ~ "Urease activity",
                                    rownames(ko.both.path) %in% c("K01066","K00532") ~ "Unknown",
                                    rownames(ko.both.path) == "K01703" ~ "Valine leucine and isoleucine biosynthesis",
                                    rownames(ko.both.path) == "K02361" ~ "Vitamin K2 biosynthesis",
                                    .default = ko.both.path$pathway)
  
# abbreviate long pathways to shorter versions
  ko.both.path$pathway <- case_when(ko.both.path$pathway == "Alanine aspartate and glutamate metabolism" ~ "Ala-Asp-Glu metabolism",
                                    ko.both.path$pathway == "Amino sugar and nucleotide sugar metabolism" ~ "Amino sugar metabolism",
                                    ko.both.path$pathway == "Arginine and proline metabolism" ~ "Arg-Pro metabolism",
                                    ko.both.path$pathway == "Biosynthesis of unsaturated fatty acids" ~ "Fatty acid biosynthesis",
                                    ko.both.path$pathway == "Cysteine and methionine metabolism" ~ "Cys-Met metabolism",
                                    ko.both.path$pathway == "Glycine serine and threonine metabolism" ~ "Gly-Ser-Thr metabolism",
                                    ko.both.path$pathway == "Histidine metabolism" ~ "His metabolism",
                                    ko.both.path$pathway == "Lipopolysaccharide biosynthesis" ~ "LPS biosynthesis",
                                    ko.both.path$pathway == "Lysine biosynthesis" ~ "Lys metabolism",
                                    ko.both.path$pathway == "Lysine degradation" ~ "Lys metabolism",
                                    ko.both.path$pathway == "Lysine metabolism" ~ "Lys metabolism",
                                    ko.both.path$pathway == "Nicotinate and nicotinamide metabolism" ~ "Vitamin B3 metabolism",
                                    ko.both.path$pathway == "Phenylalanine tyrosine and tryptophan biosynthesis" ~ "Phe-Tyr-Trp metabolism",
                                    ko.both.path$pathway == "Phosphotransferase system (PTS)" ~ "Phosphotransferase system",
                                    ko.both.path$pathway == "Porphyrin and chlorophyll metabolism" ~ "Porphyrin metabolism",
                                    ko.both.path$pathway == "Terpenoid backbone biosynthesis" ~ "Terpenoid biosynthesis",
                                    ko.both.path$pathway == "Valine leucine and isoleucine biosynthesis" ~ "Val-Leu-Ile metabolism",
                                    .default = ko.both.path$pathway)
  
# assign formerly 'Unknown' KOs a function based on current KEGG database and
# literature (see end of script for references)
  ko.both.path$pathway <- case_when(rownames(ko.both.path) == "K03281" ~ "Acid survival",
                                    rownames(ko.both.path) == "K03320" ~ "Ammonium transport",
                                    rownames(ko.both.path) == "K00561" ~ "Antimicrobial resistance",
                                    rownames(ko.both.path) == "K03402" ~ "Arg-Pro metabolism",
                                    rownames(ko.both.path) == "K03758" ~ "Arginine uptake",
                                    rownames(ko.both.path) %in% c("K01419","K03695","K03696","K03697") ~ "ATP-dependent proteolysis",
                                    rownames(ko.both.path) %in% c("K02477","K02478") ~ "Bacteriocin production",
                                    rownames(ko.both.path) == "K01448" ~ "CAMP resistance",
                                    rownames(ko.both.path) == "K03798" ~ "Cell division",
                                    rownames(ko.both.path) %in% c("K07012","K09952") ~ "CRISPR/Cas system",
                                    rownames(ko.both.path) %in% c("K01361","K01364") ~ "Cytokine degradation",
                                    rownames(ko.both.path) %in% c("K06189","K11068") ~ "Cytolysis",
                                    rownames(ko.both.path) %in% c("K03631","K04485","K01356") ~ "DNA damage repair",
                                    rownames(ko.both.path) %in% c("K02346","K02469","K02621","K02622","K03168",
                                                                  "K03169","K03346") ~ "DNA replication",
                                    rownames(ko.both.path) == "K11936" ~ "Exopolysaccharide biosynthesis",
                                    rownames(ko.both.path) %in% c("K02237","K02238","K02240","K02242","K02243",
                                                                  "K02244","K02248") ~ "Extracellular DNA uptake",
                                    rownames(ko.both.path) %in% c("K00777", "K02395","K02481","K06603","K13626") ~ "Flagellar assembly",
                                    rownames(ko.both.path) %in% c("K01795","K02444","K03436") ~ "Fructose and mannose metabolism",
                                    rownames(ko.both.path) %in% c("K02429","K02431","K07046") ~ "Fucose metabolism",
                                    rownames(ko.both.path) %in% c("K01854","K02082","K02529") ~ "Galactose metabolism",
                                    rownames(ko.both.path) == "K02445" ~ "Glycerol uptake",
                                    rownames(ko.both.path) == "K07029" ~ "Glycerolipid metabolism",
                                    rownames(ko.both.path) == "K02437" ~ "Gly-Ser-Thr metabolism",
                                    rownames(ko.both.path) %in% c("K02012","K02217","K07243") ~ "Iron acquisition",
                                    rownames(ko.both.path) == "K03303" ~ "Lactate metabolism",
                                    rownames(ko.both.path) %in% c("K04744","K06041") ~ "LPS biosynthesis",
                                    rownames(ko.both.path) %in% c("K03284","K06213","K07507") ~ "Magnesium transport",
                                    rownames(ko.both.path) %in% c("K03322","K12950") ~ "Manganese transport",
                                    rownames(ko.both.path) %in% c("K03686","K04078") ~ "Molecular chaperones",
                                    rownames(ko.both.path) == "K03297" ~ "Multidrug resistance - efflux pump",
                                    rownames(ko.both.path) %in% c("K03282","K03442") ~ "Osmotic stress adaptation",
                                    rownames(ko.both.path) %in% c("K02279","K02282","K02283","K02652") ~ "Pilus assembly",
                                    rownames(ko.both.path) == "K01991" ~ "Polysaccharide transport",
                                    rownames(ko.both.path) == "K01529" ~ "Purine metabolism",
                                    rownames(ko.both.path) == "K03445" ~ "Purine efflux",
                                    rownames(ko.both.path) == "K03756" ~ "Putrescine transport",
                                    rownames(ko.both.path) == "K02824" ~ "Pyrimidine metabolism",
                                    rownames(ko.both.path) %in% c("K01153","K01154","K01155","K01156","K03427") ~ "Restriction endonuclease",
                                    rownames(ko.both.path) %in% c("K01200", "K01215","K02438","K05992","K06896","K00701") ~ "Starch and sucrose metabolism",
                                    rownames(ko.both.path) == "K05946" ~ "Teichoic acid biosynthesis",
                                    rownames(ko.both.path) %in% c("K02483","K02484") ~ "Two-component system",
                                    rownames(ko.both.path) %in% c("K03187","K03188") ~ "Urease activity",
                                    .default = ko.both.path$pathway)
  
# finally, edit 'and' to ampersands ('&') in pathway names
  ko.both.path$pathway <- gsub("and", "&", ko.both.path$pathway)
  
# save as .Rda
  save(ko.both.path, 
       file = paste(path.to.github, "Rdata/ko.both.path.Rda",sep = ""))
}

################### health vs. BV: differential abundance #####################

# if CLR-transformation w/ scale simulation has already been done, load the .Rda
# object from 'metatranscriptome/Rdata/'. Otherwise, run ALDEx2 w/ scale sim on
# the filtered london & europe dataset which is aggregated by KO number, then
# calculate effect sizes and t-test statistics before combining the data frames
# and saving as .Rda
if(file.exists(paste(path.to.github, "Rdata/ko.both.clr.all.Rda", sep = ""))){
  load(paste(path.to.github, "Rdata/ko.both.clr.all.Rda", sep = ""))
} else{
  ko.conds <- c(rep('H',8), rep('B',12), rep('B',14), rep('H', 8)) 
  
  set.seed(2023)
  mu.matrix <- aldex.makeScaleMatrix(gamma = 0.5, mu = c(1.15,1),
                                     conditions = ko.conds)
  
  ko.both.clr <- aldex.clr(reads = ko.both, conds = ko.conds,
                           gamma = mu.matrix, verbose = TRUE)
  
  ko.both.clr.e <- aldex.effect(ko.both.clr, verbose = TRUE, 
                                include.sample.summary = TRUE)
  
  ko.both.clr.t <- aldex.ttest(ko.both.clr, verbose = TRUE)
  
  ko.both.clr.all <- cbind(ko.both.clr.e, ko.both.clr.t)
  
  save(ko.both.clr.all,
       file = paste(path.to.github, "Rdata/ko.both.clr.all.Rda", sep = ""))
}

# check the distribution of median diff.btw values (health = +ve, BV = -ve)
plot(density(ko.both.clr.all$diff.btw))
abline(v=0, lty=2, lwd=2, col='red')

# get indices for genes assigned to 'Ribosome', 'Aminoacyl-tRNA biosynthesis' 
# and 'Glycolysis / Gluconeogenesis'
path.ribo <- which(ko.both.path$pathway == "Ribosome")
path.trna <- which(ko.both.path$pathway == "Aminoacyl-tRNA biosynthesis")
path.glyc <- which(ko.both.path$pathway == "Glycolysis / Gluconeogenesis")

# plot difference vs. dispersion of median CLR values for health vs. BV and
# add coloured points for K0s with effect size >1, K0s assigned to 'Ribosome',
# 'Aminoacyl-tRNA biosynthesis' and 'Glycolysis / Gluconeogenesis' 
# Note: All purple points are also red in the un-altered aldex.plot output (i.e.
#       P < 0.01); we want to emphasise that relying on P-values isn't
#       necessarily the best idea!

# png(filename = paste(path.to.github,
#                      "figs_for_paper/lon_eur_healthBV_MW.png", sep = ""),
#     units = "in", height = 7, width = 14, res = 400)

aldex.plot(ko.both.clr.all, xlim=c(0.346,9), cutoff.pval=0.01)
title('ALDEx2 w/ ScaleSim (London/Europe filtered K0): mu = 1: 1.15, gamma = 0.75, p <0.01',
      adj=0, line= 0.8)

points(ko.both.clr.all[path.ribo ,"diff.win"],
       ko.both.clr.all[path.ribo ,"diff.btw"],
       pch= 19, col= "blue3", cex= 0.5)

points(ko.both.clr.all[path.trna,"diff.win"],
       ko.both.clr.all[path.trna,"diff.btw"],
       pch= 19, col= "royalblue1", cex= 0.5)

points(ko.both.clr.all[path.glyc,"diff.win"],
       ko.both.clr.all[path.glyc,"diff.btw"],
       pch= 19, col= "lightskyblue1", cex= 0.5)

points(ko.both.clr.all[which(abs(ko.both.clr.all$effect) >1),"diff.win"],
       ko.both.clr.all[which(abs(ko.both.clr.all$effect) >1),"diff.btw"],
       pch= 19, col= "purple3", cex= 0.5)

legend("bottomleft", pch = 19, cex = 0.75,
       col = c( "black", "purple3", "blue3", "royalblue1", "lightskyblue1"),
       legend = c("Non-significant", "Absolute effect size >1",
                  "Ribosome", "Aminoacyl-tRNA biosynthesis", "Glycolysis / Gluconeogenesis"))

# dev.off()

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

# =================== INVESTIGATING KOs OF INTEREST (START) ===================

 # make vector of KO terms to pull (2012,2217,7243 = iron, 8303 = collagenase,
# 3191 = ureI, 1429/1430 =  UreA / UreB)
ir.co.ur.ko <- c("K02012","K02217","K07243",
                 "K08303",
                 "K03191", "K01429", "K01430")

#check if all KOs are in virginia dataset 
ir.co.ur.ko %in% rownames(ko.both.all) # all TRUE

# filter vnum -> KO table for these KOs, then for vNums in database
ir.co.ur.ko.le <- KO %>% 
  filter(V2 %in% ir.co.ur.ko)

ir.co.ur.ko.le <- ir.co.ur.ko.le %>% 
  filter(rownames(ir.co.ur.ko.le) %in% rownames(new.both.data))

# make summary list object of vNumbers, their corresponding KOs and counts in
# the london/europe dataset
ir.co.ur.tax.le <- list()
for(i in ir.co.ur.ko){
  df <- ir.co.ur.ko.le %>% 
    filter(V2 == i)
  taxa <- tax.table[rownames(df), 2]
  names(taxa) <- rownames(df)
  ir.co.ur.ko.sums <- rowSums(new.both.data[rownames(df),])
  ir.co.ur.tax.le[[i]] <- data.frame(vnum= rownames(df),
                                     ko = rep(i, nrow(df)),
                                     tax.id = taxa,
                                     total.reads =ir.co.ur.ko.sums,
                                     row.names = NULL)
}

# collapse list of genes in london/europe dataset and their counts to dataframe
# and set rownames to numbers
ir.co.ur.tax.le.df <- do.call(rbind, ir.co.ur.tax.le)
rownames(ir.co.ur.tax.le.df) <- NULL

# group by KO and then sum reads
ir.co.ur.tax.le.sum <- ir.co.ur.tax.le.df %>% 
  group_by(ko,tax.id) %>% 
  summarise(sum=sum(total.reads))

# urea channel limited to sneathia amnii/sanguinigens
# non-heme ferritin = prevotella amnii / bivia / timonensis
# periplasmic iron bind = atopobium / megasphaera
# pH iron channel = gardnerella

# ==================== INVESTIGATING KOs OF INTEREST (END) ====================

##################### health vs. BV: Z-scores & heatmap #######################

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

# draw heatmap of differential pathways for lon/eur health vs. BV (non-scaled)
# png(paste(path.to.github,
#           "figs_for_paper/lon_eur_healthBV_diffAbund_K0.png", sep = ""),
#     units = "in", height = 10, width = 12.25, res = 400)

pheatmap(ko.both.clr.path.z, cutree_cols = 2, cutree_rows = 4,
         annotation_col = hm.metadata[,-3], annotation_colors = hm.column.cols,
         show_colnames = FALSE, treeheight_row = 0, treeheight_col = 10,
         clustering_method = "ward.D2", border_color = rgb(0,0,0,0.2),
         fontsize_row = 7.5, cellwidth = 15, cellheight = 6)

# dev.off()

# get range  of Z-scores
max(ko.both.clr.path.z) #  3.61
min(ko.both.clr.path.z) # -5.00

# draw heatmap of differential pathways for lon/eur health vs. BV (scaled)
# png(paste(path.to.github,
#           "figs_for_paper/lon_eur_healthBV_diffAbund_K0_scale.png", sep = ""),
#     units = "in", height = 10, width = 12.25, res = 400)

pheatmap(ko.both.clr.path.z, cutree_cols = 2, cutree_rows = 4,
         annotation_col = hm.metadata[,-3], annotation_colors = hm.column.cols,
         show_colnames = FALSE, treeheight_row = 0, treeheight_col = 10,
         clustering_method = "ward.D2", border_color = rgb(0,0,0,0.2),
         fontsize_row = 7.5, cellwidth = 15, cellheight = 6, breaks = seq(-5,5,0.1))

# dev.off()

# looks like two main BV subgroups exist in the data from these heatmaps, and
# splitting the dendrogram into 5 groups gives a vector of healthy & BV subgroup
# membership. Check to see if the BV subgroup vector exists; if not, make it
if(!file.exists(paste(path.to.github,"Rdata/lon.eur.bv.groups.Rda",sep = ""))){
  
# draw sample dendrogram as in heatmap (euclidean distance & Ward's method)
  hm.hclust <- hclust(dist(t(ko.both.clr.path.z), method = "euclidean"),
                      method = "ward.D2")
  
# plot dendrogram
  plot(hm.hclust, cex = 0.7, hang = -1)
  
# k = 5 splits the dendrogram into 2 main BV subgroups and 1 smaller subgroup of
# 3 BV samples (NOTE: one 'healthy' sample is included in BV1; dominated by
# P. bivia, but no BV symptoms!)

  rect.hclust(hm.hclust, k = 5)
  
# split tree into 3 groups (health/BV1/BV2), edit group names and save as .Rda
  lon.eur.bv.groups <- data.frame(group = cutree(hm.hclust, k = 5))
  lon.eur.bv.groups$group <- case_when(lon.eur.bv.groups$group == 1 ~ "Healthy",
                                       lon.eur.bv.groups$group == 2 ~ "BV1",
                                       lon.eur.bv.groups$group == 3 ~ "BV3",
                                       lon.eur.bv.groups$group == 4 ~ "BV2",
                                       lon.eur.bv.groups$group == 5 ~ "BV1")
  
# save as .Rda object
  save(lon.eur.bv.groups,
       file = paste(path.to.github, "Rdata/lon.eur.bv.groups.Rda", sep = ""))
}

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
ko.both.clr.path.z.sub[7,] <- ko.both.clr.path.z["Iron acquisition",1:42]
ko.both.clr.path.z.sub[9,] <- ko.both.clr.path.z["Fatty acid metabolism",1:42]
ko.both.clr.path.z.sub[10,] <- ko.both.clr.path.z["ATP-dependent proteolysis",1:42]
rownames(ko.both.clr.path.z.sub)[c(7,9,10)] <- c("Iron acquisition","Fatty acid metabolism","ATP-dependent proteolysis")

# draw scaled heatmap with top 10 differential pathways (for main paper)
# draw heatmap of differential pathways for lon/eur health vs. BV (scaled)
# png(paste(path.to.github,
#           "figs_for_paper/lon_eur_healthBV_diffAbund_K0_sub.png", sep = ""),
#     units = "in", height = 4.25, width = 9.75, res = 400)

pheatmap(ko.both.clr.path.z.sub, cutree_cols = 2, cutree_rows = 2,
         annotation_col = hm.metadata[,-3], annotation_colors = hm.column.cols,
         show_colnames = FALSE, treeheight_row = 0, treeheight_col = 50,
         clustering_method = "ward.D2", border_color = rgb(0,0,0,0.2), 
         breaks = seq(-4.5,4.5,0.09), fontsize = 10, cellwidth = 10, cellheight = 10)

# dev.off()

############### assigning Unknown KO terms: literature evidence ################

# K00561: ermC - 23S rRNA (adenine-N6)-dimethyltransferase; Antimicrobial resistance (https://card.mcmaster.ca/ontology/36389)
# K00701: cgt - cyclomaltodextrin glucanotransferase; Starch and sucrose metabolism (https://doi.org/10.1007/s10295-012-1185-y)
# K00777: degS - NarL family, sensor histidine kinase; Flagellar assembly ( https://doi.org/10.1111/j.1365-2958.2007.05923.x)
# K01005: tagT - polyisoprenyl-teichoic acid--peptidoglycan teichoic acid transferase; Teichoic acid biosynthesis (https://doi.org/10.1038/emboj.2011.358)
# K01056: pth - peptidyl-tRNA hydrolase, PTH1 family; Translation (https://doi.org/10.1093/nar/gkh924)
# K01153: hsdR - type I restriction enzyme, R subunit; Restriction endonuclease (https://doi.org/10.1128/mmbr.64.2.412-434.2000)
# K01154: hsdS - type I restriction enzyme, S subunit; Restriction endonuclease (https://doi.org/10.1128/mmbr.64.2.412-434.2000)
# K01155: E3.1.21.4 - type II restriction enzyme; Restriction endonuclease (https://doi.org/10.1016/0092-8674(86)90698-7)
# K01156: res - type III restriction endonuclease; Restriction endonuclease (https://doi.org/10.1016/0378-1119(93)90623-B)
# K01200: pulA - Pullulanase; Starch and sucrose metabolism (https://doi.org/10.1111/j.1365-2958.1990.tb02016.x)
# K01215: dexB - glucan 1,6-alpha-glucosidase; Starch and sucrose (https://doi.org/10.1099/00221287-139-9-2019)
# K01356: lexA - repressor lexA; DNA damage repair (https://doi.org/10.1073/pnas.78.7.4204) 
# K01361: E3.4.21.96 - lactocepin; Cytokine degradation (https://doi.org/10.1016/j.chom.2012.02.006)
# K01364: speB - streptopain; Cytokine degradation (https://doi.org/10.1371/journal.pone.0004769)
# K01419: clpQ - ATP-dependent HslUV protease; ATP-dependent proteolysis (https://doi.org/10.1038/nsmb898)
# K01448: amiABC - N-acetylmuramoyl-L-alanine amidase; CAMP resistance (https://doi.org/10.1053/j.gastro.2021.02.013)
# K01529: yjjX - inosine/xanthosine triphosphatase; Purine metabolism (https://doi.org/10.1016/j.str.2005.07.007)
# K01795: algG - mannuronan 5-epimerase; Fructose and mannose metabolism (https://doi.org/10.1128/JB.176.7.1821-1830.1994)
# K01854: glf - UDP-galactopyranose mutase; Galactose metabolism (https://doi.org/10.1016/j.jmb.2009.05.081)
# K01991: gcfE - polysaccharide biosynthesis/export protein; Polysaccharide transport (https://doi.org/10.1371/journal.pone.0259900)
# K02082: agaS - D-galactosamine 6-phosphate deaminase/isomerase; Galactose metabolism (https://doi.org/10.1186/1471-2180-13-94)
# K02237: comE - competence protein ComEC; Extracellular DNA uptake (https://doi.org/10.1128/jb.180.1.41-45.1998)
# K02238: comC - competence protein ComEC; Extracellular DNA uptake (https://doi.org/10.1128/jb.180.1.41-45.1998)
# K02240: comFA - competence protein ComFA; Extracellular DNA uptake (https://doi.org/10.1128/jb.180.1.41-45.1998)
# K02242: comFC - competence protein ComFC; Extracellular DNA uptake (https://doi.org/10.1128/jb.180.1.41-45.1998)
# K02243: comGA - competence protein ComGA; Extracellular DNA uptake (https://doi.org/10.1128/jb.180.1.41-45.1998)
# K02244: comGB - competence protein ComGB; Extracellular DNA uptake (https://doi.org/10.1128/jb.180.1.41-45.1998)
# K02248: comGF - competence protein ComGF; Extracellular DNA uptake (https://doi.org/10.1128/jb.180.1.41-45.1998)
# K02279: cpaB - pilus assembly protein CpaB; Pilus assembly (https://doi.org/10.1073/pnas.182411999)
# K02282: cpaE - pilus assembly protein CpaE; Pilus assembly (https://doi.org/10.1073/pnas.182411999)
# K02283: cpaF - pilus assembly protein CpaF; Pilus assembly (https://doi.org/10.1073/pnas.182411999)
# K02346: dinB - DNA polymerase IV; DNA replication (https://doi.org/10.1093/emboj/19.22.6259)
# K02395: flgJ - peptidoglycan hydrolase FlgJ; Flagellar assembly (https://doi.org/10.3389/fcimb.2020.00178)
# K02429: fucP - L-fucose permease; Fucose metabolism (https://doi.org/10.1073/pnas.1213445109)
# K02431: fucU - alpha-L-fucopyranose 1-epimerase (https://doi.org/10.1074/jbc.M402016200)
# K02437: gcvH - glycine cleavage system H protein; Gly-Ser-Thr metabolism (https://doi.org/10.1128%2FJB.01480-13)
# K02438: glgX - glycogen debranching enzyme; Starch and sucrose metabolism (https://doi.org/10.1128%2FJB.187.4.1465-1473.2005)
# K02444: glpR - glycerol-3-phosphate regulon repressor; Fructose and mannose metabolism (https://doi.org/10.1128/jb.00244-18)
# K02445: glpT - glycerol-3-phosphate transporter; Glycerol uptake (https://doi.org/10.1128%2FJB.00748-09)
# K02469: gyrA - DNA gyrase subunit A; DNA replication (https://doi.org/10.1016/j.biochi.2007.02.012)
# K02470: gyrB - DNA gyrase subunit B; DNA replication (https://doi.org/10.1016/j.biochi.2007.02.012)
# K02477: no_symbol - two-component system, LytTR family, response regulator; Bacteriocin production (https://doi.org/10.1371%2Fjournal.pgen.1007709)
# K02478: no_symbol - two-component system, LytTR family, sensor kinase; Bacteriocin production (https://doi.org/10.1371%2Fjournal.pgen.1007709)
# K02481: flgR - two-component system, NtrC family, response regulator; Flagellar assembly (https://doi.org/10.1074/jbc.m400357200)
# K02483: no_symbol - two-component system, OmpR family, response regulator; Two-component system 
# K02484: no_symbol - two-component system, OmpR family, sensor kinase; Two-component system
# K02529: galR - galactose operon repressor; Galactose metabolism (https://doi.org/10.1111/j.1365-2958.1991.tb00728.x)
# K02530: lacR - lactose phosphotransferase system repressor; Lactose metabolism (https://doi.org/10.1128%2FAEM.01370-14)
# K02621: parC - topoisomerase IV subunit A; DNA replication (https://doi.org/10.1021/bi026352v)
# K02622: parE - topoisomerase IV subunit B; DNA replication (https://doi.org/10.1093/nar/gkr018)
# K02652: pilB - type IV pilus assembly protein; Pilus assembly (https://doi.org/10.1128%2Fjb.172.6.2911-2919.1990)
# K02824: uraA - uracil permease; Pyrimidine metabolism (https://doi.org/10.1128/JB.177.8.2008-2013.1995) 
# K03168: topA - DNA topoisomerase I; DNA replication (https://doi.org/10.1016/S0167-4781(98)00125-0)
# K03169: topB - DNA topoisomerase III; DNA replication (https://doi.org/10.1016/S0167-4781(98)00125-0)
# K03187: ureE - urease accessory protein; Urease activity (https://doi.org/10.1128/JB.01265-06)
# K03188: ureF - urease accessory protein; Urease activity (https://doi.org/10.1128/JB.01265-06)
# K03281: clcA - chloride channel protein, CIC family; Acid survival (https://doi.org/10.4014%2Fjmb.2303.03009)
# K03282: mcsL - large conductance mechanosensitive channel' Osmotic stress adaptation (https://doi.org/10.1093/emboj/18.7.1730)
# K03284: corA - magnesium transporter; Magnesium transport (https://doi.org/10.1080/09687680701441883)
# K03297: emrE - small multidrug resistance pump; Multidrug resistance - efflux pump (https://doi.org/10.1074/jbc.M408187200)
# K03303: lutP - lactate permease; Lactate metabolism (https://doi.org/10.1128/JB.01464-08)
# K03320: amt - ammonium transporter, Amt family; Ammonium transport (https://doi.org/10.1080/09687680701388423)
# K03322: mntH - manganese transport protein; Manganese transport (https://doi.org/10.1046/j.1365-2958.2000.01922.x)
# K03346: dnaB - replication initiation and membrane attachment protein; DNA replication (https://doi.org/10.1111/j.1365-2958.2004.04451.x)
# K03402: argR - transcriptional regulator of arginine metabolism; Arg-Pro metabolism (https://doi.org/10.1016/0022-2836(92)91022-h)
# K03427: hsdM - type I restriction enzyme M protein; Restriction endonuclease (https://doi.org/10.1128/MMBR.64.2.412-434.2000)
# K03436: fruR - fructose operon transcriptional repressor; Fructose and mannose metabolism (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0274005)
# K03442: mcsS - small conductance mechanosensitive channel' Osmotic stress adaptation (https://doi.org/10.1093/emboj/18.7.1730)
# K03445: nepI - purine ribonucleoside efflux pump; Purine efflux (https://doi.org/10.1016/j.femsle.2005.06.051)
# K03631: recN - DNA repair protein RecN; DNA damage repair (https://doi.org/10.1038/ncomms15282)
# K03686: dnaJ - molecular chaperone dnaJ; Molecular chaperone (https://doi.org/10.1007/s00018-006-6192-6)
# K03695: clpB - ATP-dependent Clp protease ATP-binding subunit ClpB; ATP-dependent proteolysis (https://doi.org/10.1074/jbc.M209686200)
# K03695: clpC - ATP-dependent Clp protease ATP-binding subunit ClpC; ATP-dependent proteolysis (https://doi.org/10.1074/jbc.M209686200)
# K03697: clpE - ATP-dependent Clp protease ATP-binding subunit ClpE; ATP-dependent proteolysis (https://doi.org/10.1074/jbc.M209686200)
# K03756: potE - putrescine:ornithine antiporter; Putrescine transport (https://doi.org/10.1074/jbc.272.10.6318)
# K03758: arcD - arginine:ornithine antiporter; Arginine uptake (https://doi.org/10.1128%2Fjb.174.5.1568-1573.1992)
# K03798: FstH - cell division protease FtsH; Cell division (https://doi.org/10.1128/MMBR.67.1.52-65.2003)
# K04078: groES - chaperonin GroES; Molecular chaperones (https://doi.org/10.1074/jbc.M112.398628)
# K04485: radA - A repair protein RadA/Sms; DNA damage repair (https://doi.org/10.1074/jbc.M113.464883)
# K04744: lptD - LPS-assembly protein; LPS biosynthesis (https://doi.org/10.1128/JB.00270-08)
# K05946: tagA - N-acetyl-beta-D-mannosaminyltransferase; Teichoic acid biosynthesis (https://doi.org/10.1099/00221287-148-3-815)
# K05992: amyM - maltogenic alpha-amylase; Starch and sucrose metabolism (https://doi.org/10.1111/j.1574-6968.1988.tb03149.x)
# K06041: kdsD - arabinose-5-phosphate isomerase; LPS biosynthesis (https://doi.org/10.1016/j.bbrc.2009.07.154)
# K06189: tlyC - hemolysin (HlyC) family protein; Cytolysis (https://doi.org/10.1006/mpat.1994.1028)
# K06213: mgtE - magnesium transporter; Magnesium transport (https://doi.org/10.1038/emboj.2009.288)
# K06603: flaG - flagellar protein FlaG; Flagellar assembly (https://doi.org/10.1111/cmi.12886)
# K06896: mapP - maltose 6'-phosphate phosphatase; Starch and sucrose metabolism (https://doi.org/10.1111/mmi.12183)
# K07012: cas3 - CRISPR-associated endonuclease/helicase Cas3; CRISPR/Cas system (https://doi.org/10.1016/j.molcel.2012.03.018)
# K07029: dagK - diacylglycerol kinase; Glycerolipid metabolism (https://doi.org/10.1111/j.1365-2958.2006.05174.x)
# K07046: E3.1.1.120 - L-fucono-1,5-lactonase; Fucose metabolism (https://doi.org/10.1021/bi3015554)
# K07243: efeU - high-affinity iron transporter; Iron acquisition (https://doi.org/10.1146/annurev.nutr.18.1.441)
# K07507: mgtC - Mg2+ transporter-C (MgtC) family protein; Magnesium transport (https://doi.org/10.1128/IAI.66.8.3802-3809.1998)
# K09952: csn1 - CRISPR-associated endonuclease Csn1; CRISPR/Cas system (https://doi.org/10.1038/nature09886)
# K11068: hlyIII - hemolysin III; Cytolysis (https://doi.org/10.1007/s00284-004-4288-5)
# K11936: pgaC - poly-beta-1,6-N-acetyl-D-glucosamine synthase: Exopolysaccharide biosynthesis (https://doi.org/10.1128/JB.01920-07)
# K12950: ctpC - manganese-transporting P-type ATPase C; Manganese transport (https://doi.org/10.1074/jbc.M112.448175)
# K13626: fliW - flagellar assembly factor FliW; Flagellar assembly (https://doi.org/10.1128/JB.00820-06)
# K12960: mtaD - 5-methylthioadenosine/S-adenosylhomocysteine deaminase; Cys-Met metabolism (https://doi.org/10.1128/JB.01308-13)