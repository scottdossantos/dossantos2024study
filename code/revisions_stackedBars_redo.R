# remake of BV stacked bar charts for lon/eur and virginia data - new and old
# filtering

#################################### setup ####################################

library(dplyr) # for data manipulation
library(CoDaSeq) # for filtering
library(tidyr) # for data reshaping
library(ggplot2) # for bar plot
library(forcats) # for reordering factors in plot
library(ggtext) # for editing ggplot labels
library(glue) # for interpreting expressions as R code

# gh path
setwd("~/Documents/GitHub/dossantos2024study/")
github <- "~/Documents/GitHub/dossantos2024study/"

# load in original filtered and non-filtered data for lon/eur and virginia data
load(paste0(github, "Rdata/new.both.data.Rda")) # batch-corrected
load(paste0(github, "Rdata/virginia.data.Rda"))

load(paste0(github, "Rdata/new.both.filt.Rda")) # batch-corrected
load(paste0(github, "Rdata/virginia.filt.Rda"))

# species lookup table for VIRGO
lookup.spp <- read.table(paste0(github, "1_VIRGO/1.taxon.tbl.txt"),
                         header = F, sep = "\t", quote = "", row.names = 2)

# species colour lookup table 
lookup.col <- read.table(paste0(github, "Rdata/species_colors_updated.txt"),
                         header = F, sep = "\t", quote = "", row.names = 1)

# filter using less-stringent criteria, encompassing reviewer example
filt.ls.loneur <- codaSeq.filter(new.both.data, samples.by.row=F,
                                 min.occurrence=0.10, min.prop=0.00005)

filt.ls.virg <- codaSeq.filter(virginia.data, samples.by.row=F,
                               min.occurrence=0.10, min.prop=0.00005)

# load in taxa vectors for new-filtered and old-filtered datasets
load(paste0(github,"Rdata/tax.vec.loneur.ls.Rda"))
load(paste0(github,"Rdata/tax.vec.virg.ls.Rda"))

load(paste0(github,"Rdata/new.tax.vec.Rda"))
load(paste0(github,"Rdata/virginia.tax.vec.Rda"))


# remove janky BV samples from lon/eur
new.both.filt <- new.both.filt[,-c(9,22)]
filt.ls.loneur <- filt.ls.loneur[,-c(9,22)]

# load in names of samples in BV subgroup analyses
load(paste0(github,"Rdata/bvSubgroup_names.Rda"))

# load in hclust object from KO-based clustering of BV subgroups
load(paste0(github, "Rdata/virginia.bv.groups.hclust.Rda"))

# load in virginia BV grouping table
load(paste0(github, "Rdata/virginia.bv.groups.Rda"))

############################## London/Europe data ##############################

# calculate total number of reads per species across all samples
# output is a list of named integer vectors; read sum x samples, one per species
reads.loneur.total <- list()
for(i in levels(factor(new.tax.vec))){
  reads.loneur.total[[i]] <- colSums(new.both.filt[which(new.tax.vec == i),])
}

# add info for all taxa with no known taxonomy in VIRGO
reads.loneur.total[["Unknown"]]<-colSums(new.both.filt[which(is.na(new.tax.vec)),])

# convert to data frame
reads.loneur.df <- as.data.frame(do.call(rbind, reads.loneur.total))

# transpose feature table, retaining unknown taxa rows
bars.loneur.df <- as.data.frame(t(reads.loneur.df))
rownames(bars.loneur.df) <- NULL

# get the total number of genes in each species
n.genes.loneur <- vector()
for(i in levels(factor(new.tax.vec))){
  n.genes.loneur[i] <- length(which(new.tax.vec == i))
}

# get total number of genes in species with no known taxonomy in VIRGO
n.genes.loneur["Unknown"]<-length(which(is.na(new.tax.vec)))

# get taxa represented by >75 genes (including genes w/ no taxonomy in VIRGO)
tax.order.75.le <- names(which(n.genes.loneur > 75))

# calculate sum of reads in each sample corresponding to taxa with <75 genes
bars.loneur.df$Other <- rowSums(bars.loneur.df[, -which(colnames(bars.loneur.df) %in% tax.order.75.le)])

# remove taxa not in >75 gene list (now collapsed to 'Other') and move 'Other'
# column before 'Unknown'
bars.loneur.df <- bars.loneur.df[, which(colnames(bars.loneur.df) %in% c(tax.order.75.le,"Other"))]
bars.loneur.df <- bars.loneur.df %>% 
  relocate(Other, .before = Unknown)

# check total read numbers are the same in collapsed df as original 
rowSums(bars.loneur.df) == colSums(reads.loneur.df) # all TRUE

# convert to proportions
bars.loneur.prop.df <- as.data.frame(t(apply(bars.loneur.df, 1, function(x){
  x/sum(x)
})))

# add sample names to new column and relocate to front of data frame
bars.loneur.prop.df$sample <- colnames(reads.loneur.df)
bars.loneur.prop.df <- bars.loneur.prop.df %>% 
  relocate(sample, .before = Atopobium_vaginae)

# pivot dataframe into tidy format
bars.loneur.plot <- bars.loneur.prop.df %>% 
  pivot_longer(-c(sample), names_to = "species", values_to = "rel.abund")

# convert relative abundance to percentage (i.e. x100)
bars.loneur.plot$rel.abund <- bars.loneur.plot$rel.abund*100

# add number for numerical ordering of taxa, with 'Unknown' last
bars.loneur.plot$order.species <- NA
num <- 1
for(i in unique(bars.loneur.plot$species)){
  bars.loneur.plot$order.species <- case_when(bars.loneur.plot$species == i ~ as.integer(num),
                                          .default = bars.loneur.plot$order.species)
  num <- num +1
}

# gsub underscores to spaces and edit species names to include markdown 
# formatting for a call to ggtext's element_markdown() function
bars.loneur.plot$species <- gsub("_", " ", bars.loneur.plot$species)
bars.loneur.plot$species <- case_when(bars.loneur.plot$species == "BVAB1" ~ "BVAB1",
                                      bars.loneur.plot$species == "Megasphaera sp." ~ glue("<i>Megasphaera</i> sp."),
                                      bars.loneur.plot$species == "Megasphaera genomosp." ~ glue("<i>Megasphaera</i> genomosp."),
                                      bars.loneur.plot$species == "Prevotella sp." ~ glue("<i>Prevotella</i> sp."),
                                      bars.loneur.plot$species == "Other" ~ "Other",
                                      bars.loneur.plot$species == "Unknown" ~ "Unknown",
                                      .default = glue("<i>{bars.loneur.plot$species}</i>"))

# make copy of taxa colours df, edit row names to match those of barplot df and
# add a grey colours for 'Other' and 'Unknown'
bars.loneur.cols <- lookup.col[colnames(bars.loneur.prop.df[,-1]),1, drop = FALSE]
bars.loneur.cols[17,1] <- "grey80"
bars.loneur.cols[18,1] <- "grey35"
rownames(bars.loneur.cols)[17] <- "Other"
rownames(bars.loneur.cols)[18] <- "Unknown"
rownames(bars.loneur.cols) <- gsub("_", " ", rownames(bars.loneur.cols))
rownames(bars.loneur.cols) <- case_when(rownames(bars.loneur.cols) == "BVAB1" ~ "BVAB1",
                                 rownames(bars.loneur.cols) == "Megasphaera sp." ~ glue("<i>Megasphaera</i> sp."),
                                 rownames(bars.loneur.cols) == "Megasphaera genomosp." ~ glue("<i>Megasphaera</i> genomosp."),
                                 rownames(bars.loneur.cols) == "Prevotella sp." ~ glue("<i>Prevotella</i> sp."),
                                 rownames(bars.loneur.cols) == "Other" ~ "Other",
                                 rownames(bars.loneur.cols) == "Unknown" ~ "Unknown",
                                 .default = glue("<i>{rownames(bars.loneur.cols)}</i>"))

# sanity check that taxa color rownames and unique plot df species are equal
unique(bars.loneur.plot$species) == rownames(bars.loneur.cols) # all TRUE

# create column for re-ordering samples, based on order of samples in the 
# BV subgroup functional heatmap
bars.loneur.plot$order.samples <- case_when(bars.loneur.plot$sample == "v.010A" ~ 1,
                                        bars.loneur.plot$sample == "v.006A" ~ 2,
                                        bars.loneur.plot$sample == "v.009A" ~ 3,
                                        bars.loneur.plot$sample == "v.a6_1" ~ 4,
                                        bars.loneur.plot$sample == "v.012A" ~ 5,
                                        bars.loneur.plot$sample == "v.a6_6" ~ 6,
                                        bars.loneur.plot$sample == "v.a13_22" ~ 7,
                                        bars.loneur.plot$sample == "v.a17_4" ~ 8,
                                        bars.loneur.plot$sample == "v.a5_1" ~ 9,
                                        bars.loneur.plot$sample == "v.a15_6" ~ 10,
                                        bars.loneur.plot$sample == "v.012B" ~ 11,
                                        bars.loneur.plot$sample == "v.a8_6" ~ 12,
                                        bars.loneur.plot$sample == "v.a13_19" ~ 13,
                                        bars.loneur.plot$sample == "v.a5_12" ~ 14,
                                        bars.loneur.plot$sample == "v.a7_1" ~ 15,
                                        bars.loneur.plot$sample == "v.018B" ~ 16,
                                        bars.loneur.plot$sample == "v.017B" ~ 17)

# filter observations for samples in BV subgroup analyses
bars.loneur.plot <- bars.loneur.plot %>% 
  filter(sample %in% rownames(lon.eur.bv.groups))

# plot microbiome composition from metatranscriptome reads of known taxonomy
# for taxa represented by >75 genes
# png("~/BVSubgroups_lon_eur_relAbund_redo.png",
#     units = "in", height = 4, width = 8, res = 400)

bars.loneur.plot %>% 
  ggplot(aes(x = fct_reorder(sample, order.samples, sum),
             y = rel.abund, 
             fill = rev(fct_reorder(species, order.species, sum))))+
  geom_bar(position = "stack", stat = "identity",
           colour = "black", linewidth=0.3)+
  scale_fill_manual(values = rev(bars.loneur.cols$V2),
                    labels=rev(rownames(bars.loneur.cols)))+
  scale_y_continuous(limits = c(0,101), expand=c(0,0))+
  xlab("Samples")+ylab("Relative abundance (%)")+ labs(fill="<b>Taxon</b>")+
  theme_bw()+
  theme(legend.text = element_markdown(colour = "black", size = 10),
        legend.title = element_markdown(colour = "black", size = 11),
        legend.key.size = unit(0.4,"cm"), legend.box.spacing = unit(0.05, "cm"),
        legend.key.spacing.y = unit(0.085,"cm"),
        axis.title.x = element_text(colour = "black", size = 12),
        axis.title.y = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 11),
        axis.text.x = element_blank())+
  guides(fill = guide_legend(byrow = TRUE))

# dev.off()

################################ Virginia data ################################

# calculate total number of reads per species across all samples
# output is a list of named integer vectors; read sum x samples, one per species
reads.virg.total <- list()
for(i in levels(factor(tax.vec.virginia))){
  reads.virg.total[[i]] <- colSums(virginia.filt[which(tax.vec.virginia == i),])
}

# add info for all taxa with no known taxonomy in VIRGO
reads.virg.total[["Unknown"]]<-colSums(virginia.filt[which(is.na(tax.vec.virginia)),])

# collapse list to data frame
reads.virg.df <- as.data.frame(do.call(rbind, reads.virg.total))

# transpose feature table and remove sample names
bars.virg.df <- as.data.frame(t(reads.virg.df))
rownames(bars.virg.df) <- NULL

# get the total number of genes in each species
n.genes.virg <- vector()
for(i in levels(factor(tax.vec.virginia))){
  n.genes.virg[i] <- length(which(tax.vec.virginia == i))
}

# get total number of genes in species with no known taxonomy in VIRGO
n.genes.virg["Unknown"]<-length(which(is.na(tax.vec.virginia)))

# get taxa represented by >75 genes (including genes w/ no taxonomy in VIRGO)
tax.order.75.virg <- names(which(n.genes.virg > 75))

# calculate sum of reads in each sample corresponding to taxa with <75 genes
bars.virg.df$Other <- rowSums(bars.virg.df[, -which(colnames(bars.virg.df) %in% tax.order.75.virg)])

# remove taxa not in >75 gene list (now collapsed to 'Other')
bars.virg.df <- bars.virg.df[, which(colnames(bars.virg.df) %in% c(tax.order.75.virg,"Other"))]

# relocate 'Other' column before 'Unknown'
bars.virg.df <- bars.virg.df %>% 
  relocate(Other, .before = Unknown)

# check total read numbers are the same in collapsed df as original 
length(which(rowSums(bars.virg.df) == colSums(reads.virg.df))) # all TRUE

# convert to proportions
bars.virg.prop.df <- as.data.frame(t(apply(bars.virg.df, 1, function(x){
  x/sum(x)
})))

# add sample names to new column and relocate to front of data frame
bars.virg.prop.df$sample <- colnames(reads.virg.df)
bars.virg.prop.df <- bars.virg.prop.df %>% 
  relocate(sample, .before = Atopobium_vaginae)

# pivot dataframe into tidy format
bars.virg.plot <- bars.virg.prop.df %>% 
  pivot_longer(-c(sample), names_to = "species", values_to = "rel.abund")

# convert relative abundance to percentage (i.e. x100)
bars.virg.plot$rel.abund <- bars.virg.plot$rel.abund*100

# add number for numerical ordering of taxa, with 'Unknown' last
bars.virg.plot$order.species <- NA
num <- 1
for(i in unique(bars.virg.plot$species)){
  bars.virg.plot$order.species <- case_when(bars.virg.plot$species == i ~ as.integer(num),
                                          .default = bars.virg.plot$order.species)
  num <- num +1
}

# gsub underscores to spaces and edit species names to include markdown 
# formatting for a call to ggtext's element_markdown() function
bars.virg.plot$species <- gsub("_", " ", bars.virg.plot$species)
bars.virg.plot$species <- case_when(bars.virg.plot$species == "BVAB1" ~ "BVAB1",
                                  bars.virg.plot$species == "Megasphaera genomosp." ~ glue("<i>Megasphaera</i> genomosp."),
                                  bars.virg.plot$species == "Other" ~ "Other",
                                  bars.virg.plot$species == "Unknown" ~ "Unknown",
                                  .default = glue("<i>{bars.virg.plot$species}</i>"))

# make copy of taxa colours df, edit row names to match those of barplot df and
# add a grey colour for 'Other'
bars.virg.cols <- lookup.col[colnames(bars.virg.prop.df[,-1]),1, drop = FALSE]

rownames(bars.virg.cols)[c(3,5,11,18,19)] <- c("Finegoldia_magna",
                                               "Lactobacillus_coleohominis",
                                               "Peptoniphilus_harei",
                                               "Other", "Unknown")

bars.virg.cols[c(3,5,11,18,19),1] <- c("forestgreen", "blue", "chocolate4",
                                       "grey80", "grey35")

rownames(bars.virg.cols) <- gsub("_", " ", rownames(bars.virg.cols))
rownames(bars.virg.cols) <- case_when(rownames(bars.virg.cols) == "BVAB1" ~ "BVAB1",
                                 rownames(bars.virg.cols) == "Megasphaera sp." ~ glue("<i>Megasphaera</i> sp."),
                                 rownames(bars.virg.cols) == "Megasphaera genomosp." ~ glue("<i>Megasphaera</i> genomosp."),
                                 rownames(bars.virg.cols) == "Prevotella sp." ~ glue("<i>Prevotella</i> sp."),
                                 rownames(bars.virg.cols) == "Other" ~ "Other",
                                 rownames(bars.virg.cols) == "Unknown" ~ "Unknown",
                                 .default = glue("<i>{rownames(bars.virg.cols)}</i>"))

# sanity check that taxa color rownames and unique plot df species are equal
unique(bars.virg.plot$species) == rownames(bars.virg.cols) # all TRUE

# filter observations for samples in BV subgroup analyses
bars.virg.plot <- bars.virg.plot %>% 
  filter(sample %in% rownames(virginia.bv.groups))

# assign ordering integer to observations from each sample based on their order
# in the heatmap (taken from hclust object)
bars.virg.plot$order.samples <- NA
for(i in 1:nrow(virginia.bv.groups)){
  index.in.sample.vec <- virginia.bv.groups.hclust$order[i]
  bars.virg.plot$order.samples <- case_when(bars.virg.plot$sample == unique(bars.virg.plot$sample)[index.in.sample.vec] ~ i,
                                          .default = bars.virg.plot$order.samples)
}

# plot microbiome composition from metatranscriptome reads of known taxonomy
# for taxa represented by >75 genes
# png("~/BVSubgroups_virginia_relAbund.png",
#     units = "in", height = 4, width = 8, res = 400)

bars.virg.plot %>% 
  ggplot(aes(x = fct_reorder(sample, order.samples, sum),
             y = rel.abund, 
             fill = rev(fct_reorder(species, order.species, sum))))+
  geom_bar(position = "stack", stat = "identity",
           colour = "black", linewidth=0.1)+
  scale_fill_manual(values = rev(bars.virg.cols$V2),
                    labels=rev(rownames(bars.virg.cols)))+
  scale_y_continuous(limits = c(0,101), expand=c(0,0))+
  xlab("Samples")+ylab("Relative abundance (%)")+ labs(fill="<b>Taxon</b>")+
  theme_bw()+
  theme(legend.text = element_markdown(colour = "black", size = 10),
        legend.title = element_markdown(colour = "black", size = 11),
        legend.key.size = unit(0.4,"cm"), legend.box.spacing = unit(0.05, "cm"),
        legend.key.spacing.y = unit(0.1,"cm"),
        axis.title.x = element_text(colour = "black", size = 12),
        axis.title.y = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 11),
        axis.text.x = element_blank())+
  guides(fill = guide_legend(byrow = TRUE))

# dev.off()
