# virginia dataset: species-level summary - BV subgroups only

#################################### setup ####################################

library(dplyr) # for data manipulation
library(tidyr) # for data reshaping
library(ggplot2) # for bar plot
library(forcats) # for reordering factors in plot
library(ggtext) # for editing ggplot labels
library(glue) # for interpreting expressions as R code

# set path to the github repository
user <- "sds/cc"
path.to.github <- switch(user,
                         "sds/cc" = "~/Documents/GitHub/metatranscriptome/",
                         "gg" = "~/Documents/0_git/projects/metatranscriptome/")

# load in filtered virginia feature table, virginia taxa vector, virginia
# bv subgroup vector and taxa colours df
load(paste0(path.to.github,"Rdata/virginia.filt.Rda"))
load(paste0(path.to.github,"Rdata/tax.vec.virginia.Rda"))
load(paste0(path.to.github,"Rdata/virginia.bv.groups.Rda"))
load(paste0(path.to.github,"Rdata/virginia.bv.groups.hclust.Rda"))
tax.colors <- read.table(paste0(path.to.github,"code/species_colors.txt"),
                         sep="\t", header=T, row.names=1, stringsAsFactors=F)

####################### grouping of vNumbers by species #######################

# calculate total number of reads per species across all samples
# output is a list of named integer vectors; read sum x samples, one per species
no.reads.total <- list()
for(i in levels(factor(tax.vec.virginia))){
  no.reads.total[[i]] <- colSums(virginia.filt[which(tax.vec.virginia == i),])
}

# collapse list to data frame
no.reads.df <- as.data.frame(do.call(rbind, no.reads.total))

# transpose feature table and remove sample names
bars.df <- as.data.frame(t(no.reads.df))
rownames(bars.df) <- NULL

# get the total number of genes in each species
n.genes <- vector()
for(i in levels(factor(tax.vec.virginia))){
  n.genes[i] <- length(which(tax.vec.virginia == i))
}

# get taxa represented by >75 genes (ignoring genes w/ no taxonomy in VIRGO)
tax.order.75 <- names(which(n.genes > 75))

# calculate sum of reads in each sample corresponding to taxa with <75 genes
bars.df$Other <- rowSums(bars.df[, -which(colnames(bars.df) %in% tax.order.75)])

# remove taxa not in >75 gene list (now collapsed to 'Other')
bars.df <- bars.df[, which(colnames(bars.df) %in% c(tax.order.75,"Other"))]

# check total read numbers are the same in collapsed df as original 
length(which(rowSums(bars.df) == colSums(no.reads.df))) # all TRUE

# convert to proportions
bars.df.prop <- as.data.frame(t(apply(bars.df, 1, function(x){
  x/sum(x)
})))

# add sample names to new column and relocate to front of data frame
bars.df.prop$sample <- colnames(no.reads.df)
bars.df.prop <- bars.df.prop %>% 
  relocate(sample, .before = Atopobium_vaginae)

# pivot dataframe into tidy format
bars.df.plot <- bars.df.prop %>% 
  pivot_longer(-c(sample), names_to = "species", values_to = "rel.abund")

# convert relative abundance to percentage (i.e. x100)
bars.df.plot$rel.abund <- bars.df.plot$rel.abund*100

# add number for numerical ordering of taxa, with 'Other' last
bars.df.plot$order.species <- NA
num <- 1
for(i in unique(bars.df.plot$species)){
  bars.df.plot$order.species <- case_when(bars.df.plot$species == i ~ as.integer(num),
                                          .default = bars.df.plot$order.species)
  num <- num +1
}

# gsub underscores to spaces and edit species names to include markdown 
# formatting for a call to ggtext's element_markdown() function
bars.df.plot$species <- gsub("_", " ", bars.df.plot$species)
bars.df.plot$species <- case_when(bars.df.plot$species == "BVAB1" ~ "BVAB1",
                                  bars.df.plot$species == "Megasphaera genomosp." ~ glue("<i>Megasphaera</i> genomosp."),
                                  bars.df.plot$species == "Other" ~ "Other",
                                  .default = glue("<i>{bars.df.plot$species}</i>"))

# make copy of taxa colours df, edit row names to match those of barplot df and
# add a grey colour for 'Other'
bars.cols <- tax.colors[colnames(bars.df.prop[,-1]),1, drop = FALSE]
bars.cols[17,1] <- "grey80"
rownames(bars.cols)[17] <- "Other"
rownames(bars.cols) <- gsub("_", " ", rownames(bars.cols))
rownames(bars.cols) <- case_when(rownames(bars.cols) == "BVAB1" ~ "BVAB1",
                                 rownames(bars.cols) == "Megasphaera sp." ~ glue("<i>Megasphaera</i> sp."),
                                 rownames(bars.cols) == "Megasphaera genomosp." ~ glue("<i>Megasphaera</i> genomosp."),
                                 rownames(bars.cols) == "Prevotella sp." ~ glue("<i>Prevotella</i> sp."),
                                 rownames(bars.cols) == "Other" ~ "Other",
                                 .default = glue("<i>{rownames(bars.cols)}</i>"))

# make copy of taxa colours df, edit row names to match those of barplot df and
# add a grey colour for 'Other'
bars.cols <- tax.colors[colnames(bars.df.prop[,-1]),1, drop = FALSE]
rownames(bars.cols)[c(3,5,11,18)] <- c("Finegoldia_magna",
                                       "Lactobacillus_coleohominis",
                                       "Peptoniphilus_harei",
                                       "Other")

bars.cols[c(3,5,11,18),1] <- c("yellow3",
                               "steelblue",
                               "chocolate4",
                               "grey80")

rownames(bars.cols) <- gsub("_", " ", rownames(bars.cols))
rownames(bars.cols) <- case_when(rownames(bars.cols) == "BVAB1" ~ "BVAB1",
                                 rownames(bars.cols) == "Megasphaera genomosp." ~ glue("<i>Megasphaera</i> genomosp."),
                                 rownames(bars.cols) == "Other" ~ "Other",
                                 .default = glue("<i>{rownames(bars.cols)}</i>"))


# sanity check that taxa color rownames and unique plot df species are equal
unique(bars.df.plot$species) == rownames(bars.cols) # all TRUE

# filter observations for samples in BV subgroup analyses
bars.df.plot <- bars.df.plot %>% 
  filter(sample %in% rownames(virginia.bv.groups))

# assign ordering integer to observations from each sample based on their order
# in the heatmap (taken from hclust object)
bars.df.plot$order.samples <- NA
for(i in 1:nrow(virginia.bv.groups)){
  index.in.sample.vec <- virginia.bv.groups.hclust$order[i]
  bars.df.plot$order.samples <- case_when(bars.df.plot$sample == unique(bars.df.plot$sample)[index.in.sample.vec] ~ i,
                                          .default = bars.df.plot$order.samples)
}

# plot microbiome composition from metatranscriptome reads of known taxonomy
# for taxa represented by >75 genes
# png(paste(path.to.github,
#           "/figs_for_paper/virginia_BVSubgroups_relAbund.png", sep = ""),
#     units = "in", height = 4, width = 8, res = 400)
bars.df.plot %>% 
  ggplot(aes(x = fct_reorder(sample, order.samples, sum),
             y = rel.abund, 
             fill = rev(fct_reorder(species, order.species, sum))))+
  geom_bar(position = "stack", stat = "identity",
           colour = "black", linewidth=0.1)+
  scale_fill_manual(values = rev(bars.cols$speciescolor),
                    labels=rev(rownames(bars.cols)))+
  scale_y_continuous(limits = c(0,101), expand=c(0,0))+
  xlab("Samples")+ylab("Relative abundance (%)")+ labs(fill="<b>Species</b>")+
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
