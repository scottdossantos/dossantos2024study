# london/europe dataset: species-level summary - BV subgroups only

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
                         "sds/cc" = "~/Documents/GitHub/dossantos2024study/",
                         "gg" = "~/Documents/0_git/projects/dossantos2024study/")

# run setup code to get batch-corrected feature tables and taxa vector
source(paste(path.to.github, "code/setup.R", sep=""))

# remove two janky London BV samples from the filtered dataset:
# v.001A - close to 100% L. gasseri with practically no BV organisms
# v.019A - around 80 % iners, ~ 20% crispatus and a tiny bit of Gardnerella
new.both.filt <- new.both.filt[,-c(9,22)]

# NOTE: these samples were very out of place in initial analyses

# load in names of samples in BV subgroup analyses
load(paste(path.to.github,"Rdata/bvSubgroup_names.Rda",sep = ""))

####################### grouping of vNumbers by species #######################

# calculate total number of reads per species across all samples
# output is a list of named integer vectors; read sum x samples, one per species
no.reads.total <- list()
for(i in levels(factor(new.tax.vec))){
  no.reads.total[[i]] <- colSums(new.both.filt[which(new.tax.vec == i),])
}

# add info for all taxa with no known taxonomy in VIRGO
no.reads.total[["No_tax_info"]]<-colSums(new.both.filt[which(is.na(new.tax.vec)),])

# convert to data frame
no.reads.df <- as.data.frame(do.call(rbind, no.reads.total))

# transpose feature table and remove unknown taxa rows
bars.df <- as.data.frame(t(no.reads.df[-65,]))
rownames(bars.df) <- NULL

# get the total number of genes in each species
n.genes <- vector()
for(i in levels(factor(new.tax.vec))){
  n.genes[i] <- length(which(new.tax.vec == i))
}

# get total number of genes in species with no known taxonomy in VIRGO
n.genes["No_tax_info"]<-length(which(is.na(new.tax.vec)))

# get taxa represented by >75 genes (ignoring genes w/ no taxonomy in VIRGO)
tax.order.75 <- names(which(n.genes[-65] > 75))

# calculate sum of reads in each sample corresponding to taxa with <75 genes
bars.df$Other <- rowSums(bars.df[, -which(colnames(bars.df) %in% tax.order.75)])

# remove taxa not in >75 gene list (now collapsed to 'Other')
bars.df <- bars.df[, which(colnames(bars.df) %in% c(tax.order.75,"Other"))]

# check total read numbers are the same in collapsed df as original 
rowSums(bars.df) == colSums(no.reads.df[-65,]) # all TRUE

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
                                  bars.df.plot$species == "Megasphaera sp." ~ glue("<i>Megasphaera</i> sp."),
                                  bars.df.plot$species == "Megasphaera genomosp." ~ glue("<i>Megasphaera</i> genomosp."),
                                  bars.df.plot$species == "Prevotella sp." ~ glue("<i>Prevotella</i> sp."),
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

# sanity check that taxa color rownames and unique plot df species are equal
unique(bars.df.plot$species) == rownames(bars.cols) # all TRUE

# create column for re-ordering samples, based on order of samples in the 
# BV subgroup functional heatmap
bars.df.plot$order.samples <- case_when(bars.df.plot$sample == "v.010A" ~ 1,
                                        bars.df.plot$sample == "v.006A" ~ 2,
                                        bars.df.plot$sample == "v.009A" ~ 3,
                                        bars.df.plot$sample == "v.a6_1" ~ 4,
                                        bars.df.plot$sample == "v.012A" ~ 5,
                                        bars.df.plot$sample == "v.a6_6" ~ 6,
                                        bars.df.plot$sample == "v.a13_22" ~ 7,
                                        bars.df.plot$sample == "v.a17_4" ~ 8,
                                        bars.df.plot$sample == "v.a5_1" ~ 9,
                                        bars.df.plot$sample == "v.a15_6" ~ 10,
                                        bars.df.plot$sample == "v.012B" ~ 11,
                                        bars.df.plot$sample == "v.a8_6" ~ 12,
                                        bars.df.plot$sample == "v.a13_19" ~ 13,
                                        bars.df.plot$sample == "v.a5_12" ~ 14,
                                        bars.df.plot$sample == "v.a7_1" ~ 15,
                                        bars.df.plot$sample == "v.018B" ~ 16,
                                        bars.df.plot$sample == "v.017B" ~ 17)

# filter observations for samples in BV subgroup analyses
bars.df.plot <- bars.df.plot %>% 
  filter(sample %in% rownames(lon.eur.bv.groups))

# plot microbiome composition from metatranscriptome reads of known taxonomy
# for taxa represented by >75 genes
# png(paste(path.to.github,
#           "/figs_for_paper/lon_eur_BVSubgroups_relAbund.png", sep = ""),
#     units = "in", height = 4, width = 8, res = 400)
bars.df.plot %>% 
  ggplot(aes(x = fct_reorder(sample, order.samples, sum),
             y = rel.abund, 
             fill = rev(fct_reorder(species, order.species, sum))))+
  geom_bar(position = "stack", stat = "identity",
           colour = "black", linewidth=0.3)+
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