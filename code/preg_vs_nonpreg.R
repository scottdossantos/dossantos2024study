# analysis of pregnant vs. non-pregnant vaginal metatranscriptomes

#################################### setup ####################################

library(ALDEx2) # for CLR transformation (see notes in script!)
library(dplyr) # for data manipulation
library(tidyr) # for data manipulation
library(ggplot2) # for data visualisation
library(forcats) # for reordering factors
library(ggtext) # for editing ggplot labels
library(glue) # for interpreting expressions as R code
library(ggh4x) # for ggplot customisations (facet background colours)
library(CoDaSeq) # for pca plots


# set user and path to github
user <- "sds/cc"
path.to.github <- switch(user,
                         "sds/cc" = "~/Documents/GitHub/dossantos2024study/",
                         "gg" = "~/Documents/0_git/projects/dossantos2024study/")

# load in filtered merged dataset (lon/eur/vir) and corresponding taxa vector,
# health/BV group info for virginia samples, and taxa colour table
load(paste(path.to.github, "Rdata/preg.np.filt.bc.Rda", sep = ""))
load(paste(path.to.github, "Rdata/preg.np.tax.vec.Rda", sep = ""))
load(paste(path.to.github, "Rdata/virginia.groups.Rda", sep = ""))
tax.colors <- read.table(paste(path.to.github,"code/species_colors.txt", sep=""),
                         sep="\t", header=T, row.names=1, stringsAsFactors=F)

####################### CLR summary: data transformation #######################

# make metadata object for merged datasets
preg.np.meta <- data.frame(dataset = rep(c("London", "Europe", "Virginia"),
                                         c(22, 22, 297)),
                           health = c(rep(c("Healthy","BV","BV","Healthy"),
                                          c(8,14,14,8)),
                                      virginia.groups$groups.2),
                           row.names = colnames(preg.np.filt.bc))

# CLR transformation of merged lon/eur/vir batch-corrected data with scale-sim
# NOTE: This command is EXTREMELY memory intensive and threw a 'vector memory
#       exhausted' error when run on a 2023 MacBook Pro (M2 chip, 16 GB RAM).
#       As a workaround, one can follow instructions by SO user 'Purrsia' on
#       this SO thread: https://stackoverflow.com/questions/51295402/ to  
#       increase the memory available to RStudio by calling 'edit_R_environ()'
#       from the 'usethis' package, adding "R_MAX_VSIZE=100Gb" to the first line
#       of .Renviron as suggested and restarting RStudio. The resulting CLR
#       object is 23.7 GB! Naturally, the .Rda file is too large to be stored in
#       this study's github repository, but we include the code used to generate
#       it on our local machine for transparency and reproducibility, should 
#       anyone want to replicate the  analysis. We have, however, included the
#       median CLR sample summary produced by a truncated version of the 
#       aldex.effect() function from the large aldex.clr object (see below and 
#       'code/trunc.aldex.effect.R' for more information).
#       - Scott Dos Santos, 1st Dec 2023.

# if the effect size & median clr sample summary object already exists in the
# 'Rdata/' directory, load it in. Otherwise, if you want to replicate the
# analysis, set the seed for RNG, make a matrix of scale values, then run
# ALDEx2 on all 3 merged datasets (THIS ASSUMES YOU HAVE READ THE NOTE ABOVE!)
if(file.exists(paste(path.to.github, "Rdata/preg.np.clr.e.Rda", sep = ""))){
  load(paste(path.to.github, "Rdata/preg.np.clr.e.Rda", sep = ""))
} else{
  
  # set seed for RNG
  set.seed(2023)
  
  # make matrix of scale values from log normal distribution (stdev = 0.5,
  # difference = 15%) for ALDEx2
  mu.matrix <- aldex.makeScaleMatrix(gamma = 0.5, mu = c(1,1.15),
                                     conditions = preg.np.meta$health)
  
  # clr transform data, implementing scale simulation
  preg.np.clr <- aldex.clr(reads = preg.np.filt.bc, conds = preg.np.meta$health,
                           gamma = mu.matrix, verbose = TRUE)

}

# calculate summary of median CLR values across all 128 MC instances for CLR 
# transformed data
# NOTE: This command causes RStudio to throw a fatal error after ~5 minutes on 
#       a 2023 MacBook Pro (M2 chip, 16 GB RAM), likely due to running out of 
#       memory. The error occurs AFTER calculating median CLR values for all
#       genes within and between groups (and, if specified, the median CLR
#       values for all genes in each sample), during effect size calculation.
#       As a workaround, we truncated the aldex.effect() code accordingly to be
#       used below as part of this script and saved a copy in this study's 
#       github repository under 'code/trunc.aldex.effect.R'. See this file for
#       more information. Below, we load in the output obtained when setting
#       a seed, loading in the large aldex.clr object produced above on a local
#       machine, and running trunc.aldex.effect() on the latter, but show the
#       code used to generate it for transparency.

if(file.exists(paste(path.to.github, "Rdata/preg.np.clr.e.Rda", sep = ""))){
  load(paste(path.to.github, "Rdata/preg.np.clr.e.Rda", sep = ""))
} else{
  
  # load in truncated aldex.effect() function from 'code/'
  source(paste(path.to.github, "code/trunc_aldex_effect.R", sep = ""))
  
  # set seed for RNG
  set.seed(2023)
  
  # run truncated aldex.effect() code to generate median CLR for all features
  # (this assumes you have generated the aldex.clr object as above)
  preg.np.clr.e <- trunc.aldex.effect(preg.np.clr, verbose = TRUE,
                                      include.sample.summary = TRUE)
  # save to .Rda
  save(preg.np.clr.e,
       file = paste(path.to.github, "Rdata/preg.np.clr.e.Rda", sep = ""))
}

# extract sample summary table and remove 'rab.sample.' from column names
preg.np.clr <- preg.np.clr.e[, grep("rab.sample.", colnames(preg.np.clr.e))]
colnames(preg.np.clr) <- gsub("rab.sample.", "", colnames(preg.np.clr))

################### CLR summary: data cleaning and plotting ###################

# create columns corresponding to vNumber and taxonomy in sample summary table,
# then move columns to positions 1 and 2 of data frame and remove any rows
# where species is NA (i.e. we don't know the species assignment of these genes)
preg.np.clr$vnum <- rownames(preg.np.clr)
preg.np.clr$species <- preg.np.tax.vec
preg.np.clr <- preg.np.clr %>% 
  relocate(c(vnum, species), .before = h.001B) %>% 
  filter(!is.na(species))

# convert to long format for ggplot, specifying each observation's vnumber and
# species
preg.np.plot <- preg.np.clr %>% 
  pivot_longer(-c(vnum, species),
               names_to = "sample",
               values_to = "median.clr")

# add dataset information for all observations 
# NOTE: A full stop/period ('.') is a regex for any character, while a backslash
#       ('\') indicates the start of an escape sequence ('\n' for new line or 
#       '\t' for tab). R will try to check what '\.' signifies; however it  
#       doesn't mean anything and will not be recognised (type '?Quotes' for 
#       more info). Therefore, one needs to escape the escape character using a
#       double backslash, allowing '.' to be interpreted literally
preg.np.plot$dataset <- case_when(preg.np.plot$sample %in% grep("h\\.0", colnames(preg.np.filt.bc), value = T) ~ "London",
                                  preg.np.plot$sample %in% grep("v\\.0", colnames(preg.np.filt.bc), value = T) ~ "London",
                                  preg.np.plot$sample %in% grep("h\\.b", colnames(preg.np.filt.bc), value = T) ~ "Europe",
                                  preg.np.plot$sample %in% grep("v\\.a", colnames(preg.np.filt.bc), value = T) ~ "Europe",
                                  .default = "Virginia")

# get total number of genes per species in merged dataset
n.genes.all <- vector()
for(i in levels(factor(preg.np.tax.vec))){
  n.genes.all[i] <- length(which(preg.np.tax.vec == i))
}

# get species represented by at least 75 genes
n.genes.75 <- n.genes.all[which(n.genes.all >=75)]

# filter plotting dataframe for species not represented by at least 75 genes
# (goes from 7,369,010 observations to 7,151,793 observations)
preg.np.plot <- preg.np.plot %>% 
  filter(species %in% names(n.genes.75))

# get colours associated with species represented by >75 genes, and add colours
# for Finegoldia magna & Peptoniphilus harei
tax.cols.75 <- tax.colors[names(n.genes.75),1]
tax.cols.75[c(4,12)] <- c("yellow3", "chocolate4")
names(tax.cols.75) <- names(n.genes.75)

# add colour information to observations (takes a few seconds!)
preg.np.plot$colour <- NA
for(i in names(tax.cols.75)){
  preg.np.plot$colour <- case_when(preg.np.plot$species == i & preg.np.plot$dataset %in% c("London","Europe") ~ tax.cols.75[i],
                                   preg.np.plot$dataset == "Virginia" ~ rgb(0,0,0,0.025),
                                   .default = preg.np.plot$colour)
}

# add size information to observations
preg.np.plot$size <- case_when(preg.np.plot$dataset == "London" ~ 0.5,
                               preg.np.plot$dataset == "Europe" ~ 0.5,
                               .default = 0.1)

# add shape information to observations
shape.healthy <- grep("h\\.", colnames(preg.np.filt.bc), value = TRUE)
shape.bv <- grep("v\\.", colnames(preg.np.filt.bc), value = TRUE)

# squares are healthy, triangles are bv
preg.np.plot$shape <- case_when(preg.np.plot$sample %in% shape.healthy ~ 15,
                                preg.np.plot$sample %in% shape.bv ~ 17,
                                .default = 19)

# edit the species names to replace underscores with spaces 
preg.np.plot$species <- gsub("_", " ", preg.np.plot$species)

# create markdown expressions to italicise species names where appropriate (i.e.
# not "BVAB1" or "genomosp."). Here, glue() allows evaluation of expressions in
# braces as R code
preg.np.plot$species <- case_when(preg.np.plot$species == "BVAB1" ~ "BVAB1",
                                  preg.np.plot$species == "Megasphaera sp." ~ glue("<i>Megasphaera</i> sp."),
                                  preg.np.plot$species == "Megasphaera genomosp." ~ glue("<i>Megasphaera</i> genomosp."),
                                  preg.np.plot$species == "Prevotella sp." ~ glue("<i>Prevotella</i> sp."),
                                  .default = glue("<i>{preg.np.plot$species}</i>"))

# make stripchart of median CLR values for species (takes a while to render all
# seven million points!)
# png(paste(path.to.github,
#           "figs_for_paper/preg_np_medianCLR_species_strip.png", sep = ""),
#     units = "in", height = 6, width = 12, res = 400)

preg.np.plot %>% 
  ggplot(aes(x = median.clr, y = species))+
  geom_point(shape = preg.np.plot$shape, color = preg.np.plot$colour,
             size = preg.np.plot$size, position = "jitter")+
  xlab("Median CLR")+
  theme_bw()+
  theme(axis.text.x = element_text(colour = "black", size = 9),
        axis.text.y = element_text(colour = "black"),
        axis.title.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  theme(axis.text.y = element_markdown())

# dev.off()

# the number and density of data points on the stripchart makes it difficult to
# gauge the true distribution. Violin plots would make for a better figure, but
# requires further transformation

# edit names attribute of taxa colours so they match the edited species names
# in markdown format (necessary for a call to 'scale_fill_manual()' later)
names(tax.cols.75) <- gsub("_", " ", names(tax.cols.75))

names(tax.cols.75) <- case_when(names(tax.cols.75) == "BVAB1" ~ "BVAB1",
                                names(tax.cols.75) == "Megasphaera sp." ~ glue("<i>Megasphaera</i> sp."),
                                names(tax.cols.75) == "Megasphaera genomosp." ~ glue("<i>Megasphaera</i> genomosp."),
                                names(tax.cols.75) == "Prevotella sp." ~ glue("<i>Prevotella</i> sp."),
                                .default = glue("<i>{names(tax.cols.75)}</i>"))

# edit dataset column to include bold & italic text simultaneously, AND specify
# text colour in markdown format
preg.np.plot$dataset <- case_when(preg.np.plot$dataset == "Europe" ~ glue("<b style = 'color: white'>Europe (<i>n</i> = 44)<b/>"),
                                  preg.np.plot$dataset == "London" ~ glue("<b style = 'color: black'>London (<i>n</i> = 44)<b/>"),
                                  .default = glue("<b style: = 'color: black'>Virginia (<i>n</i> = 297)<b/>"))

# pull sample names of healthy vs. bv samples across all 3 datasets, then add a
# column to the un-summarised plot dataframe corresponding to health status (to
# facet across a second dimension: health vs. bv)
group.health <- c(grep("h", rownames(preg.np.meta), value = TRUE),
                  rownames(virginia.groups)[which(virginia.groups$groups.2 == "Healthy")])

group.bv <- c(grep("v", rownames(preg.np.meta), value = TRUE),
              rownames(virginia.groups)[which(virginia.groups$groups.2 == "BV")])

preg.np.plot$group <- case_when(preg.np.plot$sample %in% group.health ~ "Healthy",
                                preg.np.plot$sample %in% group.bv ~ "BV")

# calculate median CLR values by dataset and species and summarise median, min
# and max median CLR values per species per dataset (NO LONGER USED!)
preg.np.plot.violin <- preg.np.plot %>% 
  group_by(species, dataset, group) %>% 
  summarise(median = median(median.clr),
            min = min(median.clr),
            max = max(median.clr))

# add a column to the summarised plot dataframe corresponding to violin fill 
# colours (for adding text labels of median CLR with 'geom_text()')
preg.np.plot.violin$col <- NA
for(i in preg.np.plot.violin$species){
  preg.np.plot.violin$col <- case_when(preg.np.plot.violin$species == i ~ tax.cols.75[i],
                                    .default = preg.np.plot.violin$col)
}

# the colour "ivory2" is not easily distinguished from the white background, so
# obtain the rgb code, reduce all values by 10 and replace "ivory2" with this
# colour in the plotting dataframe
col2rgb("ivory2") # equivalent to 'rgb(238, 238, 224)'
rgb(218, 218, 204, maxColorValue = 255) # returns "#DADACC"

preg.np.plot.violin$col <- gsub("ivory2", "#DADACC", preg.np.plot.violin$col)

# add a column to the summarised plot dataframe for setting the 'nudge_x' value
# in the 'geom_text()' call (-1.125 is fine for all except Europe, Healthy,
# A. tradius and Europe, Healthy P. buccalis), then edit accordingly
preg.np.plot.violin$nudge.x <- case_when(preg.np.plot.violin$species =="<i>Anaerococcus tetradius</i>" &
                                           preg.np.plot.violin$dataset == "<b>Europe (<i>n</i> = 44)<b/>"&
                                           preg.np.plot.violin$group == "BV" ~ 0.5, 
                                         preg.np.plot.violin$species =="<i>Prevotella buccalis</i>" &
                                           preg.np.plot.violin$dataset == "<b>Europe (<i>n</i> = 44)<b/>"&
                                           preg.np.plot.violin$group == "Healthy" ~ 1,
                                         .default = -1.125)

# create a call to a 'strip_' function for ggh4x to interpret in order to set  
# facet strip background colours for both x and y dimensions while retaining
# markdown formatting of bold/italics
preg.np.strip <- strip_themed(background_x = elem_list_rect(fill = c("snow", "black","grey")),
                              text_x = element_markdown(face = "bold", size = 10.5),
                              background_y = elem_list_rect(fill = c("goldenrod2","skyblue1")),
                              text_y = element_markdown(face = "bold", size = 10.5))

# visualise CLR distributions per species as violin plots faceted by dataset
# (takes a few seconds!). Commented out 'geom_text()' call for plotting median 
# CLR values as it distracts from the message of the figure 
# png(paste(path.to.github,
#           "figs_for_paper/preg_np_medianCLR_species_violin.png", sep = ""),
#     units = "in", height = 8.5, width = 12, res = 400)

preg.np.plot %>% 
  ggplot(aes(x = median.clr, 
             y = fct_relevel(species, "BVAB1", after = 2,), 
             fill= species))+
  geom_violin(draw_quantiles = 0.5, linewidth = 0.35,
              scale = "width", show.legend = FALSE)+
  # geom_text(data = preg.np.plot.violin,
  #           aes(x = preg.np.plot.violin$max,
  #               y = fct_relevel(preg.np.plot.violin$species,
  #                               "BVAB1", after = 2)),
  #           label = format(signif(preg.np.plot.violin$median, 3),
  #                          nsmall = 1),
  #           colour = preg.np.plot.violin$col, size = 3,
  #           nudge_x = preg.np.plot.violin$nudge.x,
  #           nudge_y = 0.325, show.legend = FALSE)+
  scale_y_discrete(limits = rev)+
  scale_fill_manual(values = tax.cols.75)+
  xlab("Median CLR abundance") +
  theme_bw()+
  theme(axis.text.y = element_markdown(colour = "black", size = 10),
        axis.text.x = element_text(colour = "black", size = 10),
        axis.title.x = element_text(colour = "black"),
        axis.title.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  facet_grid2(group ~ dataset, strip = preg.np.strip)

# dev.off()

################### biplot: pregnant vs. non-pregnant - all ###################

# extract sample summary table and remove 'rab.sample.' from column names again
preg.np.clr <- preg.np.clr.e[, grep("rab.sample.", colnames(preg.np.clr.e))]
colnames(preg.np.clr) <- gsub("rab.sample.", "", colnames(preg.np.clr))

# PCA of CLR-transformed data (expects samples in rows)
pca.preg.np <- prcomp(t(preg.np.clr))

# get list of indices of each species' genes in the taxa vector/dataframe
pca.ind.load <- list()
for(i in levels(factor(names(n.genes.75)))){
  pca.ind.load[[i]] <- which(preg.np.tax.vec == i)
}

# add all other genes to "Other"
pca.ind.load[["Other"]] <- setdiff(1:length(preg.np.tax.vec),
                               unlist(pca.ind.load))

# add a transparent grey to the taxa colour vector to represent 'Other'
pca.load.cols <- c(tax.cols.75, rgb(0,0,0,0.05))

# make list of group indices
pca.ind.grp <- list(Pregnant = 45:341,
                `Non-pregnant` = 1:44)

# make vector of grouping colours
pca.grp.cols <- c("grey70", "black")

# make biplot of pregnant vs. non-pregnant vaginal metatranscriptomes
# png(paste(path.to.github, "figs_for_paper/preg_np_biplot_species.png", sep = ""),
#     units = "in", height = 6, width = 12, res = 400)

codaSeq.PCAplot(pca.preg.np, plot.groups = TRUE, plot.loadings = TRUE,
                plot.ellipses = NULL, plot.density = "loadings", PC = c(1,2),
                grp = pca.ind.grp, grp.col = pca.grp.cols, grp.sym = "text",
                grp.cex = 0.475, load.grp = pca.ind.load, load.sym = rep(19,21),
                load.col = pca.load.cols, load.cex = 0.2, plot.legend = "groups",
                leg.position = "topleft", leg.cex = 0.8, leg.columns = 1,
                title = "All datasets - Pregnant vs. Non-pregnant")

legend(x = 0.8025, y = 1.08, legend = names(n.genes.75), box.lwd = 0,
       pch = 19, col = tax.cols.75, cex = 0.6, ncol = 2)

# dev.off()
