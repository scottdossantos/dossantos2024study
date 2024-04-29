# virginia dataset: named gardnerella species biplots

#################################### setup ####################################

# load required packages 
library(dplyr)
library(CoDaSeq)

# set path to the github repository
user <- "sds/cc"
path.to.github <- switch(user,
                         "sds/cc" = "~/Documents/GitHub/dossantos2024study/",
                         "gg" = "~/Documents/0_git/projects/dossantos2024study/")

# load in following .Rda objects (all created in 'lon_eur_species_biplots.R'):
#    - virginia.filt (BV-only vNumber feature table for virginia dataset)
#    - tax.vec.virginia (BV-only vNumber taxa vector for virginia dataset)
#    - virginia.groups (health vs. BV metadata for virginia dataset)

load(paste(path.to.github, "Rdata/virginia.filt.Rda", sep = ""))
load(paste(path.to.github, "Rdata/tax.vec.virginia.Rda", sep = ""))
load(paste(path.to.github, "Rdata/virginia.groups.Rda", sep = ""))

# load pre-computed taxa table (specifying column headers)
tax.table <- read.table(paste0(path.to.github,'1_VIRGO/1.taxon.tbl.txt'),
                        col.names = c("cluster", "vnumber", "taxon", "length"),
                        row.names = 2, header = F)

# load in BLAST results for genes unique to each named gardnerella spp. against
# database of V number genes present in virginia dataset (all 297)
load(paste(path.to.github, "Rdata/blast_98pc_virginia_leopoldii.Rda", sep = ""))
load(paste(path.to.github, "Rdata/blast_98pc_virginia_piotii.Rda", sep = ""))
load(paste(path.to.github, "Rdata/blast_98pc_virginia_swidsinskii.Rda", sep = ""))
load(paste(path.to.github, "Rdata/blast_98pc_virginia_vaginalis.Rda", sep = ""))

###################### renaming species in taxonomy table ######################

# get subset of tax.table with only vNumbers corresponding to gardnerella
# i.e. all gardnerella genes in virgo
tax.gardnerella<-tax.table[grep("Gardnerella",tax.table$taxon),]

# for each species dataframe, pull gene lengths and percent query coverage for 
# each blast hit, add them to the corresponding species blast dataframes, then
# filter blast hits using conservative E value and query coverage thresholds
blast.dfs<- list()
for(i in ls(pattern = "blast.98.virginia.")){
  blast.dfs[[i]] <- get(i)
  blast.dfs[[i]]$vnum.length <- tax.gardnerella$length[match(get(i)$sseqid, rownames(tax.gardnerella))]
  blast.dfs[[i]]$p.coverage <- blast.dfs[[i]]$length / blast.dfs[[i]]$vnum.length * 100
  blast.dfs[[i]] <- blast.dfs[[i]][c(names(blast.dfs[[i]])[1:4], "vnum.length",
                                     "p.coverage", names(blast.dfs[[i]][5:12]))]
  assign(i, blast.dfs[[i]])
  assign(i, get(i) %>% filter(evalue < 1e-20 & p.coverage >90))
}

# rename genes in virgo taxonomy table based on blast hits of genes unique to
# each named gardnerella species (takes ~10 seconds)
for(i in 1:nrow(tax.gardnerella)){
  tax.gardnerella$taxon[i] <- case_when(rownames(tax.gardnerella)[i] %in% blast.98.virginia.leopoldii$sseqid ~ "Gardnerella_leopoldii",
                                        rownames(tax.gardnerella)[i] %in% blast.98.virginia.piotii$sseqid ~ "Gardnerella_piotii",
                                        rownames(tax.gardnerella)[i] %in% blast.98.virginia.swidsinskii$sseqid ~ "Gardnerella_swidsinskii",
                                        rownames(tax.gardnerella)[i] %in% blast.98.virginia.vaginalis$sseqid ~ "Gardnerella_vaginalis",
                                        .default = "Gardnerella")
}

# check to see how many V numbers have been renamed- should be 153
length(grep("Gardnerella_", tax.gardnerella$taxon)) # 153: no duplicate issues

######################### data transformation for PCA #########################

# subset the virginia dataset for only BV samples, then check for rows with
# zero sums
virginia.filt.bv <- virginia.filt[, which(virginia.groups$groups.2 == "BV")]
any(rowSums(virginia.filt.bv) == 0) # FALSE

# set seed for RNG and perform zero replacement
set.seed(2023)
virginia.filt.bv.zr <- cmultRepl(t(virginia.filt.bv), label = 0,
                                 method = "CZM", z.warning = 0.99)

# CLR transformation and PCA
virginia.filt.bv.clr <- apply(t(virginia.filt.bv.zr), 2, function(x) log(x) - mean(log(x)))
pca.virginia.bv <- prcomp(t(virginia.filt.bv.clr))

######################### virginia gardnerella biplot #########################

# set colours for new gardnerella biplot
tax.cols.gard <-c ("lightgrey", "red", "orange",
                   "lightpink2", "indianred4", rgb(0,0,0,0))

names(tax.cols.gard)<- c("Gardnerella",
                         "Gardnerella_leopoldii",
                         "Gardnerella_piotii",
                         "Gardnerella_swidsinskii",
                         "Gardnerella_vaginalis",
                         "Other")

# make a gardnerella-only taxa vector
tax.vec.virginia.gard <- tax.vec.virginia[grep("Gardnerella", tax.vec.virginia)]

# pull vNumbers and renamed (rn) taxonomic IDs from renamed gardnerella df
tax.gardnerella.rn <- tax.gardnerella[names(tax.vec.virginia.gard),]

# replace old Gardnerella vaginalis labels in virginia taxa vector with new ones
for (i in 1:length(tax.vec.virginia)) {
  if(names(tax.vec.virginia)[i] %in% rownames(tax.gardnerella.rn)){
    vnum <- names(tax.vec.virginia)[i]
    vnum.tax <- grep(vnum, rownames(tax.gardnerella.rn))
    tax.vec.virginia[i] <- tax.gardnerella.rn$taxon[vnum.tax]
  }
}

# rename object containing edited virginia taxa vector and save as .Rda
tax.vec.virginia.gv <- tax.vec.virginia
save(tax.vec.virginia.gv,
     file = paste0(path.to.github, "Rdata/tax.vec.virginia.gv.Rda"))

# pull the indices of the V numbers corresponding to the Gardnerella genus and
# each named species in the virginia dataset (for colouring species)
ind.gard <- list()
for (i in levels(factor(tax.gardnerella.rn$taxon))) {
  ind.gard[[i]] <- which(tax.vec.virginia == i)
}

# pull indices of all other vNumbers (those NOT found in the above index list)
ind.gard[["Other taxa"]] <- setdiff(c(1:nrow(virginia.filt.bv)),
                                    unlist(ind.gard))

# initiate graphics device, plot the BV-only biplot, emphasising only the 
# vNumbers corresponding to Gardnerella and save as .png
# png(filename = paste(path.to.github,
#                      "figs_for_paper/suppl_virginia_biplot_gardnerella.png", sep = ""),
#     units = "in", height = 9, width = 12, res = 400)

codaSeq.PCAplot(pca.virginia.bv, plot.groups = FALSE, plot.loadings = TRUE,
                plot.ellipses = NULL, plot.density = "loadings", grp.cex = 0.3,
                load.grp = ind.gard, load.cex = 0.5, load.col = tax.cols.gard,
                load.sym = 19, plot.legend = "loadings", leg.position = "topleft",
                leg.cex = 0.9, leg.columns = 1, PC = c(1,3),
                title= "Virginia: BV - Gardnerella spp.")

# dev.off()
