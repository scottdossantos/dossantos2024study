# all datasets: named gardnerella species biplots

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
#    - preg.np.preg.filt.bc (batch-correct feature table for all 3 datasets)
#    - preg.np.tax.vec (vNumber taxa vector for all 3 datasets)
#    - virginia.groups (health vs. BV metadata for all 3 datasets)

load(paste(path.to.github, "Rdata/preg.np.filt.bc.Rda", sep = ""))
load(paste(path.to.github, "Rdata/preg.np.tax.vec.Rda", sep = ""))
load(paste(path.to.github, "Rdata/virginia.groups.Rda", sep = ""))

# make metadata object for merged datasets
preg.np.meta <- data.frame(dataset = rep(c("London", "Europe", "Virginia"),
                                         c(22, 22, 297)),
                           health = c(rep(c("Healthy","BV","BV","Healthy"),
                                          c(8,14,14,8)),
                                      virginia.groups$groups.2),
                           row.names = colnames(preg.np.filt.bc))

# load pre-computed taxa table (specifying column headers)
tax.table <- read.table(paste0(path.to.github,'1_VIRGO/1.taxon.tbl.txt'),
                        col.names = c("cluster", "vnumber", "taxon", "length"),
                        row.names = 2, header = F)

# load in BLAST results for genes unique to each named gardnerella spp. against
# database of V number genes present in virginia dataset (all 297)
load(paste(path.to.github, "Rdata/blast_98pc_all3_leopoldii.Rda", sep = ""))
load(paste(path.to.github, "Rdata/blast_98pc_all3_piotii.Rda", sep = ""))
load(paste(path.to.github, "Rdata/blast_98pc_all3_swidsinskii.Rda", sep = ""))
load(paste(path.to.github, "Rdata/blast_98pc_all3_vaginalis.Rda", sep = ""))

###################### renaming species in taxonomy table ######################

# get subset of tax.table with only vNumbers corresponding to gardnerella
# i.e. all gardnerella genes in virgo
tax.gardnerella<-tax.table[grep("Gardnerella",tax.table$taxon),]

# for each species dataframe, pull gene lengths and percent query coverage for 
# each blast hit, add them to the corresponding species blast dataframes, then
# filter blast hits using conservative E value and query coverage thresholds
blast.dfs<- list()
for(i in ls(pattern = "blast.98.all3.")){
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
  tax.gardnerella$taxon[i] <- case_when(rownames(tax.gardnerella)[i] %in% blast.98.all3.leopoldii$sseqid ~ "Gardnerella_leopoldii",
                                        rownames(tax.gardnerella)[i] %in% blast.98.all3.piotii$sseqid ~ "Gardnerella_piotii",
                                        rownames(tax.gardnerella)[i] %in% blast.98.all3.swidsinskii$sseqid ~ "Gardnerella_swidsinskii",
                                        rownames(tax.gardnerella)[i] %in% blast.98.all3.vaginalis$sseqid ~ "Gardnerella_vaginalis",
                                        .default = "Gardnerella")
}

# check to see how many V numbers have been renamed- should be 198
length(grep("Gardnerella_", tax.gardnerella$taxon)) # 198: no duplicate issues

######################### data transformation for PCA #########################

# subset the combined dataset for only BV samples, then check for rows with
# zero sums
preg.np.bv <- preg.np.filt.bc[, which(preg.np.meta$health == "BV")]
any(rowSums(preg.np.bv) == 0) # FALSE

# set seed for RNG and perform zero replacement
set.seed(2023)
preg.np.bv.zr <- cmultRepl(t(preg.np.bv), label = 0,
                                 method = "CZM", z.warning = 0.99)

# CLR transformation and PCA
preg.np.bv.clr <- apply(t(preg.np.bv.zr), 2, function(x) log(x) - mean(log(x)))
pca.all3.bv <- prcomp(t(preg.np.bv.clr))

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
preg.np.tax.vec.gard <- preg.np.tax.vec[grep("Gardnerella", preg.np.tax.vec)]

# pull vNumbers and renamed (rn) taxonomic IDs from renamed gardnerella df
tax.gardnerella.rn <- tax.gardnerella[names(preg.np.tax.vec.gard),]

# replace old Gardnerella vaginalis labels in virginia taxa vector with new ones
for (i in 1:length(preg.np.tax.vec)) {
  if(names(preg.np.tax.vec)[i] %in% rownames(tax.gardnerella.rn)){
    vnum <- names(preg.np.tax.vec)[i]
    vnum.tax <- grep(vnum, rownames(tax.gardnerella.rn))
    preg.np.tax.vec[i] <- tax.gardnerella.rn$taxon[vnum.tax]
  }
}

# pull the indices of the V numbers corresponding to the Gardnerella genus and
# each named species in the virginia dataset (for colouring species)
ind.gard <- list()
for (i in levels(factor(tax.gardnerella.rn$taxon))) {
  ind.gard[[i]] <- which(preg.np.tax.vec == i)
}

# pull indices of all other vNumbers (those NOT found in the above index list)
ind.gard[["Other taxa"]] <- setdiff(c(1:nrow(preg.np.bv)),
                                    unlist(ind.gard))

# initiate graphics device, plot the BV-only biplot, emphasising only the 
# vNumbers corresponding to Gardnerella and save as .png
# png(paste(path.to.github,
#           "figs_for_paper/suppl_preg_np_biplot_gardnerella.png", sep = ""),
#     units = "in", height = 9, width = 12, res = 400)

codaSeq.PCAplot(pca.all3.bv, plot.groups = FALSE, plot.loadings = TRUE,
                plot.ellipses = NULL, plot.density = "loadings", grp.cex = 0.2,
                load.grp = ind.gard, load.cex = 0.5, load.col = tax.cols.gard,
                load.sym = 19, plot.legend = "loadings", leg.position = "topright",
                leg.cex = 1.15, leg.columns = 1, PC = c(1,3),
                title= "All datasets: BV - Gardnerella spp.")

# dev.off()
