# london/europe dataset: named gardnerella species biplots

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
#    - new.both.filt.bv (BV-only vNumber feature table, without entry 17852)
#    - new.tax.vec.bv (BV-only vNumber taxa vector, without entry 17852)
#    - both.bv.pca (PCA object of CLR-transformed, BV-only dataset)

load(paste(path.to.github, "Rdata/new.both.filt.bv.Rda", sep = ""))
load(paste(path.to.github, "Rdata/new.tax.vec.bv.Rda", sep = ""))
load(paste(path.to.github, "Rdata/both.bv.pca.Rda", sep = ""))

# load pre-computed taxa table (specifying column headers)
tax.table <- read.table(paste0(path.to.github,'1_VIRGO/1.taxon.tbl.txt'),
                        col.names = c("cluster", "vnumber", "taxon", "length"),
                        row.names = 2, header = F)

# load in BLAST results for genes unique to each named gardnerella spp. against
# database of V number genes present in london & europe samples (all 44)
load(paste(path.to.github, "Rdata/blast_98pc_leopoldii.Rda", sep = ""))
load(paste(path.to.github, "Rdata/blast_98pc_piotii.Rda", sep = ""))
load(paste(path.to.github, "Rdata/blast_98pc_swidsinskii.Rda", sep = ""))
load(paste(path.to.github, "Rdata/blast_98pc_vaginalis.Rda", sep = ""))

###################### renaming species in taxonomy table ######################

# get subset of tax.table with only vNumbers corresponding to gardnerella
# i.e. all gardnerella genes in virgo
tax.gardnerella<-tax.table[grep("Gardnerella",tax.table$taxon),]

# for each species dataframe, pull gene lengths and percent query coverage for 
# each blast hit, add them to the corresponding species blast dataframes, then
# filter blast hits using conservative E value and query coverage thresholds
blast.dfs<- list()
for(i in ls(pattern = "blast.98.")){
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
  tax.gardnerella$taxon[i] <- case_when(rownames(tax.gardnerella)[i] %in% blast.98.leopoldii$sseqid ~ "Gardnerella_leopoldii",
                                        rownames(tax.gardnerella)[i] %in% blast.98.piotii$sseqid ~ "Gardnerella_piotii",
                                        rownames(tax.gardnerella)[i] %in% blast.98.swidsinskii$sseqid ~ "Gardnerella_swidsinskii",
                                        rownames(tax.gardnerella)[i] %in% blast.98.vaginalis$sseqid ~ "Gardnerella_vaginalis",
                                        .default = "Gardnerella")
}

# check to see how many V numbers have been renamed- should be 302
length(grep("Gardnerella_", tax.gardnerella$taxon)) # 302: no duplicate issues

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ COMMENT @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# NOTE: this was an issue in prior versions of this script where no filtering
#       by E value / query coverage was done- several instances where multiple
#       blast hits to the same vNumber were causing duplicate issues. Kept this
#       code for future me, in case this happens with the virginia dataset!

# bind all blast rows together and check for duplicate V numbers
# blast.98.all<-rbind(blast.98.leopoldii,blast.98.piotii,blast.98.swidsinskii,blast.98.vaginalis)

# grep(TRUE,duplicated(blast.98.all$sseqid)) 

# dups<-blast.98.all %>%      # 10 x duplicate and 1x triplicate
#   group_by(sseqid)%>%
#   summarise(counts=n())%>%
#   arrange(desc(counts))

# arrange by number of vnumber (sseqid) occurrences and pull out this subset

# dups<-dups%>%
#   arrange(desc(counts))
# dups<-dups[1:11,]

# use the above vnumbers to pull out blast hits aligning to duplicated vnumbers
# and arrange by v numbers

# vnum.grep<-dups$sseqid
# blast.dups<-blast.98.all[grep(paste(vnum.grep,collapse = "|"),blast.98.all$sseqid),]%>%
#   arrange(sseqid)

# manual inspection of these shows no instances where genes from different species
# align to the same v number. 9 cases where a single gene from the same species 
# aligns to the same v number twice (small alignment lengths & different parts 
# of the gene) and 3 cases where two different genes from the same species align
# to the same v number.

# regardless, this doesn't complicate the assigning of new taxonomic IDs to
# V numbers!

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END COMMENT @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# set colours for new gardnerella biplot
tax.cols.gard <-c ("lightgrey", "red", "orange",
                   "lightpink2", "indianred4", rgb(0,0,0,0))

names(tax.cols.gard)<- c("Gardnerella",
                         "Gardnerella_leopoldii",
                         "Gardnerella_piotii",
                         "Gardnerella_swidsinskii",
                         "Gardnerella_vaginalis",
                         "Other")

# load in KEGG pathway & module tables (to assign to Gardnerella V numbers)
path.table <- read.table(paste(path.to.github,
                               "1_VIRGO/8.C.kegg.pathway.copy.txt",
                               sep = ""),
                         sep= "\t", header = TRUE, row.names = 1, fill = TRUE)

module.table <- read.table(paste(path.to.github,
                                 "1_VIRGO/8.B.kegg.module.txt",
                                 sep = ""),
                           sep= "\t", header = FALSE, row.names = 1, fill = TRUE)

# pull out unique V numbers for each gardnerella spp.
vnum.gard.funct <- tax.gardnerella[grep("Gardnerella_",tax.gardnerella$taxon),] %>% 
  arrange(taxon)

colnames(vnum.gard.funct)<- c("cluster", "taxonomic_id", "gene_length")

# add KO pathways & modules to unique V numbers for gardnerella spp.
vnum.gard.funct$pathway <- path.table[rownames(vnum.gard.funct),2]
vnum.gard.funct$module <- module.table[rownames(vnum.gard.funct),2]

# filter genes for only known functions
vnum.gard.funct <- vnum.gard.funct[!is.na(vnum.gard.funct$pathway) |
                                    !is.na(vnum.gard.funct$module),]

# pull gardnerella vNumbers present in london & europe bv datasets
new.tax.vec.gard <- new.tax.vec.bv[grep("Gardnerella", new.tax.vec.bv)]

# pull vNumbers and renamed (rn) taxonomic IDs from renamed gardnerella df
tax.gardnerella.rn <- tax.gardnerella[names(new.tax.vec.gard),]

# replace old Gardnerella vaginalis labels in lon/eur taxa vector with new ones
for (i in 1:length(new.tax.vec.bv)) {
  if(names(new.tax.vec.bv)[i] %in% rownames(tax.gardnerella.rn)){
    vnum <- names(new.tax.vec.bv)[i]
    vnum.tax <- grep(vnum, rownames(tax.gardnerella.rn))
    new.tax.vec.bv[i] <- tax.gardnerella.rn$taxon[vnum.tax]
  }
}

# rename and save renamed taxa vector
loneur.tax.vec.gvag <- new.tax.vec.bv
save(loneur.tax.vec.gvag, file = paste0(path.to.github, "Rdata/loneur.tax.vec.gvag.Rda"))

# pull the indices of the V numbers corresponding to the Gardnerella genus and
# each named species in the london/europe dataset (for colouring species)
ind.gard <- list()
for (i in levels(factor(tax.gardnerella.rn$taxon))) {
  ind.gard[[i]] <- which(new.tax.vec.bv == i)
}

# pull indices of all other vNumbers (those NOT found in the above index list)
ind.gard[["Other taxa"]] <- setdiff(c(1:nrow(new.both.filt.bv)),
                                  unlist(ind.gard))

# initiate graphics device, plot the BV-only biplot, emphasising only the 
# vNumbers corresponding to Gardnerella and save as .png
# png(filename = paste(path.to.github,
#                      "figs_for_paper/suppl_lon_eur_biplot_gardnerella.png", sep = ""),
#     units = "in", height = 9, width = 12, res = 400)

codaSeq.PCAplot(both.bv.pca, plot.groups = FALSE, plot.loadings = TRUE,
                plot.ellipses = NULL, plot.density = "loadings",
                load.grp = ind.gard, load.cex = 0.5, load.col = tax.cols.gard,
                load.sym = 19, plot.legend = "loadings", leg.position = "topright",
                leg.cex = 0.9, leg.columns = 1, PC = c(1,2),
                title= "London & Europe: BV - Gardnerella spp.")

# dev.off()
