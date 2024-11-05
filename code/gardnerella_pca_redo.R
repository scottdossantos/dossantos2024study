# redo of gardnerella plots for reviewer

#################################### setup ####################################

# set github directory
github <- "~/Documents/GitHub/dossantos2024study/"

# load required libraries
library(CoDaSeq) # for plotting PCAs

# need to load in several .Rda objects containing updated version of Lon/Eur
# and Virginia datasets mapped to updated version of VIRGO with all Gardnerella
# spp.
#   - PCAs for updated BV only ordinations (all datasets)
#   - Gene loading lists for updated BV only ordinations (all datasets)
#   - Sample loading list for updated BV only ordinations (all datasets)
#   - Colour lookup table
load(paste0(github,"Rdata/gardRedo.pca.le.bv.Rda"))
load(paste0(github,"Rdata/gardRedo.pca.virg.bv.Rda"))
load(paste0(github,"Rdata/gardRedo.pca.pnp.bv.Rda"))
load(paste0(github,"Rdata/gardRedo.load.le.bv.Rda"))
load(paste0(github,"Rdata/gardRedo.load.virg.bv.Rda"))
load(paste0(github,"Rdata/gardRedo.load.pnp.bv.Rda"))
load(paste0(github,"Rdata/gardRedo.grp.le.bv.Rda"))
load(paste0(github,"Rdata/gardRedo.grp.virg.bv.Rda"))
load(paste0(github,"Rdata/gardRedo.grp.pnp.bv.Rda"))
lookup.col <- read.table(paste0(github,"Rdata/species_colors_reviewers.txt"),
                         sep = "\t", header = F, quote = "", row.names = 1)

################################### lon/eur ###################################

# make colour vector for gardnerella (virginia dataset)
cols.gard.le <- c("grey80","red","orange3",
                  "palevioletred1","indianred4",
                  rgb(0,0,0,0))

names(cols.gard.le) <- c("Gardnerella","Gardnerella_leopoldii",
                         "Gardnerella_pickettii","Gardnerella_swidsinskii",
                         "Gardnerella_vaginalis","Other taxa")

# subset loadings list for Gardnerella (virginia dataset), re-order the species
# alphabetically, then pull loading indices that aren't already in there for
# 'Other taxa'
load.gard.le <- load.le.bv[grep("Gardnerella", names(load.le.bv))]
load.gard.le <- load.gard.le[sort(names(load.gard.le))]
load.gard.le[["Other taxa"]] <- setdiff(1:max(unlist(load.gard.le)),
                                        unlist(load.gard.le))

# sanity check for duplicates and names -> colours
any(duplicated(unlist(load.gard.le))) # returns FALSE
names(cols.gard.le) == names(load.gard.le) # all TRUE

# png(paste0(github, "figs_for_paper/R_FigsLondonEurope/suppl_redo_lon_eur_biplot_gardnerella.png"),
#     units = "in", height = 6, width = 10, res = 400)

codaSeq.PCAplot(pca.le.bv, plot.groups = T, plot.loadings = T,
                plot.density = "loadings", plot.legend = "loadings",
                grp = grp.le.bv, grp.cex = 0.7, grp.col = "navy",
                PC = c(1,2), load.grp = load.gard.le, load.sym = 19,
                load.col = cols.gard.le, load.cex = 0.4, leg.columns = 2,
                leg.position = "topleft", leg.cex = 0.7,
                title = "London/Europe dataset: BV only - Gardnerella")

# dev.off()

################################### virginia ###################################

# make colour vector for gardnerella (virginia dataset)
cols.gard.virg <- c("grey80","burlywood4","red","orange3",
                    "orange","palevioletred1","indianred4",
                    rgb(0,0,0,0))

names(cols.gard.virg) <- c(grep("Gardnerella", rownames(lookup.col), value = T),
                           "Other taxa")

# subset loadings list for Gardnerella (virginia dataset), re-order the species
# alphabetically, then pull loading indices that aren't already in there for
# 'Other taxa'
load.gard.virg <- load.virg.bv[grep("Gardnerella", names(load.virg.bv))]
load.gard.virg <- load.gard.virg[sort(names(load.gard.virg))]
load.gard.virg[["Other taxa"]] <- setdiff(1:max(unlist(load.gard.virg)),
                                          unlist(load.gard.virg))

# sanity check for duplicates and names -> colours
any(duplicated(unlist(load.gard.virg))) # returns FALSE
names(cols.gard.virg) == names(load.gard.virg) # all TRUE

# plot BV only PCA for virginia dataset with Gardnerella spp. highlighted
# png(paste0(github, "figs_for_paper/R_FigsVirginia/suppl_redo_virginia_biplot_gardnerella.png"),
#            units = "in", height = 6, width = 10, res = 400)

codaSeq.PCAplot(pca.virg.bv, plot.groups = T, plot.loadings = T,
                plot.density = "loadings", plot.legend = "loadings",
                grp = grp.virg.bv, grp.cex = 0.25, grp.col = "navy",
                PC = c(1,2), load.grp = load.gard.virg, load.sym = 19,
                load.col = cols.gard.virg, load.cex = 0.4, leg.columns = 2,
                leg.position = "topleft", leg.cex = 0.7,
                title = "Virginia dataset: BV only - Gardnerella")

# dev.off()

############################### preg vs. nonpreg ###############################

# make colour vector for gardnerella (combined dataset)
cols.gard.pnp <- c("grey80","burlywood4","red","orange3",
                   "orange","palevioletred1","indianred4",
                   rgb(0,0,0,0))

names(cols.gard.pnp) <- c(grep("Gardnerella", rownames(lookup.col), value = T),
                           "Other taxa")

# subset loadings list for Gardnerella (combined dataset), re-order the species
# alphabetically, then pull loading indices that aren't already in there for
# 'Other taxa'
load.gard.pnp <- load.pnp.bv[grep("Gardnerella", names(load.pnp.bv))]
load.gard.pnp <- load.gard.pnp[sort(names(load.gard.pnp))]
load.gard.pnp[["Other taxa"]] <- setdiff(1:max(unlist(load.gard.pnp)),
                                          unlist(load.gard.pnp))

# sanity check for duplicates and names -> colours
any(duplicated(unlist(load.gard.pnp))) # returns FALSE
names(cols.gard.pnp) == names(load.gard.pnp) # all TRUE

# plot BV only PCA for combined dataset with Gardnerella spp. highlighted
# png(paste0(github, "figs_for_paper/R_FigsPregNonPreg/suppl_redo_preg_np_biplot_gardnerella.png"),
#            units = "in", height = 6, width = 10, res = 400)

codaSeq.PCAplot(pca.pnp.bv, plot.groups = T, plot.loadings = T,
                plot.density = "loadings", plot.legend = "loadings",
                grp = grp.pnp.bv, grp.cex = 0.25, grp.col = "navy",
                PC = c(1,2), load.grp = load.gard.pnp, load.sym = 19,
                load.col = cols.gard.pnp, load.cex = 0.4, leg.columns = 2,
                leg.position = "bottomright", leg.cex = 0.7)

# dev.off()
