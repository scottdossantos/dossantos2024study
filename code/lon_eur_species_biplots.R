# london/europe dataset: species-level summary biplots

#################################### setup ####################################

# set path to the github repository
user <- "sds/cc"
path.to.github <- switch(user,
                         "sds/cc" = "~/Documents/GitHub/metatranscriptome/",
                         "gg" = "~/Documents/0_git/projects/metatranscriptome/")

# run setup code to get batch-corrected feature tables and taxa vector
source(paste(path.to.github, "code/setup.R", sep=""))


# remove two janky London BV samples from the filtered dataset:
#   - v.001A - close to 100% L. gasseri with practically no BV organisms
#   - v.019A - around 80 % iners, ~ 20% crispatus and a tiny bit of Gardnerella
new.both.filt <- new.both.filt[,-c(9,22)]

# NOTE: these samples looked very out of place in previous heatmaps and in
# differential abundance analyses

################################ PCA & biplots ################################

# set seed for RNG 
set.seed(2023)

# Bayesian multiplicative replacement of zeroes (transposed as samples expected 
# in rows)
clr.input <- cmultRepl(t(new.both.filt), method = "CZM", 
                       label = "0", z.warning = 0.99)

# CLR transformation of batch-corrected, zero-replaced feature table
# (expects samples as columns, so have to transpose again)
both.clr <- apply(t(clr.input) , 2, function(x) log(x) - mean(log(x)))

# run PCA on CLR-transformed data (expects samples as rows, so transpose again)
both.pca <- prcomp(t(both.clr))

# get the total number of genes in each species
n.genes <- vector()
for(i in levels(factor(new.tax.vec))){
  n.genes[i] <- length(which(new.tax.vec == i))
}

# get taxa represented by >75 genes
tax.order.75 <- names(which(n.genes > 75))
tax.colours <- c(tax.colors[tax.order.75,], rgb(0,0,0,0.05))

# make a list of row indices corresponding to the rows of the feature table
loadings.group <- list()
for(i in 1:length(tax.order.75)){
    loadings.group[[tax.order.75[i]]] <- which(new.tax.vec == tax.order.75[i])
}

# Get indicies of all taxa represented by <75 genes
loadings.group[["Other"]] <- setdiff(c(1:nrow(both.clr)),
                                     unlist(loadings.group))

# draw biplot with legend and density plots for loadings, then save as .png

png(paste(path.to.github,
          "figs_for_paper/lon_eur_summarySpecies_biplot.png", sep = ""),
    units = "in", height=9, width=12, res = 400)

codaSeq.PCAplot(both.pca, plot.groups=FALSE, plot.loadings=TRUE,
                plot.ellipses = NULL, plot.density = "loadings",  
                load.grp=loadings.group, load.col=tax.colours, PC = c(1,2),
                plot.legend = "loadings", leg.position = "bottomleft",
                leg.cex = 0.65, leg.columns = 3, 
                title = "London & Europe - Healthy vs. BV")

dev.off()

############################## top species biplot ##############################

# get the total number of reads (batch-corrected) of each species in each sample
n.reads <- list()
for(i in levels(factor(new.tax.vec))){
  n.reads[[i]] <- colSums(new.both.filt[which(new.tax.vec == i),])
}

# collapse list into a data frame by calling rbind on all list elements
n.reads.df <- as.data.frame(do.call(rbind, n.reads))
rownames(n.reads.df) <- names(n.reads)

# summarise reads per sample / spp.
reads.sample <- colSums(n.reads.df) # total mapped reads per sample
reads.species <- rowSums(n.reads.df) # total mapped reads per species
length(reads.sample) <- length(reads.species)

reads.summary <- cbind(name.samples = names(reads.sample),
                       reads.samples = reads.sample, 
                       name.species = names(reads.species),
                       reads.species = reads.species)

rownames(reads.summary)<-NULL

# generate matrix for identifying good cutoff value for defining 'main' species
top.sp.mat <- matrix(data=NA, ncol=ncol(n.reads.df), nrow=7)
colnames(top.sp.mat)<-colnames(n.reads.df)
rownames(top.sp.mat)<-c("total",paste0(rep("10^",5),c(3:7)), ">75_genes")

# calculate total read counts for all species
top.sp.mat[1,] <- colSums(n.reads.df)

# get total read counts using various read cutoff values: 10^3 - 10^7
top.sp.mat[2,] <- colSums(n.reads.df[which(rowSums(n.reads.df) > 1e3),])
top.sp.mat[3,] <- colSums(n.reads.df[which(rowSums(n.reads.df) > 1e4),])
top.sp.mat[4,] <- colSums(n.reads.df[which(rowSums(n.reads.df) > 1e5),])
top.sp.mat[5,] <- colSums(n.reads.df[which(rowSums(n.reads.df) > 1e6),])
top.sp.mat[6,] <- colSums(n.reads.df[which(rowSums(n.reads.df) > 1e7),])
top.sp.mat[7,] <- colSums(n.reads.df[tax.order.75,])

# calculate proportion of reads retained when using various cutoff values
top.sp.cut<-t(t(top.sp.mat)/top.sp.mat[1,])

# using >75 genes retains >94% of reads across all samples
summary(top.sp.cut[7,])

# filter dataframe to exclude taxa represented by fewer than 75 genes
top.species <- n.reads.df[tax.order.75,]

# CLR transformation (running cMultRepl() before doesn't make much difference to
# the overall plot interpretation at all!)
top.clr <- apply(top.species+0.5, 2, function(x) log(x) - mean(log(x)))

# run PCA: this is a species representation, which should go in supplement
top.pca <- prcomp(t(top.clr))

# make list of health vs. BV status
top.group.list <- list(Healthy=grep("h.", colnames(top.species)),
                     BV=grep("v.", colnames(top.species)))

# make list of top species
top.load.list<-list()
for(i in levels(factor(rownames(top.species)))){
  top.load.list[[i]] <- which(rownames(top.species) == i)
}

# draw biplot with density plots for groups and save as .png
# png(paste(path.to.github, 
#           "figs_for_paper/suppl_lon_eur_biplot_topSpecies.png", sep = ""),
#     units = "in", height=9, width=12, res = 400)

codaSeq.PCAplot(top.pca, plot.groups = TRUE, plot.loadings = TRUE, 
                plot.ellipses = NULL, plot.density = "groups",
                grp = top.group.list, grp.col = c("skyblue2","goldenrod2"),
                grp.sym = "text", grp.cex = 1.1, load.grp = top.load.list,
                load.col = rep("blue",16), load.sym = "text", load.cex = 0.7,
                PC = c(1,2), plot.legend = "groups", leg.cex = 0.75,
                leg.position = "bottomleft",
                title = "London & Europe: Healthy vs. BV - Top Species")

# dev.off()

############################ BV samples only biplot ############################

# pull out only BV samples from batch-corrected dataset
new.both.filt.bv <- new.both.filt[, grep("v.",colnames(new.both.filt))]

any(rowSums(new.both.filt.bv)==0) # shows TRUE
which(rowSums(new.both.filt.bv)==0) # culprit is row 17852 (V1614838)

# removing two samples prior to batch correction means row 17852 (V1614838) now
# has 0 counts for all BV samples, so we have to remove it from the data frame
# and also `tax.vec` (will also have to remove this for the Gardnerella biplots)

# remove row from batch-corrected BV dataframe and taxa vector, and save as .Rda
new.both.filt.bv <- new.both.filt.bv[-17852,]

# save(new.both.filt.bv, file = paste(path.to.github,
#                              "Rdata/new.both.filt.bv.Rda", sep = ""))

new.tax.vec[17852] # confirm V number = V1614838 (taxonomic ID is unknown)
new.tax.vec.bv <- new.tax.vec[-17852] 

# save(new.tax.vec.bv, file = paste(path.to.github,
#                                   "Rdata/new.tax.vec.bv.Rda", sep = ""))

any(rowSums(new.both.filt.bv)==0) # now shows FALSE


# set seed for RNG and perform zero replacement
set.seed(2023)
bv.clr.input <- cmultRepl(t(new.both.filt.bv), method = "CZM", 
                          label = "0", z.warning = 0.99)

# CLR transformation
bv.clr <- apply(t(bv.clr.input), 2, function(x) log(x) - mean(log(x)))

# run PCA on CLR-transformed, zero-replaced, BV-only data

bv.pca <- prcomp(t(bv.clr))

# save pca object to Rdata/
both.bv.pca<-bv.pca

# save(both.bv.pca,file = paste(path.to.github,
#                               "Rdata/both.bv.pca.Rda", sep = ""))

# get the total number of genes in each species
n.genes.bv <- vector()
for(i in levels(factor(new.tax.vec.bv))){
  n.genes.bv[i] <- length(which(new.tax.vec.bv == i))
}

# get list of taxa represented by >75 genes
bv.taxa.75<- names(which(n.genes.bv >75))

# get taxa colours for above taxa
bv.taxa.cols<-c(tax.colors[bv.taxa.75,], rgb(0,0,0,0.05))

# get loading indices
bv.loadings<-list()
for (j in levels(factor(bv.taxa.75))) {
  bv.loadings[[j]]<-which(new.tax.vec.bv==j)
}

bv.loadings[["Other"]]<-setdiff(c(1:nrow(new.both.filt.bv)),
                                 unlist(bv.loadings))

# draw PCA plot of BV only data with density plots and save as .png
# png(filename = paste(path.to.github,
#                      "figs_for_paper/suppl_lon_eur_biplot_bvOnly.png",sep = ""),
#     units = "in", height=9, width=12,res = 400)

codaSeq.PCAplot(bv.pca, plot.groups = FALSE, plot.loadings = TRUE, 
                plot.ellipses = NULL, plot.density = "loadings",
                load.grp = bv.loadings, load.col = bv.taxa.cols, load.sym = 19,
                load.cex = 0.5, PC = c(1,2), plot.legend = "loadings",
                leg.position = "bottomright", leg.columns = 3, leg.cex = 0.7,
                title = "London & Europe: BV Only")

# dev.off()
