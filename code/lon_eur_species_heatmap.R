# london/europe dataset: species-level summary heatmap

#################################### setup ####################################

library(dplyr) # for data manipulation
library(pheatmap) # for heatmap
library(viridisLite) # for colourblind-friendly palettes
library(vegan) # for rearranging dendrogram

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

# get total & median number of reads in each species for all / healthy / bv
n.reads<-rowSums(no.reads.df)
n.reads.h<-rowSums(no.reads.df[,grep("h.", colnames(no.reads.df))])
n.reads.bv<-rowSums(no.reads.df[,grep("v.", colnames(no.reads.df))])

med.reads<-apply(no.reads.df, 1, median)
med.reads.h<-apply(no.reads.df[,grep("h.", colnames(no.reads.df))], 1, median)
med.reads.bv<-apply(no.reads.df[,grep("v.", colnames(no.reads.df))], 1, median)

# get the total number of genes in each species
n.genes <- vector()
for(i in levels(factor(new.tax.vec))){
  n.genes[i] <- length(which(new.tax.vec == i))
}

# get total number of genes in species with no known taxonomy in VIRGO
n.genes["No_tax_info"]<-length(which(is.na(new.tax.vec)))

# make summary data frame
no.reads.summ<-as.data.frame(cbind(No.Genes = n.genes,
                                   No.Reads.All = n.reads,
                                   No.Reads.Healthy = n.reads.h,
                                   No.Reads.BV = n.reads.bv,
                                   Med.Reads.All = as.integer(med.reads),
                                   Med.Reads.Healthy = as.integer(med.reads.h),
                                   Med.Reads.BV = as.integer(med.reads.bv)))

# get taxa represented by >75 genes (ignoring genes w/ no taxonomy in VIRGO)
tax.order.75 <- names(which(n.genes[-65] > 75))

# calculate total number of reads for species w/ >75 genes, across all samples
# output is a list of named integer vectors; read sum x samples, one per species
no.reads.total.75 <- list()
for(i in levels(factor(tax.order.75))){
  no.reads.total.75[[i]] <- colSums(new.both.filt[which(new.tax.vec == i),])
}

# collapse to data frame
no.reads.total.75.df <- as.data.frame(do.call(rbind, no.reads.total.75))

# generate matrix for checking what % of data captured by >75 genes
top.sp.mat <- matrix(data=NA, ncol=ncol(no.reads.df), nrow=7)
colnames(top.sp.mat)<-colnames(no.reads.total.75.df)
rownames(top.sp.mat)<-c("total",paste0(rep("10^",5),c(3:7)), ">75 genes")

# calculate total read counts for all species (excl. spp with unknown taxonomy)
top.sp.mat[1,] <- colSums(no.reads.df[-65,])

# get total read counts using various read cutoff values: 10^3 - 10^7
top.sp.mat[2,] <- colSums(no.reads.df[which(rowSums(no.reads.df[-65,]) > 1e3),])
top.sp.mat[3,] <- colSums(no.reads.df[which(rowSums(no.reads.df[-65,]) > 1e4),])
top.sp.mat[4,] <- colSums(no.reads.df[which(rowSums(no.reads.df[-65,]) > 1e5),])
top.sp.mat[5,] <- colSums(no.reads.df[which(rowSums(no.reads.df[-65,]) > 1e6),])
top.sp.mat[6,] <- colSums(no.reads.df[which(rowSums(no.reads.df[-65,]) > 1e7),])
top.sp.mat[7,] <- colSums(no.reads.total.75.df)

# calculate proportion of reads retained when using various cutoff values
top.sp.cut<-t(t(top.sp.mat)/top.sp.mat[1,])

# NOTE: have to transpose input matrix so samples are in rows as R will divide
# the input matrix values by the supplied vector values in order BY COLUMN. Also
# have to then transpose the output to get samples in columns again

# using >75 genes retains >94% of reads across all samples
summary(top.sp.cut[7,])

############################### heatmap plotting ###############################

# convert to relative abundance
hm.input <- as.matrix(t(t(no.reads.total.75.df)/colSums(no.reads.total.75.df)))

# create metadata object for london/europe samples
hm.metadata <- data.frame(Dataset = rep(c("London","Europe"), c(20,22)),
                          `BV status` = rep(c("Healthy","BV","Healthy"), c(8,26,8)),
                          row.names = colnames(hm.input), check.names = FALSE)

# generate dendrogram from distance matrix calculation on batch-corrected props
hm.dendrogram<-hclust(dist(t(hm.input), method = "euclidean"),
                      method = "ward.D2")

plot(hm.dendrogram)

# get order of samples in dendrogram as a named integer vector
dend.labs<-c(1:length(rownames(hm.metadata)))
for(i in 1:length(rownames(hm.metadata))){
  ind<-hm.dendrogram$order[i]
  names(dend.labs)[i]<-hm.dendrogram$labels[ind]
}

# edit order as desired and sort vector according to new order
# made it so all CSTs are unseparated and so health/BV split is easier to see
# NOTE: the only reason I know which indices are which is because I made the
#       heatmap once already with a non-reordered dendrogram

dend.labs[13:18] <- c(8:13)
dend.labs[8:12] <- c(14:18)
dend.labs[14] <- 13
dend.labs[15:18] <- c(9:12)
dend.labs[31:42] <- c(24:35)
dend.labs[24:30] <- c(36:42)
dend.labs<-sort(dend.labs)

# make vector of new orders that corresponds to order of samples in metadata
dend.labs.new<-vector(length = length(rownames(hm.metadata)))
for(i in 1:length(rownames(hm.metadata))){
  ind.in.meta<-which(names(dend.labs)[i]==rownames(hm.metadata))
  dend.labs.new[ind.in.meta]<-i
}

# add new order to metadata as new column
hm.metadata$Order <- dend.labs.new

# re-ordering dendrogram so h.b13_22 next to other dark blue samples
hm.dendrogram.new <- vegan:::reorder.hclust(hm.dendrogram, hm.metadata$Order)

# plot new dendrogram
plot(hm.dendrogram.new)

# add cluster membership to metadata (assigned by cutting tree into six groups)
hm.metadata$CST<-cutree(hm.dendrogram, k=6)

# this gives the following 'CSTs':
#         1 = L. crispatus dominant (CST I)
#         2 = L. iners dominant (CST III)
#         3 = Various anaerobes (CST IV)
#         4 = L. jensenii dominant (CST V)
#         5 = Mix of L. iners and BV organisms (CST III / CST IV)
#         6 = Gardnerella dominant (CST IV)

# change grouping for h.b13_22 & h.b6_10 (L. iners dominant- should be CST III)
hm.metadata$CST[42] <- 2
hm.metadata$CST[38] <- 2

# replace cutree integers with heatmap labels of canonical CSTs
#    1 =  crispatus          2 = iners
#    3 =  BV (various)       4 = jensenii
#    5 =  iners / BV         6 = BV (gardnerella)

hm.metadata$CST <- case_when(hm.metadata$CST == 1 ~ "I",
                             hm.metadata$CST == 2 ~ "III",
                             hm.metadata$CST == 3 ~ "IV",
                             hm.metadata$CST == 4 ~ "V",
                             hm.metadata$CST == 5 ~ "III / IV",
                             hm.metadata$CST == 6 ~ "IV")

# save metadata and CST vector as a .Rda object for use in making a heatmap of
# differential functions between health and BV in London/Europe datasets
# save(hm.metadata, 
#      file=paste(path.to.github, "Rdata/hm.metadata.Rda", sep = ""))

# make list for column colour bars
hm.column.cols<-list(Dataset=c(Europe="grey10", London="snow"),
                     `BV status`=c(BV="goldenrod2", Healthy="skyblue1"),
                     CST=c(I=viridis(5)[1],
                           III=viridis(5)[2],
                           `III / IV`=viridis(4)[3],
                           IV=viridis(5)[4],
                           V=viridis(5)[5]))

# change colour for CST IV in accordance with reviewer comment
hm.column.cols[["CST"]][4] <- "#9AC869FF"

# save(hm.column.cols,
#      file = paste(path.to.github, "Rdata/hm.column.cols.Rda", sep = ""))

# make list of row names to be italicised
hm.row.italics<-lapply(gsub("_", " ", rownames(hm.input)),
                          function(x){
                            bquote(italic(.(x)))
                          })

# edit taxa to undo italics where needed
hm.row.italics[2]<-"BVAB1"
hm.row.italics[8]<-expression(italic("Megasphaera")~"genomosp.")
hm.row.italics[9]<-expression(italic("Megasphaera")~"sp.")
hm.row.italics[13]<-expression(italic("Prevotella")~"sp.")

# make vector of colourblind-friendly heatmap colours
hm.colours<-rev(magma(100))

# draw heatmap and save as .png (cheesed the legend break as 1.0 wouldn't show)
# png(paste(path.to.github,
#           "figs_for_paper/lon_eur_summarySpecies_heatmap.png", sep = ""),
#     units = "in", height = 6, width = 10, res = 400)

pheatmap(hm.input, color = hm.colours, border_color = rgb(0,0,0,0.1),
         cluster_rows = FALSE, cluster_cols = hm.dendrogram.new,
         show_colnames = FALSE, labels_row = as.expression(hm.row.italics),
         annotation_col = hm.metadata[,-3], annotation_colors = hm.column.cols,
         legend_breaks = c(0, 0.2, 0.4, 0.6, 0.8, 0.9955),
         legend_labels = c("0","0.2","0.4","0.6","0.8","1.0"))

# dev.off()
