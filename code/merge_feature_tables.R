# create merged feature table from T2T-mapped VIRGO outputs

############################ Input data description ############################

# Europe reads downloaded from ENA
# London reads stored on Gloor lab server (human/rRNA-depleted files only)
# Virginia (MOMS-PI) reads downloaded from dbGaP

# note: the London reads in this location have been mapped to HG38 and SILVA,
# but NOT T2T (complete human genome produced by long read sequencing). Final
# feature table of London reads HAS been mapped to T2T genome assembly.

# Input for VIRGO was london/europe/virginia reads that were trimmed for quality 
# (minimum Q20 over sliding window of 4 bp) and cropped to 75 bp, and that did
# not align to HG38 (human), T2T (complete human) or SILVA (rRNA)

# VIRGO step 1 produced '.out' files containing all V numbers present in a given
# sample, the frequency of reads aligning to that V number, and the length of
# the gene represented by a given V number

# VIRGO step 2 was run in three batches: one run corresponding to each dataset,
# and produced summary tables- most important being "summary.NR.abundance.txt".
# This is a feature table of all non-redundant V numbers and their frequency 
# across the supplied samples.

############################# setup and load data #############################

library(CoDaSeq)

# set user and path to github
user <- "sds/cc"
path.to.github <- switch(user,
                         "sds/cc" = "~/Documents/GitHub/dossantos2024study/",
                         "gg" = "~/Documents/0_git/projects/dossantos2024study/")

# load feature tables from VIRGO step 2 (have to use unz() on the virginia table
# as it is a .zip archive)
initial.london <- read.table(paste(path.to.github,
                                   "Rdata/t2t_summary.london.NR.abundance.txt",
                                   sep = ""),
                             header = TRUE, row.names = 1, check.names = FALSE)

initial.europe <- read.table(paste(path.to.github,
                                   "Rdata/t2t_summary.europe.NR.abundance.txt",
                                   sep = ""),
                             header = TRUE, row.names = 1, check.names = FALSE)

initial.virginia <- read.table(unz(description = paste(path.to.github,
                                                       "Rdata/t2t_summary.virginia.NR.abundance.txt.zip",
                                                       sep = ""),
                                   filename = "t2t_summary.virginia.NR.abundance.txt"),
                               header = T, row.names = 1, quote = '', sep = "\t")

# load taxa table from VIRGO
tax.table <- read.table(paste(path.to.github,'1_VIRGO/1.taxon.tbl.txt', sep=""),
                        header=F, row.names=2)

# load in KEGG pathway table from VIRGO
path.table <- read.table(paste(path.to.github,
                               '1_VIRGO/8.C.kegg.pathway.copy.txt',
                               sep=""), 
                         sep="\t", header=T, row.names=1, fill=TRUE)

# load in health vs. bv groupings for virginia dataset
load(paste(path.to.github, "Rdata/virginia.groups.Rda", sep = ""))

############################# London feature table #############################

# pull out healthy / BV sample names
# note:008B and 016B are 'intermediate BV' based on Nugent scoring and are
# excluded from downstream analysis (see 'Rdata/m.24.Rda')
h.london <- c("001B", "002B", "004B", "006B", "009B", "010B", "015B", "020B")

bv.london <- c("001A", "003A", "006A", "008A", "009A", "010A", "012B",
               "013B", "014B", "018B", "017B", "012A", "013A", "019A")

# subset london VIRGO output by above sample names
renamed.london <- data.frame(initial.london[,h.london], 
                             initial.london[,bv.london],
                             check.names = FALSE)

# check for and remove any rows with 0 read sums
any(rowSums(renamed.london) ==0)
renamed.london <- renamed.london[-c(which(rowSums(renamed.london)==0)),]

# edit column names to include health / BV status
colnames(renamed.london) <- c(paste('h.',h.london, sep=''),
                              paste('v.', bv.london, sep=''))

############################# Europe feature table #############################

# pull current names of healthy & BV samples
# note 1: bv europe samples are BV prior to metronidazole treatment and
#         healthy europe samples are 'normal' after metronidazole treatment
# note 2: samples excluded from downstream analysis are those whose 
#         reaction to metronidazole was recorded as 'na'

bv.europe <- c("ERR2014354", "ERR2014359", "ERR2014361", "ERR2014363", 
            "ERR2014365", "ERR2014370", "ERR2014372", "ERR2014374", 
            "ERR2014376", "ERR2014378", "ERR2014383", "ERR2014388", 
            "ERR2014390", "ERR2014392")

bv.europe <- paste(bv.europe,"_1",sep = "")

h.europe <- c("ERR2014355", "ERR2014362", "ERR2014366", "ERR2014373", 
              "ERR2014377", "ERR2014379", "ERR2014384", "ERR2014389")

h.europe <- paste(h.europe,"_1",sep = "")

# make new names of healthy and bv samples (corresponds with '/code/setup.R')
bv.europe.names <- paste('v.',c("a4_1", "a5_1", "a5_12", "a6_1", "a6_4", "a6_6",
                                "a6_10", "a7_1", "a8_3", "a8_6", "a13_19",
                                "a13_22", "a15_6", "a17_4"), sep='')

h.europe.names <- paste('h.',c("b4_1", "b5_12", "b6_4", "b6_10", "b8_3", "b8_6",
                               "b13_19", "b13_22"), sep='')

# subset europe VIRGO output by above sample names
renamed.europe <- data.frame(initial.europe[,bv.europe],
                             initial.europe[,h.europe])

# edit column names to include health / BV status
colnames(renamed.europe) <- c(bv.europe.names, h.europe.names)

# check for and remove any rows with 0 read sums
any(rowSums(renamed.europe) ==0)
renamed.europe <- renamed.europe[-c(which(rowSums(renamed.europe)==0)),]

############################# Merge feature tables #############################

# get vector of non-redundant v numbers in europe & london datasets
nr.vnumbers.all <- (unique(c(rownames(renamed.europe),
                             rownames(renamed.london))))

# make dataframe of zeros, with dimensions equal to combined number of samples
# from both datasets by number of non-redundant v numbers found in both datasets
# 44 cols x 330,123 rows
new.both.data <-data.frame(matrix(data=0,
                                  ncol=ncol(renamed.london) + ncol(renamed.europe),
                                  nrow=length(nr.vnumbers.all)))

# set column and row names as per both datasets
colnames(new.both.data) <- c(colnames(renamed.london),
                             colnames(renamed.europe))

rownames(new.both.data) <- nr.vnumbers.all

# combining of dataframes is done in two steps to ensure that v numbers found in
# both datasets are merged into the same row

# fill in counts for all rows found in the london dataset across the first 22
# columns (corresponds to the 22 london samples)
new.both.data[rownames(renamed.london), 1:22] <- renamed.london

# now fill in counts for all rows found in europe dataset across the next 22 
# columns (corresponds to 22 europe samples)
new.both.data[rownames(renamed.europe), 23:44] <- renamed.europe

# batch-correction of dataset using negative binomial regression with
# ComBat_seq (from sva package), including biological condition of interest
# (health vs. BV) in the model
conds <- c(rep('h',8), rep('b',14), rep('b',14), rep('h',8)) # health vs. BV
batch <- c(rep(1,22), rep(2,22)) # london batch vs. europe batch

new.both.data <- as.data.frame(sva::ComBat_seq(as.matrix(new.both.data),
                                               batch=batch,
                                               group=conds,
                                               full_mod=T))

# save dataframe as .Rda object
# save(new.both.data, file=paste(path.to.github,
#                                'Rdata/new.both.data.Rda',
#                                sep = ""))

# filter dataset with codaSeq.filter(), using same parameters as '/code/setup.R'
# minimum prevalence: gene found in >=30 % of samples
# minimum prop: gene has relative abundance of >=0.005 % in 1 or more samples
# returns 19,964 genes across 44 samples
new.both.filt <- codaSeq.filter(new.both.data, samples.by.row=F,
                                min.occurrence=0.30, min.prop=0.00005)

# save filtered dataframe as .Rda object
# save(new.both.filt, file = paste(path.to.github,
#                                  "Rdata/new.both.filt.Rda",
#                                  sep = ""))

# make taxa vector from filtered dataset and save to .Rda object
new.tax.vec <- tax.table[rownames(new.both.filt),2]
names(new.tax.vec) <- rownames(new.both.filt)

# save(new.tax.vec, file = paste(path.to.github,
#                                "Rdata/new.tax.vec.Rda",
#                                sep = ""))

#################### Merge feature tables: preg & non-preg ####################

# remove "_1" from all sample headers and rename this for consistency with the
# code for london & europe datasets
colnames(initial.virginia) <- gsub("_1", "", colnames(initial.virginia))
renamed.virginia <- initial.virginia

# get vector of non-redundant v numbers in europe, london and virginia datasets 
nr.vnumbers.all.v <- unique(c(rownames(renamed.europe),
                              rownames(renamed.london),
                              rownames(initial.virginia)))

# make dataframe of zeros, with dimensions equal to combined number of samples
# from both datasets by number of non-redundant v numbers found in both datasets
# 341 cols x 537,322 rows
preg.np.data <-data.frame(matrix(data=0,
                                 ncol=ncol(renamed.london) + ncol(renamed.europe) + ncol(renamed.virginia),
                                 nrow=length(nr.vnumbers.all.v)))

# set column and row names as per both datasets
colnames(preg.np.data) <- c(colnames(renamed.london),
                            colnames(renamed.europe),
                            colnames(renamed.virginia))

rownames(preg.np.data) <- nr.vnumbers.all.v

# combining of dataframes is done in three steps to ensure that v numbers shared 
# among datasets are merged into the correct row

# fill in counts for all rows found in the london dataset across the first 22
# columns (corresponds to the 22 london samples)
preg.np.data[rownames(renamed.london), 1:22] <- renamed.london

# now fill in counts for all rows found in europe dataset across the next 22 
# columns (corresponds to 22 europe samples)
preg.np.data[rownames(renamed.europe), 23:44] <- renamed.europe

# now fill in counts for all rows found in london dataset across the remaining 
# columns (corresponds to 297 virginia samples)
preg.np.data[rownames(renamed.virginia), 45:341] <- renamed.virginia

# batch correction is probably not feasible for this dataset: the discrepancy
# between colSums() for the london & european samples and colSums() for the 
# virginia samples might lead to the l/e samples being 'drowned out' ? Also,
# ComBat_seq was running for >6 hours and still hadn't finished.

# filter the merged 3 datasets using parameters used in initial analyses (gives
# 31,654 genes across 341 samples)
preg.np.filt <- codaSeq.filter(preg.np.data, samples.by.row=F,
                                min.occurrence=0.30, min.prop=0.00005)

# save unfiltered and filtered merged dataframes of all datasets
# save(preg.np.data,
#      file = paste(path.to.github, "Rdata/preg.np.data.Rda", sep = ""))
# save(preg.np.filt,
#      file = paste(path.to.github, "Rdata/preg.np.filt.Rda", sep = ""))

# make taxa vector from filtered preg vs. non-preg dataset and save to .Rda
preg.np.tax.vec <- tax.table[rownames(preg.np.filt),2]
names(preg.np.tax.vec) <- rownames(preg.np.filt)
# save(preg.np.tax.vec,
#      file = paste(path.to.github, "Rdata/preg.np.tax.vec.Rda", sep = ""))

# create vectors for batches (1/2/3) and condition (health/bv)
dataset <- rep(c(1, 2, 3), c(22,22,297))
health <- c(rep(c("Healthy","BV","BV","Healthy"), c(8,14,14,8)),
            virginia.groups$groups.2)

# try batch correction on the filtered pregnant vs. non-pregnant dataset
# (takes ~5 mins rather than >6 hours!)
preg.np.filt.bc <- sva::ComBat_seq(counts = as.matrix(preg.np.filt),
                                   batch = dataset, group = health,
                                   full_mod = TRUE)

# check that the rownames of the batch-corrected dataset are equal to the names
# of the filtered taxa vector (should be 31,654 TRUE)
length(which(rownames(preg.np.filt.bc) == names(preg.np.tax.vec))) # yep!

# save batch-corrected, filtered london/europe/virginia dataset as .Rda
# save(preg.np.filt.bc,
#      file = paste(path.to.github, "Rdata/preg.np.filt.bc.Rda", sep = ""))

####################### Build EggNOG/EC/KO .Rda objects #######################

# read in eggNOG / EC / KO tables

egg <- read.table(paste(path.to.github,"1_VIRGO/3.eggnog.NOG.txt", sep=""),
                  header=F, row.names=2, sep="\t", stringsAsFactors=F)

EC <- read.table(paste(path.to.github,"1_VIRGO/5.EC.txt", sep=""),
                 header=F, row.names=2, sep="\t", stringsAsFactors=F, quote='')

KO <- read.table(paste(path.to.github,"1_VIRGO/8.A.kegg.ortholog.txt", sep=""),
                 header=F, row.names=1, sep="\t", stringsAsFactors=F, quote='')


# eggNOG - filtered data
egg.both.vec <- egg[rownames(new.both.filt),6]
names(egg.both.vec) <- rownames(new.both.filt)
egg.both <- aggregate(new.both.filt, by=list(egg.both.vec), FUN=sum)
rownames(egg.both) <- egg.both$Group.1
egg.both$Group.1 <- NULL
# save(egg.both, file = paste(path.to.github, 'Rdata/egg.both.Rda', sep = ""))

# eggNOG - all data
egg.both.vec <- egg[rownames(new.both.data),6]
names(egg.both.vec) <- rownames(new.both.data)
egg.both.all <- aggregate(new.both.data, by=list(egg.both.vec), FUN=sum)
rownames(egg.both.all) <- egg.both.all$Group.1
egg.both.all$Group.1 <- NULL
# save(egg.both.all, file = paste(path.to.github,
#                                 'Rdata/egg.both.all.Rda', sep = ""))



# EC - filtered data
ec.both.vec <- EC[rownames(new.both.filt),2]
names(ec.both.vec) <- rownames(new.both.filt)
ec.both <- aggregate(new.both.filt, by=list(ec.both.vec), FUN=sum)
rownames(ec.both) <- ec.both$Group.1
ec.both$Group.1 <- NULL
# save(ec.both, file = paste(path.to.github, 'Rdata/ec.both.Rda', sep = ""))

# EC - all data
ec.both.vec <- EC[rownames(new.both.data),2]
names(ec.both.vec) <- rownames(new.both.data)
ec.both.all <- aggregate(new.both.data, by=list(ec.both.vec), FUN=sum)
rownames(ec.both.all) <- ec.both.all$Group.1
ec.both.all$Group.1 <- NULL
# save(ec.both.all, file = paste(path.to.github,
#                                'Rdata/ec.both.all.Rda', sep = ""))



# KO - filtered data
ko.both.vec <- KO[rownames(new.both.filt),1]
names(ko.both.vec) <- rownames(new.both.filt)
ko.both <- aggregate(new.both.filt, by=list(ko.both.vec), FUN=sum)
rownames(ko.both) <- ko.both$Group.1
ko.both$Group.1 <- NULL
# save(ko.both, file = paste(path.to.github,'Rdata/ko.both.Rda', sep = ""))

# KO - all data
ko.both.vec <- KO[rownames(new.both.data),1]
names(ko.both.vec) <- rownames(new.both.data)
ko.both.all <- aggregate(new.both.data, by=list(ko.both.vec), FUN=sum)
rownames(ko.both.all) <- ko.both.all$Group.1
ko.both.all$Group.1 <- NULL
# save(ko.both.all, file = paste(path.to.github,
#                                'Rdata/ko.both.all.Rda', sep = ""))

################## Describe genes filtered out of new dataset ##################

# load in original 'both.data' object
load(paste(path.to.github, "Rdata/both.data.Rda", sep = ""))

# pull v numbers found in original 'both.data' that are no longer present in
# 'new.both.data (returns 17483 genes)

dropped.genes<-setdiff(rownames(both.data),rownames(new.both.data))

# assign taxonomy to genes and (in case where V numbers didn't have a known
# taxonomy, add the V number label back in and edit column 3 header)

dropped.genes.tax<-tax.table[dropped.genes,]
rownames(dropped.genes.tax)<-dropped.genes
colnames(dropped.genes.tax)[2]<- "taxonomic.id"

# replace NA values with text

dropped.genes.tax$taxonomic.id[is.na(dropped.genes.tax$taxonomic.id)]<-"No_taxa_info"

# count occurrences of unique taxonomic IDs removed from new dataset

unique.drop.tax <- vector()

for(tax in levels(factor(dropped.genes.tax$taxonomic.id))){
  unique.drop.tax[tax]<-length(which(dropped.genes.tax$taxonomic.id==tax))
}

# assign KEGG functions to dropped V numbers, if function is known

dropped.genes.tax$pathway <- path.table[rownames(dropped.genes.tax),2]

# replace NA values with text

dropped.genes.tax$pathway[is.na(dropped.genes.tax$pathway)]<-"No_pathway_info"

# count occurrences of unique taxonomic IDs removed from new dataset

unique.drop.path <- vector()

for(path in levels(factor(dropped.genes.tax$pathway))){
  unique.drop.path[path]<-length(which(dropped.genes.tax$pathway==path))
}

# bind above two vectors to data frame

length(unique.drop.path)<-11
dropped.uniques<-data.frame(cbind(tax.name = names(unique.drop.tax),
                                 tax.dropped = unique.drop.tax,
                                 path.name = names(unique.drop.path),
                                 path.dropped = unique.drop.path),
                           row.names = NULL)

# get subset of both.data containing counts of dropped genes across all samples

dropped.both <- both.data[dropped.genes,]

# sum reads aligning to each 'taxon' for all dropped genes, across all samples
no.reads.dropped<-list()

for(i in levels(factor(dropped.genes.tax$taxonomic.id))){
  no.reads.dropped[[i]] <- colSums(dropped.both[which(dropped.genes.tax$taxonomic.id == i),])
}

no.reads.dropped.df <- as.data.frame(do.call(rbind, no.reads.dropped))

# most reads dropped are of unknown taxonomy - no surprise
# max no reads dropped from single sample = 23,643

summary(t(no.reads.dropped.df)[,10])
