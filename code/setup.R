# load the minimum analysis tools

library(zCompositions)
library(ALDEx2)
library(CoDaSeq)

# set the user when you run this code
# gg sets path for greg
# cc sets path for scott/clara

user <- 'cc'
if(user=='gg'){
	locn <- "~/Documents/0_git/projects/dossantos2024study/"
}else if(user=='cc'){
	locn <- "~/Documents/GitHub/dossantos2024study/"
}
# load in the lookup tables
# data tables are loaded in merge_table.R

# standard color table
tax.colors <- read.table(paste(locn,"code/species_colors.txt", sep=""), sep="\t",
    header=T, row.names=1, stringsAsFactors=F)

# VIRGO taxon table
tax.table <- read.table(paste(locn,'1_VIRGO/1.taxon.tbl.txt', sep=""),
    header=F, row.names=2)

# load eggnog / EC / KO lookup tables
egg <- read.table(paste(locn,"1_VIRGO/3.eggnog.NOG.txt", sep=""),
                  header=F, row.names=2, sep="\t", stringsAsFactors=F)

egg.copy <- read.table(paste(locn,"1_VIRGO/3.eggnog.NOG copy.txt", sep=""),
                       header=F, row.names=2, sep="\t", stringsAsFactors=F)

EC <- read.table(paste(locn,"1_VIRGO/5.EC.txt", sep=""),
                 header=F, row.names=2, sep="\t", stringsAsFactors=F, quote='')

KO <- read.table(paste(locn,"1_VIRGO/8.A.kegg.ortholog.txt", sep=""),
                 header=F, row.names=1, sep="\t", stringsAsFactors=F, quote='')

# get a standard filtered set of data
# merge the groups together each group

# load in data tables
# make sure these are all available in one directory
# make the merged sample set of all twntyfr data and the pre and
# post metronidazole samples from ENA dataset

# London data: d.24
# London metadata: m.24
# London taxonomy table: d.24.tax

# ENA data: nr.ena
# ENA metadata: m.ena
# ENA taxonomy table:  nr.ena.tax

# combined metadata: all.meta (age, pH, Nugent scores x 44 samples)
# combined data: new.both.data (44 samples, 330,123 genes)
# combined filtered data: new.both.filt (44 samples, 22,606 genes: >0.3 prev > 0.00005 abund)
# taxonomy lookup table: new.tax.vec (22,606 genes with taxonomy)

# EC lookup table: EC (111,599 enzyme classifications)
# EC collapsed function table: ec.both (603 enzyme classifications)
# EC collapsed all gene table: ec.both.all (1,050 enzyme classifications)

# eggNOG lookup table (314,735 functions)
# eggnog collapsed function table: egg.both (2,097 functions)
# eggnog collapsed all gene table: egg.both.all (8,480 functions)

# KO lookup table: KO (415,247 ontology terms)
# KO collapsed function table: ko.both (1,708 ontology terms)
# KO collapsed all gene table: ko.both.all (3,666 ontology terms)


# controls if data regenerated from scratch 
# default is to load precomputed data

scratch <- FALSE 
# keep this false as GG/SDS will generate the .Rda files


if(file.exists(paste(locn,'Rdata/d.24.Rda',sep = ""))){
  if(!file.exists(paste(locn,'Rdata/d.24.Rda',sep = ""))) stop("place the data files in the Rdata directory")
  if(!file.exists(paste(locn,'Rdata/m.24.Rda',sep = ""))) stop("place the data files in the Rdata directory")
  if(!file.exists(paste(locn,'Rdata/d.24.tax.Rda',sep = ""))) stop("place the data files in the Rdata directory")
  print('scratch set to FALSE: loading previously generated data (T2T-mapped!)')
  
  # London data
  print('loading London data (T2T-mapped!)') 
  load(paste(locn,'Rdata/d.24.Rda',sep = ""))       # haven't updated object as we don't use it
  load(paste(locn,'Rdata/m.24.Rda',sep = ""))       # haven't updated object as we don't use it
  load(paste(locn,'Rdata/d.24.tax.Rda',sep = ""))   # haven't updated object as we don't use it

  # Europe data
  print('loading Europe (ENA) data (T2T-mapped!)')
  load(paste(locn,'Rdata/nr.ena.Rda',sep = ""))     # haven't updated object as we don't use it
  load(paste(locn,'Rdata/m.ena.Rda',sep = ""))      # haven't updated object as we don't use it
  load(paste(locn,'Rdata/nr.ena.tax.Rda',sep = "")) # haven't updated object as we don't use it

  # Merged data
  print('loading new merged data (T2T-mapped!)')
  load(paste(locn,'Rdata/all.meta.Rda',sep = ""))        # metadata (unchanged)
  load(paste(locn,'Rdata/new.both.data.Rda',sep = ""))   # new 'all' dataset
  load(paste(locn,'Rdata/new.both.filt.Rda',sep = ""))   # new 'filt' dataset
  load(paste(locn,'Rdata/new.tax.vec.Rda',sep = ""))     # new taxa object
  
  # eggNOG / EC / KO tables
  print('loading new eggNOG / EC / KO tables (T2T-mapped!)')
  load(paste(locn,'Rdata/coglist.Rda',sep = ""))         # ---- these
  load(paste(locn,'Rdata/egg.both.Rda',sep = ""))        # ---- are
  load(paste(locn,'Rdata/egg.both.all.Rda',sep = ""))    # ---- all
  load(paste(locn,'Rdata/ec.both.all.Rda',sep = ""))     # ---- new
  load(paste(locn,'Rdata/ec.both.Rda',sep = ""))         # ---- T2T-mapped
  load(paste(locn,'Rdata/ko.both.Rda',sep = ""))         # ---- .Rda
  load(paste(locn,'Rdata/ko.both.all.Rda',sep = ""))     # ---- objects
} 

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# THE CODE BELOW IS FOR GENERATING THE OLD FEATURE TABLES PRIOR TO T2T MAPPING
# TO GENERATE NEW DATA FROM SCRATCH, RUN '/code/merge_feature_tables.R' INSTEAD!

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

if(scratch == TRUE){ #get the data from source
  locn <- '~/Documents/0_git/twntyfr/VIRGO/'
  d.24.init <- 
  read.table(paste(locn,'v-twntyfr/summary.NR.abundance.txt', sep=''), 
  header=T, row.names=1, check.names=F)
  # starting with 2492333 entries
  m.24 <- read.table("~/Documents/0_git/twntyfr/data/metadata_full.txt",
	header=T, row.names=1, sep="\t", stringsAsFactors=F)

  # pull out the H and BV samples
  h.24 <- c("001B", "002B", "004B", "006B", "009B", "010B", "015B", "020B")

  bv.24 <- c("001A", "003A", "006A", "008A", "009A", "010A", "012B",
	  "013B", "014B", "018B", "017B", "012A", "013A", "019A")

  d.24 <- data.frame(d.24.init[,h.24], d.24.init[,bv.24], check.names=F)

  d.24 <- d.24[rowSums(d.24) > 0,]
  colnames(d.24) <- c(paste('h.',h.24, sep=''), paste('v.', bv.24, sep=''))

  # get list of taxa in d.24
  tax_in_d.24 <- rownames(tax.table) %in% rownames(d.24)
  rows.tax <- rownames(tax.table)[tax_in_d.24]
  d.24.tax <- d.24[rows.tax,]
  
  save(d.24, file='Rdata/d.24.Rda')
  save(m.24, file='Rdata/m.24.Rda')
  save(d.24.tax, file='Rdata/d.24.tax.Rda')

  d.ena <- read.table('~/Documents/0_git/twntyfr/VIRGO/v-ENA/summary.NR.abundance.txt', header=T, row.names=1)
  # need to get this from backup
  m.ena <- read.table("~/Documents/0_git/Log-Ratio-Publication/data/pub.metadata.txt", header=T, row.names=1, sep="\t", comment.char="", quote="")

  # this pulls out the time 1 BV samples
  ena.bv <- c("ERR2014354", "ERR2014359", "ERR2014361", "ERR2014363", "ERR2014365", "ERR2014370", "ERR2014372", "ERR2014374", "ERR2014376", "ERR2014378", "ERR2014383", "ERR2014388", "ERR2014390", "ERR2014392")
  ena.bv.names <- paste('v.',c("a4_1", "a5_1", "a5_12", "a6_1", "a6_4", "a6_6", "a6_10", "a7_1", "a8_3", "a8_6", "a13_19", "a13_22", "a15_6", "a17_4"), sep='')

  # this pulls out the time 2 N samples after metronidazole treatment
  ena.n <- c("ERR2014355", "ERR2014362", "ERR2014366", "ERR2014373", "ERR2014377", "ERR2014379", "ERR2014384", "ERR2014389")
  ena.n.names <- paste('h.',c("b4_1", "b5_12", "b6_4", "b6_10", "b8_3", "b8_6", "b13_19", "b13_22"), sep='')

  nr.ena <- data.frame(d.ena[,ena.bv], d.ena[,ena.n])
  colnames(nr.ena) <- c(ena.bv.names, ena.n.names)
  nr.ena <- nr.ena[rowSums(nr.ena) > 0,]
  # get list of taxa in ena
  tax_in_nr.ena <- rownames(tax.table) %in% rownames(nr.ena)
  rows.tax <- rownames(tax.table)[tax_in_nr.ena]
  nr.ena.tax <- nr.ena[rows.tax,]
  
  save(nr.ena, file='Rdata/nr.ena.Rda')
  save(m.ena, file='Rdata/m.ena.Rda')
  save(nr.ena.tax, file='Rdata/nr.ena.tax.Rda')


names.all <- (unique(c(rownames(nr.ena), rownames(d.24))))

# NOTE: merge() runs out of memory
#both <- merge(d.24, nr.ena, incomparables=0)

#this works and is simple to understand
both.data <- data.frame(matrix(data=0,
    ncol=ncol(d.24) + ncol(nr.ena), nrow=length(names.all)))

colnames(both.data) <- c(colnames(d.24), colnames(nr.ena))
rownames(both.data) <- names.all

both.data[rownames(d.24), 1:22] <- d.24
both.data[rownames(nr.ena), 23:44] <- nr.ena
save(both.data, file='Rdata/both.data.Rda')

# merge metadata tables
meta.names.ldn <- c(h.24, bv.24)
meta.names.ena <- c(ena.bv, ena.n)

m.alt.ena <- m.ena
rownames(m.alt.ena) <- m.ena$Run.Assembly.accession

age <- c(m.24[meta.names.ldn,"age"], m.alt.ena[meta.names.ena,'Age'])
ph <- c(m.24[meta.names.ldn,"ph"], m.alt.ena[meta.names.ena,'pH'])
nugent <- c(m.24[meta.names.ldn,"nugent_score"], m.alt.ena[meta.names.ena,'Nugent.score'])

all.meta <- data.frame(cbind(age,ph,nugent))
rownames(all.meta) <- colnames(both.data)
save(all.meta, file=('Rdata/all.meta.Rda'))

# filter to a reasonable number of features
# reduces from 347540 to 22601 features
# gene counts between 20 and 800 depending on sample read depth
# filtering to 0.000005 gives 50457 features

# must occur in 30% of samples, with a minimal proportional abundance of 
both.filt <- codaSeq.filter(both.data, min.occurrence=0.30, min.prop=0.00005, samples.by.row=F)
save(both.filt, file='Rdata/both.filt.Rda')

# assign a taxonomy to all rows
tax.vec <- tax.table[rownames(both.filt),2]
names(tax.vec) <- rownames(both.filt)
save(tax.vec, file='Rdata/tax.vec.Rda')

# group genes by eggnog classification
coglist <- read.table("~/Desktop/coglist.txt", header=T, stringsAsFactors=F, sep="\t",  quote="", comment.char="", row.names=2, fill=T)
save(coglist, file='Rdata/coglist.Rda')

both.vec <- egg[rownames(both.filt),6]
names(both.vec) <- rownames(both.filt)
egg.both <- aggregate(both.filt, by=list(both.vec), FUN=sum)
rownames(egg.both) <- egg.both$Group.1
egg.both$Group.1 <- NULL
save(egg.both, file='Rdata/egg.both.Rda')

both.vec <- egg[rownames(both.data),6]
names(both.vec) <- rownames(both.data)
egg.both.all <- aggregate(both.data, by=list(both.vec), FUN=sum)
rownames(egg.both.all) <- egg.both.all$Group.1
egg.both.all$Group.1 <- NULL
save(egg.both.all, file='Rdata/egg.both.all.Rda')

# EC
ec.both.vec <- EC[rownames(both.filt),2]
names(ec.both.vec) <- rownames(both.filt)
ec.both <- aggregate(both.filt, by=list(ec.both.vec), FUN=sum)
rownames(ec.both) <- ec.both$Group.1
ec.both$Group.1 <- NULL
save(ec.both, file='Rdata/ec.both.Rda')

both.vec <- EC[rownames(both.data),2]
names(both.vec) <- rownames(both.data)
ec.both.all <- aggregate(both.data, by=list(both.vec), FUN=sum)
rownames(ec.both.all) <- ec.both.all$Group.1
ec.both.all$Group.1 <- NULL
save(ec.both.all, file='Rdata/ec.both.all.Rda')

# KO numbers
ko.both.vec <- KO[rownames(both.filt),1]
names(ko.both.vec) <- rownames(both.filt)
ko.both <- aggregate(both.filt, by=list(ko.both.vec), FUN=sum)
rownames(ko.both) <- ko.both$Group.1
ko.both$Group.1 <- NULL
save(ko.both, file='Rdata/ko.both.Rda')

both.vec <- KO[rownames(both.data),1]
names(both.vec) <- rownames(both.data)
ko.both.all <- aggregate(both.data, by=list(both.vec), FUN=sum)
rownames(ko.both.all) <- ko.both.all$Group.1
ko.both.all$Group.1 <- NULL
save(ko.both.all, file='Rdata/ko.both.all.Rda')

}
