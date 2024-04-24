# association of pre-term birth with vaginal metatranscriptome

# The data from the MOMS-PI study on pre-term birth is considered too personal
# to be hostd on dbGaP and will not be included in an authorised access request.
# However, the MOMS-PI study authors are more than happy to share the data.
# If you have an authorised access request from dbGaP approved and would like 
# access to this data, please contact Dr. Gregory Buck from the MOMS-PI study
# team.

# We have included the code below for transparency on how we analysed this data.
# The file containing the pre-term birth data is not uploaded to the study's
# github repository and can only be accessed by requesting it from the MOMS-PI
# authors.

#################################### setup ####################################

library(dplyr) # data manipulation
library(stringr) # string manipulation
library(ALDEx2) # for CLR transformations and correlation

user <- "sds/cc"
path.to.github <- switch(user,
                         "sds/cc" = "~/Documents/GitHub/metatranscriptome/",
                         "gg" = "~/Documents/0_git/projects/metatranscriptome/")

# load in virginia metadata, grouping object, KO-aggregated and vNumber feature
# tables, virginia taxa vector, and VIRGO taxa table
load(paste(path.to.github, "Rdata/virginia.meta.Rda", sep = ""))
load(paste(path.to.github, "Rdata/virginia.groups.Rda", sep = ""))
load(paste(path.to.github, "Rdata/ko.virginia.filt.Rda", sep = ""))
load(paste(path.to.github, "Rdata/virginia.filt.Rda", sep = ""))
load(paste(path.to.github, "Rdata/tax.vec.virginia.Rda", sep = ""))

# load in taxa table from github
tax.table <- read.table(paste0(path.to.github, "1_VIRGO/1.taxon.tbl.txt"),
                        header = FALSE, sep = "\t", quote = "", row.names = 2)

# load in PTB data from local machine
virginia.ptb <- read.table("~/Documents/PostDoc_Western/Data/VaginalMetatranscriptome/R_work/virginia_PTBdata.txt",
                           header = TRUE, sep = "\t", quote = "")

################################ assigning PTB ################################

# split gestational age at delivery column into numerical weeks and days columns
virginia.ptb$weeks <- as.integer(str_sub(virginia.ptb$ga_atDelivery, 1, 2))
virginia.ptb$days <- as.integer(str_sub(virginia.ptb$ga_atDelivery, -2, -2))

# filter out any rows where GA at delivery is unknown
virginia.ptb <- virginia.ptb %>% 
  filter(!is.na(weeks))

# convert days as a proportion of a week
virginia.ptb$days.dec <- virginia.ptb$days/7

# add weeks and days columns to get decimal gestational age
virginia.ptb$ga.dec <- virginia.ptb$weeks + virginia.ptb$days.dec

# create new variable indicating if baby was born pre-term or not (PTB usually
# defined as birth prior to 37 weeks gestation)
virginia.ptb$ptb <- case_when(virginia.ptb$ga.dec< 37 ~ "PTB",
                              .default = "Non-PTB")

# subset the ptb data for only samples present in metatranscriptome dataset
virginia.ptb <- virginia.ptb %>% 
  filter(subject_ID %in% virginia.meta$subject_id)

# filter metadata for samples with known PTB data
virginia.meta <- virginia.meta %>% 
  filter(subject_id %in% virginia.ptb$subject_ID)

# create ptb vector
ptb <- vector(length = length(virginia.meta))
for(i in 1:length(virginia.meta$subject_id)){
  index <- grep(virginia.meta$subject_id[i], virginia.ptb$subject_ID)
  ptb[i] <- virginia.ptb$ptb[index]
  }

# create ga at del vector
ga <- vector(length = length(virginia.meta))
for(i in 1:length(virginia.meta$subject_id)){
  index <- grep(virginia.meta$subject_id[i], virginia.ptb$subject_ID)
  ga[i] <- virginia.ptb$ga.dec[index]
}

ptb.meta <- virginia.meta
ptb.meta$ga.del <- ga
ptb.meta$ptb <- ptb
ptb.meta <- ptb.meta %>% 
  relocate(c(ga.del, ptb), .after = sra_id)

############################# PTB correlation: KO #############################

# remove five samples lacking PTB data from virginia grouping object and 
# KO-aggregated feature table
ptb.virginia.groups <- virginia.groups[virginia.meta$sra_id,]
ptb.ko.virginia.filt <- ko.virginia.filt[,virginia.meta$sra_id]

# check for rows with zero sums after removal of above five samples
any(rowSums(ptb.ko.virginia.filt) == 0) # returns FALSE

# set seed for RNG
set.seed(2023)

# make a scale matrix for running aldex w/ scale simulation
mu.matrix <- aldex.makeScaleMatrix(gamma = 0.5, mu = c(1, 1.15),
                                   ptb.virginia.groups$groups.2)

# CLR transformation with scale-reliant inference
ptb.ko.virginia.filt.clr <- aldex.clr(reads = ptb.ko.virginia.filt,
                                      conds = ptb.virginia.groups$groups.2,
                                      mc.samples = 128, denom = "all",
                                      gamma = mu.matrix, verbose = TRUE)

# correlate gestational age at delivery and KOs present in filtered virginia 
# dataset
ptb.corr.ko <- aldex.corr(ptb.ko.virginia.filt.clr, ptb.meta$ga.del)

# only a single KO is significantly correlated with GA across 3 different 
# measures of correlation when using an FDR threshold of 1% as in previous
# analyses: 

#KO3771 (surA, peptidyl-prolyl cis-trans isomerase)
#             corr        ePval         ePvalBH
# pearson:    -0.2301185  0.0002021394  0.009993878
# spearman:   -0.2210061  0.0002919605  0.01273319
# kendall:    -0.15239437 0.0002844615  0.01288342

########################### PTB correlation: species ###########################

# calculate read sums for each species in the dataset, one list element per
# species
reads.spp <- list()
for(i in levels(factor(tax.vec.virginia))){
  reads.spp[[i]] <- colSums(virginia.filt[which(tax.vec.virginia == i),])
}

# collapse list to dataframe and remove five samples lacking PTB data
reads.spp.df <- as.data.frame(do.call(rbind, reads.spp))
reads.spp.df <- reads.spp.df[,ptb.meta$sra_id]

# remove samples lacking PTB data from virginia grouping table
ptb.virginia.groups <- virginia.groups[virginia.meta$sra_id,]

# set seed for RNG
set.seed(2023)

# make a scale matrix for running aldex w/ scale simulation
mu.matrix <- aldex.makeScaleMatrix(gamma = 0.5, mu = c(1, 1.15),
                                   ptb.virginia.groups$groups.2)

# CLR transformation with scale-reliant inference
ptb.spp.virginia.filt.clr <- aldex.clr(reads = reads.spp.df,
                                       conds = ptb.virginia.groups$groups.2,
                                       mc.samples = 128, denom = "all",
                                       gamma = mu.matrix, verbose = TRUE)

# correlate gestational age at delivery and species present in filtered virginia 
# dataset
ptb.corr.spp <- aldex.corr(ptb.spp.virginia.filt.clr, ptb.meta$ga.del)

# several taxa correlated significantly (FDR <1%) with GA at delivery across
# all 3 measures:
#         - L. crispatus positively corr
#         - BVAB1/P.amnii/P.timonensis/S.amnii negatively corr

####################### PTB correlation: KO BV subgroups #######################

# load in vector with virginia BV subgroup samples
load(paste0(path.to.github, "Rdata/virginia.bv.groups.Rda"))

# subset ptb metadata for only bv subgroup samples
sub.ptb.meta <- ptb.meta
rownames(sub.ptb.meta) <- sub.ptb.meta$sra_id
sub.ptb.meta <- sub.ptb.meta[colnames(sub.reads.spp.df),]

# subset KO-aggregated feature table for only bv subgroup samples
sub.ko.filt <- ko.virginia.filt[,rownames(virginia.bv.groups)]
any(rowSums(sub.ko.filt) == 0) # returns FALSE

# CLR-transform w/ scale
set.seed(2023)

sub.mu.matrix <- aldex.makeScaleMatrix(gamma = 0.5, mu = c(1,1.15),
                                       virginia.bv.groups$group)

sub.ptb.ko.clr <- aldex.clr(reads = sub.ko.filt,
                            conds = virginia.bv.groups$group,
                            mc.samples = 128, denom = "all",
                            verbose = TRUE, gamma = sub.mu.matrix)

# correlate gestational age at delivery and KOs present among bv subgroup
# samples in filtered virginia dataset
sub.ptb.corr.ko <- aldex.corr(sub.ptb.ko.clr, sub.ptb.meta$ga.del)

# no KO terms correlated with gestational age at delivery across any measure

#################### PTB correlation: species BV subgroups ####################

# subset species-aggregated data frame by bv subgroups
sub.reads.spp.df <- reads.spp.df[,rownames(virginia.bv.groups)]
any(rowSums(sub.reads.spp.df) == 0) # returns FALSE

# CLR-transform w/ scale
set.seed(2023)

sub.mu.matrix <- aldex.makeScaleMatrix(gamma = 0.5, mu = c(1,1.15),
                                       virginia.bv.groups$group)

sub.ptb.spp.clr <- aldex.clr(reads = sub.reads.spp.df,
                             conds = virginia.bv.groups$group,
                             mc.samples = 128, denom = "all",
                             verbose = TRUE, gamma = sub.mu.matrix)

# correlate gestational age at delivery and species present among bv subgroup
# samples in filtered virginia dataset
sub.ptb.corr.spp <- aldex.corr(sub.ptb.spp.clr, sub.ptb.meta$ga.del)

# nothing correlated at all: all eBH p vals >0.01

# welch's t-test of ga at delivery for BV1 vs. BV2
t.test(sub.ptb.meta$ga.del~virginia.bv.groups$group, var.equal = FALSE)

# mean BV1 = 36.52632
# mean BV2 = 38.70663
# t = -3.0902    df = 97.218     p-value = 0.002609
