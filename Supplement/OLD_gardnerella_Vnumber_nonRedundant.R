#identifying non-redundant V numbers in vaginal metatranscriptome datasets

# This script pulls out the V numbers present in the London & Europe datasets
# and converts them to a list of FASTA headers, as found in the VIRGO database,
# located on agrajag at: /Volumes/data/twntyfr_2018/VIRGO/0_db/NT.fasta

#################################### setup ####################################

library(dplyr) # for data manipulation

# set path to github directory
# set path to github directory (edit 'path.to.github' to reflect your machine!)
user <- "sds/cc"
path.to.github <- switch(user,
                         "sds/cc" = "~/Documents/GitHub/dossantos2024study/",
                         "gg" = "~/Documents/0_git/projects/dossantos2024study/")

# load filtered feature tables for lon/eur, virginia and all 3 datasets
load(paste0(path.to.github, "Rdata/new.both.filt.Rda")) # london/europe
load(paste0(path.to.github, "Rdata/virginia.filt.Rda")) # virginia
load(paste0(path.to.github, "Rdata/preg.np.filt.bc.Rda")) # all 3

# rename for clarity
filt.lon.eur <- new.both.filt 
filt.virginia <- virginia.filt 
filt.all.3 <- as.data.frame(preg.np.filt.bc)

# remove original dataframes
rm(new.both.filt)
rm(virginia.filt)
rm(preg.np.filt.bc)

######################## london & europe datasets ########################

# sanity check for rows with zero reads across datasets (FALSE in all cases)
any(rowSums(filt.lon.eur)==0)
any(rowSums(filt.virginia)==0)
any(rowSums(filt.all.3)==0)

# pull unique V numbers (to use for BLAST/mapping against gardnerella genes)
vnum.lon.eur <-rownames(filt.lon.eur)
vnum.virginia <-rownames(filt.virginia)
vnum.all.3 <-rownames(filt.all.3)

# write headers to .txt file (note: output for lon/eur follows different naming
# convention as the gardnerella analysis was performed on this dataset first,
# and the companion 'walkthrough_GardnerellaPangenome.txt' file uses this 
# name)
# write(vnum.lon.eur, "vnumbers_to_search_lon_eur.txt")
# write(vnum.virginia, "vnumbers_to_search_virginia.txt")
# write(vnum.all.3, "vnumbers_to_search_all3.txt")
