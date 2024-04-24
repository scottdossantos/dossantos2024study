# filtering gardnerella blast results: all datasets

# this script was used for filtering Gardnerella BLAST results against the
# VIRGO database and generating the .Rda objects used in the analyses in
# '[dataset]_gardnerella_biplots.R' scripts.

#################################### setup ####################################

library(dplyr) # for data manipulation

# set path to github directory
user <- "sds/cc"
path.to.github <- switch(user,
                         "sds/cc" = "~/Documents/GitHub/metatranscriptome/",
                         "gg" = "~/Documents/0_git/projects/metatranscriptome/")

# define column headers (for more info on what these are exactly, see
# https://www.metagenomics.wiki/tools/blast/blastn-output-format-6)
headers.blast<- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                  "qstart", "qend", "sstart", "send", "evalue", "bitscore")

###################### london & europe: filter BLAST hits ######################

# load BLAST results for london / europe datasets (and add in column headers)
blast.leopoldii <- read.csv(paste0(path.to.github, "Rdata/blast_g_leopoldii.tsv"),
                          sep = "\t", header = FALSE, col.names = headers.blast)

blast.piotii <- read.csv(paste0(path.to.github, "Rdata/blast_g_piotii.tsv"), 
                       sep = "\t", header = FALSE, col.names = headers.blast)

blast.swidsinskii <- read.csv(paste0(path.to.github, "Rdata/blast_g_swidsinskii.tsv"),
                            sep = "\t", header = FALSE, col.names = headers.blast)

blast.vaginalis <- read.csv(paste0(path.to.github, "Rdata/blast_g_vaginalis.tsv"),
                          sep = "\t", header = FALSE, col.names = headers.blast)

# filter dataframes to retain only hits with >98% sequence similarity
outputs.blast <- list()

for(i in ls(pattern = "blast.")){
  outputs.blast[[i]] <- get(i)
  outputs.blast[[i]] <- outputs.blast[[i]] %>% 
    filter(pident > 98)
  
  assign(gsub("blast", "blast.98", i), outputs.blast[[i]])
}

# write london/europe .Rda objects
# save(blast.98.leopoldii, file = paste0(path.to.github,
#                                        "Rdata/blast_98pc_leopoldii.Rda"))
# 
# save(blast.98.piotii, file = paste0(path.to.github, 
#                                     "Rdata/blast_98pc_piotii.Rda"))
# 
# save(blast.98.swidsinskii, file = paste0(path.to.github, 
#                                          "Rdata/blast_98pc_swidsinskii.Rda"))
# 
# save(blast.98.vaginalis, file = paste0(path.to.github, 
#                                        "Rdata/blast_98pc_vaginalis.Rda"))

######################### virginia: filter BLAST hits #########################

# load BLAST results for virginia dataset (and add in column headers)
blast.virginia.leopoldii <- read.csv(paste0(path.to.github, "Rdata/blast_virginia_g_leopoldii.tsv"),
                            sep = "\t", header = FALSE, col.names = headers.blast)

blast.virginia.piotii <- read.csv(paste0(path.to.github, "Rdata/blast_virginia_g_piotii.tsv"), 
                         sep = "\t", header = FALSE, col.names = headers.blast)

blast.virginia.swidsinskii <- read.csv(paste0(path.to.github, "Rdata/blast_virginia_g_swidsinskii.tsv"),
                              sep = "\t", header = FALSE, col.names = headers.blast)

blast.virginia.vaginalis <- read.csv(paste0(path.to.github, "Rdata/blast_virginia_g_vaginalis.tsv"),
                            sep = "\t", header = FALSE, col.names = headers.blast)

# filter dataframes to retain only hits with >98% sequence similarity
outputs.blast <- list()

for(i in ls(pattern = "blast.virginia.")){
  outputs.blast[[i]] <- get(i)
  outputs.blast[[i]] <- outputs.blast[[i]] %>% 
    filter(pident > 98)
  
  assign(gsub("blast", "blast.98", i), outputs.blast[[i]])
}

# write virginia .Rda objects
# save(blast.98.virginia.leopoldii, file = paste0(path.to.github,
#                                                 "Rdata/blast_98pc_virginia_leopoldii.Rda"))
# 
# save(blast.98.virginia.piotii, file = paste0(path.to.github,
#                                              "Rdata/blast_98pc_virginia_piotii.Rda"))
# 
# save(blast.98.virginia.swidsinskii, file = paste0(path.to.github,
#                                                   "Rdata/blast_98pc_virginia_swidsinskii.Rda"))
# 
# save(blast.98.virginia.vaginalis, file = paste0(path.to.github,
#                                                 "Rdata/blast_98pc_virginia_vaginalis.Rda"))

###################### all 3 datasets: filter BLAST hits ######################

# load BLAST results for all 3 combined datasets (and add in column headers)
blast.all3.leopoldii <- read.csv(paste0(path.to.github, "Rdata/blast_all3_g_leopoldii.tsv"),
                            sep = "\t", header = FALSE, col.names = headers.blast)

blast.all3.piotii <- read.csv(paste0(path.to.github, "Rdata/blast_all3_g_piotii.tsv"), 
                         sep = "\t", header = FALSE, col.names = headers.blast)

blast.all3.swidsinskii <- read.csv(paste0(path.to.github, "Rdata/blast_all3_g_swidsinskii.tsv"),
                              sep = "\t", header = FALSE, col.names = headers.blast)

blast.all3.vaginalis <- read.csv(paste0(path.to.github, "Rdata/blast_all3_g_vaginalis.tsv"),
                            sep = "\t", header = FALSE, col.names = headers.blast)

# filter dataframes to retain only hits with >98% sequence similarity
outputs.blast <- list()

for(i in ls(pattern = "blast.all3.")){
  outputs.blast[[i]] <- get(i)
  outputs.blast[[i]] <- outputs.blast[[i]] %>% 
    filter(pident > 98)
  
  assign(gsub("blast", "blast.98", i), outputs.blast[[i]])
}

# write london/europe .Rda objects
# save(blast.98.all3.leopoldii, file = paste0(path.to.github,
#                                        "Rdata/blast_98pc_all3_leopoldii.Rda"))
# 
# save(blast.98.all3.piotii, file = paste0(path.to.github,
#                                     "Rdata/blast_98pc_all3_piotii.Rda"))
# 
# save(blast.98.all3.swidsinskii, file = paste0(path.to.github,
#                                          "Rdata/blast_98pc_all3_swidsinskii.Rda"))
# 
# save(blast.98.all3.vaginalis, file = paste0(path.to.github,
#                                        "Rdata/blast_98pc_all3_vaginalis.Rda"))
