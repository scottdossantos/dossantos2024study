################################################################################
# identifying marker genes for each gardnerella species
# Scott Dos Santos, 2024-09-10

# pulled 44 gardnerella genomes from NCBI; all <100 contigs, aside from four
# G. pickettii genomes which are >100 but <200 contigs

# -------------------------------------------------
# GENOME FILE NAME				                        N
# -------------------------------------------------
# Gardnerella_greenwoodii_00703Dmash.fasta	     11
# Gardnerella_greenwoodii_101.fasta		           43
# Gardnerella_greenwoodii_UMB1686.fasta		        7
# Gardnerella_greenwoodii_c31Ua_26.fasta	       11
# Gardnerella_leopoldii_UMB0682.fasta		          6
# Gardnerella_leopoldii_UMB0742.fasta		          4
# Gardnerella_leopoldii_UMB0913.fasta		          4
# Gardnerella_leopoldii_UMB1489.fasta		          4
# Gardnerella_leopoldii_Ugent_06.41.fasta	        1
# Gardnerella_pickettii_00703Bmash.fasta	       16
# Gardnerella_pickettii_00703C2mash.fasta	       22
# Gardnerella_pickettii_GED7275B.fasta		       18
# Gardnerella_pickettii_JCP7659.fasta		        181
# Gardnerella_pickettii_JCP7719.fasta		        145
# Gardnerella_pickettii_JCP8017A.fasta		      163
# Gardnerella_pickettii_JCP8017B.fasta		      158
# Gardnerella_pickettii_UMB0830.fasta		          4
# Gardnerella_pickettii_UMB0833.fasta		         12
# Gardnerella_pickettii_UMB9259.fasta		         21
# Gardnerella_pickettii_c17Ua_112.fasta		        3
# Gardnerella_piotii_GH007.fasta		              2
# Gardnerella_piotii_GH020.fasta		              2
# Gardnerella_piotii_JNFY15.fasta		              1
# Gardnerella_piotii_UGent_18.01.fasta		        5
# Gardnerella_swidsinskii_GS_9838-1.fasta	        9
# Gardnerella_swidsinskii_GV37.fasta		          1
# Gardnerella_swidsinskii_JNFY3.fasta		          1
# Gardnerella_swidsinskii_UMB0170.fasta		        4
# Gardnerella_swidsinskii_UMB0264.fasta		       13
# Gardnerella_swidsinskii_UMB0411.fasta		        5
# Gardnerella_swidsinskii_UMB1698.fasta		        5
# Gardnerella_vaginalis_0288E.fasta		           17
# Gardnerella_vaginalis_1400E.fasta              28
# Gardnerella_vaginalis_1500E.fasta		           27
# Gardnerella_vaginalis_284V.fasta		           16
# Gardnerella_vaginalis_315-A.fasta		           13
# Gardnerella_vaginalis_55152.fasta		           25
# Gardnerella_vaginalis_6119V5.fasta		         12
# Gardnerella_vaginalis_75712.fasta		            3
# Gardnerella_vaginalis_ATCC_14018.fasta	        1
# Gardnerella_vaginalis_ATCC_14018-2.fasta	      2
# Gardnerella_vaginalis_ATCC_14018-3.fasta	     14
# Gardnerella_vaginalis_ATCC_14019.fasta	        1
# Gardnerella_vaginalis_HMP9231.fasta		          1

# annotated with bakta and defined pangenome from resulting GFF files with 
# panaroo. Now want to parse the panaroo presence/absence file and identify
# genes which are unique to each species. These will be considered 'marker'
# genes for each Gardnerella species. The VIRGO taxonomy will be updated
# accordingly and all vNumbers previously assigned to 'Gardnerella_vaginalis'
# will be renamed 'Gardnerella'

#################################### setup ####################################

library(dplyr)

setwd("~/Documents/GitHub/dossantos2024study")
wd <- paste0(getwd(),"/")

# load in presence/absence table
pa <- read.csv(paste0(wd,"Rdata/panaroo_gene_presence_absence.csv"), header = T,
               quote = "", row.names = 1, na.strings = "")

# retain only binary rows
pa.bin <- pa %>% 
  select(-c(Non.unique.Gene.name,Annotation))

# convert all values in pres/abs df to binary 1/0
for(col in 1:ncol(pa.bin)){
  pa.bin[,col] <- case_when(is.na(pa.bin[,col]) ~ 0, .default = 1)
}

# replace colnames with shortened version of species & number
colnames(pa.bin) <-c(paste0("g_gw",1:4),paste0("g_lp",1:5),paste0("g_pk",1:11),
                     paste0("g_po",1:4),paste0("g_sw",1:7),paste0("g_va",1:13))

# column references: greenwoodii = 1:4, leopoldii = 5:9, pickettii = 10:20,
#                    piotii = 21:24, swidsinskii = 25:31, vaginalis = 32:44

# identify marker genes for each species
spp <- vector(length = nrow(pa.bin))
for(row in 1:nrow(pa.bin)){
  to.scan <- pa.bin[row,]
  if(sum(to.scan) == 44){
    spp[row] <- "Gardnerella"
  } else{
    spp[row] <- case_when(any(to.scan[1:4] == 1) && sum(to.scan[5:44]) == 0 ~ "Gardnerella_greenwoodii",
                          any(to.scan[5:9] == 1) && sum(to.scan[c(1:4, 10:44)]) == 0 ~ "Gardnerella_leopoldii",
                          any(to.scan[10:20] == 1) && sum(to.scan[c(1:9, 21:44)]) == 0 ~ "Gardnerella_pickettii",
                          any(to.scan[21:24] == 1) && sum(to.scan[c(1:20, 25:44)]) == 0 ~ "Gardnerella_piotii",
                          any(to.scan[25:31] == 1) && sum(to.scan[c(1:24, 32:44)]) == 0 ~ "Gardnerella_swidsinskii",
                          any(to.scan[32:44] == 1) && sum(to.scan[c(1:31)]) == 0 ~ "Gardnerella_vaginalis",
                          .default = "Gardnerella")
  }
}

pa.bin.df <- data.frame(spp)
rownames(pa.bin.df) <- rownames(pa)

# count how many marker genes per species
pa.bin.df %>% 
  group_by(spp) %>% 
  count()

#  spp                         n
# 
# Gardnerella              2,217
# Gardnerella_greenwoodii     62
# Gardnerella_leopoldii      101
# Gardnerella_pickettii      175
# Gardnerella_piotii          60
# Gardnerella_swidsinskii    273
# Gardnerella_vaginalis      432

# read in pangenome headers fro 'pan_genome_reference.fa'
headers <- read.table(paste0(wd,"Rdata/panaroo_headers.txt"),
                      header = F, sep = "\t", quote = "")

# identify which headers are present in the pres/abs table but not in the
# gene FASTA file produced by panaroo
pa.bin.df.notinfasta <- pa.bin.df[which(rownames(pa) %in% headers$V1), , drop = F]

# recount marker genes per taxon (minus suspected paralogs)
pa.bin.df.notinfasta %>% 
  group_by(spp) %>% 
  count()

#  spp                         n
# 
# Gardnerella              2,176
# Gardnerella_greenwoodii     62
# Gardnerella_leopoldii      101
# Gardnerella_pickettii      175
# Gardnerella_piotii          60
# Gardnerella_swidsinskii    272
# Gardnerella_vaginalis      430

# only lost 3 marker genes, the remaining 41 were ambiguously 'Gardnerella'

# now load in the taxonomy lookup table and pull the vNumbers that are assigned
# to 'Gardnerella_vaginalis'
tax <- read.table(paste0(wd,"1_VIRGO/1.taxon.tbl.txt"),
                  header = F, sep = "\t", quote = "", row.names = 2)

# filter for G.vaginalis
tax.gard <- tax %>% 
  filter(V3 == "Gardnerella_vaginalis")

# these are the vNumbers to extract from the VIRGO NT.fasta file; write to file
write(rownames(tax.gard),
      file = paste0(wd,"Rdata/gardnerella_vnums.txt"),
      ncolumns = 1)

# follow steps in 'walkthrough_GardnerellaPangenome_REVISED.txt' to make BLAST
# db of the above gardnerella genes in VIRGO and align all genes in
# 'pan_genome_reference.fa' against it

# now read in the BLAST aligment results and gene lengths for pangenome genes
pan.length <- read.table(paste0(wd, "Rdata/panaroo_reference_length.txt"),
                         header = F, sep = "\t", quote = "", row.names = 1,
                         col.names =c("gene","length"))

blast <- read.table(paste0(wd, "Rdata/panaroo_blast.tsv"), row.names = 1,
                    header = F, sep = "\t", quote = "")

colnames(blast) <- c("sseqid", "pident", "length", "mismatch", "gapopen",
                     "qstart,", "qend", "sstart", "send", "evalue", "bitscore")

# add length of both query genes (pangenome) and subject hits (VIRGO genes)
blast$slength <- tax[blast$sseqid,"V4"]
blast$qlenth <- pan.length[rownames(blast),"length"]

# calculate alignment overlap and move columns to beside pident
# NOTE: overlap = length of aligned sequence expressed as proportion of the 
#       shorter of query and subject sequences. Where gap introduction means 
#       that alignment length > length of shorter sequence, overlap = 100 %
overlap <- vector(length = nrow(blast))
for(i in 1:nrow(blast)){
  if(blast$length[i] <= min(c(blast$slength[i], blast$qlenth[i]))){
    overlap[i] <- (blast$length[i] / min(c(blast$slength[i], blast$qlenth[i])))*100
  } else{
    overlap[i] <- 100
  }
}

blast$overlap <- round(overlap, digits = 3)
blast <- blast %>% 
  relocate(overlap, .before = pident)

# check for query genes that align to the same vNumber
any(duplicated(blast$sseqid)) # TRUE

# group blast hits by duplication, filter for duplicates only, then sort by 
# number of hits
blast.dup <- blast %>% 
  group_by(sseqid) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  filter(n >1)

# 50 vNumbers with duplicate hits, 2 vNumbers with triplicate hits. Filter by
# duplicates
blast.dup.only <- blast %>% 
  filter(sseqid %in% blast.dup$sseqid) %>% 
  arrange(sseqid)

# add the taxonomy for these duplicated genes
blast.dup.only$tax <- pa.bin.df[rownames(blast.dup.only),"spp"]
blast.dup.only <- blast.dup.only %>% 
  relocate(tax, .before = mismatch)

# sort duplicate hit list in descending order, first by overlap, then alignment
# length, and finally % identity, then take the top hit and use that to classify
# the vNumber's taxonomy

dup.tokeep <- data.frame(matrix(data = NA, ncol = ncol(blast.dup.only),
                                nrow = nrow(blast.dup)))

colnames(dup.tokeep) <- colnames(blast.dup.only)

for(i in 1:nrow(dup.tokeep)){
  dup.df <- blast.dup.only %>% 
    filter(sseqid == levels(factor(blast.dup.only$sseqid))[i]) %>% 
    arrange(desc(overlap),desc(length),desc(pident))
  
  # in the case of blast hits to the same vNumber where the taxonomy of the
  # pangenome gene is different (i.e. potentially NOT a species marker gene), 
  # only change the taxonomy to indicate a non-marker gene if the difference in
  # % identity and overlap are both <5 %, otherwise leave as the top blast hit
  if(all(identical(dup.df$tax[1],dup.df$tax[2])==FALSE,
         abs(dup.df$pident[1]-dup.df$pident[2])<5,
         abs(dup.df$overlap[1]-dup.df$overlap[2])<5)){
    dup.df$tax <- "Gardnerella"
  }
  
  dup.tokeep[i,] <- dup.df[1,]
}

# manually edit taxonomy for V1870020 as this is a triplicate result and the two
# 'top' results are both G.vag, while the third is G. pick and all are pretty 
# similar by length and pident
dup.tokeep[51,"tax"] <- "Gardnerella"

# sanity check for duplicates
any(duplicated(dup.tokeep$sseqid)) # returns FALSE

# add taxonomy to the blast hits and then filter results based on presence of 
# vNumber in blast.dup 
blast$tax <- pa.bin.df[rownames(blast),"spp"]

blast.uniq <- blast %>% 
  relocate(tax, .before = mismatch) %>% 
  filter(!sseqid %in% blast.dup$sseqid)

# sanity check for no duplicates and identical colnames to dup.tokeep
any(duplicated(blast.uniq$sseqid)) # returns FALSE
all(colnames(blast.uniq) == colnames(dup.tokeep)) # returns TRUE

# combine dfs containing unique Gard genes -> vNumber blast results
blast.nodups <- rbind(blast.uniq, dup.tokeep)

# count how many genes there are per species
blast.nodups %>% 
  group_by(tax) %>% 
  count()

# tax                          n    vs. duplicates not removed
#
# Gardnerella               2135                         (-41)
# Gardnerella_greenwoodii     61                          (-1)
# Gardnerella_leopoldii       93                          (-8)
# Gardnerella_pickettii      161                         (-14)
# Gardnerella_piotii          53                          (-7)
# Gardnerella_swidsinskii    234                         (-38)
# Gardnerella_vaginalis      411                         (-19)

# now filter for overlap >90 and pident >95 and count again
# NOTE: filters out ~190 genes
blast.nodups.filt <- blast.nodups %>% 
  filter(overlap >90 & pident >95)

blast.nodups.filt %>% 
  group_by(tax) %>% 
  count()

# tax                         n
#
# Gardnerella              2057
# Gardnerella_greenwoodii    52
# Gardnerella_leopoldii      78
# Gardnerella_pickettii     141
# Gardnerella_piotii         32
# Gardnerella_swidsinskii   204
# Gardnerella_vaginalis     395


# finally, replace all instances of 'Gardnerella_vaginalis' with 'Gardnerella'
# and replace the taxonomy of all Gardnerella genes according to the filtered,
# non-duplicate BLAST results
tax$V3 <- gsub("Gardnerella_vaginalis", "Gardnerella", tax$V3)

for(i in 1:nrow(blast.nodups.filt)){
  tax$V3 <- case_when(rownames(tax) == blast.nodups.filt$sseqid[i] ~ blast.nodups.filt$tax[i],
                      .default = tax$V3)
}

# this taxonomy table was used to colour the PCA plots of BV-only datasets for
# the plots in Figure 4B and Supplementary Figure S5. For code, see  the script
# 'gardnerella_pca_redo.R'
