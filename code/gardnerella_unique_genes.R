########## Identifying genes unique to each named Gardnerella species ##########

# This script is run as part of the process to identify genes unique to each
# gardnerella species, using two publicly available genomes per species. For
# detailed info about the genomes used and their processing, see:
# /[GitHubPath]/dossantos2024study/code/gardnerella_unique_genes_what_i_did.txt

#################################### setup ####################################

# load dplyr
library(dplyr)

# set path to github directory
user <- "sds/cc"
if(user=="sds/cc"){
  path.to.github<-"~/Documents/GitHub/dossantos2024study/"
}else if(user=="gg"){
  path.to.github<-"~/Documents/0_git/projects/dossantos2024study/"
}

# load panaroo output (see 'gardnerella_unique_genes_what_i_did.txt' for info)
gene.pres.abs<-read.csv(paste(path.to.github,
                              "Rdata/gene_presence_absence.csv", sep = ""))

################################ data clean-up ################################

# rename columns for clarity
new.col.names<-c("leopoldii_type","leopoldii_oth",
                 "piotii_oth","piotii_type",
                 "swidsinskii_type","swidsinskii_oth",
                 "vaginalis_type","vaginalis_oth")

colnames(gene.pres.abs)[4:11]<-new.col.names

# move piotii type strain column to the left for consistency
gene.pres.abs<-gene.pres.abs%>%
  relocate(piotii_oth,
           .after = piotii_type)

# copy of table for editing (paranoia)
bin.gene.pres.abs<-gene.pres.abs

# replace internal prokka gene locus tags with binary pres/abs (1/0) values
for (i in c(4:11)) {
  vec<-bin.gene.pres.abs[,i]
  vec.one<-replace(vec,vec!="",1)
  vec.one.zero<-replace(vec.one,vec.one=="",0)
  bin.gene.pres.abs[,i]<-vec.one.zero
}

############################## pull unique genes ##############################

# get genes unique to leopoldii
unique.leopoldii<-bin.gene.pres.abs%>%
  filter(leopoldii_type==1,
         leopoldii_oth==1,
         piotii_type==0,
         piotii_oth==0,
         swidsinskii_type==0,
         swidsinskii_oth==0,
         vaginalis_type==0,
         vaginalis_oth==0)

# get genes unique to piotii
unique.piotii<-bin.gene.pres.abs%>%
  filter(leopoldii_type==0,
         leopoldii_oth==0,
         piotii_type==1,
         piotii_oth==1,
         swidsinskii_type==0,
         swidsinskii_oth==0,
         vaginalis_type==0,
         vaginalis_oth==0)

# get genes unique to swidsinskii
unique.swidsinskii<-bin.gene.pres.abs%>%
  filter(leopoldii_type==0,
         leopoldii_oth==0,
         piotii_type==0,
         piotii_oth==0,
         swidsinskii_type==1,
         swidsinskii_oth==1,
         vaginalis_type==0,
         vaginalis_oth==0)

# get genes unique to vaginalis
unique.vaginalis<-bin.gene.pres.abs%>%
  filter(leopoldii_type==0,
         leopoldii_oth==0,
         piotii_type==0,
         piotii_oth==0,
         swidsinskii_type==0,
         swidsinskii_oth==0,
         vaginalis_type==1,
         vaginalis_oth==1)

# write text file for each species --> to be used for BLAST/mapping
write(unique.leopoldii$Gene,"unique_genes_g_leopoldii.txt")
write(unique.piotii$Gene,"unique_genes_g_piotii.txt")
write(unique.swidsinskii$Gene,"unique_genes_g_swidsinskii.txt")
write(unique.vaginalis$Gene,"unique_genes_g_vaginalis.txt")

# # get max vector length for cbind()
# maxLength<-max(length(headers.leopoldii),
#                length(headers.piotii),
#                length(headers.swidsinskii),
#                length(headers.vaginalis))
# 
# # adjust max vector length for cbind()
# length(headers.leopoldii)<-231
# length(headers.piotii)<-231
# length(headers.swidsinskii)<-231
# 
# # combine into one dataframe
# unique.all<-as.data.frame(cbind(headers.leopoldii,
#                   headers.piotii,
#                   headers.swidsinskii,
#                   headers.vaginalis))
# 
# # write .csv files without row names or quotes around character strings
# write.csv(unique.all,
#           "unique_genes_gardnerella_spp.csv",
#           row.names = FALSE,
#           quote = FALSE)
