# reading, processing and summarising virginia FASTQ QC input

# This script contains code to create a summary of the mapping steps, showing 
# the number of reads each sample contained at each step. Because there were
# issues with running the bash script containing commands for mapping, I have
# documented EXACTLY what happened 1) for my own benefit,  to remember what I
# did in case I ever have to go back to the files on the server, and 2) so that
# anyone reading this can have a right old laugh at how daft I am while exhaling
# through their nose and muttering "idiot". It's highly unlikely that anyone 
# will ever actually read this, but just in case they do, I hope the reader
# is having a nice day. Now please enjoy a verbose retelling of what I did and 
# how dumb I am - Scott Dos Santos, 22nd April 2024.

############################# mapping explanation #############################

# So here's what happened: the first time I ran the mapping script, all was 
# going well, but I stopped the process at SRR6743975_1 because I was zipping
# each FASTQ file as part of the extraction loop (which was taking ages) and
# using .gz files in the master script would have been spending too long zipping
# and unzipping each file. By the time I realised this, the extraction script
# had ~ 20 samples left to go and ~20 samples had been mapped. After stopping
# the mapping script, I also stopped the extraction script and edited the 
# extraction script to do fasterq-dump and remove reverse read files but NOT zip 
# the remaining forward read file. I then restarted the extraction script.

# I then realised I would have to re-make the symbolic links in the virginia
# directory as the input files were now '.fastq' rather than '.fastq.gz', so
# I made new symbolic links which ended in 'fastq.gz' but acutally pointed to 
# '.fastq' files. I thought I was being clever, but this came back to bite me in
# the arse. Thinking all would be well, I edited the 'basenames_virginia' file
# to remove all samples prior to SRR6743975_1, waited for the extraction script
# to finish and re-started the mapping, outputting to 'output_take2_...txt'.

# All WAS well until I hit sample SRR6744094_1, which was the first sample to
# get stuck in a seemingly-infinite loop at the end of the VIRGO step 1 script.
# I stopped the mapping and tried to restart it, outputting to 
# 'output_take3_...txt'. This got stuck too, leaving me with only one set of
# output text in 'output_take3_mapping_virginia.txt'; for a sample that failed,
# too! After trying and failing to run the script manually on the rRNA-mapped
# FASTQ file for this sample, I decided to skip it and try a smaller file next,
# 'SRR6744101_1'. This was because the next sample in line, 'SRR6744095_1', was
# of a similar size to 'SRR6744094_1', and I wondered if size played into why 
# the loop got stuck (it didn't, btw). This is why SRR6744095_1 comes after
# SRR6744101_1 in the output of this R script below. Again, I thought all would
# be well and restarted the mapping, outputting to 'output_take4_...txt'.

# It was not. There were a further six instances where the loop got stuck and
# I summarily re-started the mapping, skipping over the failed samples each 
# time. These samples were: SRR6744094_1, SRR6744448_1, SRR6744460_1, 
# SRR6744507_1, SRR6744553_1, SRR6744564_1, and SRR6744918_1. After restarting
# the script following the most recent failure, I was outputting to 
# 'output_take10_...txt'. At one point, I can't remember which, I upped the
# threads to 90 on Greg's advice.

# During the loop which output to 'output_take10_...txt', I was bitten in the
# arse by my own stupidity when I renamed those symbolic links earlier. When
# the loop encountered the first of these symbolic links, which was named
# 'SRR6744928_1.fastq.gz' but actually pointed to 'SRR6744928_1.fastq',
# Trimmomatic's automatic input format detection became unhappy, quite correctly
# reporting that the input file was NOT in GZIP format. This caused Trimmomatic,
# and every command after it, to fail.

# When I discovered this, I renamed the offending symbolic links such that they
# now correctly ended in '.fastq' and then restarted the mapping, outputting to
# 'output_take11_...txt'. However, because I am not smart, I forgot to alter the
# last line of the script where I remove all files ending in '.fastq'. No, 
# really. That was an actual thing I did. Accordingly, the std.out log in
# 'output_take11_...txt' shows that the first sample, SRR6744928_1, completes
# all the steps properly, then, as per the script, everything ending in '.fastq'
# is deleted, including the new symbolic links. All commands after this point
# fail because the corresponding symbolic links are gone.

# I realised this the next morning and edited the rm command to instead be:
# rm $i'.unalign.*', which was daft, because the shell will interpret the 
# wildcard character literally when it is enclosed in single quotes, therefore
# the unaligned read files would not be deleted, instead, throwing an error
# saying so. However, this didn't actually affect the loop and all of the 
# remaining files completed the mapping successfully.

# Because I can't leave well enough alone, I also wanted to see if I could
# get the 'failed' samples done by running virgo step 1 on my local machine.
# Spoilers: I can. I copied and edited the mapping script so that it only runs
# on the failed symbolic links and properly deletes the unneeded FASTQ files, 
# but stops after mapping to SILVA. As each file completed, I scp'd them to
# my local machine and ran step 1 of VIRGO. It worked, and .out files are being
# produced. I will later scp these to agrajag and re-run step 2 of virgo.

# So, now that you're all caught up, here's what I have done to automate the 
# counting of the read numbers at each step:
#    - scp'd all of the output files to my local machine
#    - told this R script to remove 'output_take3_...txt' from the data frame of
#      files to read in
#    - opened 'output_take10_...txt' and 'output_take11_...txt' in TextEdit to
#      manually remove all output after the last successful sample, ending in
#      'Reported x alignments', as expected
#    - opened 'output_take12_...txt' in TextEdit and deleted the last line which
#      contained the error about not being able to 'rm' all '.unalign' files
#    - opened 'output_mapping_virginia.txt' and deleted the information for
#      SRR6743975_1 so it wouldn't be duplicated in this script
#    - opened 'output_failed_mapping.txt' and removed lines starting with
#      'Setting the index via positional argument will....' as this throws
#      off the index for 'failed' samples by 1
#    - Wrote this R script to parse all output files and calculate the reads
#      at each step for the 290 virginia samples which didn't get stuck in the
#      loop

#################################### setup #####################################

library(dplyr) # dataframe manipulation
library(stringr) # string manipulation

# set path to github repository
user <- "sds/cc"
path.to.github <- switch(user,
                         "sds/cc" = "~/Documents/GitHub/metatranscriptome/",
                         "gg" = "~/Documents/0_git/projects/metatranscriptome/")

# get list of files in directory
output.files <- list.files(path = paste(path.to.github,
                                        "Rdata/virginia_output/",
                                        sep = ""),
                           pattern = "output")

# add a vector of the order the files should be read
output.files <- as.data.frame(cbind(files=output.files, order= c(13,1,10,11,12,2:9)))
output.files$order <- as.numeric(output.files$order)

# rearrange based on 'order' and remove 'output_take3_mapping_virginia.txt',
# which contains info for only 1 failed sample
output.files <- output.files %>%
  arrange(order) %>% 
  slice(1:2, 4:13)

# total number of samples is 297; generate vectors of this length for sample id
# and all read values
vec.1.sample.id <- vector(length = 297)
vec.2.raw <- vector(length = 297)
vec.3.trim <- vector(length = 297)
vec.4.hg38 <- vector(length = 297)
vec.5.t2t <- vector(length = 297)
vec.6.rrna <- vector(length = 297)
vec.7.virgo <- vector(length = 297)

############################ count virginia reads #############################

# initialise main loop
for(file in 1:nrow(output.files)){
  # read in output text file, making sure to disable interpretation of comment
  # characters as some read info lines begin with '#'
  std.out<-read.table(paste(path.to.github, "Rdata/virginia_output/",
                            output.files$files[file], sep = ""),
                      header = FALSE, sep = "\t", quote = "",
                      check.names = FALSE, comment.char = "")

  # get number of samples processed in this file
  num.samples <- length(grep("TrimmomaticSE: Started with arguments:", std.out$V1))
  
  # check for output file that completed normally
  if("Reported"  %in% str_split_i(std.out$V1[nrow(std.out)], " ", 1)){
    
    # initialise sub-loop for all samples in output file
    for(samp in 1:num.samples){
    
      # get index of SRR sample ID row
      index.1.sample.id <- grep("-threads", std.out$V1)[samp]

      # get sample ID and remove '.fastq.gz' from it (all 'SRR' IDs are same length)
      item.1.sample.id <- str_split_i(std.out$V1[index.1.sample.id], " ", 4)
      item.1.sample.id <- substr(item.1.sample.id,start = 1, stop = 12)

      # get index of raw and trimmed reads row 
      index.2.reads.raw.trim<- grep("Input", std.out$V1)[samp]

      # get number of raw reads and number of trimmed reads
      item.2.reads.raw <- as.integer(str_split_i(std.out$V1[index.2.reads.raw.trim], " ", 3))
      item.3.reads.trim <- as.integer(str_split_i(std.out$V1[index.2.reads.raw.trim], " ", 5))

      # get the index of the "TrimmomaticSE: Completed successfully" row and use this
      # as a base index to calculate the indices of each mapping step
      index.3.mapping.base <- grep("TrimmomaticSE: Completed successfully",std.out$V1)[samp]

      # calculate indices of HG38-mapped, T2T-mapped, SILVA-mapped and VIRGO-mapped 
      # read information
      index.4.mapping.hg38 <- as.integer(index.3.mapping.base+3)
      index.5.mapping.t2t <- as.integer(index.3.mapping.base+9)
      index.6.mapping.rrna <- as.integer(index.3.mapping.base+15)
      index.7.mapping.virgo <- as.integer(index.3.mapping.base+24)

      # get number of reads that don't align to HG38, T2T and SILVA, as well as the
      # number of reads that map to VIRGO
      item.4.mapping.hg38 <- as.integer(str_split_i(std.out$V1[index.4.mapping.hg38], " ", 5))
      item.5.mapping.t2t <- as.integer(str_split_i(std.out$V1[index.5.mapping.t2t], " ", 5))
      item.6.mapping.rrna <- as.integer(str_split_i(std.out$V1[index.6.mapping.rrna], " ", 5))
      item.7.mapping.virgo <- as.integer(str_split_i(std.out$V1[index.7.mapping.virgo], " ", 9))
      
      # get index of the position to add the items
      position.chr <- as.integer(min(which(vec.1.sample.id==FALSE)))
      position.int <- as.integer(min(which(vec.7.virgo==0)))
      
      # add all items as the position'th element of reads.x objects
      if(is.logical(vec.7.virgo)){
        vec.1.sample.id[position.chr] <- item.1.sample.id
        vec.2.raw[position.chr] <- item.2.reads.raw
        vec.3.trim[position.chr] <- item.3.reads.trim
        vec.4.hg38[position.chr] <- item.4.mapping.hg38
        vec.5.t2t[position.chr] <- item.5.mapping.t2t
        vec.6.rrna[position.chr] <- item.6.mapping.rrna
        vec.7.virgo[position.chr] <- item.7.mapping.virgo
      }else if(is.integer(vec.2.raw)){
        vec.1.sample.id[position.chr] <- item.1.sample.id
        vec.2.raw[position.int] <- item.2.reads.raw
        vec.3.trim[position.int] <- item.3.reads.trim
        vec.4.hg38[position.int] <- item.4.mapping.hg38
        vec.5.t2t[position.int] <- item.5.mapping.t2t
        vec.6.rrna[position.int] <- item.6.mapping.rrna
        vec.7.virgo[position.int] <- item.7.mapping.virgo  
        }
      }
    }else if(std.out$V1[nrow(std.out)] == "VIRGO path   = /Volumes/data/twntyfr_2018/VIRGO/"){
    
      for(samp in 1:(num.samples-1)){
    
        # get index of SRR sample ID row
        index.1.sample.id <- grep("-threads", std.out$V1)[samp]
    
        # get sample ID and remove '.fastq.gz' from it (all 'SRR' IDs are same length)
        item.1.sample.id <- str_split_i(std.out$V1[index.1.sample.id], " ", 4)
        item.1.sample.id <- substr(item.1.sample.id,start = 1, stop = 12)
    
        # get index of raw and trimmed reads row 
        index.2.reads.raw.trim<- grep("Input", std.out$V1)[samp]
    
        # get number of raw reads and number of trimmed reads
        item.2.reads.raw <- as.integer(str_split_i(std.out$V1[index.2.reads.raw.trim], " ", 3))
        item.3.reads.trim <- as.integer(str_split_i(std.out$V1[index.2.reads.raw.trim], " ", 5))
    
        # get the index of the "TrimmomaticSE: Completed successfully" row and use this
        # as a base index to calculate the indices of each mapping step
        index.3.mapping.base <- grep("TrimmomaticSE: Completed successfully",std.out$V1)[samp]
    
        # calculate indices of HG38-mapped, T2T-mapped, SILVA-mapped and VIRGO-mapped 
        # read information
        index.4.mapping.hg38 <- as.integer(index.3.mapping.base+3)
        index.5.mapping.t2t <- as.integer(index.3.mapping.base+9)
        index.6.mapping.rrna <- as.integer(index.3.mapping.base+15)
        index.7.mapping.virgo <- as.integer(index.3.mapping.base+24)
    
        # get number of reads that don't align to HG38, T2T and SILVA, as well as the
        # number of reads that map to VIRGO
        item.4.mapping.hg38 <- as.integer(str_split_i(std.out$V1[index.4.mapping.hg38], " ", 5))
        item.5.mapping.t2t <- as.integer(str_split_i(std.out$V1[index.5.mapping.t2t], " ", 5))
        item.6.mapping.rrna <- as.integer(str_split_i(std.out$V1[index.6.mapping.rrna], " ", 5))
        item.7.mapping.virgo <- as.integer(str_split_i(std.out$V1[index.7.mapping.virgo], " ", 9))
        
        # get index of the position to add the items
        position.chr <- as.integer(min(which(vec.1.sample.id==FALSE)))
        position.int <- as.integer(min(which(vec.7.virgo==0)))
        
        # add all items as the position'th element of reads.x objects
        if(is.logical(vec.7.virgo)){
          vec.1.sample.id[position.chr] <- item.1.sample.id
          vec.2.raw[position.chr] <- item.2.reads.raw
          vec.3.trim[position.chr] <- item.3.reads.trim
          vec.4.hg38[position.chr] <- item.4.mapping.hg38
          vec.5.t2t[position.chr] <- item.5.mapping.t2t
          vec.6.rrna[position.chr] <- item.6.mapping.rrna
          vec.7.virgo[position.chr] <- item.7.mapping.virgo
        }else if(is.integer(vec.2.raw)){
          vec.1.sample.id[position.chr] <- item.1.sample.id
          vec.2.raw[position.int] <- item.2.reads.raw
          vec.3.trim[position.int] <- item.3.reads.trim
          vec.4.hg38[position.int] <- item.4.mapping.hg38
          vec.5.t2t[position.int] <- item.5.mapping.t2t
          vec.6.rrna[position.int] <- item.6.mapping.rrna
          vec.7.virgo[position.int] <- item.7.mapping.virgo  
        }
      }
    }
}

final.read.summary <- data.frame(sample.id = vec.1.sample.id,
                                 nreads.raw = as.numeric(vec.2.raw),
                                 reads.trim = as.numeric(vec.3.trim),
                                 reads.hg38 = as.numeric(vec.4.hg38),
                                 reads.t2t = as.numeric(vec.5.t2t),
                                 reads.rRNA = as.numeric(vec.6.rrna),
                                 reads.virgo = as.numeric(vec.7.virgo))

final.read.summary <- final.read.summary %>% 
  arrange(sample.id) %>% 
  tibble::column_to_rownames("sample.id")

final.mean <- apply(final.read.summary, 2, mean)
final.median <- apply(final.read.summary, 2, median)

final.summary <- data.frame(mean = final.mean, 
                            median = final.median,
                            row.names = colnames(final.read.summary))

# write.table(final.read.summary,
#             file = paste(path.to.github, 
#                          "Rdata/t2t_read_counts_summary_Virginia.txt",
#                          sep= ""),
#             quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

############################ compare europe reads #############################

# read in summary reads for europe dataset (calculated during mapping steps for
# europe data)

europe.summary <- read.table(paste(path.to.github,
                                   "Rdata/t2t_read_counts_summary_London_Europe.txt",
                                   sep = ""),
                             header = TRUE, sep = "\t", quote = "", row.names = 1)

# retain only rows corresponding to london/europe samples which WERE used
# (sample names taken from 'merge_feature_tables.R')
names.london.bv <- c("001A", "003A", "006A", "008A", "009A", "010A", "012B",
                     "013B", "014B", "018B", "017B", "012A", "013A", "019A")

names.london.h <- c("001B", "002B", "004B", "006B", "009B", "010B", "015B", "020B")

names.europe.bv <- paste(c("ERR2014354", "ERR2014359", "ERR2014361", "ERR2014363", 
                           "ERR2014365", "ERR2014370", "ERR2014372", "ERR2014374", 
                           "ERR2014376", "ERR2014378", "ERR2014383", "ERR2014388", 
                           "ERR2014390", "ERR2014392"), "_1", sep = "")

names.europe.h <- paste(c("ERR2014355", "ERR2014362", "ERR2014366", "ERR2014373", 
                          "ERR2014377", "ERR2014379", "ERR2014384", "ERR2014389"),
                        "_1", sep = "")

europe.summary <- europe.summary %>% 
  filter(rownames(europe.summary) %in% c(names.london.bv, names.london.h,
                                         names.europe.bv, names.europe.h))

# calculate median raw reads for london and europe datasets
median(europe.summary$reads_start[1:22]) # london: 14,479,786 reads
median(europe.summary$reads_start[23:44]) # europe: 16,584,536 reads

# calculate median reads post-VIRGO mapping for london and europe datasets
median(europe.summary$reads_VIRGO_step1[1:22]) # london: 11,107,052 reads
median(europe.summary$reads_VIRGO_step1[23:44]) # europe: 8,096,724 reads
