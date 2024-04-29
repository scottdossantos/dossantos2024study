# processing Virginia metadata

#################################### setup ####################################

library(dplyr) # for data cleaning

# set path to github
user <- "sds/cc"
path.to.github <- switch(user,
                         "sds/cc" = "~/Documents/GitHub/dossantos2024study/",
                         "gg" = "~/Documents/0_git/projects/dossantos2024study/")

# NOTE: The files imported into R as part of this script are not provided on the
#       study's github repo as they contain data hosted on the NCBI Database of
#       Genotypes and Phenotypes. For transparency, we share the following code
#       describing how these data were used. The corresponding .Rda objects are
#       also intentionally omitted from the study's github repo for the same 
#       reason. The only information shown below are several column names for
#       multiple metadata files, and some possible values containing within the
#       columns, all of which are also shown on the MOMS-PI study dbGaP page, 
#       under the 'Variables' tab.

########################## merge: clinical & SRA data ##########################

# read in clinical sample metadata
clinical<-read.csv("[LOCAL-PATH-TO]/meta_1_clinicalPhenotype_phs001523.v1.pht007510.v1.p1.c1.phenotype_clinical.DS-PREG-COM-IRB-PUB-MDS.csv")

# convert 'Visit_ID' and 'vaginal_pH' to factor
clinical$visit_num <- as.numeric(gsub("Visit_", "", clinical$visit_num))
clinical$vaginal_pH <- gsub(0, "Unknown", clinical$vaginal_pH)

# make sure only values in 'reported_ga' are 1st/2nd/3rd trimester or NA
unique(clinical$reported_ga) # above is correct

# rename column and replace 1st/2nd/3rd with 1/2/3 as numerics
colnames(clinical)[6] <- "visit_trimester"
clinical$visit_trimester <- gsub("First Trimester", as.numeric(1), clinical$visit_trimester)
clinical$visit_trimester <- gsub("Second Trimester", as.numeric(2), clinical$visit_trimester)
clinical$visit_trimester <- gsub("Third Trimester", as.numeric(3), clinical$visit_trimester)
clinical$visit_trimester <- as.numeric(clinical$visit_trimester)

# select only relevant column and arrange by subject ID and then visit number)
# note: have to specify select from dplyr if 'MASS' package is loaded
clinical <- clinical %>% 
  dplyr::select(-c(category, pulse, systolic, diastolic, fundal_height)) %>% 
  arrange(dbGaP_Subject_ID, visit_num)

# read in SRA sample sheet, remove unnecessary columns and re-order columns
sra <- read.csv("[LOCAL-PATH-TO]/meta_7_SraRunTable.csv")
sra.filt <- sra %>% 
  dplyr::select(-c(analyte_type:biospecimen_repository,
            body_site:study_name)) %>% 
  relocate(Run, .after = submitted_subject_id)

# rename "submitted_subject_id" for consistency with 'clinical'
colnames(sra.filt)[2] <- "SUBJECT_ID"

# remove "_MV1R" from all rows of "biospecimen_repository_sample_id" and rename 
# this column for consistency with 'clinical'
sra.filt$biospecimen_repository_sample_id <- gsub("_MV1R", "", sra.filt$biospecimen_repository_sample_id)
colnames(sra.filt)[1] <- "visit_ID"

# merge clinical data and filtered sra information
merged.virginia.1 <- left_join(sra.filt, clinical, by= "visit_ID")

##################### merge: clinical/sra & health history #####################

# read in health history data (this table has 'NA' values as text, with several
# blank cells, so need to specify 'na.strings' argument)
health <- read.csv("[LOCAL-PATH-TO]/meta_5_healthHistory_phs001523.v1.pht007512.v1.p1.c1.phenotype_healthHistorySurvey.DS-PREG-COM-IRB-PUB-MDS.csv",
                   na.strings = c("", "NA"))

# remove unnecessary columns from health history survey data
health.filt <- health %>% 
  dplyr::select(-c(visit_num:age,
            education:infertility_treatment, 
            birthcontrol:high_bp_outside_pregNAcy,
            uti_history:yeast_infection_lifetime_number,
            current_medications:medications_past_6months,
            smoker_lifetime:diet_physician_prescribed))

# rename several columns to be more informative
colnames(health.filt)[which(colnames(health.filt)=="bv")] <- "prior_diagnosis_recurrent_bv"
colnames(health.filt)[which(colnames(health.filt)=="yeast")] <- "prior_diagnosis_recurrent_yeast"
colnames(health.filt)[which(colnames(health.filt)=="vdischarge")] <- "current_vag_abnorm_discharge"
colnames(health.filt)[which(colnames(health.filt)=="vodor")] <- "current_vag_bad_odor"
colnames(health.filt)[which(colnames(health.filt)=="vitching")] <- "current_vag_itching"

# check to see how many people have reported >1 ethnicity and subset into list
list.ethn <- list()
ethn.mult <- vector()
for(i in 1:nrow(health.filt)){
  row <- as.character(health.filt[i,c(4:9)])
  ethn <- which(row == "Yes")
  list.ethn[[i]] <- ethn
  ethn.mult[i] <- length(ethn)
}

list.ethn.mult <- list.ethn[which(ethn.mult >1)]

# set list element length to 3 to add NA valeus and enable collapsing to df
for(i in 1:length(list.ethn.mult)){
  length(list.ethn.mult[[i]]) <- 3
}

# collapse to data frame
list.ethn.mult.df <- as.data.frame(do.call(rbind,list.ethn.mult))

#add column of concatenated multiple ethnicities
ethn.mult.uniq <- vector()
for(i in 1:nrow(list.ethn.mult.df)){
  ethn.mult.uniq[i] <- paste(list.ethn.mult.df[i,], collapse = "-")
}
list.ethn.mult.df$concat <- ethn.mult.uniq


# count unique combinations
list.ethn.mult.df %>% 
  group_by(concat) %>% 
  count() %>% 
  arrange(desc(n))

#   concat      n
#   <chr>   <int>
#   1-4-NA      6     African-American / Caucasian
#   3-4-NA      6     Asian / Caucasian
#   4-5-NA      5     Caucasian / Hispanic or Latino
#   1-5-NA      3     African-American / Hispanic or Latino
#   1-2-NA      2     African-American / American Indian or Alaska Native
#   2-4-NA      2     American Indian or Alaska Native / Caucasian
#   1-2-4       1     African-American / American Indian or Alaska Native / Caucasian
#   1-4-5       1     African-American / Caucasian / Hispanic or Latino

# make vector of ethnicities in dataset
ethnicities <- c("African-American", "American Indian or Alaska Native",
                 "Asian", "Caucasian", "Hispanic or Latino", "Native Hawaiian")

# collapse six ethnicity columns into one vector
health.ethnicity <- vector()
for(i in 1:nrow(health.filt)){
  row <- as.character(health.filt[i,c(4:9)])
  ethn <- which(row == "Yes")
  if("Yes"  %in%  row == FALSE){
    ethn <- NA
  }
  health.ethnicity[i] <- case_when(length(ethn) == 1L && !is.na(ethn) ~ ethnicities[ethn],
                                   is.na(ethn) ~ as.character("Unknown"),
                                   length(ethn) >1L ~ paste(ethnicities[ethn],
                                                       collapse= " / "))[1]
}

# add vector to filtered health survey data, re-position it at the start of the
# old ethnicity columns and remove the old ethnicity columns

health.filt <- health.filt %>% 
  mutate(ethnicity= health.ethnicity) %>% 
  relocate(ethnicity, .before = african_american) %>% 
  dplyr::select(-c(african_american:native_hawaiian))

# pull indices of any duplicate patient IDs in initial health history survey
which(duplicated(health.filt$SUBJECT_ID)) # 40, 255, 538, 568

# one of these patients is in the 297 samples we care about: row 40 (row 130 in
# the sra/clinical merged data frame). Can merge on subject ID for getting 
# ethnicity data, but will lose info on discharge/odor/itching at visit 1
# (unless samples we care about ARE visit 1 samples --> fix later)
health.filt$SUBJECT_ID.x[which(duplicated(health.filt$SUBJECT_ID))] %in% merged.virginia.1$SUBJECT_ID.x

# rename subject ID in health history for consistency with sra/clinical data
colnames(health.filt)[2] <- "SUBJECT_ID.x"

# merge health survey data with sra/clinical data based on subject ID; this will
# throw two warnings: 
#  row 130 in sra/clinical matches multiple rows in health history
#   - extra row created in merged table as a result
#  row 202 in health history matches multiple rows in sra/clinical
#   - corresponding rows in merged table are identical (no extra rows created)
merged.virginia.2 <- left_join(merged.virginia.1, health.filt, by= "SUBJECT_ID.x")


# delete the extra row introduced in the merge (row 131 = EP081373_K40 ; most of
# the data in usable columns in the second row is NA anyway!)

merged.virginia.2 <- merged.virginia.2[-131,]

# sanity check to confirm that dbGaP subject ID is the same (unique to a PATIENT
# but shared between samples from the same patient): 296 = TRUE, EP175456 is not
# present in 'health.filt' and all corresponding rows are empty in the merged
# table

merged.virginia.2$dbGaP_Subject_ID.x == merged.virginia.2$dbGaP_Subject_ID.y

# remove visit.ID.y from sra/clinical/health data as it doesn't correspond to 
# the rows anymore, as well as the duplicate 'SUBJECT_ID.y' and 
# 'dbGaP_Subject_ID.y'

merged.virginia.2 <- merged.virginia.2 %>% 
  dplyr::select(-c(visit_ID.y, SUBJECT_ID.y, dbGaP_Subject_ID.y))

# clean up column names and convert 'vaginal_ph' to factor
colnames(merged.virginia.2)[1] <- "visit_id"
colnames(merged.virginia.2)[2] <- "subject_id"
colnames(merged.virginia.2)[3] <- "sra_id"
colnames(merged.virginia.2)[4] <- "dbgap_id"
colnames(merged.virginia.2)[5] <- "visit_number"
colnames(merged.virginia.2)[7] <- "vaginal_ph"

##################### merge: clinical/sra/health & reviews #####################

# read in review of systems data (this table has 'NA' values as text, with 
# several blank cells, so need to specify 'na.strings' argument)
review <- read.csv("[LOCAL-PATH-TO]/meta_4_reviewOfSystems_phs001523.v1.pht007511.v1.p1.c1.phenotype_reviewOfSystemsSurvey.DS-PREG-COM-IRB-PUB-MDS.csv",
                   na.strings = c("", "NA"))

# remove unnecessary columns from review of systems survey data
review.filt <- review %>% 
  dplyr::select(-c(visit_num:survey_ID, 
            bed_rest_since_last_visit, 
            cerclage:contractions_3rdtri, 
            diarrhea_1sttri:diet_vegetarian,
            ear_sinus_since_last_visit:oral_steroids_since_last_visit,
            plan_to_breastfeed:respiratory_3rdtri,
            vbleeding_1sttri:vbleeding_3rdtri,
            vdischarge_1sttri:vdischarge_3rdtri,
            vitching_1sttri:vitching_3rdtri,
            vodor_1sttri:vodor_3rdtri))

# none of this data needs processing: rename 'visit_ID' for consistency
colnames(review.filt)[3] <- "visit_id"

# merge sra/clinical/health history data with review data based on 'visit_id'
# (this throws no warnings; happy days)
merged.virginia.3 <- left_join(merged.virginia.2, review.filt, by= "visit_id")

# get number of samples (that we care about) which didn't have any info in 
# 'review.filt': 72 samples

length(which(is.na(merged.virginia.3$SUBJECT_ID)))

# get breakdown of samples with no info in 'review.filt' (most missing samples
# should be visit 1)
merged.virginia.3[which(is.na(merged.virginia.3$SUBJECT_ID)),] %>% 
  group_by(visit_number) %>% 
  count()

# A tibble: 5 × 2
# Groups:   visit_number [5]
# visit_number     n
#  <dbl> <int>
#            1    56
#            2     4
#            3     7
#            4     3
#            6     2

# fill in 'vdischarge', 'vitching', 'vodor', 'doctor_yeast_meds', and 
# 'otc_yeast_meds' info for visit_1 samples which have no entries in 
# 'review.filt' (this info is found in 'current_vag_' and '_yeast_meds' fields
# in 'merged.virginia.3')

merged.virginia.3$vdischarge[which(merged.virginia.3$visit_number==1 & is.na(merged.virginia.3$SUBJECT_ID))] <-
  merged.virginia.3$current_vag_abnorm_discharge[which(merged.virginia.3$visit_number==1 & is.na(merged.virginia.3$SUBJECT_ID))]

merged.virginia.3$vitching[which(merged.virginia.3$visit_number==1 & is.na(merged.virginia.3$SUBJECT_ID))] <-
  merged.virginia.3$current_vag_itching[which(merged.virginia.3$visit_number==1 & is.na(merged.virginia.3$SUBJECT_ID))]

merged.virginia.3$vodor[which(merged.virginia.3$visit_number==1 & is.na(merged.virginia.3$SUBJECT_ID))] <-
  merged.virginia.3$current_vag_bad_odor[which(merged.virginia.3$visit_number==1 & is.na(merged.virginia.3$SUBJECT_ID))]

merged.virginia.3$doctor_yeast_meds_since_last_visit[which(merged.virginia.3$visit_number==1 & is.na(merged.virginia.3$SUBJECT_ID))] <-
  merged.virginia.3$doctor_yeast_meds[which(merged.virginia.3$visit_number==1 & is.na(merged.virginia.3$SUBJECT_ID))]

merged.virginia.3$otc_yeast_meds_since_last_visit[which(merged.virginia.3$visit_number==1 & is.na(merged.virginia.3$SUBJECT_ID))] <-
  merged.virginia.3$otc_yeast_meds[which(merged.virginia.3$visit_number==1 & is.na(merged.virginia.3$SUBJECT_ID))]

# remove 'dbGaP_Subject_ID' and 'SUBJECT_ID fields (many are blank), and the
# 'current_vag_' and '_yeast_meds' fields (these only applied to visit 1 samples
# and we have merged these entries with the 'vdischarge', 'vitching', and 
# 'vodor' fields from 'review.filt')
merged.virginia.3 <- merged.virginia.3 %>% 
  dplyr::select(-c(dbGaP_Subject_ID,
            SUBJECT_ID, 
            current_vag_abnorm_discharge,
            current_vag_bad_odor,
            current_vag_itching,
            doctor_yeast_meds,
            otc_yeast_meds))

# rename current vaginal symtoms columns from 'review.filt'
colnames(merged.virginia.3)[which(colnames(merged.virginia.3)=="vdischarge")] <- "current_vag_abnorm_discharge"
colnames(merged.virginia.3)[which(colnames(merged.virginia.3)=="vodor")] <- "current_vag_bad_odor"
colnames(merged.virginia.3)[which(colnames(merged.virginia.3)=="vitching")] <- "current_vag_itching"

# re-position 'current_vag_' symptom columns to be just after 'vaginal_ph'
merged.virginia.3 <- merged.virginia.3 %>%
  relocate(current_vag_itching, .after = current_vag_bad_odor) %>% 
  relocate(current_vag_abnorm_discharge:current_vag_itching,
           .after = vaginal_ph)

# 'current_vag_abnorm_discharge' and 'current_vag_bad_odor' have "Not_sure" and
# "Not_Sure" responses, so gsub these to be all "Not_Sure" for consistency with
# other columns

merged.virginia.3$current_vag_abnorm_discharge <- gsub("Not_sure", "Not_Sure",
                                                       merged.virginia.3$current_vag_abnorm_discharge)

merged.virginia.3$current_vag_bad_odor <- gsub("Not_sure", "Not_Sure",
                                                       merged.virginia.3$current_vag_bad_odor)

# 'current_vag_itching' and 'diabetes' has "Not_sure" and responses, so gsub
# these to be all "Not_Sure" for consistency with other columns

merged.virginia.3$current_vag_itching <- gsub("Not_sure", "Not_Sure",
                                              merged.virginia.3$current_vag_itching)

merged.virginia.3$diabetes <- gsub("Not_sure", "Not_Sure",
                                   merged.virginia.3$diabetes)

# 'doctor/otc_yeast_meds_since_last_visit' have "Uncertain" responses so gsub 
# these to be all "Not_Sure" for consistency with other columns

merged.virginia.3$doctor_yeast_meds_since_last_visit <- gsub("Uncertain",
                                                             "Not_Sure",
                                                             merged.virginia.3$doctor_yeast_meds_since_last_visit)

merged.virginia.3$otc_yeast_meds_since_last_visit <- gsub("Uncertain",
                                                             "Not_Sure",
                                                             merged.virginia.3$otc_yeast_meds_since_last_visit)

# rename merged metadata file and save as .Rda object
# virginia.meta <- merged.virginia.3
# save(virginia.meta, file = "[LOCAL-PATH-TO]/virginia.meta.Rda")
