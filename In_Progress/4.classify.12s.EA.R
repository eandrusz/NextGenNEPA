## NGN Classifying ASVs: MiFish/MiMammal all runs after run 1 
# Author: Eily Allan - modified from Erin D'Agnese and Ramon Gallego 
# Person running: Eily
# Last modified: 3/3/22 by Eily
# Date of run: 3/3/22 by Eily 

# Overview 
# The general overview is that ASVs from dada2 will be read in, we will find hashes that were already previously annotated from other runs 
#and pull them out (because we do not need to re-assign taxonomy). For ASVs that have not yet been classified, we classify using insect 
#(found here - https://cran.r-project.org/web/packages/insect/vignettes/insect-vignette.html). 

#Then we will add the newly classified hashes via insect to the previous effort to keep growing our list of classified things.


###### NOTE you will run this twice per run - MiFish and MiMammal

####################################################################
# Set up
####################################################################

## Load packages
library(tidyverse)
library(insect)
library(seqinr)
library(here)
library(taxonomizr)

## REQUIRES HUMAN INPUT 
#marker <- "MiFish"
marker <- "MiMammal"
run.num <- "6"

# manually change file path because of the date in the folder for the dada2 output 
Hash_key <- here("Output", "dada2_output", "run6_20220127", marker, "hash_key.csv")
ASVs <- here("Output","dada2_output","run6_20220127",marker,"ASV_table.csv")

## REQUIRES HUMAN INPUT AND HARD CODE BECAUSE MUST BE STORED LOCALLY
# note mifish and mimammal use the same classifier right now
classifier <- "/Users/elizabethandruszkiewicz/GoogleDrive/UW/GitHub/NextGenNEPA_LOCAL/Input/classifiers/classifier_12S_v1.rds"

## Don't need to change 
previous_good_effort <- paste0(here("Output","classification_output"),"/",marker,".all.good.previous.hashes.annotated.rds")
previous_effort<- paste0(here("Output","classification_output"),"/",marker,".all.previous.hashes.annotated.rds")

run_output_folder <- paste0(here("Output","classification_output"),"/run",run.num) 
dir.create(path = run_output_folder)
run_marker_output_folder <- paste0(run_output_folder,"/",marker)
dir.create(path = run_marker_output_folder)

## Read in files 
Hash_key <- read_csv(Hash_key)
Hash <- Hash_key %>% 
  select(Hash, Sequence) %>% 
  distinct()
ALL.ASVs <- read_csv(ASVs)
tree <- read_rds(classifier)
previous.good.effort <- read_rds(previous_good_effort)
previous.effort <- read_rds(previous_effort)


####################################################################
# Classify
## Note that the classifier is maintained by insect. We are using version 1, which was last updated on 20181111. 
## Check here to see if there is a newer version: https://cran.r-project.org/web/packages/insect/vignettes/insect-vignette.html
####################################################################

## Prepare and run insect

## FOR 12S REMOVE ANYTHING TOO BIG (aka bacteria)
hash.length <- nchar(Hash$Sequence)
hash.keep <- hash.length < 200
Hash <- Hash[hash.keep,]

## here we are going to separate out the ASVs that are NOT in our previous effort file (i.e., only classify things that did not get a good classification previously)
# we probably only want to use the good classifications from previous - but you could change this to params$previous_effort
new.set <- anti_join(Hash, previous.good.effort, by = c("Hash" = "representative")) # remove anything previously classified
new.hashes.insect <- char2dna(new.set$Sequence)
names(new.hashes.insect) <- new.set$Hash

# just check it out 
new.hashes.insect

# do the classification! this takes a while
clasif.hashes <- classify (x = new.hashes.insect, tree = tree, cores = 4)

# rename columns to be useful
names(clasif.hashes) <- c('representative', 'taxID', 'taxon', 'rank', 'score', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')

## Check them out

# look at family/genus/species classifications
clasif.hashes %>% 
  unite (family, genus, species, sep = "|", col = "taxa")

# see how many were not assigned a rank
clasif.hashes %>% dplyr::count (rank) %>% arrange(desc(n))

# How many have a valid family but no phylum info
clasif.hashes %>% 
  filter(family!= "" & phylum == "") %>% 
  distinct(class) 

# Add new phylum info
clasif.hashes %>% 
  mutate(phylum = case_when(phylum != "" ~ phylum,
                            TRUE   ~ class))


####################################################################
# Save all new results
## Note that we will save two versions - all classifications and just the "good" ones - both as RDS and fasta
####################################################################

# save it in a subfolder of Output/classification_output/ and the run - so we know what was classified only from run 1
# that is our "run_marker_output_folder" file path that we defined above
# here we will save it as both an rds and a csv 
saveRDS(clasif.hashes, paste0(run_marker_output_folder,"/new.hashes.annotated.rds"))
clasif.hashes <- readRDS(paste0(run_marker_output_folder,"/new.hashes.annotated.rds"))
write.csv(clasif.hashes, paste0(run_marker_output_folder,"/new.hashes.annotated.csv"))

# also write them as a tax table in the same place
source(here("functions", "tax.table.R"))
taxtable <- tax.table(clasif.hashes)
write.csv(taxtable,file=paste0(run_marker_output_folder,"/tax.table.csv"))

####################################################################
# Save only good new results
## For MiFish/MiMammal this is to species level
####################################################################

# but let's only save the "good ones" (SPECIES for 12S) to use as a previous effort 
good.clasif.hashes <- clasif.hashes %>% filter(rank == "species")

# first, let's save it in a subfolder of Output/classification_output/ and the run - so we know what was classified only from run 1
# that is our "run_marker_output_folder" file path that we defined above
# here we will save it as both an rds and a csv 
saveRDS(good.clasif.hashes, paste0(run_marker_output_folder,"/new.good.hashes.annotated.rds"))
clasif.hashes <- readRDS(paste0(run_marker_output_folder,"/new.good.hashes.annotated.rds"))
write.csv(clasif.hashes, paste0(run_marker_output_folder,"/new.good.hashes.annotated.csv"))

# also write them as a tax table in the same place
source(here("functions", "tax.table.R"))
taxtable <- tax.table(good.clasif.hashes)
write.csv(taxtable,file=paste0(run_marker_output_folder,"/good.tax.table.csv"))


####################################################################
# Combine new results with results from previous runs and re-write the previous effor files
# first all classifications
####################################################################

# make and save combined rds file with the previous effort and the new one 
combined <- rbind(previous.effort, clasif.hashes)
combined.taxtable <- tax.table(combined)

# we want to put the combined / new previous hashes in the general Output/classification_output folder that we can re-write over and over again as we add new classifications each time that we add another run 

# so here we are going to change the file path - and the file name to say "all.previous.hashes" - and let's also add the date on the end of the file name so we know when it was last updated - this will be slightly annoying later when we have to change the file name when we read in the classifications, but probably worth it so we keep track of things
# save it again as an rds and then also as a csv - and tax table
saveRDS(combined, file=paste0(here("Output","classification_output"),"/", marker,".all.previous.hashes.annotated.rds"))
combined <- readRDS(file=paste0(here("Output","classification_output"),"/", marker,".all.previous.hashes.annotated.rds"))
write.csv(combined, file=paste0(here("Output","classification_output"),"/", marker,".all.previous.hashes.annotated.csv"))
write.csv(combined.taxtable,file=paste0(here("Output","classification_output"),"/", ".", marker,"combined.tax.table.csv"))

####################################################################
# Combine new results with results from previous runs and re-write the previous effor files
# then just the good classifications
####################################################################

combined.good <- rbind(previous.good.effort,good.clasif.hashes)
combined.good.taxtable <- tax.table(combined.good)

# we want to put the combined / new previous hashes in the general Output/classification_output folder that we can re-write over and over again as we add new classifications each time that we add another run 

# so here we are going to change the file path - and the file name to say "all.previous.hashes" 
# save it again as an rds and then also as a csv - and tax table
saveRDS(combined.good, file=paste0(here("Output","classification_output"),"/", marker,".all.good.previous.hashes.annotated.rds"))
combined.good <- readRDS(file=paste0(here("Output","classification_output"),"/", marker,".all.good.previous.hashes.annotated.rds"))
write.csv(combined.good, file=paste0(here("Output","classification_output"),"/", marker,".all.good.previous.hashes.annotated.csv"))
write.csv(combined.good.taxtable,file=paste0(here("Output","classification_output"),"/", marker,"combined.good.tax.table.csv"))

