## NGN Classifying ASVs: MiFish/MiMammal Run 1 
# Author: Eily Allan - modified from Erin D'Agnese and Ramon Gallego 
# Person running: Eily
# Last modified: 3/3/22 by Eily
# Date of run: 3/3/22 by Eily 

# Overview 
# This code is meant to take output from dada2 and assign taxonomy. This is adapted from Ramon Gallego's "insect.all.Rmd" script found here 
#(https://github.com/ramongallego/eDNA.and.Ocean.Acidification.Gallego.et.al.2020/tree/master/Scripts). The general overview is that ASVs from 
#dada2 will be read in and classified via insect, using the classifier publicly available from the creators of insect (found here - 
#https://cran.r-project.org/web/packages/insect/vignettes/insect-vignette.html). We will save the classifications for future runs 
#(so we don't need to keep classifying the same hash). 

#We will save two versions - the first if a hash is given any taxonomic rank by insect, and the second if a hash is given a "good" classification 
#by insect. We have different thresholds for what is a "good" classification for COI and 12S so we want to set different parameters - we want 12S 
#to get all the way to species level but with COI we expect less resolved taxonomic resolution so we will classify anything to order level as 
# "good". After we save these from run 1, in future runs, we will first check to see if the hashes were already classified here and then we will 
#only classify *new* hashes. Then we will add those new ones to our list of previously classified for the next run, and so on...

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
marker <- "MiFish"
#marker <- "MiMammal"
run.num <- "1"

# manually change file path because of the date in the folder for the dada2 output 
Hash_key <- here("Output", "dada2_output", "run1_20211117", marker, "hash_key.csv")
ASVs <- here("Output","dada2_output","run1_20211117",marker,"ASV_table.csv")

## REQUIRES HUMAN INPUT AND HARD CODE BECAUSE MUST BE STORED LOCALLY
# note mifish and mimammal use the same classifier right now
classifier <- "/Users/elizabethandruszkiewicz/GoogleDrive/UW/GitHub/NextGenNEPA_LOCAL/Input/classifiers/classifier_12S_v1.rds"

## Don't need to change 
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

# reformat to be as insect likes it
all.hashes.insect <- char2dna(Hash$Sequence)
names (all.hashes.insect) <- Hash$Hash
# check it
all.hashes.insect

# do the classification! this takes a while
clasif.hashes <- classify (x = all.hashes.insect, tree = tree, cores = 4)

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
# Save all results
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

# now, let's also save a copy of them in the general Output/classification_output folder that we can re-write over and over again as we add new classifications each time that we add another run 

# so here we are going to change the file path - and the file name to say "all.previous.hashes" 
# save it again as an rds and then also as a csv - and tax table
saveRDS(clasif.hashes, file=paste0(here("Output","classification_output"),"/", ".", marker,".all.previous.hashes.annotated.rds"))
clasif.hashes <- readRDS(file=paste0(here("Output","classification_output"),"/", ".", marker,".all.previous.hashes.annotated.rds"))
write.csv(clasif.hashes, file=paste0(here("Output","classification_output"),"/", ".", marker,".all.previous.hashes.annotated.csv"))
write.csv(taxtable,file=paste0(here("Output","classification_output"),"/", ".", marker,"previous.tax.table.csv"))


####################################################################
# Save only good results
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

# now, let's also save a copy of them in the general Output/classification_output folder that we can re-write over and over again as we add new classifications each time that we add another run 

# so here we are going to change the file path - and the file name to say "all.previous.hashes" 
# save it again as an rds and then also as a csv - and tax table
saveRDS(clasif.hashes, file=paste0(here("Output","classification_output"),"/", ".", marker,".all.good.previous.hashes.annotated.rds"))
clasif.hashes <- readRDS(file=paste0(here("Output","classification_output"),"/", ".", marker,".all.good.previous.hashes.annotated.rds"))
write.csv(clasif.hashes, file=paste0(here("Output","classification_output"),"/", ".", marker,".all.good.previous.hashes.annotated.csv"))
write.csv(taxtable,file=paste0(here("Output","classification_output"),"/", ".", marker,"previous.good.tax.table.csv"))

