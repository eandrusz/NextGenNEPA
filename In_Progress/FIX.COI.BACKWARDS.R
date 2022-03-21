library(here)
library(tidyverse)
library(vegan)
library(reshape2)


# Read the merged ASV tables
COI.ASV.table <- read_csv(paste0(here("Output", "dada2_output", "20220316.combined.COI.ASV.table.csv")))
COI.hash.key <- read_csv(paste0(here("Output", "dada2_output", "20220316.combined.COI.hash.key.csv")))

# Read in the classification output for each run 
# the files written to the general folder (not run specific) should have everything classified from all runs because they are our "databases" to use for future runs - we want to use the "DATE.MARKER.all.good.previous.hashes.annotated.rds" file 
COI.annotations <- readRDS(file=paste0(here("Output","classification_output"),"/COI.all.previous.hashes.annotated.rds"))

all.metadata <- read.csv(here("Input", "sequencing_metadata_files", "master_sequencing_datasheet_20211026.csv"))


# first map annotations onto hash key - try to see if things are still backwards 
roots <- COI.annotations %>% 
  filter(taxon == "root") %>% 
  rename(Hash = representative) %>% 
  left_join(COI.ASV.table, by= "Hash") %>% 
  rename(Sample_ID = Sample_name) %>% 
  left_join(all.metadata, by= "Sample_ID")
  

notroots <- COI.annotations %>% 
  filter(taxon != "root")
