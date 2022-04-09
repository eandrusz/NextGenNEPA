# Mock community sleuthing 


library(tidyverse)
library(here)

mf.taxa <- read_csv(paste0(here("Output","classification_output","rs"),"/MiFish/hashes.annotated.csv"))
mm.taxa <- read_csv(paste0(here("Output","classification_output","rs"),"/MiMammal/hashes.annotated.csv"))
coi.taxa <- read_csv(paste0(here("Output","classification_output","rs"),"/COI/hashes.annotated.csv"))

mf.species <- unique(mf.taxa$species) # 53
mm.species <- unique(mm.taxa$species) # 38
coi.species <- unique(coi.taxa$species) # 75

# how many unique total across all markers
all.species <- unique(c(mf.species, mm.species, coi.species)) # 99 total 

# how many are unique to each marker 
mf.only <- mf.species[!(mf.species %in% union(coi.species,mm.species))] # 12 
mm.only <- mm.species[!(mm.species %in% union(coi.species,mf.species))] # 2
coi.only <- coi.species[!(coi.species %in% union(mf.species,mm.species))] # 44

# how many are shared 
all.shared <- intersect(coi.species, intersect(mf.species, mm.species)) # 26
mf.mm.shared <- intersect(mf.species, mm.species) # subtract out the all.shared so 36 - 26 = 10 
mf.coi.shared <- intersect(mf.species, coi.species) # subtract out the all.shared so 31 - 26 = 5 
mm.coi.shared <- intersect(mm.species, coi.species) # subtract out the all.shared so 26 - 26 = 0

# what about when we add in blast 

mf.blast <- paste0(here("Output","classification_output"),"/rs/MiFish/Rosetta_mifish.txt")
mf.blast <- read_delim(mf.blast, col_names = c("qseqid", "sseqid", "sacc", "pident", "length", "mismatch", "gapopen", "qcovus", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "staxid", "qlen", "sscinames", "sseq"), delim = "\t" )

coi.blast <- paste0(here("Output","classification_output"),"/rs/COI/Rosetta_COI.txt")
coi.blast <- read_delim(coi.blast, col_names = c("qseqid", "sseqid", "sacc", "pident", "length", "mismatch", "gapopen", "qcovus", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "staxid", "qlen", "sscinames", "sseq"), delim = "\t" )

coi.blast.notinsect <- setdiff(coi.blast$sscinames, coi.species)


## look at blast results for cottus
mf.blast.cottus <- mf.blast %>% 
  filter(str_detect(sscinames,"Cottus"))

mf.blast.salvelinus <- mf.blast %>% 
  filter(str_detect(sscinames,"Salvelinus"))


coi.blast.cottus <- coi.blast %>% 
  filter(str_detect(sscinames,"Cottus"))

coi.blast.salvelinus <- coi.blast %>% 
  filter(str_detect(sscinames,"Salvelinus"))

