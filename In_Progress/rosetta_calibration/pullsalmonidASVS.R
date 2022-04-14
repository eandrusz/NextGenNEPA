## Salmonids manual annotation 
## EAA 
## 4/12/22 


library(here)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(DECIPHER)
select <- dplyr::select

## Start with MiFish data 
marker <- "MiFish"

# Read in MiSeq data -- ASV table, annotations, and metadata
# HARDCODE IN FILEPATHS FOR NOW BC NOT SYNCING TO GOOGLE DRIVE - sorry this is annoying
ASV.table <- read_csv(paste0(here("Output","dada2_output","rs_20220404"),"/",marker,"/ASV_table.csv"))
hash.key <- read_csv(paste0(here("Output","dada2_output","rs_20220404"),"/",marker,"/Hash_key.csv"))
taxonomy <- read_csv(paste0(here("Output","classification_output","rs"),"/",marker,"/hashes.annotated.csv"))
input.metadata <- read_csv(paste0(here("Input","metadata"),"/rosetta_stone_start.csv"))                      
sequencing.metadata <- read_csv(paste0(here("Input","sequencing_metadata_files"),"/metadata-input-rs.csv"))
blast.file <- paste0(here("Output","classification_output"),"/rs/", marker,"/Rosetta_",marker,".txt")
blast.results <- read_delim(blast.file, col_names = c("qseqid", "sseqid", "sacc", "pident", "length", "mismatch", "gapopen", "qcovus", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "staxid", "qlen", "sscinames", "sseq"), delim = "\t" )

# take only what we need from the sequencing metadata file
sequencing.metadata <- sequencing.metadata %>% 
  select(c(Sample_name, Community, Tech_rep))

# take only salmonids from the insect taxonomy file
salmonidae.insect.hashes <- taxonomy %>% 
  filter(family == "Salmonidae")
  
# take only salmonids from the blast file
salmonidae.blast.hashes <- blast.results %>% 
  filter(str_detect(sscinames,"Salmon") | str_detect(sscinames,"Salvelinus") | str_detect(sscinames,"Oncorhynchus"))

# take hashes for some outgroups - just use insect for now 
#outgroup.insect.hashes <- taxonomy %>% 
  #filter(family == "Cottidae" | family == "Percidae")

# grab all hash ids of these groups
#all.hashes.keep <- unique(c(salmonidae.insect.hashes$representative, salmonidae.blast.hashes$qseqid, outgroup.insect.hashes$representative))
all.hashes.keep <- unique(c(salmonidae.insect.hashes$representative, salmonidae.blast.hashes$qseqid))


ASVs.for.tree <- ASV.table %>% 
  filter(Hash %in% all.hashes.keep)

hash.key.for.tree <- hash.key %>% 
  filter(Hash %in% all.hashes.keep)

DNAStringSet(hash.key.for.tree$Sequence) -> seqs.for.fasta
names(seqs.for.fasta) <- hash.key.for.tree$Hash

writeXStringSet(seqs.for.fasta, file=here("In_Progress","rosetta_calibration","output","salmon_seqs_MiMammal.fasta"))


#### compare insect vs. blast for salmonids
ASVs.for.tree2 <- ASVs.for.tree %>% 
  group_by(Hash) %>% 
  summarize(totalReads = sum(nReads))






