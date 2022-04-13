## Check environmental samples for salmonids
## EAA 
## 4/12/22


library(here)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(DECIPHER)
select <- dplyr::select

# for mifish 
# take the salmonids from the mock community and see if the same ASVs are found in environmental samples

salmonids <- readDNAStringSet(here("In_Progress","rosetta_calibration","output","salmon_seqs_MiMammal.fasta"))
salmon.hashes <- names(salmonids)

rename <- dplyr::rename

## look at environmental hash key 
all.enviro.asv.table <- read_csv(here("Output","dada2_output", "20220314.combined.MiMammal.ASV.table.csv"))
intersect(salmon.hashes, all.enviro.asv.table$Hash)

## things that are annotated
all.enviro.hash.key <- readRDS(here("Output","classification_output", "MiMammal.all.previous.hashes.annotated.rds"))
intersect(salmon.hashes, all.enviro.hash.key$representative)

# Ok first pull all salmonid data out from the hash key 
enviro.salmon.hashes <- all.enviro.hash.key %>% 
  filter(family == "Salmonidae")
intersect(salmon.hashes, enviro.salmon.hashes$representative)

enviro.salmon.asv.table <- all.enviro.asv.table %>% 
  filter(Hash %in% enviro.salmon.hashes$representative) %>% 
  dplyr::rename(representative = Hash) %>% 
  left_join(enviro.salmon.hashes, by = "representative")

enviro.all.salmon.reads <- sum(enviro.salmon.asv.table$nReads)
# 2702513

enviro.salmon.MChashes <- enviro.salmon.asv.table %>% 
  filter(representative %in% salmon.hashes)
enviro.salmon.MChashes.reads <- sum(enviro.salmon.MChashes$nReads)
# 2454614

# WOW THAT IS GREAT- 
percent = 2454614/2702513*100
# 91%

enviro.all.reads <- sum(all.enviro.asv.table$nReads)
# 16364771

# annotated hashes sum of reads
enviro.annotated.reads <- all.enviro.asv.table %>% 
  filter(Hash %in% all.enviro.hash.key$representative) 
enviro.annotated.reads <- sum(enviro.annotated.reads$nReads)
# 14229799

percentannotated <- 14229799/16364771*100
# 87% not bad

#### write hash key and dada2 output for only salmonids with only hashes found in mock community 
taxonomy.to.write <- enviro.salmon.asv.table %>% 
  filter(representative %in% salmon.hashes) %>% 
  dplyr::select(c(representative, taxon)) %>% 
  distinct()

asv.table.to.write <- all.enviro.asv.table %>% 
  filter(Hash %in% salmon.hashes) 

tax.table.to.write <- asv.table.to.write %>% 
  dplyr::rename(representative = Hash) %>% 
  left_join(taxonomy.to.write, by = "representative") %>% 
  select(-Locus) %>% 
  group_by(Sample_name, taxon) %>% 
  summarize(Nreads = sum(nReads)) %>% 
  filter(!is.na(taxon)) %>% 
  rename(species = taxon) %>% 
  filter(!str_detect(Sample_name, "Kangaroo")) %>% 
  mutate(., Sample_name = case_when(Sample_name == "MiMammal.0421.3Chk.Up.1TR2" ~ "MiMammal.0421.3Chk.Up.1.TR2",
                             Sample_name == "MiMammal.0421.3Chk.Up.1TR3" ~ "MiMammal.0421.3Chk.Up.1.TR3",
                             Sample_name == "MiMammal.0821.4Pad.Up5.2TR2" ~ "MiMammal.0821.4Pad.Up5.2.TR2",
                             Sample_name == "MiMammal.0821.4Pad.Up5.2TR3" ~ "MiMammal.0821.4Pad.Up5.2.TR3", 
                             TRUE ~ Sample_name)) %>% 
  separate(Sample_name, into=c("marker", "time", "creek", "site", "biol", "tech")) %>%
  select(-marker) %>% 
  mutate(., tech = case_when(tech == "TR2" ~ 2,
                             tech == "TR3" ~ 3,
                             is.na(tech) ~ 1)) 
  

write_rds(tax.table.to.write, file=paste0(here("In_Progress","rosetta_calibration","data"),"/MiMammal.salmonidonly.envirodata.RDS"))
