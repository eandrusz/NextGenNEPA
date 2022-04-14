## Check environmental samples for salmonids
## EAA 
## 4/12/22


library(here)
library(tidyverse)
library(ggplot2)
library(gridExtra)
select <- dplyr::select

# for mifish 

## look at environmental hash key 
all.enviro.asv.table <- read_csv(here("Output","dada2_output", "20220314.combined.MiFish.ASV.table.csv"))
all.enviro.hash.key <- readRDS(here("Output","classification_output", "MiFish.all.previous.hashes.annotated.rds"))

enviro.all.reads <- sum(all.enviro.asv.table$nReads)
# 15150742

# annotated hashes sum of reads
enviro.annotated.reads <- all.enviro.asv.table %>% 
  filter(Hash %in% all.enviro.hash.key$representative) 
enviro.annotated.reads <- sum(enviro.annotated.reads$nReads)
# 11925664

percentannotated <- 11925664/15150742*100
# 78% not bad

annotated.asv.table <- all.enviro.asv.table %>% 
  rename(representative = Hash) %>% 
  left_join(all.enviro.hash.key, by = "representative")



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
  mutate(., Sample_name = case_when(Sample_name == "MiFish.0321.1Prt.Dn.1TR3" ~ "MiFish.0321.1Prt.Dn.1.TR3",
                             Sample_name == "MiFish.0421.3Chk.Dn.3TR2" ~ "MiFish.0421.3Chk.Dn.3.TR2",
                             Sample_name == "MiFish.0421.3Chk.Dn.3TR3" ~ "MiFish.0421.3Chk.Dn.3.TR3",
                             Sample_name == "MiFish.0821.2Brn.Dn.2TR2" ~ "MiFish.0821.2Brn.Dn.2.TR2", 
                             Sample_name == "MiFish.0821.2Brn.Dn.2TR3" ~ "MiFish.0821.2Brn.Dn.2.TR3",
                             TRUE ~ Sample_name)) %>% 
  separate(Sample_name, into=c("marker", "time", "creek", "site", "biol", "tech")) %>%
  select(-marker) %>% 
  mutate(., tech = case_when(tech == "TR2" ~ 2,
                             tech == "TR3" ~ 3,
                             is.na(tech) ~ 1)) 
  

write_rds(tax.table.to.write, file=paste0(here("In_Progress","rosetta_calibration","data"),"/MiFish.salmonidonly.envirodata.RDS"))


tax.table.to.write %>% 
  filter(time == "0321") %>%  
  ggplot(aes(x = species, y = Nreads)) +
  geom_point() +  
  facet_grid(~site ~ creek, scales="free_y") +
  theme_bw() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  ggtitle(label="MiFish - 0321")

