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

all.enviro.asv.table <- read_csv(here("Output","dada2_output", "20220314.combined.MiMammal.ASV.table.csv"))
all.enviro.hash.key <- readRDS(here("Output","classification_output", "MiMammal.all.previous.hashes.annotated.rds"))

all.enviro.asv.table <- read_csv(here("Output","dada2_output", "20220401.combined.COI.ASV.table.csv"))
all.enviro.hash.key <- readRDS(here("Output","classification_output", "COI.all.previous.hashes.annotated.rds"))


enviro.all.reads <- sum(all.enviro.asv.table$nReads)
# MiFish 15150742
# MiMammal 16364770
# COI 20515157

enviro.annotated.reads <- all.enviro.asv.table %>% 
  filter(Hash %in% all.enviro.hash.key$representative) 
enviro.annotated.reads <- sum(enviro.annotated.reads$nReads)
# MiFish 11925664
# MiMammal 14229799
# COI 9475749

percentannotated <- 11925664/15150742*100
# MiFish 78% not bad
percentannotated <- 14229799/16364770*100
# MiMammal 86% not bad
percentannotated <- 9475749/20515157*100
# COI 46% not bad


simple.taxonomy <- all.enviro.hash.key %>% 
  select(c(representative,species))

tax.table.to.write <- all.enviro.asv.table %>% 
  rename(representative = Hash) %>% 
  left_join(simple.taxonomy, by = "representative") %>% 
  select(-Locus) %>% 
  group_by(Sample_name, species) %>% 
  summarize(Nreads = sum(nReads)) %>% 
  filter(!is.na(species)) %>% 
  rename(species = species) %>% 
  filter(!str_detect(Sample_name, "Kangaroo")) %>% 
  # mutate(., Sample_name = case_when(Sample_name == "MiFish.0321.1Prt.Dn.1TR2" ~ "MiFish.0321.1Prt.Dn.1.TR2",
  #                            Sample_name == "MiFish.0321.1Prt.Dn.1TR3" ~ "MiFish.0321.1Prt.Dn.1.TR3",
  #                            Sample_name == "MiFish.0421.3Chk.Dn.3TR2" ~ "MiFish.0421.3Chk.Dn.3.TR2",
  #                            Sample_name == "MiFish.0421.3Chk.Dn.3TR3" ~ "MiFish.0421.3Chk.Dn.3.TR3",
  #                            Sample_name == "MiFish.0821.2Brn.Dn.2TR2" ~ "MiFish.0821.2Brn.Dn.2.TR2", 
  #                            Sample_name == "MiFish.0821.2Brn.Dn.2TR3" ~ "MiFish.0821.2Brn.Dn.2.TR3",
  #                            TRUE ~ Sample_name)) %>% 
  # mutate(., Sample_name = case_when(Sample_name == "MiMammal.0321.2Brn.Up.1TR2" ~ "MiMammal.0321.2Brn.Up.1.TR2",
  #                                   Sample_name == "MiMammal.0321.2Brn.Up.1TR3" ~ "MiMammal.0321.2Brn.Up.1.TR3",
  #                                   Sample_name == "MiMammal.0421.3Chk.Up.1TR2" ~ "MiMammal.0421.3Chk.Up.1.TR2",
  #                                   Sample_name == "MiMammal.0421.3Chk.Up.1TR3" ~ "MiMammal.0421.3Chk.Up.1.TR3",
  #                                   Sample_name == "MiMammal.0821.4Pad.Up5.2TR2" ~ "MiMammal.0821.4Pad.Up5.2.TR2", 
  #                                   Sample_name == "MiMammal.0821.4Pad.Up5.2TR3" ~ "MiMammal.0821.4Pad.Up5.2.TR3",
  #                                   TRUE ~ Sample_name)) %>% 
filter(!str_detect(Sample_name, "OLD|MBT|HIP|Delta")) %>%
  mutate(., Sample_name = case_when(Sample_name == "COI.0321.2Brn.Up.2TR2" ~ "COI.0321.2Brn.Up.2.TR2",
                                  Sample_name == "COI.0321.2Brn.Up.2TR3" ~ "COI.0321.2Brn.Up.2.TR3",
                                  Sample_name == "COI.0421.4Pad.Up11.2TR2" ~ "COI.0421.4Pad.Up11.2.TR2",
                                  Sample_name == "COI.0421.4Pad.Up11.2TR3" ~ "COI.0421.4Pad.Up11.2.TR3",
                                  Sample_name == "COI.0821.5Sqm.Dn.1TR2" ~ "COI.0821.5Sqm.Dn.1.TR2", 
                                  Sample_name == "COI.0821.5Sqm.Dn.1TR3" ~ "COI.0821.5Sqm.Dn.1.TR3",
                                  TRUE ~ Sample_name)) %>% 
  separate(Sample_name, into=c("marker", "time", "creek", "site", "biol", "tech")) %>%
  select(-marker) %>% 
  mutate(., tech = case_when(tech == "TR2" ~ 2,
                             tech == "TR3" ~ 3,
                             is.na(tech) ~ 1)) 
  

write_rds(tax.table.to.write, file=paste0(here("In_Progress","rosetta_calibration","data"),"/MiFish.all.envirodata.RDS"))
write_rds(tax.table.to.write, file=paste0(here("In_Progress","rosetta_calibration","data"),"/MiMammal.all.envirodata.RDS"))
write_rds(tax.table.to.write, file=paste0(here("In_Progress","rosetta_calibration","data"),"/COI.incomplete.all.envirodata.RDS"))


tax.table.to.write %>% 
  filter(time == "0321") %>%  
  ggplot(aes(x = species, y = Nreads)) +
  geom_point() +  
  facet_grid(~site ~ creek, scales="free_y") +
  theme_bw() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  #ggtitle(label="MiFish - 0321")
  #ggtitle(label="MiMammal - 0321")
  ggtitle(label="COI - 0321")

