## NGN - all technical replicates 
## EAA 
## 3/17/22


### note COI is fucked - come back to this 

library(here)
library(tidyverse)
library(patchwork)

date.dada <- "20220316"
marker <- "COI"

ASV.table <- read_csv(paste0(here("Output", "dada2_output"), "/",date.dada, ".combined.", marker,".ASV.table.csv"))
annotations <- readRDS(file=paste0(here("Output","classification_output"),"/",marker, ".all.good.previous.hashes.annotated.rds"))
all.metadata <- read.csv(here("Input", "sequencing_metadata_files", "master_sequencing_datasheet_20220317.csv"))

## add annotations to asvs

df <- ASV.table %>% 
  rename(representative = Hash) %>% 
  left_join(annotations, by = "representative") %>% 
  rename(Sample_ID = Sample_name) %>% 
  left_join(all.metadata,by = "Sample_ID") %>% 
  group_by(Sample_ID) %>% 
  mutate(sample_read_depth = sum(nReads))

## separate into kangaroo and not kangaroo 

kang <- df %>% 
  filter(str_detect(Sample_ID,"Kangaroo")) %>%  # select only the kangaroo samples 
  select(c(Sample_ID, nReads, taxon, sample_read_depth)) %>% 
  group_by(Sample_ID, taxon) %>% 
  summarize(nReads = sum(nReads)) %>% 
  mutate(totReads = sum(nReads)) %>% 
  mutate(propReads = nReads/totReads) %>% 
  #mutate(Sample_ID = recode(Sample_ID, "MiFish.Kangaroo.Run.2" = "Run 2", "MiFish.Kangaroo.Run.3" = "Run 3", "MiFish.Kangaroo.Run.4" = "Run 4", "MiFish.Kangaroo.Run.5" = "Run 5", "MiFish.Kangaroo.Run.6" = "Run 6", "MiFish.Kangaroo.Run.7.NA" = "Run7"))
  mutate(Sample_ID = recode(Sample_ID, "COI.0321.Kangaroo.Kangaroo.1" = "Run 1", "COI.Kangaroo.Run.2" = "Run 2", "COI.Kangaroo.Run.3" = "Run 3", "COI.Kangaroo.Run.4" = "Run 4", "MiFish.Kangaroo.Run.5" = "Run 5", "MiFish.Kangaroo.Run.6" = "Run 6", "MiFish.Kangaroo.Run.7.NA" = "Run7"))


notkang <- df %>% 
  filter(!str_detect(Sample_ID,"Kangaroo")) %>%  # select only real environmental samples 
  filter(genus=="Osphranter")   # keep only samples that have at least one ASV that is kangaroo 

## plot kangaroo sample compositions 

kang.raw.plot <- ggplot(kang, aes(x = Sample_ID, y = nReads, fill = taxon)) +
  geom_col() +  theme_bw() +
  labs(x="Sequencing Run", y="Total Number of Reads", title = "Kangaroo Samples") 
  
kang.prop.plot <- ggplot(kang, aes(x = Sample_ID, y = propReads, fill = taxon)) +
  geom_col() +  theme_bw() +
  labs(x="Sequencing Run", y="Proportion of Reads") 


kang.raw.plot + kang.prop.plot + plot_layout(ncol = 1, nrow = 2)
ggsave(paste0(here("Output","figures"), marker,"_kangaroo_samples_taxa.png"))
  
         