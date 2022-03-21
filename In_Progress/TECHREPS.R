## NGN technical replicates 
## EAA 
## 3/17/21

library(here)
library(tidyverse)
library(patchwork)

date.dada <- "20220314"
marker <- "MiMammal"

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

## pull out only the technical replicates 

tech.reps <- df %>% 
  filter(Type == "tech_rep") %>% 
  #filter(Month.year == "821") %>% 
  group_by(Sample_ID, taxon) %>% 
  summarize(nReads = sum(nReads)) %>% 
  mutate(totReads = sum(nReads)) %>% 
  mutate(propReads = nReads/totReads) %>% 
  separate(Sample_ID, into = c("marker", "month", "creek", "station", "biorep"), remove=FALSE) %>% 
  unite(sample, c("month", "creek", "station"), remove=FALSE)
  
  
ggplot(tech.reps, aes(x = Sample_ID, y = propReads, fill = taxon)) +
  geom_col() +  theme_bw() +
  facet_wrap(~sample, scales = "free_x") +
  labs(x="Sequencing Run", y="Proportion of Reads") +
  #ggtitle("MiFish Technical Replicates")
  ggtitle("MiMammal Technical Replicates")

ggsave(paste0(here("Output","figures"), "/",marker,"_techreps.png"))


