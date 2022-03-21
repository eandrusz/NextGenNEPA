## NGN Run to run variability with 16S 
# EAA 
# 3/17/22


library(here)
library(tidyverse)
library(patchwork)

asvs <- list.files(path = here("Output","dada2_output"), pattern = "^ASV_table.csv", recursive = T, full.names = T)
Ac16S.asvs <- str_subset(asvs, "Ac16S")
all.metadata <- read.csv(here("Input", "sequencing_metadata_files", "master_sequencing_datasheet_20220317.csv"))

# merge the asv tables for both runs 

run1 <- read_csv(Ac16S.asvs[1])
run1$Sample_name <- as.character(run1$Sample_name)
run1 <- run1 %>% 
  left_join(all.metadata, by = "Sample_name") %>% 
  select(Sample_ID, Hash, nReads, Sequencing.run) %>% 
  rename("Sample_name" = "Sample_ID") 

run5 <- read_csv(Ac16S.asvs[2])
run5.meta <- all.metadata %>% 
  filter(Sequencing.run == 5)
run5 <- run5 %>% 
  rename("Sample_ID" = "Sample_name") %>% 
  left_join(run5.meta, by = "Sample_ID") %>% 
  select(Sample_ID, Hash, nReads, Sequencing.run) %>% 
  rename("Sample_name" = "Sample_ID") 

run6 <- read_csv(Ac16S.asvs[3])
run6$Sample_name <- as.character(run6$Sample_name)
run6 <- run6 %>% 
  left_join(all.metadata, by = "Sample_name") %>% 
  select(Sample_ID, Hash,nReads, Sequencing.run) %>% 
  rename("Sample_name" = "Sample_ID")

run7 <- read_csv(Ac16S.asvs[4])
run7.meta <- all.metadata %>% 
  filter(Sequencing.run == 7)
run7 <- run7 %>% 
  rename("Sample_ID" = "Sample_name") %>% 
  left_join(run7.meta, by = "Sample_ID") %>% 
  select(Sample_ID, Hash, nReads) %>% 
  mutate(Sequencing.run = 7) %>% 
  rename("Sample_name" = "Sample_ID") %>% 
  separate(Sample_name, into = c("marker", "creekstn", "biorep", "month"), remove=FALSE) %>% 
  unite(Sample_name, c("marker", "creekstn", "biorep", "month"), sep=".",na.rm = TRUE,remove=TRUE)

# combine
all.runs <- rbind(run1,run5,run6,run7)
all.runs <- all.runs %>% 
  unite(sample.run, c(Sample_name, Sequencing.run), remove=FALSE)

# take top 25 ASVs and then plot side by side stacked bar charts with relative abundance 
topasvs <- all.runs %>% 
  group_by(Hash) %>% 
  mutate(totReadsASV= sum(nReads)) %>% 
  arrange(desc(totReadsASV)) %>% 
  select(Hash, totReadsASV) %>% 
  distinct() 

topasvs <- topasvs[1:10,]
  
# check 
plottopasvs <- all.runs %>%
  filter(Hash %in% topasvs$Hash) %>% 
  group_by(sample.run) %>% 
  mutate(totReads = sum(nReads)) %>% 
  mutate(propReads = nReads/totReads) %>% 
  separate(Sample_name, into = c("marker", "Creek.Station", "Month.Year"), extra="drop", remove=FALSE)

ggplot(plottopasvs, aes(x = sample.run, y = propReads, fill=Hash)) +
  geom_col() +  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~Sample_name, scales="free_x") +
  #facet_grid(Creek.Station ~ Month.Year, scales = "free") +
  labs(x="Sample / Run", y="Proportion of Reads") +
  #ggtitle("Run to run variability with 16S samples - All ASVs") +
  ggtitle("Run to run variability with 16S samples - Top 10 ASVs")




