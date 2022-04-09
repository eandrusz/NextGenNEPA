# NGN mock community data prep 
# EAA 
# 4/5/22

library(here)
library(tidyverse)
library(ggplot2)
library(gridExtra)
select <- dplyr::select

## Start with MiFish data 
marker <- "MiMammal"

# Read in MiSeq data -- ASV table, annotations, and metadata
# HARDCODE IN FILEPATHS FOR NOW BC NOT SYNCING TO GOOGLE DRIVE - sorry this is annoying
ASV.table <- read_csv(paste0(here("Output","dada2_output","rs_20220404"),"/",marker,"/ASV_table.csv"))
taxonomy <- read_csv(paste0(here("Output","classification_output","rs"),"/",marker,"/hashes.annotated.csv"))
#blast.taxonomy <- read_csv(paste0(here("Output","classification_output","rs"),"/",marker,"/hashes.annotated.blast.csv"))
input.metadata <- read_csv(paste0(here("Input","metadata"),"/rosetta_stone_start.csv"))                      
sequencing.metadata <- read_csv(paste0(here("Input","sequencing_metadata_files"),"/metadata-input-rs.csv"))


# take only what we need from the sequencing metadata file
sequencing.metadata <- sequencing.metadata %>% 
  select(c(Sample_name, Community, Tech_rep))

# take only what we need from the taxonomy file
taxonomy <- taxonomy %>% 
  rename(Hash = representative) %>% 
  select(c(Hash, species))


## MANUAL FIXING THINGS 

taxonomy <- taxonomy %>% 
  add_row(Hash="7818939be318eeca036d8b3906645d17e81c5012", species="Random_gBlock") %>% 
  add_row(Hash="c2b79afb7f9c2c9714de591bf2d93d3cfdd1749e", species="Random_gBlock") %>% 
  add_row(Hash="2b8b0c7200e210117e842d5622e31d69571e398d", species="Random_gBlock") %>% 
  add_row(Hash="be22865680c34f64c6ee65c21164a45929dddca6", species="Random_gBlock") %>% 
  add_row(Hash="85f094efc653f5e3928552bee2d3594c53e9d9a4", species="Mammoth") %>% 
  add_row(Hash="20c7f9dbdf79b3ad51debc2428d4a2d27daf3069", species="Mammoth") %>% 
  add_row(Hash="c2c8136d70447488b26ff4951c4f2f3c65122da0", species="Mammoth") %>%
  add_row(Hash="db00dff518d6ea45063e9365b3244dfc93d0191d", species="Moa") %>% 
  add_row(Hash="dc1779e8fa382387a1bdee2fff5959c7611b90b3", species="Moa")

# for now, take only mock communities A, B, and C -- and then add on Zack's metadata 
if(marker == "MiFish" | marker == "COI") {
  start.metadata <- rbind(input.metadata, zack) 
} else if ( marker == "MiMammal") {
  start.metadata <- input.metadata
} 
start.metadata <- expand_grid(start.metadata,Tech_rep = 1:3)

seqs.tax <- ASV.table %>% 
  left_join(sequencing.metadata, by ="Sample_name") %>% 
  select(-Locus) %>% 
  left_join(taxonomy, by = "Hash") %>% 
  select(c(Community, species, nReads, Tech_rep)) %>% 
  filter(!is.na(species)) %>% # for now, just remove anything that did not get annotated to species level (go back later?)
  rename(Species = species) %>% 
  filter(!str_detect(Community,"gB")) %>%   # for now, just take MCA, MCB, MCC (no gblocks yet)
  group_by(Community, Species, Tech_rep) %>% 
  summarize(nReads = sum(nReads)) %>% 
  group_by(Community, Tech_rep) %>% 
  mutate(totReads = sum(nReads)) %>% 
  mutate(propReads = nReads/totReads) %>% 
  select(-totReads)

gblocks <- ASV.table %>% 
  left_join(sequencing.metadata, by ="Sample_name") %>% 
  select(-Locus) %>% 
  left_join(taxonomy, by = "Hash") %>% 
  select(c(Community, species, nReads, Tech_rep)) %>% 
  filter(!is.na(species)) %>% # for now, just remove anything that did not get annotated to species level (go back later?)
  rename(Species = species) %>% 
  filter(str_detect(Community,"gB")) %>%   # for now, just take MCA, MCB, MCC (no gblocks yet)
  group_by(Community, Species, Tech_rep) %>% 
  summarize(nReads = sum(nReads)) %>% 
  group_by(Community, Tech_rep) %>% 
  mutate(totReads = sum(nReads)) %>% 
  mutate(propReads = nReads/totReads) %>% 
  select(-totReads)

# issue with piranga ludoviciana not matching... hard code it for now
#seqs.tax[103:105,2] <- start.metadata[121:123,2]

mock.data <- start.metadata %>% 
  full_join(seqs.tax) %>% 
  mutate(Mass_Input = if_else(is.na(Mass_Input), 0, Mass_Input)) %>% 
  mutate(b_proportion = if_else(is.na(Prop_Input), 0, Prop_Input)) %>% 
  mutate(Nreads = if_else(is.na(nReads), 0, nReads)) %>% 
  mutate(proportion_reads = if_else(is.na(propReads), 0, propReads)) %>% 
  mutate(N_pcr_mock = 43) %>%  # 35 cycles of regular PCR and then 8 of indexing (Zack just used the first PCR for this column...)
  rename(species = Species) %>% 
  rename(tech = Tech_rep) %>% 
  rename(site = Community) %>% 
  select(-c(Mass_Input, nReads, Prop_Input, propReads))
  
# write output to file 
write_rds(mock.data, file=paste0(here("In_Progress","rosetta_calibration","data"),"/",marker,"mockdata.RDS"))


# what the hell does it look like? 

plot(mock.data$b_proportion, mock.data$proportion_reads)

ggplot(mock.data, aes(x = tech, y = proportion_reads, fill = species)) +
  geom_col() +  
  guides(fill = "none") +
  facet_wrap(~site, scales="free_x") +
  theme_bw() +
  labs(x="Technical Replicate", y="Proportion of Input DNA", title = "Expected Proportions") +
  ggtitle(label=paste0(marker))

ggsave(paste0(here("In_Progress","rosetta_calibration","output"),"/",marker,"mockdata.png"))


gbstart <- data.frame(Community=rep("gB", times=9), Species=rep(c("Moa","Mammoth","Random_gBlock"), each=3),Tech_rep=rep(1:3, times=3),b_proportion=rep(c(195/729,267/729,267/729), each=3))

onlygb <- gblocks %>% 
  filter(Community == "gB") %>% 
  left_join(gbstart, by=c("Community", "Species","Tech_rep")) %>% 
  #filter(Tech_rep != 1) %>% 
  ungroup() %>% 
  add_row(Community= "gB Input", Species = "Moa", Tech_rep = 0, nReads = 0, propReads = 0.2674897, b_proportion=0) %>% 
  add_row(Community= "gB Input", Species = "Mammoth", Tech_rep = 0, nReads = 0, propReads = 0.3662551, b_proportion=0) %>% 
  add_row(Community= "gB Input", Species = "Random_gBlock", Tech_rep = 0, nReads = 0, propReads = 0.3662551, b_proportion=0)

ggplot(onlygb, aes(x = Tech_rep, y = propReads, fill = Species)) +
  geom_col() +  
  #guides(fill = "none") +
  #facet_wrap(~Community, scales="free_x") +
  theme_bw() +
  labs(x="Technical Replicate", y="Proportion of Input DNA", title = "Expected Proportions") +
  ggtitle(label=paste0(marker)) + 
  scale_fill_brewer(palette="Dark2")

ggsave(paste0(here("In_Progress","rosetta_calibration","output"),"/",marker,"gblocks.png"))

