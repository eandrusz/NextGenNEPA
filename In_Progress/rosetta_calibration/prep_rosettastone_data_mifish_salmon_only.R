# NGN mock community data prep 
# EAA 
# 4/5/22

library(here)
library(tidyverse)
library(ggplot2)
library(gridExtra)
select <- dplyr::select

## Start with MiFish data 
marker <- "MiFish"

# Read in MiSeq data -- ASV table, annotations, and metadata
# HARDCODE IN FILEPATHS FOR NOW BC NOT SYNCING TO GOOGLE DRIVE - sorry this is annoying
ASV.table <- read_csv(paste0(here("Output","dada2_output","rs_20220404"),"/",marker,"/ASV_table.csv"))
taxonomy <- read_csv(paste0(here("Output","classification_output","rs"),"/",marker,"/hashes.annotated.csv"))
#blast.taxonomy <- read_csv(paste0(here("Output","classification_output","rs"),"/",marker,"/manual.hashes.annotated.blast.csv"))
input.metadata <- read_csv(paste0(here("Input","metadata"),"/rosetta_stone_start.csv"))                      
sequencing.metadata <- read_csv(paste0(here("Input","sequencing_metadata_files"),"/metadata-input-rs.csv"))

# Get Zack's starting metadata - probably should double check this with him - pulling from the github here
zack_all <- readRDS(here("In_Progress","rosetta_calibration","data","mifish_mock_community_data.RDS"))
zack <- zack_all %>% 
  filter(community == "North_Even") %>% 
  filter(tech_rep == 1) %>% 
  filter(Cycles == 39) %>% 
  filter(start_conc_ng != 0) %>% 
  select(c(community,ID_mifish, start_conc_ng)) %>% 
  rename(Mass_Input = start_conc_ng) %>% 
  mutate(Prop_Input = Mass_Input / sum(Mass_Input)) %>% 
  rename(Species = ID_mifish) %>% 
  rename(Community = community) %>% 
  mutate(Community = "MCZ")

# take only what we need from the sequencing metadata file
sequencing.metadata <- sequencing.metadata %>% 
  select(c(Sample_name, Community, Tech_rep))

# take only what we need from the taxonomy file
taxonomy <- taxonomy %>% 
  rename(Hash = representative) %>% 
  filter(family == "Salmonidae") %>% 
  select(c(Hash, taxon))
  
## MANUAL FIXING THINGS 
taxonomy <- taxonomy %>% 
  filter(Hash != "f3b1b944f738536215cd355d76bc8c2bffd66624") %>% 
  add_row(Hash="f3b1b944f738536215cd355d76bc8c2bffd66624", taxon="Oncorhynchus mykiss") %>% 
  filter(Hash != "631402d8ed10c8b197e1d67c208b29b9ba254a64") %>% 
  add_row(Hash="631402d8ed10c8b197e1d67c208b29b9ba254a64", taxon="Salvelinus malma") %>% 
  filter(Hash != "2c79a691f61d8877b4c2f3953447d6e033680d33") %>% 
  add_row(Hash="2c79a691f61d8877b4c2f3953447d6e033680d33", taxon="Salvelinus confluentus") %>% 
  filter(Hash != "1a49f11fcb4d876f9b04aebbdae6602f8586ecc9") %>% 
  rename(species = taxon)


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

mock.data <- start.metadata %>% 
  filter(Species %in% taxonomy$species) %>% 
  full_join(seqs.tax) %>% 
  mutate(Mass_Input = if_else(is.na(Mass_Input), 0, Mass_Input)) %>% 
  mutate(b_proportion = if_else(is.na(Prop_Input), 0, Prop_Input)) %>% 
  mutate(Nreads = if_else(is.na(nReads), 0, nReads)) %>% 
  mutate(proportion_reads = if_else(is.na(propReads), 0, propReads)) %>% 
  mutate(N_pcr_mock = 43) %>%  # 35 cycles of regular PCR and then 8 of indexing (Zack just used the first PCR for this column...)
  rename(species = Species) %>% 
  rename(tech = Tech_rep) %>% 
  rename(site = Community) %>% 
  select(-c(Mass_Input, nReads, Prop_Input, propReads)) %>% 
  unite(com_rep, c(site,tech), remove=FALSE) %>% 
  filter(com_rep != "MCZ_1") %>%  
  dplyr::select(-com_rep)
  
# write output to file 
write_rds(mock.data, file=paste0(here("In_Progress","rosetta_calibration","data"),"/",marker,".salmonidonly.mockdata.RDS"))


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
