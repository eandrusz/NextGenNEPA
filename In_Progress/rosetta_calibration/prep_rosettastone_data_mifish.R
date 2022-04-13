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
  dplyr::select(c(Hash, species))

goodtaxonomy <- taxonomy %>% 
  filter(!is.na(species))

blast.taxonomy <- blast.taxonomy %>% 
  rename(Hash = qseqid, species = sscinames) %>% 
  select(c(Hash, species))

taxonomytoadd <- taxonomy %>% 
  filter(is.na(species)) %>% 
  left_join(blast.taxonomy, by= "Hash") %>% 
  dplyr::select(-species.x) %>% 
  rename(species = species.y) %>% 
  filter(!is.na(species))

taxonomy <- rbind(goodtaxonomy, taxonomytoadd)

## MANUAL FIXING THINGS 
# "Anaxyrus exsul" -> "Anaxyrus boreas"
#  Ardea purpurea -> Ardea herodias
# Micropterus salmoides salmoides -> "Micropterus salmoides"
# Neogale vison -> Neovison vison
# Oncorhynchus clarkii lewisi -> Oncorhynchus clarkii 
# Oreochromis niloticus x Oreochromis aureus ->  Tilapia sparrmanii
# philander opossum -> Didelphis virginiana
# Salmo trutta trutta -> Salmo trutta 

taxonomy <- taxonomy %>% 
  mutate(., species = case_when(species == "Anaxyrus exsul" ~ "Anaxyrus boreas",
                                species == "Ardea purpurea" ~ "Ardea herodias",
                                #species == "Micropterus salmoides salmoides" ~ "Micropterus salmoides",
                                #species == "Neogale vison" ~ "Neovison vison",
                                #species == "Oncorhynchus clarkii lewisi" ~ "Oncorhynchus clarkii",
                                species == "Oreochromis niloticus x Oreochromis aureus" ~ "Tilapia sparrmanii",
                                species == "Philander opossum" ~ "Didelphis virginiana",
                                #species == "Salmo trutta trutta" ~ "Salmo trutta",
                                species == "Lithobates catesbeianus" ~ "Rana catesbeiana",
                                species == "Macropus cf. giganteus TL-2021" ~ "Osphranter rufus",
                                species == "Macropus fuliginosus" ~ "Osphranter rufus",
                                species == "Macropus giganteus" ~ "Osphranter rufus",
                                species == "Piranga leucoptera" ~ "Piranga ludoviciana",
                                #species == "Rhinoraja longicauda;Bathyraja trachouros" ~ "Bathyraja abyssicola",
                                #species == "Puma concolor" ~ "Felis catus", 
                                species == "Castor fiber" ~ "Castor canadensis", 
                                species == "Petrogale xanthopus" ~ "Osphranter rufus", 
                                species == "Setophaga kirtlandii" ~ "Cardellina pusilla",
         TRUE ~ species)) %>% 
  add_row(Hash="7818939be318eeca036d8b3906645d17e81c5012", species="Random_gBlock") %>% 
  add_row(Hash="c2b79afb7f9c2c9714de591bf2d93d3cfdd1749e", species="Random_gBlock")


taxonomy <- taxonomy %>% 
  filter(representative != "2c79a691f61d8877b4c2f3953447d6e033680d33") %>% 
  add_row(representative="2c79a691f61d8877b4c2f3953447d6e033680d33", sscinames="Salvelinus confluentus") %>% 
  filter(representative != "f3b1b944f738536215cd355d76bc8c2bffd66624") %>% 
  add_row(representative="f3b1b944f738536215cd355d76bc8c2bffd66624", sscinames="Oncorhynchus mykiss") %>% 
  filter(representative != "631402d8ed10c8b197e1d67c208b29b9ba254a64") %>% 
  add_row(representative="631402d8ed10c8b197e1d67c208b29b9ba254a64", sscinames="Salvelinus malma") %>% 
  add_row(representative="7818939be318eeca036d8b3906645d17e81c5012", species="Random_gBlock") %>% 
  add_row(representative="c2b79afb7f9c2c9714de591bf2d93d3cfdd1749e", species="Random_gBlock")

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
seqs.tax[111:113,2] <- start.metadata[121:123,2]

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
write_rds(mock.data, file=paste0(here("In_Progress","rosetta_calibration","data"),"/",marker,"mockdatablast.RDS"))


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
  mutate(., Species = case_when(Species == "Dinornis robustus" ~ "Moa",
                                Species == "Mammuthus primigenius" ~ "Mammoth",
                                TRUE ~ Species)) %>% 
  left_join(gbstart, by=c("Community", "Species","Tech_rep")) %>% 
  filter(Tech_rep != 1) %>% 
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

