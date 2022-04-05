## Processing mock communities 
# EA 
# 4/4/22

## Read in files 
metadata <- read.csv(here("Input", "sequencing_metadata_files","metadata-input-rs.csv"))
MiFish.ASV.table <- read_csv(paste0(here("Output", "dada2_output", "rs_20220404","MiFish","ASV_table.csv")))
MiFish.annotations <- readRDS(file=paste0(here("Output","classification_output","rs","MiFish","hashes.annotated.rds")))


## start joining 
MiFish.ASV.taxon <- MiFish.ASV.table %>% 
  rename(representative = Hash) %>% 
  left_join(MiFish.annotations, by = "representative") %>% 
  select(-Locus)

# and same - keep things that do have an annotation (reminder to species or genus level so "good" annotation)
MiFish.ASV.yes.taxon <- MiFish.ASV.taxon[! is.na(MiFish.ASV.taxon$taxon),]

# we can also do this at the species level 
MiFish.by.species <- MiFish.ASV.yes.taxon %>% 
  select(-representative) %>% # we don't need the hash identifier anymore
  filter(species != "") %>% 
  group_by(Sample_name, class, species) %>% # for each sample that has multiple asvs that assign to the same taxa...
  summarise(tot = sum(nReads)) %>% 
  left_join(metadata, by = "Sample_name")

ggplot(MiFish.by.species, aes(x = Sample_name, y = tot, fill = species)) +
  geom_col() +  
  facet_wrap(~Community, scales="free_x") +
  theme_bw() +
  labs(x="Sample", y="Total Number of Reads", title = "MiFish Mock Communities") 

