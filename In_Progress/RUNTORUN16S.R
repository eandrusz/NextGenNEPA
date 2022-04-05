## NGN Run to run variability with 16S 
# EAA 
# 3/17/22


library(here)
library(tidyverse)
library(patchwork)

asvs <- list.files(path = here("Output","dada2_output"), pattern = "^ASV_table.csv", recursive = T, full.names = T)
Ac16S.asvs <- str_subset(asvs, "Ac16S")
all.metadata <- read.csv(here("Input", "sequencing_metadata_files", "master_sequencing_datasheet_20220324_FIXED.csv"))

blast.file <- "/Users/elizabethandruszkiewicz/GoogleDrive/UW/GitHub/NextGenNEPA/Output/classification_output/run1/Ac16S/Ac16Sruns1567toblast.output.txt"
blast.results <- read_delim(blast.file, col_names = c("qseqid", "sseqid", "sacc", "pident", "length", "mismatch", "gapopen", "qcovus", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "staxid", "qlen", "sscinames", "sseq"), delim = "\t")
blast.results.simple <- blast.results %>% 
  select(c(qseqid, sscinames)) %>% 
  rename("Hash" = "qseqid") %>% 
  separate(sscinames, into="species", sep=";", extra="drop") %>% 
  separate(species, into=c("genus", "species"), extra="drop") %>% 
  distinct() %>% 
  group_by(Hash) %>% 
  top_n(n=1, genus) %>% 
  top_n(n=1, Hash)

# merge the asv tables for both runs 

run1 <- read_csv(Ac16S.asvs[1])
run1$Sample_name <- as.character(run1$Sample_name)
run1 <- run1 %>% 
  left_join(all.metadata, by = "Sample_name") %>% 
  select(Sample_ID, Hash, nReads, Sequencing.run) %>% 
  rename("Sample_name" = "Sample_ID") 

run5 <- read_csv(Ac16S.asvs[2])
run5 <- run5 %>% 
  select(-Locus) %>%
  mutate(Sequencing.run = 5) 

run6 <- read_csv(Ac16S.asvs[3])
run6 <- run6 %>% 
  left_join(all.metadata, by = "Sample_name") %>% 
  select(Sample_ID, Hash,nReads, Sequencing.run) %>% 
  rename("Sample_name" = "Sample_ID")

run7 <- read_csv(Ac16S.asvs[4])
run7 <- run7 %>% 
  separate(Sample_name, into = c("marker", "creekstn", "biorep", "month"), remove=TRUE) %>% 
  unite(Sample_name, c("marker", "creekstn", "biorep", "month"), sep=".",na.rm = TRUE,remove=TRUE) %>% 
  select(-Locus) %>% 
  mutate(Sequencing.run = 7)
  
# combine
all.runs <- rbind(run1,run5,run6,run7)
all.runs <- rbind(run1,run5,run7)

all.runs <- all.runs %>% 
  unite(sample.run, c(Sample_name, Sequencing.run), remove=FALSE)

all.runs.taxa <- all.runs %>% 
  left_join(blast.results.simple, by="Hash") %>% 
  select(-Hash) %>% 
  group_by(sample.run, genus) %>% 
  summarize(totreads = sum(nReads)) %>% 
  left_join(all.runs, by= "sample.run")
  

# take top 25 ASVs and then plot side by side stacked bar charts with relative abundance 
topasvs <- all.runs %>% 
  group_by(Hash) %>% 
  mutate(totReadsASV= sum(nReads)) %>% 
  arrange(desc(totReadsASV)) %>% 
  select(Hash, totReadsASV) %>% 
  distinct() 
topasvs <- topasvs[1:20,]
  
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

plottaxa <- all.runs.taxa %>%
  group_by(sample.run) %>% 
  mutate(totReads = sum(totreads)) %>% 
  mutate(propReads = totreads/totReads) %>% 
  separate(Sample_name, into = c("marker", "Creek.Station", "Month.Year"), extra="drop", remove=FALSE)

ggplot(plottaxa, aes(x = sample.run, y = propReads, fill=genus)) +
  geom_col() +  theme_bw() +
  #theme(legend.position = "none") +
  facet_wrap(~Sample_name, scales="free_x") +
  #facet_grid(Creek.Station ~ Month.Year, scales = "free") +
  labs(x="Sample / Run", y="Proportion of Reads") +
  #ggtitle("Run to run variability with 16S samples - All ASVs") +
  ggtitle("Run to run variability with 16S samples - Genus level")

## try plotting as nmds
for.nmds <- plottopasvs %>% 
  select(sample.run, Hash, nReads) %>% 
  pivot_wider(names_from = sample.run, values_from = nReads) %>% 
  column_to_rownames(var="Hash")
for.nmds[is.na(for.nmds)] = 0

for.nmds <- all.runs.taxa %>% 
  select(sample.run, genus, nReads) %>% 
  group_by(sample.run, genus) %>% 
  mutate(allreads = sum(nReads)) %>% 
  pivot_wider(names_from = sample.run, values_from = nReads) %>% 
  column_to_rownames(var="genus")
for.nmds[is.na(for.nmds)] = 0

bc.nmds <- metaMDS(t(for.nmds), distance = "bray")
bc.MDS1 = bc.nmds$points[,1] #store nmds values
bc.MDS2 = bc.nmds$points[,2] #store nmds values

meta.for.nmds.pad <- plottopasvs %>% 
  select(sample.run) %>% 
  filter(str_detect(sample.run,"Pad")) %>% 
  separate(sample.run, into = c("marker", "Creek.Station", "biorep", "Month.Year","Run"), extra="drop", remove=FALSE) %>% 
  unite(Creek.Station.Month.Year, c("Creek.Station", "Month.Year"), sep=".",na.rm = TRUE,remove=FALSE) %>%
  unite(Creek.Station.Month.Year.BioRep, c("Creek.Station.Month.Year","biorep"), sep=".",na.rm = TRUE,remove=FALSE) %>% 
  distinct()
meta.for.nmds.chk <- plottopasvs %>% 
  select(sample.run) %>% 
  filter(str_detect(sample.run,"Chk")) %>% 
  separate(sample.run, into = c("marker", "Creek.Station.biorep", "Month.Year","Run"), extra="drop", remove=FALSE) %>% 
  extract(Creek.Station.biorep, into=c("Creek.Station", "biorep"), "(.{5})(.{1})", remove=TRUE) %>% 
  unite(Creek.Station.Month.Year, c("Creek.Station", "Month.Year"), sep=".",na.rm = TRUE,remove=FALSE) %>% 
  unite(Creek.Station.Month.Year.BioRep, c("Creek.Station.Month.Year","biorep"), sep=".",na.rm = TRUE,remove=FALSE) %>% 
  distinct()
meta.for.nmds <- rbind(meta.for.nmds.pad, meta.for.nmds.chk)

pca.to.plot <- cbind(meta.for.nmds, bc.MDS1, bc.MDS2)

ggplot(pca.to.plot, aes(x=bc.MDS1, y=bc.MDS2)) +
  geom_point(size=4, aes(color=factor(Creek.Station.Month.Year), shape=factor(Run))) +
  theme_bw() +
  scale_shape_manual(values = c(15,16,17,18)) +
  labs(x="PC1",y="PC2", color="Creek", shape = "Sequencing Run") +
  ggtitle('Ac16S Run to Run Variability - Bray Curtis')

#ggsave(file="/Users/erdag/github/NextGenNEPA-main/Output/MiFish.allruns.salmonids.index.creeks.bysite.braycurtisPCA.png", dpi = "retina")

#doesn't have much of a discernible pattern, may need to remove tech reps to do this?


