## NGN technical replicates 
## EAA 
## 3/17/21

library(here)
library(tidyverse)
library(patchwork)

#date.dada <- "20220314"
marker <- "MiFish"

asvs <- list.files(path = here("Output","dada2_output"), pattern = "^ASV_table.csv", recursive = T, full.names = T)
marker.asvs <- str_subset(asvs, marker)
all.metadata <- read.csv(here("Input", "sequencing_metadata_files", "master_sequencing_datasheet_20220324_FIXED.csv"))
annotations <- readRDS(file=paste0(here("Output","classification_output"),"/",marker, ".all.good.previous.hashes.annotated.rds"))

##### only look at runs 1, 5, and 7 right now
 
run1 <- read_csv(marker.asvs[1])
run1$Sample_name <- as.character(run1$Sample_name)
run1 <- run1 %>% 
  left_join(all.metadata, by = "Sample_name") %>% 
  select(Sample_ID, Hash, nReads, Sequencing.run) %>% 
  rename("Sample_name" = "Sample_ID") %>% 
  separate(Sample_name, into = c("marker", "month", "creek", "stn", "biorep"), extra= "drop", remove=TRUE) %>% 
  unite(Sample_name, c("marker", "month", "creek", "stn", "biorep"), sep=".", na.rm = TRUE,remove=TRUE)


run5 <- read_csv(marker.asvs[5])
run5$Sequencing.run <- 5
run5 <- run5 %>% 
  select(-Locus)

run7 <- read_csv(marker.asvs[7])
run7 <- run7 %>% 
  separate(Sample_name, into = c("marker", "month", "creek", "stn", "biorep","techrep"), remove=TRUE) %>% 
  unite(Sample_name, c("marker", "month", "creek", "stn", "biorep","techrep"), sep=".",na.rm = TRUE,remove=TRUE) %>% 
  mutate(Sequencing.run = 7) %>% 
  select(-Locus)


# combine

all.runs <- rbind(run1,run5,run7)

all.runs <- all.runs %>% 
  unite(sample.run, c(Sample_name, Sequencing.run), remove=FALSE)

# plot by proportions of ASVs

plotasvs <- all.runs %>%
  group_by(sample.run) %>% 
  mutate(totReads = sum(nReads)) %>% 
  mutate(propReads = nReads/totReads) %>% 
  separate(Sample_name, into = c("marker", "month", "creek", "stn", "biorep","techrep"), extra="drop", remove=FALSE) %>% 
  unite(creekstn, c("creek","stn"), remove = FALSE) %>% 
  unite(biotech, c("biorep","techrep"),na.rm = TRUE, remove = FALSE)

ggplot(plotasvs, aes(x = as.factor(biotech), y = propReads, fill=Hash)) +
  geom_col() +  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~creekstn ~ month, scales="free_x") +
  #facet_grid(Creek.Station ~ Month.Year, scales = "free") +
  labs(x="Sample / Run", y="Proportion of Reads") +
  #scale_x_discrete(limits = c("1","5","7")) +
  #ggtitle("Run to run variability with 16S samples - All ASVs") +
  ggtitle("Run to run + technical + biological variability with MiFish samples")

ggsave(filename="/Users/elizabethandruszkiewicz/Desktop/mifish_variability.png")

## try plotting as nmds

for.nmds <- plotasvs %>% 
  filter(!str_detect(sample.run, "Kangaroo")) %>% 
  select(sample.run, Hash, nReads) %>% 
  pivot_wider(names_from = sample.run, values_from = nReads) %>% 
  column_to_rownames(var="Hash")
for.nmds[is.na(for.nmds)] = 0

bc.nmds <- metaMDS(t(for.nmds), distance = "bray")
bc.MDS1 = bc.nmds$points[,1] #store nmds values
bc.MDS2 = bc.nmds$points[,2] #store nmds values

meta.for.nmds <- plotasvs %>% 
  filter(!str_detect(sample.run, "Kangaroo")) %>%
  select(sample.run,creekstn,month,Sequencing.run) %>% 
  unite(Creek.Station.Month.Year, c("creekstn","month"), sep=".", remove=FALSE) %>% 
  distinct()

pca.to.plot <- cbind(meta.for.nmds, bc.MDS1, bc.MDS2)

ggplot(pca.to.plot, aes(x=bc.MDS1, y=bc.MDS2)) +
  geom_point(size=4, aes(color=factor(creekstn), shape=factor(month))) +
  theme_bw() +
  scale_shape_manual(values = c(15,16,17,18,19)) +
  labs(x="PC1",y="PC2", color="Creek-Station", shape = "Month") +
  ggtitle('MiFish Run to Run Variability - Bray Curtis')

ggsave(filename="/Users/elizabethandruszkiewicz/Desktop/mifish_ndms.png")


## add annotations to asvs

df <- plotasvs %>% 
  rename(representative = Hash) %>% 
  left_join(annotations, by = "representative") #%>% 
  #rename(Sample_ID = Sample_name) %>% 
  #left_join(all.metadata,by = "Sample_ID") %>% 
  #group_by(Sample_name) %>% 
  #mutate(sample_read_depth = sum(nReads))

## pull out only the technical replicates 

df.species <- df %>% 
  #filter(Type == "tech_rep") %>% 
  #filter(Month.year == "821") %>% 
  group_by(Sample_name, taxon) %>% 
  summarize(nReads = sum(nReads)) %>% 
  mutate(totReads = sum(nReads)) %>% 
  mutate(propReads = nReads/totReads) %>% 
  separate(Sample_name, into = c("marker", "month", "creek", "stn", "biorep","techrep"), extra="drop", remove=FALSE) %>% 
  unite(creekstn, c("creek","stn"), remove = FALSE) %>% 
  unite(biotech, c("biorep","techrep"),na.rm = TRUE, remove = FALSE)
  
  
ggplot(df.species, aes(x = biotech, y = propReads, fill = taxon)) +
  geom_col() +  theme_bw() +
  facet_wrap(~creekstn ~ month, scales = "free_x") +
  labs(x="Replicate", y="Proportion of Reads") +
  #ggtitle("MiFish Technical Replicates")
  ggtitle("MiFish by species")

ggsave(paste0(here("Output","figures"), "/",marker,"_techreps.png"))


## try plotting as nmds

for.nmds <- df.species %>% 
  filter(!str_detect(Sample_name, "Kangaroo")) %>% 
  select(Sample_name, taxon, nReads) %>% 
  filter(!is.na(taxon)) %>% 
  pivot_wider(names_from = Sample_name, values_from = nReads) %>% 
  column_to_rownames(var="taxon")
for.nmds[is.na(for.nmds)] = 0

bc.nmds <- metaMDS(t(for.nmds), distance = "bray")
bc.MDS1 = bc.nmds$points[,1] #store nmds values
bc.MDS2 = bc.nmds$points[,2] #store nmds values

meta.for.nmds <- plotasvs %>% 
  filter(!str_detect(sample.run, "Kangaroo")) %>%
  select(Sample_name,creekstn,month,Sequencing.run) %>% 
  unite(Creek.Station.Month.Year, c("creekstn","month"), sep=".", remove=FALSE) %>% 
  distinct() %>% 
  filter(Sample_name %in% colnames(for.nmds))

pca.to.plot <- cbind(meta.for.nmds, bc.MDS1, bc.MDS2)

ggplot(pca.to.plot, aes(x=bc.MDS1, y=bc.MDS2)) +
  geom_point(size=4, aes(color=factor(creekstn), shape=factor(month))) +
  theme_bw() +
  scale_shape_manual(values = c(15,16,17,18,19)) +
  labs(x="PC1",y="PC2", color="Creek-Station", shape = "Month") +
  ggtitle('MiFish Run to Run Variability - Bray Curtis')

ggsave(filename="/Users/elizabethandruszkiewicz/Desktop/mifish_ndms.png")

