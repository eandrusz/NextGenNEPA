## NGN Run to run variability with 16S 
# EAA 
# 3/17/22


library(here)
library(tidyverse)
library(patchwork)

asvs <- list.files(path = "/Users/elizabethandruszkiewicz/Desktop/Ac16S_runtorun/output/dada2", pattern = "^ASV_table.csv", recursive = T, full.names = T)
Ac16S.asvs <- str_subset(asvs, "Ac16S")
metadata <- list.files(path ="/Users/elizabethandruszkiewicz/Desktop/Ac16S_runtorun/metadata/forcutadaptnew", pattern = "*fix.csv", recursive = T, full.names = T)


# merge the asv tables for runs 

run1 <- read_csv(Ac16S.asvs[1])
run1$Run <- 1

run5 <- read_csv(Ac16S.asvs[2])
run5$Run <- 5

run7 <- read_csv(Ac16S.asvs[3])
run7$Run <- 7

# combine

all.runs <- rbind(run1,run5,run7)

all.runs <- all.runs %>% 
  unite(sample.run, c(Sample_name, Run), remove=FALSE)
  
# plot by proportions of ASVs

plotasvs <- all.runs %>%
  group_by(sample.run) %>% 
  mutate(totReads = sum(nReads)) %>% 
  mutate(propReads = nReads/totReads) %>% 
  separate(Sample_name, into = c("marker", "Creek.Station", "Bio.Rep","Month.Year"), extra="drop", remove=FALSE)

ggplot(plotasvs, aes(x = as.factor(Run), y = propReads, fill=Hash)) +
  geom_col() +  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~Sample_name, scales="free_x") +
  #facet_grid(Creek.Station ~ Month.Year, scales = "free") +
  labs(x="Sample / Run", y="Proportion of Reads") +
  scale_x_discrete(limits = c("1","5","7")) +
  #ggtitle("Run to run variability with 16S samples - All ASVs") +
  ggtitle("Run to run variability with 16S samples")

ggsave(filename="/Users/elizabethandruszkiewicz/Desktop/Ac16S_runtorun/output/plots/stackedbars.png")

## try plotting as nmds

for.nmds <- plotasvs %>% 
  select(sample.run, Hash, nReads) %>% 
  pivot_wider(names_from = sample.run, values_from = nReads) %>% 
  column_to_rownames(var="Hash")
for.nmds[is.na(for.nmds)] = 0

bc.nmds <- metaMDS(t(for.nmds), distance = "bray")
bc.MDS1 = bc.nmds$points[,1] #store nmds values
bc.MDS2 = bc.nmds$points[,2] #store nmds values

meta.for.nmds <- plotasvs %>% 
  select(sample.run,Creek.Station,Month.Year,Run) %>% 
  unite(Creek.Station.Month.Year, c("Creek.Station","Month.Year"), sep=".") %>% 
  distinct()

pca.to.plot <- cbind(meta.for.nmds, bc.MDS1, bc.MDS2)

ggplot(pca.to.plot, aes(x=bc.MDS1, y=bc.MDS2)) +
  geom_point(size=4, aes(color=factor(Creek.Station.Month.Year), shape=factor(Run))) +
  theme_bw() +
  scale_shape_manual(values = c(15,16,17,18)) +
  labs(x="PC1",y="PC2", color="Creek", shape = "Sequencing Run") +
  ggtitle('Ac16S Run to Run Variability - Bray Curtis')

ggsave(filename="/Users/elizabethandruszkiewicz/Desktop/Ac16S_runtorun/output/plots/ndms.png")

