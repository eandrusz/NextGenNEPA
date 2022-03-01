## DADA2 for NGN sequenicing runs 
# Author: Eily via Moncho
# Person running: Eily
# Last modified: 2/28/22 by Eily
# Date of run: 2/28/22 by Eily 

# Overview
# After removing primers by cutadapt and splitting primers into subfolder, we are going to apply a denoising algorithm `DADA2`, https://github.com/benjjneb/dada2 to generate ASV tables and hash keys. On the parameters section you have chosen a trimming length for each of the two files, whether to use Hashes or not, and a folder where to drop the output files. 
# The trimming length depends on the marker and the run. In general, for a paired end 2x300 run, we use trimming lengths of:
#MiFish/MiMammal: 150/150
#COI: 250/200
#This is because our 12S (MF/MM) fragment is only ~170 bp long and our COI (Leray XT) fragment is ~313 bp long. These should not be blindly used, but instead, after plotting the Q score plots for a subset (or aggregate) samples they should be revisited. The Q score plots vary for each run and if the run is poor you should adjust to have these be shorter and if it is a great run you can extend the trimming lengths for the COI marker to get more overlap. The rule of thumb is to trim when Q scores drop below 30. So go back to the parameters after running the chunk of plotting Q scores. 

####################################################################
# Set up
####################################################################

## Load packages 
library (here)
library (tidyverse)
library (dada2)
library (digest)
library (seqinr)
library (kableExtra)

## Set parameters - run number, marker, etc. - NEEDS HUMAN INPUT
### CHANGE EACH TIME - run number and marker 
run.num = 5
marker = "COI"

### CHANGE EACH TIME: THE NO PRIMER FOLDER MUST BE HARDCODED TO LOCAL HARD DRIVE BECAUSE IT IS TOO BIG FOR GITHUB
# also note that if you run this on a different day than when you did cutadapt change to the actual date rather than sys.date()
#noprimerfolder = paste0("/Users/elizabethandruszkiewicz/GoogleDrive/UW/GitHub/NextGenNEPA_LOCAL/Output/cutadapt_output/run",run.num,"_",format(Sys.Date(), "%Y%m%d"), "/noprimers/", marker)
noprimerfolder = paste0("/Users/elizabethandruszkiewicz/GoogleDrive/UW/GitHub/NextGenNEPA_LOCAL/Output/cutadapt_output/run",run.num,"_","20220226", "/noprimers/", marker)

### set parameters that don't change usually
hash = TRUE  # you rarely do NOT want to hash
keep.mid.files = FALSE  # I find I never look at these / use these and they just take up space 

## Set parameters automatically based on human input above
### set working directory to no primer folder
setwd(noprimerfolder)

### set file path for metadata based on run and marker 
metadata = paste0(noprimerfolder,"/cutadapt_output_metadata_",marker,".csv")
sample.metadata <- read_csv(metadata)

### set trimming lengths based on marker 
if(marker == "MiFish" | marker == "MiMammal") {
  trimming.length.Read1 = 150
  trimming.length.Read2 = 150 
} else if ( marker == "COI" | marker == "Ac16S") {
  trimming.length.Read1 = 250
  trimming.length.Read2 = 200
} else {
  print("Different marker -- set for trimming lengths.")
}

### set up output file paths
output.folder = paste0(here("Output","dada2_output"),"/run",run.num,"_",format(Sys.Date(), "%Y%m%d"), "/", marker)
# Check if output directory exists - if not create as a subfolder of input dir
if(!dir.exists(file.path(output.folder))){
  dir.create(path = file.path(output.folder),recursive = T)
  output.dir <- file.path(output.folder)
}else{
  output.dir <- file.path(output.folder)
}

### create a subfolder in noprimerfolder for filtered reads from dada2
filt_path <- file.path(noprimerfolder, "filtered")


####################################################################
# QC - remove any bad samples 
# *Be careful here to make sure you comment out and don't by accident drop good samples!*
####################################################################

## Plot quality score plots and change trimming lengths as needed
# It is probably fine to just plot a subset of you samples (4) and check Q score plots. You can also plot the whole aggregate R1 and R2 Q scores (this will take longer to run.) See where Q scores start dropping below 30 to determine trimming length. Change parameters if you need to!! 

### create a subset of four files to plot Q scores for
ifelse(nrow(sample.metadata)>4,
       subset <- sample.metadata %>%  sample_n(4),
       subset <- sample.metadata)

subset %>% pull(file1) %>%
  plotQualityProfile(.)
ggsave(file.path(output.dir,"R1_subset.png"))

subset %>% pull(file2) %>%
  plotQualityProfile(.)
ggsave(file.path(output.dir,"R2_subset.png"))

### plot aggregate - takes a while so comment out if you don't want
sample.metadata %>% pull(file1) %>%
  plotQualityProfile(., aggregate=TRUE)
ggsave(file.path(output.dir,"R1_aggregate.png"))

sample.metadata %>% pull(file2) %>%
  plotQualityProfile(., aggregate=TRUE)
ggsave(file.path(output.dir,"R2_aggregate.png"))


## REMOVE ANY BAD SAMPLES - if the aggregate *doesn't run*, some files are bad in there so we need to remove them
## removed any samples removed during filtering so edit metadata file to remove it
#Locus_MiFish_MiFish-0621-4Pad-Up11-3_S45_L001_R1_001.fastq)
#Locus_MiMammal_MiMammal-0721-4Pad-Up11-2_S19_L001_R1_001.fastq
#Locus_Ac16S_Ac16S-PadUp11-4-0321_S113_F1_filt.fastq.gz
#sample.metadata <- sample.metadata[-5,]


## CHANGE TRIMMING LENGTHS IF YOU NEED TO AFTER LOOKING AT Q SCORE PLOT - REQUIRES HUMAN INPUT
trimming.length.Read1 = trimming.length.Read1
trimming.length.Read2 = trimming.length.Read2 

####################################################################
# dada2 
# You should not need to change *anything* below this point - it should run perfectly. It is broken up into different chunks so you can see it as you go along, but just run each chunk sequentially. For reference, it usually takes on the order of minutes for 12S and Ac16S and on the order of hours for COI data.
####################################################################

### SOME OF COI IS BACKWARDS - check before runing

## run dada2 
sample.metadata %>%
  #filter(rc == 1) %>% # ONLY SELECT THE BACKWARDS ONES (1) OR FORWARDS ONES (0)
  separate(file1, into = "basename", sep= "_L001_R1", remove = F) %>% 
  mutate(filtF1 = file.path("filtered", paste0(basename, "_F1_filt.fastq.gz")),
         filtR1 = file.path("filtered", paste0(basename, "_R1_filt.fastq.gz"))) %>%
  select(-basename) %>% 
  mutate (outFs = pmap(.l= list (file1, filtF1, file2, filtR1),
                       .f = function(a, b, c, d) {
                         filterAndTrim(a,b,c,d,
                                       truncLen=c(trimming.length.Read1,trimming.length.Read2),
                                       maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                                       compress=TRUE, multithread=TRUE )
                       } ),
          errF1 = map(filtF1, ~ learnErrors(.x, multithread=TRUE,verbose = 0)),     # Calculate errors
          errR1 = map(filtR1, ~ learnErrors(.x, multithread=TRUE,verbose = 0)),
          derepF1 = map(filtF1, derepFastq),                   # dereplicate seqs
          derepR1 = map(filtR1, derepFastq),
          dadaF1  = map2(derepF1,errF1, ~ dada(.x, err = .y, multithread = TRUE)),  # dada2
          dadaR1  = map2(derepR1,errR1, ~ dada(.x, err = .y, multithread = TRUE)),
          mergers = pmap(.l = list(dadaF1,derepF1, dadaR1,derepR1),                 # merge things
                         .f = mergePairs )) -> output.dada2

if (keep.mid.files==TRUE){
  write_rds(output.dada2, path = "output.halfway.rds")}

seqtabF <- makeSequenceTable(output.dada2$mergers)

# check it
dim(seqtabF)
table(nchar(getSequences(seqtabF)))

## Remove chimera sequences. This is done within dada2 but you could use other ways

seqtab.nochim <- removeBimeraDenovo(seqtabF, method="consensus", multithread=TRUE)
dim(seqtab.nochim)
seqtab.nochim.df <- as.data.frame(seqtab.nochim)



## Set up output files
### Copy the metadata so it is all in one place
sample.metadata %>% write_csv(file.path(output.dir,"metadata.csv"))

### Output files
conv_file <- file.path(output.dir,"hash_key.csv")
conv_file.fasta <- file.path(output.dir,"hash_key.fasta")
ASV_file <-  file.path(output.dir,"ASV_table.csv")

## Hash or not?
# usually always hash
if (hash==TRUE)
{conv_table <- tibble( Hash = "", Sequence ="")
  map_chr (colnames(seqtab.nochim.df), ~ digest(.x, algo = "sha1", serialize = F, skip = "auto")) -> Hashes
  conv_table <- tibble (Hash = Hashes,
                        Sequence = colnames(seqtab.nochim.df))
  colnames(seqtab.nochim.df) <- Hashes

  write_csv(conv_table, conv_file) # write the table into a file
  write.fasta(sequences = as.list(conv_table$Sequence),
              names     = as.list(conv_table$Hash),
              file.out = conv_file.fasta)
  seqtab.nochim.df <- bind_cols(sample.metadata %>%
                                  select(Sample_name, Locus),
                                seqtab.nochim.df)
  seqtab.nochim.df %>%
    pivot_longer(cols = c(- Sample_name, - Locus),
                 names_to = "Hash",
                 values_to = "nReads") %>%
    filter(nReads > 0) -> current_asv
  write_csv(current_asv, ASV_file)    }else{
    #What do we do if you don't want hashes: two things - Change the header of the ASV table, write only one file
    seqtab.nochim.df %>%
      pivot_longer(cols = c(- Sample_name, - Locus),
                   names_to = "Sequence",
                   values_to = "nReads") %>%
      filter(nReads > 0) -> current_asv
    write_csv(current_asv, ASV_file)
  }



####################################################################
# Check how samples look
# Track fate of reads 
####################################################################

## Generate output summary
getN <- function(x) sum(getUniques(x))

output.dada2 %>%
  select(-file1, -file2, -filtF1, -filtR1, -errF1, -errR1, -derepF1, -derepR1) %>%
  mutate_at(.vars = c("dadaF1", "dadaR1", "mergers"),
            ~ sapply(.x,getN)) %>%
  #  pull(outFs) -> test
  mutate(input = map_dbl(outFs, ~ .x[1]),
         filtered = map_dbl(outFs, ~ .x[2]),
         tabled  = rowSums(seqtabF),
         nonchim = rowSums(seqtab.nochim)) %>%
  select(Sample_name,
         Locus,
         input,
         filtered,
         denoised_F = dadaF1,
         denoised_R = dadaR1,
         merged = mergers,
         tabled,
         nonchim) -> track
write_csv(track, file.path(output.dir,"dada2_summary.csv"))

## drop
if (keep.mid.files==FALSE){
  unlink(filt_path, recursive = T)
}


## Make output_summary table and fig
kable (track, align= "c", format = "html") %>%
  kable_styling(bootstrap_options= "striped", full_width = T, position = "center") %>%
  column_spec(2, bold=T)

track %>%
  mutate_if(is.numeric, as.integer) %>%
  pivot_longer(cols = c(-Sample_name, -Locus),
               names_to = "Step",
               values_to = "Number of Sequences") %>%
  mutate (Step = fct_relevel(Step,
                             levels = c( "input","filtered","denoised_F" ,"denoised_R" , "merged" , "tabled", "nonchim"))) %>%
  ggplot(aes(x = Step, y = `Number of Sequences`, group =  Sample_name, color = Sample_name)) +
  geom_line() +
  facet_wrap(~Sample_name) +
  guides(color = "none")

