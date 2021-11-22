#make.phyloseq.obj.R usage make.phyloseq.obj.R(asv.mat,metadata,tax.table)

make.phyloseq.obj <- function(asv.mat,metadata, tax.mat) {
 
  require(tidyverse)
  require(vegan)
  require(phyloseq)
  

  ASV = otu_table(asv.mat, taxa_are_rows = T)
  TAX = tax_table(tax.mat)
  
  physeq = phyloseq(ASV,TAX)
  
  samples = sample_data(metadata)
  sample_names(samples) <- metadata$Sample_name
  #we either need line 17 or we need to change the sample_id to Sample_name in the metadata to match the dada2 output ASV.csv
  physeq = merge_phyloseq(physeq,samples)
  
return(physeq)  
  
}
