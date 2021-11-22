# tax.table.R to make the tax.table from the annotated hashes files created from insect
# usage tax.tabe(df)

tax.table <- function(df){
  require(tidyverse)
 
  taxonomy.file <- as.data.frame(df)
  row.names(taxonomy.file) = taxonomy.file[,1]
  tax.mat = taxonomy.file[,-1]
  tax.table <- as.matrix(tax.mat)
  
  
return(tax.table)  
  
}
