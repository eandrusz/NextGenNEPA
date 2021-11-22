tibble_to_matrix <- function (tb) {
  
  require(tidyverse)
  require(vegan)
  
  tb %>% 
    group_by(Sample_name, taxon) %>% 
    summarise(tot = sum(Normalized.reads)) %>% 
    spread ( key = "taxon", value = "tot", fill = 0) -> matrix_1
  samples <- pull (matrix_1, Sample_name)
  matrix_1 %>% 
    ungroup() %>% 
    dplyr::select ( - Sample_name) -> matrix_1
  data.matrix(matrix_1) -> matrix_1
  dimnames(matrix_1)[[1]] <- samples
  vegdist(matrix_1) -> matrix_1
}