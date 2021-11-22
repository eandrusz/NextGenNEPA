custom.lca <- function (df, cutoff = 90) {
  require(taxonomizr)
  require(tidyverse)
  df %>%  # this function allows to change cutoff parameters for a specified dataframe (df)
    group_by(qseqid) %>%
    select( pident, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    nest() %>% # for each query, calculate the agreed taxonomy
    # ungroup %>% slice (1:10) %>%
    mutate(consensus = purrr::map(data,  function(.x) { 
      # If there are 100% matches - keep those and the 90s
      # If not, keep everything
      
      if(max(.x$pident > cutoff )){
        
        .x %>% 
          filter(pident > cutoff) %>% 
          select(-pident) %>% 
          condenseTaxa() %>% # agreement in Phylogeny
          paste(., collapse = "%") # Collapse all the taxa data separated by %
        
      }else{
        .x %>% 
          select(-pident) %>% 
          condenseTaxa() %>%
          paste(., collapse = "%")}
    }
    )) %>%
    select(qseqid, consensus) %>%
    unnest(consensus)} 