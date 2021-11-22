# Usage rc_rehash(asv.table, hash.key) all elements without quotes

rc_rehash<- function(asv.table, hash.key){ 
  require(dplyr)
  require(insect)
  require(digest)
  
  #reverse complement the representative sequence for each hash
  new_hash_key <- hash_key
  rc_seq <- vector(length = dim(hash_key)[1])
  # this is ridiculous and shouldn't be for loop but wtf i can't get it to work 
  for (i in 1:dim(hash_key)[1]) {
    rc_seq[i] = rc(hash_key[i,2])
  }
  new_hash_key$Sequence <- rc_seq
  
  # great - but now we need to rehash - be careful because there are different algorithms
  new_hash_key$Hash <- sapply(new_hash_key$Sequence, digest, algo = "sha1", serialize = F, skip = "auto")
  
  #now change hashes in asv table and resave it}
  old_to_new_hash_key <- data.frame(hash_key$Hash, new_hash_key$Hash)
  colnames(old_to_new_hash_key) <- c("Hash", "NewHash")
  new_asv_table <- inner_join(asv_table, old_to_new_hash_key) 
  new_asv_table <- new_asv_table[, !(colnames(new_asv_table) %in% "Hash")]
  
  new_asv_table <- new_asv_table %>% 
    rename(Hash = NewHash)

  return(new_hash_key)
  return(new_asb_table)
  }
