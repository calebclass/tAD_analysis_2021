#make acid-aro-basic features
library(tidyverse)

#load gcn4 and hsf1 datasets
gcn4_tads <- readRDS('Datasets/Datasets/hahn_tads.rds') %>% as_tibble()

hsf1_tads <- readRDS('Datasets/Datasets/random_library.rds') %>% as_tibble() %>% 
  mutate(
    AD_set = ifelse(binary_stop == 'live', 'AD_positive','AD_negative'),
    len = nchar(sequence)
  ) %>% filter(len > 4) %>% select(AD_set, sequence) %>% as_tibble()

#get positions of basic aromatic and acidic amino acids in hsf1----
##acidic info----
hsf1_tads$acid_pos <- 0

for(i in 1:nrow(hsf1_tads)){
  seq <- strsplit(hsf1_tads$sequence[i],'')[[1]]
  pos <- grep('[DE]', seq)
  avg_pos <- mean(pos)
  hsf1_tads$acid_pos[i] <- avg_pos
}

hsf1_tads$acidic_pos <- lapply(hsf1_tads$sequence,function(x){
  str_split(x,'')[[1]] %>% grep('[DE]',.) %>% mean()
}) %>% unlist()

##basic----
hsf1_tads$basic_pos <- 0

for(i in 1:nrow(hsf1_tads)){
  seq <- strsplit(hsf1_tads$sequence[i],'')[[1]]
  pos <- grep('[RHK]', seq)
  avg_pos <- mean(pos)
  hsf1_tads$basic_pos[i] <- avg_pos
}

##aromatic----
hsf1_tads$aromatic_pos <- 0

for(i in 1:nrow(hsf1_tads)){
  seq <- strsplit(hsf1_tads$sequence[i],'')[[1]]
  pos <- grep('[WYF]', seq)
  avg_pos <- mean(pos)
  hsf1_tads$aromatic_pos[i] <- avg_pos
}
################################################
#####TAKES A LONG TIME to run ---- MAKE SURE TO SAVE ----
#get positions of basic aromatic and acidic amino acids in gcn4
##acidic info----
gcn4_tads$acid_pos <- 0

for(i in 1:nrow(gcn4_tads)){
  seq <- strsplit(gcn4_tads$sequence[i],'')[[1]]
  pos <- grep('[DE]', seq)
  avg_pos <- mean(pos)
  gcn4_tads$acid_pos[i] <- avg_pos
}

gcn4_tads$acidic_pos <- lapply(gcn4_tads$sequence,function(x){
  str_split(x,'')[[1]] %>% grep('[DE]',.) %>% mean()
}) %>% unlist()

##basic----
gcn4_tads$basic_pos <- 0

for(i in 1:nrow(gcn4_tads)){
  seq <- strsplit(gcn4_tads$sequence[i],'')[[1]]
  pos <- grep('[RHK]', seq)
  avg_pos <- mean(pos)
  gcn4_tads$basic_pos[i] <- avg_pos
}

##aromatic----
gcn4_tads$aromatic_pos <- 0

for(i in 1:nrow(gcn4_tads)){
  seq <- strsplit(gcn4_tads$sequence[i],'')[[1]]
  pos <- grep('[WYF]', seq)
  avg_pos <- mean(pos)
  gcn4_tads$aromatic_pos[i] <- avg_pos
}

saveRDS(gcn4_tads,'gcn4_avg_pos.rds')
