  ##Script to build dipeptide heatmap,
  #data on live_percent
  
  #load packages
  library(tidyverse)
  library(caret)
  
  # Load in tad datasets----
  gcn4_tads <- readRDS('Datasets/Datasets/hahn_tads.rds') %>% 
    select(sequence, AD_set) %>% 
    as_tibble()
  
  hsf1_tads <- readRDS('Datasets/Datasets/random_library.rds') %>% 
    mutate(AD_set = ifelse(binary_stop == 'live', 'AD_positive', 'AD_negative')) %>% 
    select(sequence, AD_set) %>% 
    as_tibble()
  
  # Build list of possible dipeptides----
  #create aa vector
  aa_vector <- c('R','H','K','D','E',
                 'S','T','N','Q',
                 'A','V','L','I','M','F','Y','W',
                 'C','G','P')
  
  #make a list of dipeptides
  aa_combos <- expand.grid(aa_vector, aa_vector) %>% as_tibble() %>% 
    mutate(dipep = paste(Var1 %>% as.character(),
                         Var2 %>% as.character(), sep = ''))
  
  
  # Add dipep count columns to hsf1 and gcn4----
  for(i in aa_combos$dipep){
    hsf1_tads$hold <- str_count(hsf1_tads$sequence, i)
    #little hack to correct column names
    colnames(hsf1_tads)[which(colnames(hsf1_tads) == 'hold')] <- i
  }
  
  for(i in aa_combos$dipep){
    gcn4_tads$hold <- str_count(gcn4_tads$sequence, i)
    #little hack to correct column names
    colnames(gcn4_tads)[which(colnames(gcn4_tads) == 'hold')] <- i
  }

#########################################################
  ##get live_percent of each dipep  ##chane datasets gcn4 to hsf1 to run different libraries
  live <- gcn4_tads %>% filter(AD_set == 'AD_positive') %>% select(-sequence, -AD_set)
  live <- map(live, sum)
  live <- lapply(live, function(x){
    x/(sum(unlist(live)))
  })
  
  die <- gcn4_tads %>% filter(AD_set != 'AD_positive') %>% select(-sequence, -AD_set)
  die <- map(die, sum)
  die <- lapply(die, function(x){
    x/(sum(unlist(die)))
  })
  
  enrichments <- tibble(live_freq = live %>% unlist(),
                        dipep = names(live),
                        die_freq = die %>% unlist())

  enrichments$enrichment <- log2(enrichments$live_freq/enrichments$die_freq)
  
  #fix -inf value
  which(enrichments$dipep == 'MK')
  enrichments$enrichment[54] <- 0
  
  #split up dipeptide into first and second aa
  enrichments$first_aa <- gsub('.$', '', enrichments$dipep)
  enrichments$second_aa <- gsub('^.', '', enrichments$dipep)
  
  enrichments$first_aa <- factor(enrichments$first_aa, levels = aa_vector)
  enrichments$second_aa <- factor(enrichments$second_aa, levels = aa_vector)
  
  ggplot(enrichments, aes(first_aa, second_aa, fill = enrichment))+
    geom_tile(color = 'grey')+
    scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0)+
    coord_flip()+
    ggtitle('Enrichment',
            subtitle = 'Gcn4')+ #change this line to either hsf1 or gcn4
    xlab('First Amino Acid')+
    ylab('Second Amino Acid')
  ??geom_tile
  