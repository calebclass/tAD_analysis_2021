#Build motif datasets
library(tidyverse)

#recruitment motifs----
#p53 - [ALVIMWYF]..[ALVIMWYF][ALVIMWYF]
#RelA_2 - [ALVIMWYF][ALVIMWYF]..[ALVIMWYF]..[ALVIMWYF]
#CREBZF - D[VILM][VILM][RKDEQNHSTYC][RKDEQNHSTYC][VILM][VILFWYM]
#AR - F..LF
#ANAC013 - [DE].{1,2}[YF].{1,4}[DE]L
#EKLF - [RKDEQNHSTYC][ALVIMWYF][ALVIMWYF]..[ALVIMWYF]..[RKDEQNHSTYC][RKDEQNHSTYC]
#WxxLF - W..LF
#stringent_9aa - [MDENQSTYG][^KRHCGP][ILVFWM][^KRHCGP][^CGP][^KRHCGP][ILVFWM][ILVFWMAY][^KRHC]
#moderate_9aa - [MDENQSTYG][^KRHCGP][ILVFWM][^KRHCGP][^CGP][^CGP][ILVFWM][^CGP][^CGP]
#lenient_9aa - [MDENQSTYCPGA].[ILVFWMAY][^KRHCGP][^CGP][^CGP][ILVFWMAY]..

#load gcn4 and hsf1 datasets----
hahn_tads <- readRDS('Datasets/Datasets/hahn_tads.rds') %>% select(sequence, AD_set) %>% as_tibble()

#I filtered for only hsf1 tads 5 aa or longer, I also use ad_set instead of binary stop
hsf1_tads <- readRDS('Datasets/Datasets/random_library.rds') %>% select(sequence, binary_stop) %>% 
  mutate(
    AD_set = ifelse(binary_stop == 'live', 'AD_positive', 'AD_negative'),
    len = nchar(sequence)
  ) %>% filter(len > 4) %>%  select(sequence, AD_set) %>% as_tibble()

#need a motif list to loop through----
motif_list <- tibble(motif = c('p53','RelA_2','CREBZF','AR','ANAC013','EKLF',
                               'WxxLF','stringent_9aa','moderate_9aa','lenient_9aa'),
                     regex = c('[ALVIMWYF]..[ALVIMWYF][ALVIMWYF]',
                               '[ALVIMWYF][ALVIMWYF]..[ALVIMWYF]..[ALVIMWYF]',
                               'D[VILM][VILM][RKDEQNHSTYC][RKDEQNHSTYC][VILM][VILFWYM]',
                               'F..LF','[DE].{1,2}[YF].{1,4}[DE]L',
                               '[RKDEQNHSTYC][ALVIMWYF][ALVIMWYF]..[ALVIMWYF]..[RKDEQNHSTYC][RKDEQNHSTYC]',
                               'W..LF',
                               '[MDENQSTYG][^KRHCGP][ILVFWM][^KRHCGP][^CGP][^KRHCGP][ILVFWM][ILVFWMAY][^KRHC]',
                               '[MDENQSTYG][^KRHCGP][ILVFWM][^KRHCGP][^CGP][^CGP][ILVFWM][^CGP][^CGP]',
                               '[MDENQSTYCPGA].[ILVFWMAY][^KRHCGP][^CGP][^CGP][ILVFWMAY]..'))

#add motif counts to gcn4----
for(i in 1:nrow(motif_list)){
  #first step is check if motif is in a sequence
  #this cuts down on the number of sequences for the look_through regex
  hahn_tads$hold <- grepl(motif_list$regex[i], hahn_tads$sequence)
  
  #get the rows for when sequences have a motif
  #use rows for adding counts back after look_through regex
  #small_set is sequences I need to apply the look_through regex
  rows <- which(hahn_tads$hold == T)
  small_set <- hahn_tads %>% filter(hold == T)
  
  #need look_through to find overlapping matches
  look_through <- paste('(?=', motif_list$regex[i], ')', sep = '')
  
  #this will find counts of regex in each sequence
  small_set$hold <- lapply(small_set$sequence, function(x){
    gregexpr(look_through, x, perl = T) %>% 
      unlist() %>% 
      length()
  }) %>% unlist()
  
  #set all counts to 0, and then add actual counts from small_set
  hahn_tads$hold <- 0
  hahn_tads$hold[rows] <- small_set$hold
  
  #two columns one is count, one is presence or absence
  hahn_tads$hold2 <- ifelse(hahn_tads$hold > 0, 1, 0)
  
  #this peice changes column name from hold to the actual motif name
  colnames(hahn_tads)[which(colnames(hahn_tads)=='hold')] <- motif_list$motif[i]
  colnames(hahn_tads)[which(colnames(hahn_tads)=='hold2')] <- paste(motif_list$motif[i], 
                                                                    '_match', sep = '')
}

#add motif counts to HSF1----
for(i in 1:nrow(motif_list)){
  #first step is check if motif is in a sequence
  #this cuts down on the number of sequences for the look_through regex
  hsf1_tads$hold <- grepl(motif_list$regex[i], hsf1_tads$sequence)
  
  #get the rows for when sequences have a motif
  #use rows for adding counts back after look_through regex
  #small_set is sequences I need to apply the look_through regex
  rows <- which(hsf1_tads$hold == T)
  small_set <- hsf1_tads %>% filter(hold == T)
  
  #need look_through to find overlapping matches
  look_through <- paste('(?=', motif_list$regex[i], ')', sep = '')
  
  #this will find counts of regex in each sequence
  small_set$hold <- lapply(small_set$sequence, function(x){
    gregexpr(look_through, x, perl = T) %>% 
      unlist() %>% 
      length()
  }) %>% unlist()
  
  #set all counts to 0, and then add actual counts from small_set
  hsf1_tads$hold <- 0
  hsf1_tads$hold[rows] <- small_set$hold
  
  #two columns one is count, one is presence or absence
  hsf1_tads$hold2 <- ifelse(hsf1_tads$hold > 0, 1, 0)
  
  #this peice changes column name from hold to the actual motif name
  colnames(hsf1_tads)[which(colnames(hsf1_tads)=='hold')] <- motif_list$motif[i]
  colnames(hsf1_tads)[which(colnames(hsf1_tads)=='hold2')] <- paste(motif_list$motif[i], 
                                                                    '_match', sep = '')
}

#save data to use later----
#saveRDS(hahn_tads,'motifs_in_gcn4.rds')
#saveRDS(hsf1_tads, 'motifs_in_hsf1.rds')
#saveRDS(motif_list,'motif_list.rds')

