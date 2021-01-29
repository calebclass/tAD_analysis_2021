#acid-aro-basic position figure 7

library(tidyverse)

pos_lib <- readRDS('gcn4_avg_pos.rds') %>% as_tibble()

#find sequences, with Na in acid, base, or aromatic and make them 0
pos_lib$acid_pos[which(is.na(pos_lib$acid_pos))] <- 0
pos_lib$basic_pos[which(is.na(pos_lib$basic_pos))] <- 0
pos_lib$aromatic_pos[which(is.na(pos_lib$aromatic_pos))] <- 0

#look into how count of acid, aromatic, basic affects this trend
pos_lib$acid_count <- str_count(pos_lib$sequence, pattern = '[DE]')
pos_lib$basic_count <- str_count(pos_lib$sequence, pattern = '[RHK]')
pos_lib$aromatic_count <- str_count(pos_lib$sequence, pattern = '[WYF]')

######
base <- 0
acid <- 0
aro <- 0
######
#need six bins
#1 acid-aro-base, 2 aro-acid-base, 3 aro-base-acid, 4 acid-base-aro, 5 base-aro-acid, 6 base-acid-aro

find_library_live <- function(pos_lib, base = 1, aro = 1, acid = 1){
pos_lib <- pos_lib %>% filter(basic_count == base, aromatic_count == aro,
                              acid_count == acid)

group1 <- pos_lib %>% filter(acid_pos < aromatic_pos, aromatic_pos < basic_pos)
group2 <- pos_lib %>% filter(aromatic_pos < acid_pos, acid_pos < basic_pos)
group3 <- pos_lib %>% filter(aromatic_pos < basic_pos, basic_pos < acid_pos)
group4 <- pos_lib %>% filter(acid_pos < basic_pos, basic_pos < aromatic_pos)
group5 <- pos_lib %>% filter(basic_pos < aromatic_pos, aromatic_pos < acid_pos)
group6 <- pos_lib %>% filter(basic_pos < acid_pos, acid_pos < aromatic_pos)

hold <- tibble(group = paste(rep('group',6), 1:6, sep = '_'), live_percent = 0)

hold[1,2] <- table(group1$AD_set) %>% (function(x){x[2]/(x[1]+x[2])*100})
hold[2,2] <- table(group2$AD_set) %>% (function(x){x[2]/(x[1]+x[2])*100})
hold[3,2] <- table(group3$AD_set) %>% (function(x){x[2]/(x[1]+x[2])*100})
hold[4,2] <- table(group4$AD_set) %>% (function(x){x[2]/(x[1]+x[2])*100})
hold[5,2] <- table(group5$AD_set) %>% (function(x){x[2]/(x[1]+x[2])*100})
hold[6,2] <- table(group6$AD_set) %>% (function(x){x[2]/(x[1]+x[2])*100})

return(hold)
}

#base,aro,acid
find_library_live(pos_lib,7,7,7) %>% t()

table(pos_lib$aromatic_count) #0-11
table(pos_lib$acid_count) #0-15
table(pos_lib$basic_count) #0-15

#11 * 15 * 15 -- 2475

##make a big loop to capture all these live% across group and 
#1 acid-aro-base, 2 aro-acid-base
#3 aro-base-acid, 4 acid-base-aro
#5 base-aro-acid, 6 base-acid-aro
