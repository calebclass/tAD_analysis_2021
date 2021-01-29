#acid-aro-basic position figure 7

library(tidyverse)

pos_lib <- readRDS('gcn4_avg_pos.rds') %>% as_tibble()

#find sequences, with Na in acid, base, or aromatic and remove them
pos_lib <- pos_lib[-which(is.na(pos_lib$acid_pos)),]

pos_lib <- pos_lib[-which(is.na(pos_lib$basic_pos)),]

pos_lib <- pos_lib[-which(is.na(pos_lib$aromatic_pos)),]

#need six bins
  #1 acid-aro-base, 2 aro-acid-base, 3 aro-base-acid, 4 acid-base-aro, 5 base-aro-acid, 6 base-acid-aro

group1 <- pos_lib %>% filter(acid_pos < aromatic_pos, aromatic_pos < basic_pos)
group2 <- pos_lib %>% filter(aromatic_pos < acid_pos, acid_pos < basic_pos)
group3 <- pos_lib %>% filter(aromatic_pos < basic_pos, basic_pos < acid_pos)
group4 <- pos_lib %>% filter(acid_pos < basic_pos, basic_pos < aromatic_pos)
group5 <- pos_lib %>% filter(basic_pos < aromatic_pos, aromatic_pos < acid_pos)
group6 <- pos_lib %>% filter(basic_pos < acid_pos, acid_pos < aromatic_pos)

#1 - 4.3%
#2 - 4.1%
#3 - 1.4%
#4 - 2.6%
#5 - 3.3%
#6 - 5.7%

#just 1 basic.. 1 - 17.0%, 2 - 14.1%, 3 - 7.1%, 4 - 13.4%, 5 - 12.4%, 6 - 18.6%
#just 2 basic.. 1 - 9.1%, 2 - 8.4%, 3 - 4.3%, 4 - 7.8%, 5 - 7.1%, 6 - 10.6%
#just 3 basic.. 1 - 4.3%, 2 - 4.1%, 3 - 1.8%, 4 - 3.4%, 5 - 3.3%, 6 - 6.0%
#just 4 basic.. 1 - 1.7%, 2 - 1.8%, 3 - 0.7%, 4 - 1.4%, 5 - 1.5%, 6 - 2.9%
#just 5 basic..

table(group1$AD_set) %>% (function(x){
  x[2]/(x[1]+x[2])*100
})
table(group2$AD_set) %>% (function(x){
  x[2]/(x[1]+x[2])*100
})
table(group3$AD_set) %>% (function(x){
  x[2]/(x[1]+x[2])*100
})
table(group4$AD_set) %>% (function(x){
  x[2]/(x[1]+x[2])*100
})
table(group5$AD_set) %>% (function(x){
  x[2]/(x[1]+x[2])*100
})
table(group6$AD_set) %>% (function(x){
  x[2]/(x[1]+x[2])*100
})

#look into how count of acid, aromatic, basic affects this trend
pos_lib$acid_count <- str_count(pos_lib$sequence, pattern = '[DE]')
pos_lib$basic_count <- str_count(pos_lib$sequence, pattern = '[RHK]')
pos_lib$aromatic_count <- str_count(pos_lib$sequence, pattern = '[WYF]')



pos_lib <- pos_lib %>% filter(acid_count == 1, aromatic_count == 1)


#reset
hold <- pos_lib
pos_lib <- hold

