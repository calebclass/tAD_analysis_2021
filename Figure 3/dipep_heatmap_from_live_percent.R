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

###################################
#shouldn't have to change anything above, 
##get live_percent of each dipep  ##chane datasets gcn4 to hsf1 to run different libraries
#control F and swap gcn4 for hsf1 and vice versa
live <- gcn4_tads %>% filter(AD_set == 'AD_positive') %>% select(-sequence, -AD_set)
live <- map(live, sum)
live %>% names()
live_percent <- tibble(live_count = unlist(live),
                       dipep = names(live))

die <- gcn4_tads %>% filter(AD_set != 'AD_positive') %>% select(-sequence, -AD_set)
die <- map(die, sum)
live_percent$die_count <- unlist(die)

live_percent <- mutate(live_percent,
                       percent = live_count/(live_count +die_count) * 100
)

#split dipep into two columns
live_percent$first_aa <- gsub('.$', '', live_percent$dipep)
live_percent$second_aa <- gsub('^.', '', live_percent$dipep)

live_percent$first_aa <- factor(live_percent$first_aa, levels = aa_vector)
live_percent$second_aa <- factor(live_percent$second_aa, levels = aa_vector)

ggplot(live_percent, aes(first_aa, second_aa, fill = percent))+
  geom_tile(color = 'grey')+
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 3.62)+ #change midpoint to library baseline
  coord_flip()+                                                                     #hsf1 1.15%, gcn4 3.62%
  ggtitle('Live Percent', 
          subtitle = 'Gcn4')+ #change this line to either hsf1 or gcn4
  xlab('First Amino Acid')+
  ylab('Second Amino Acid')


gcn4_tads$AD_set %>% table()

