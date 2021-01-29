#Split constant test and train set
library(tidyverse)
gcn4_tads <- readRDS('datasets/hahn_tads.rds') %>% as_tibble()

set.seed(1234)
#need a set that is 80% of lives and 80% of die, and a set that is equal live and die
live_set <- gcn4_tads %>% filter(AD_set == 'AD_positive')
die_set <- gcn4_tads %>% filter(AD_set == 'AD_negative')

rows <- sample(1:nrow(live_set), round(0.80 * nrow(live_set)))
live_train <- live_set[rows,]
live_test <- live_set[-rows,]

rows <- sample(1:nrow(die_set), round(0.80 * nrow(die_set)))
die_train <- die_set[rows,]
die_test <- die_set[-rows,]

gcn4_train_set <- rbind(live_train,die_train)
gcn4_test_set <- rbind(live_test,die_test)

##make a balanced train set
rows <- sample(1:nrow(die_set), round(0.80 * nrow(live_set)))
balanced_die_train <- die_set[rows,]
balanced_die_test <- die_set[-rows,]

balanced_gcn4_train_set <- rbind(live_train,balanced_die_train)
balanced_gcn4_test_set <- rbind(live_test, balanced_die_test)

##make hsf1 train and test sets----
hsf1_tads <- readRDS('datasets/random_library.rds') %>% as_tibble() %>% 
  mutate(
    AD_set = ifelse(binary_stop == 'live', 'AD_positive','AD_negative'),
    len = nchar(sequence)
  ) %>% filter(len > 4) %>% select(-len)

#need a set that is 80% of lives and 80% of die, and a set that is equal live and die
live_set <- hsf1_tads %>% filter(AD_set == 'AD_positive')
die_set <- hsf1_tads %>% filter(AD_set == 'AD_negative')

rows <- sample(1:nrow(live_set), round(0.80 * nrow(live_set)))
live_train <- live_set[rows,]
live_test <- live_set[-rows,]

rows <- sample(1:nrow(die_set), round(0.80 * nrow(die_set)))
die_train <- die_set[rows,]
die_test <- die_set[-rows,]

hsf1_train_set <- rbind(live_train,die_train)
hsf1_test_set <- rbind(live_test,die_test)

##make a balanced train set
rows <- sample(1:nrow(die_set), round(0.80 * nrow(live_set)))
balanced_die_train <- die_set[rows,]
balanced_die_test <- die_set[-rows,]

balanced_hsf1_train_set <- rbind(live_train,balanced_die_train)
balanced_hsf1_test_set <- rbind(live_test, balanced_die_test)

#save train and test sets
#saveRDS(gcn4_train_set,'gcn4_train_set.rds')
#saveRDS(gcn4_test_set,'gcn4_test_set.rds')
#saveRDS(balanced_gcn4_train_set,'balanced_gcn4_train_set.rds')
#saveRDS(balanced_gcn4_test_set,'balanced_gcn4_test_set.rds')

#saveRDS(hsf1_train_set,'hsf1_train_set.rds')
#saveRDS(hsf1_test_set,'hsf1_test_set.rds')
#saveRDS(balanced_hsf1_train_set,'balanced_hsf1_train_set.rds')
#saveRDS(balanced_hsf1_test_set,'balanced_hsf1_test_set.rds')