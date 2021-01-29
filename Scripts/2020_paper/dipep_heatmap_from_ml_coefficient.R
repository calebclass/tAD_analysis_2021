##Script to build dipeptide heatmap, 
#data on ridge or lasso coeffecients

#load packages
library(tidyverse)
library(caret)

# Load in tad datasets----
gcn4_tads <- readRDS('train_test_split_9_19/gcn4_train_set.rds') %>% 
  select(sequence, AD_set) %>% 
  as_tibble()

hsf1_tads <- readRDS('train_test_split_9_19/balanced_hsf1_train_set.rds') %>% 
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

# Remove sequence column and train on ADset ~ .----
gcn4_train <- gcn4_tads %>% select(-sequence)
hsf1_train <- hsf1_tads %>% select(-sequence)

##fix NA column to be na
colnames(gcn4_train)['NA'] <- 'na'
colnames(hsf1_train)['NA'] <- 'na'


#####################################################

# Rerun this bottom part to generate the 4 different graphs
    #change data = ... to the dataset you want to train
    #change alpha = ... to 0 for ridge and to 1 for lasso

# Build the model----
lambda <- 10^seq(-3, 3, length = 100)

ridge <- train(
  AD_set ~., data = hsf1_train, method = "glmnet",   #change data = ... on this line
  trControl = trainControl("cv", number = 5),
  tuneGrid = expand.grid(alpha = 1, lambda = lambda)  #change alpha = ... on this line
)

# Get dipeptide coefficients from the model----
#pretty messy but it takes some reordering to get into a heatmap
df_coef <- coef(ridge$finalModel, ridge$bestTune$lambda)
df_coef <- df_coef[,1]
df_coef <- tibble(features = names(df_coef),coef = df_coef) #tibble with features, and coef

#remove unneeded top row, and fix NA column
df_coef <- df_coef[-1,]
df_coef$features[which(df_coef$features == '`NA`')] <- 'NA'

#split features into their aa parts
df_coef$first_aa = ''
df_coef$second_aa = ''
for(i in 1:nrow(df_coef)){
  split <- strsplit(df_coef$features[i] %>% as.character(), '')[[1]]
  df_coef$first_aa[i] <- split[1]
  df_coef$second_aa[i] <- split[2]
}

#need these as factors
df_coef$first_aa <- factor(df_coef$first_aa, levels = aa_vector)
df_coef$second_aa <- factor(df_coef$second_aa, levels = aa_vector)

# Heatmap----
ggplot(df_coef, aes(first_aa, second_aa, fill = coef))+
  geom_tile()+
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red')+
  coord_flip()+
  ggtitle('Ridge Feature Coefficients', #change this line to Ridge or Lasso
          subtitle = 'Hsf1')+ #change this line to either hsf1 or gcn4
  xlab('First Amino Acid')+
  ylab('Second Amino Acid')

