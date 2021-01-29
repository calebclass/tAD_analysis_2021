#this is going to invesitgate sequences which have more basic than acidic residues 
library(stringr)

Hahnhold<-hahn_tads[,c('sequence', 'AD_set')]


Hahnhold$Acidic_count<-0
Hahnhold$Acidic_count<-str_count(Hahnhold$sequence, '[DE]')

Hahnhold$Aromatic_count<-0
Hahnhold$Aromatic_count<-str_count(Hahnhold$sequence, '[WYF]')

Hahnhold$Basic_count<-0
Hahnhold$Basic_count<-str_count(Hahnhold$sequence, '[RHK]')

Hahnhold$Aliphatic_count<-0
Hahnhold$Aliphatic_count<-str_count(Hahnhold$sequence, '[AVILM]')

Hahnhold$Polar_count<-0
Hahnhold$Polar_count<-str_count(Hahnhold$sequence, '[STNQC]')

Hahnhold$Proline_count<-0
Hahnhold$Proline_count<-str_count(Hahnhold$sequence, '[P]')

Hahnhold$Glycine_count<-0
Hahnhold$Glycine_count<-str_count(Hahnhold$sequence, '[G]')

head(Hahnhold)

Hahnhold$charge<-0

Hahnhold$charge<-Hahnhold$Basic_count- Hahnhold$Acidic_count

head(Hahnhold)

#making a data set where the sequences all have a 'positive' charge (more basic than acidic)

bhold<-which(Hahnhold$charge>0)

bhold<-Hahnhold[bhold,]

head(bhold)

nrow(bhold)

bhold<-subset(bhold, select = -c(charge))

head(bhold)




X<-which(bhold$AD_set=='AD_positive')
live<-bhold[X,]

head(live)
nrow(live)

Y<-which(bhold$AD_set=='AD_negative')
die<-bhold[Y,]

head(die)
nrow(die)


#install.packages('dplyr')
library(dplyr)
#install.packages('tidyverse')
library(tidyverse)
library(ggplot2)
library(caret)
library(pROC)


#remove sequence column, add AD_set back to my datasets
live <- live %>% 
  as_tibble() %>% 
  select(-sequence) %>% 
  mutate(AD_set = 'AD_positive')

die <- die %>% 
  as_tibble() %>% 
  select(-sequence) %>% 
  mutate(AD_set = 'AD_negative')

#get train (80% of data) and test (20%) sets
#set seed for reproducible results
set.seed(123)

rows <- sample(1:nrow(live), round(nrow(live) %>% as.numeric() * 0.8))
live_train <- live[rows,]
live_test <- live[-rows,]


#changing to be the same amount of sequences from die and live pool... 
#took too long to train with 80% of data from dying sequence
rows <- sample(1:nrow(die), round(nrow(live) %>% as.numeric() * 0.8))
die_train <- die[rows,]
die_test <- die[-rows,]

train_set <- rbind(live_train, die_train)
test_set <- rbind(live_test, die_test)


#remove variables to clean up environment at this point
rm(list = c('live_train', 'live_test',
            'die_train', 'die_test', 'rows'))

#train a model
train_set$AD_set <- factor(train_set$AD_set)

#tricky tricky can't have a column called NA with the following train() function
colnames(train_set)[which(colnames(train_set) == 'NA')] <- 'nA'
colnames(test_set)[which(colnames(test_set) == 'NA')] <- 'nA'

lambda <- 10^seq(-3, 3, length = 100)


# Build the model
ridge <- train(
  AD_set ~., data = train_set, method = "glmnet",
  trControl = trainControl("cv", number = 10),
  tuneGrid = expand.grid(alpha = 0, lambda = lambda)
)


# Plot Variable importance----
df_varImp <- varImp(ridge)
plot_var_imp <- tibble(features = row.names(df_varImp$importance),
                       importance = df_varImp$importance %>% unlist())

#order features by importance
plot_var_imp$features <- factor(plot_var_imp$features,
                                levels = plot_var_imp$features[order(plot_var_imp$importance)])

groupcolors<-c(Basic_count = 'blue', Acidic_count = 'red', Aromatic_count = 'orange', Aliphatic_count = 'yellow', Glycine_count = 'green', Proline_count = 'pink', Polar_count= 'purple', charge = 'light green')


ggplot(plot_var_imp, aes(features, importance, fill = features))+
  geom_col()+
  ggtitle('Ridge Basic Sequences NO charge')+
  coord_flip()+
  scale_fill_manual(values = groupcolors)



rm(list = c('df_varImp', 'plot_var_imp'))

# Get ROC curve AUC results... graph if you dare. it takes a while ----
predicted_prob<-predict.train(ridge, test_set,type = "prob") 
myroc<-roc(test_set$AD_set, predicted_prob$AD_negative)

myroc$auc

#ROC for positive charge is 0.7695

#clean up environment
rm(list = c('predicted_prob', 'myroc'))
rm('mystats')

lasso <- train(
  AD_set ~., data = train_set, method = "glmnet",
  trControl = trainControl("cv", number = 10),
  tuneGrid = expand.grid(alpha = 1, lambda = lambda)
)


# Plot Variable importance----
df_varImp <- varImp(lasso)
plot_var_imp <- tibble(features = row.names(df_varImp$importance),
                       importance = df_varImp$importance %>% unlist())

#order features by importance
plot_var_imp$features <- factor(plot_var_imp$features,
                                levels = plot_var_imp$features[order(plot_var_imp$importance)])

groupcolors<-c(Basic_count = 'blue', Acidic_count = 'red', Aromatic_count = 'orange', Aliphatic_count = 'yellow', Glycine_count = 'green', Proline_count = 'pink', Polar_count= 'purple', charge = 'light green')

ggplot(plot_var_imp, aes(features, importance, fill = features))+
  geom_col()+
  ggtitle('Lasso Variable Importance Basic NO Charge')+
  coord_flip()+
  scale_fill_manual(values = groupcolors)

#clean up environment
rm(list = c('df_varImp', 'plot_var_imp'))

# Get ROC curve AUC results... graph if you dare. it takes a while ----
predicted_prob<-predict.train(lasso, test_set,type = "prob") 
myroc<-roc(test_set$AD_set, predicted_prob$AD_negative)

myroc$auc

#ROC: 0.7702







