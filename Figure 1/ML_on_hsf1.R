#make motif graphs based on live_percent of presence or absence of motifs
library(tidyverse)
library(caret)
library(pROC)

#load in motif datasets
hsf1 <- readRDS('Figure 1/HSF data/motifs_in_hsf1.rds')

#need to add no basic to hsf1 lib
hsf1$no_basic <- ifelse(grepl('[RK]', hsf1$sequence), 0, 1)


##ML with no_basic feature----
#split into test and train
#remove sequence column, add AD_set back to my datasets
live <- hsf1 %>% 
  select(AD_set, ends_with('match'), no_basic) %>% 
  filter(AD_set == 'AD_positive')

die <- hsf1 %>% 
  select(AD_set, ends_with('match'), no_basic) %>% 
  filter(AD_set == 'AD_negative')

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


#Run ML----
#grid of lamdas to search over
lambda <- 10^seq(-6, -3, length = 20)

# Build the model
ridge <- train(
  AD_set ~., data = train_set, method = "glmnet",
  trControl = trainControl("cv", number = 5),
  tuneGrid = expand.grid(alpha = 0, lambda = lambda)
)

#add auc to this graph

predicted_prob<-predict.train(ridge, test_set,type = "prob") 
myroc<-roc(test_set$AD_set, predicted_prob$AD_positive)
myroc$auc

#coef----
coef <- coef(ridge$finalModel, ridge$bestTune$lambda)
coef <- coef[,1]
coef <- coef[-1]

coef <- tibble(motif = names(coef),
               coef = coef)

coef$xaxis <- 2:12
coef$xaxis <- factor(coef$xaxis,
                     levels = 2:12)
coef$color <- c(rep('tan',10),'red')

ggplot(coef, aes(xaxis, coef,fill = color))+
  geom_col(show.legend = F,color='black')+
  scale_fill_manual(values = c('#ff4c4c','#dcbc67'))+
  theme_bw()+
  xlab('Motifs')+
  ylab('ML Feature Coefficient')+
  ggtitle('HSF1', subtitle = 'ROC AUC: 0.532')
ggsave("Figure 1/Figures/Fig1d_HSF-Regression.tiff",
       height = 3.0, width = 2.55)



#write.csv(coef, 'figure1_ML_hsf1_coeff_plus_basic.csv')

#remove sequence column, add AD_set back to my datasets
##without basic
live <- hsf1 %>% 
  select(AD_set, ends_with('match')) %>% 
  filter(AD_set == 'AD_positive')

die <- hsf1 %>% 
  select(AD_set, ends_with('match')) %>% 
  filter(AD_set == 'AD_negative')

#get train (80% of data) and test (20%) sets
#set seed for reproducible results
set.seed(123)

rows <- sample(1:nrow(live), round(nrow(live) %>% as.numeric() * 0.8))
live_train <- live[rows,]
live_test <- live[-rows,]

#changing to be the same amount of sequences from die and live pool... 
#took too long to train with 80% of data from dying sequence
rows <- sample(1:nrow(die), round(nrow(die) %>% as.numeric() * 0.8))
die_train <- die[rows,]
die_test <- die[-rows,]

train_set <- rbind(live_train, die_train)
test_set <- rbind(live_test, die_test)


##run some ridge or lasso
#grid of lamdas to search over
lambda <- 10^seq(-9, -6, length = 20)

# Build the model
ridge <- train(
  AD_set ~., data = train_set, method = "glmnet",
  trControl = trainControl("cv", number = 5),
  tuneGrid = expand.grid(alpha = 0, lambda = lambda)
)

#add auc to this graph

predicted_prob<-predict.train(ridge, test_set,type = "prob") 
myroc<-roc(test_set$AD_set, predicted_prob$AD_positive)
myroc$auc





coef <- coef(ridge$finalModel, ridge$bestTune$lambda)
coef <- coef[,1]
coef <- coef[-1]

coef <- tibble(motif = names(coef),
               coef = coef)

coef$xaxis <- 2:11
coef$xaxis <- factor(coef$xaxis,
                     levels = 2:11)
coef$color <- c(rep('tan',10))

ggplot(coef, aes(xaxis, coef,fill = color))+
  geom_col(show.legend = F,color='black')+
  scale_fill_manual(values = c('#dcbc67'))+
  theme_bw()+
  xlab('Motifs')+
  ylab('Feature Coefficient')+
  ggtitle('hsf1', subtitle = 'Ridge ROC AUC: 0.515')


#write.csv(coef, 'figure1_ML_hsf1_coeff_without_basic.csv')

