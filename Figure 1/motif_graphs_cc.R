#make motif graphs based on live_percent of presence or absence of motifs
library(tidyverse)

#load in motif datasets
gcn4 <- readRDS('Figure 1/motifs_in_gcn4.rds')
hsf1 <- readRDS('Figure 1/HSF data/motifs_in_hsf1.rds')
motif_list <- readRDS('Figure 1/motif_list.rds')

#get live% for presence of motifs in gcn4----
hold <- gcn4 %>% select(AD_set, ends_with('_match'))

motif_list$gcn4_live_count <- 0
motif_list$gcn4_die_count <- 0

for(i in 1:nrow(motif_list)){
  rows <- which(hold[,i+1] == 1)
  live_rows <- hold[rows,] %>% filter(AD_set == 'AD_positive')
  motif_list$gcn4_live_count[i] <- nrow(live_rows)
  motif_list$gcn4_die_count[i] <- length(rows) - nrow(live_rows)
}

motif_list$gcn4_live_percent <- motif_list$gcn4_live_count/
  (motif_list$gcn4_die_count + motif_list$gcn4_live_count) * 100

#get live% for presence of motifs in hsf1----
hold <- hsf1 %>% select(AD_set, ends_with('_match'))

motif_list$hsf1_live_count <- 0
motif_list$hsf1_die_count <- 0

for(i in 1:nrow(motif_list)){
  rows <- which(hold[,i+1] == 1)
  live_rows <- hold[rows,] %>% filter(AD_set == 'AD_positive')
  motif_list$hsf1_live_count[i] <- nrow(live_rows)
  motif_list$hsf1_die_count[i] <- length(rows) - nrow(live_rows)
}

motif_list$hsf1_live_percent <- motif_list$hsf1_live_count/
  (motif_list$hsf1_die_count + motif_list$hsf1_live_count) * 100

#add no basic and library baseline to motif_list----
motif_list[11,1] <- 'library_baseline'
motif_list[11,2] <- ''
motif_list[11,3] <- gcn4 %>% filter(AD_set == 'AD_positive') %>% nrow()
motif_list[11,4] <- gcn4 %>% filter(AD_set == 'AD_negative') %>% nrow()
motif_list[11,5] <- motif_list[11,3]/(motif_list[11,4] + motif_list[11,3]) * 100

motif_list[11,6] <- hsf1 %>% filter(AD_set == 'AD_positive') %>% nrow()
motif_list[11,7] <- hsf1 %>% filter(AD_set == 'AD_negative') %>% nrow()
motif_list[11,8] <- motif_list[11,6]/(motif_list[11,7] + motif_list[11,6]) * 100

motif_list[12,1] <- 'no_basic'
motif_list[12,2] <- '[^RHK] in sequence'

#have to find the no basic sequences
rows <- grep('[RHK]', gcn4$sequence)
hold <- gcn4[-rows,]
motif_list[12,3] <- hold %>% filter(AD_set == 'AD_positive') %>% nrow()
motif_list[12,4] <- hold %>% filter(AD_set == 'AD_negative') %>% nrow()
motif_list[12,5] <- motif_list[12,3]/(motif_list[12,4] + motif_list[12,3]) * 100

#find no basic in hsf1
rows <- grep('[RHK]', hsf1$sequence)
hold <- hsf1[-rows,]
motif_list[12,6] <- hold %>% filter(AD_set == 'AD_positive') %>% nrow()
motif_list[12,7] <- hold %>% filter(AD_set == 'AD_negative') %>% nrow()
motif_list[12,8] <- motif_list[12,6]/(motif_list[12,7] + motif_list[12,6]) * 100

#saveRDS(motif_list,'motif_list_with_live%.rds')

################11/27 work making graphs look nice----
#set the order of the motifs
#get manual color scheme -- library baseline - green, tan for motifs, red for no basic
#need a ML graph too, will need to find script for that

#set order of motifs
motif_list$motif <- factor(motif_list$motif,
                           levels = c('library_baseline', 'p53','RelA_2',
                                      'CREBZF', 'AR', 'ANAC013', 'EKLF',
                                      'WxxLF', 'stringent_9aa', 'moderate_9aa',
                                      'lenient_9aa','no_basic'))

motif_list$color <- c(rep('tan',10),'green','red')
motif_list$xaxis <- c(2:11,1,12)
motif_list$xaxis <- factor(motif_list$xaxis,
                           levels = c(1:12))

#graphing live%----
ggplot(motif_list, aes(xaxis, hsf1_live_percent,fill = color))+
  geom_col(show.legend = F,color='black')+
  scale_fill_manual(values = c('#2ba443','#ff4c4c','#dcbc67'))+
  theme_bw()+
  xlab('Motifs')+
  ylab('Percent Functional')+
  ggtitle('HSF1')
ggsave("Figures/new_Fig1b_HSF1.tiff",
       height = 2.8, width = 2.4)

ggplot(motif_list, aes(xaxis, gcn4_live_percent,fill = color))+
  geom_col(show.legend = F,color='black')+
  scale_fill_manual(values = c('#2ba443','#ff4c4c','#dcbc67'))+
  theme_bw()+
  xlab('Motifs')+
  ylab('Percent Functional')+
  ggtitle('Gcn4')
ggsave("Figures/new_Fig1b_Gcn4.tiff",
       height = 2.8, width = 2.5)

#make it similar to excel graphs

################11/27 ML under here and then graph results-----
##need to add no basic to gcn4 lib
colnames(gcn4)

gcn4$no_basic <- ifelse(grepl('[RK]', gcn4$sequence), 0, 1)

#split into test and train
######## THIS IS ON PRESENCE OF MOTIF
#remove sequence column, add AD_set back to my datasets
live <- gcn4 %>% 
  select(AD_set, ends_with('match'), no_basic) %>% 
  filter(AD_set == 'AD_positive')

die <- gcn4 %>% 
  select(AD_set, ends_with('match'), no_basic) %>% 
  filter(AD_set == 'AD_negative')

######### THIS IS ON COUNT OF MOTIF
colnames(gcn4)
live <- gcn4[,c(2,3,5,7,9,11,13,15,17,19,21,23)] %>% filter(AD_set == 'AD_positive')
die <- gcn4[,c(2,3,5,7,9,11,13,15,17,19,21,23)] %>% filter(AD_set == 'AD_negative')

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

library(caret)
##run some ridge or lasso
#grid of lamdas to search over
lambda <- 10^seq(-3, 3, length = 100)

# Build the model
ridge <- train(
  AD_set ~., data = train_set, method = "glmnet",
  trControl = trainControl("cv", number = 5),
  tuneGrid = expand.grid(alpha = 0, lambda = lambda)
)

df_varImp <- varImp(ridge)
coefs <- coef(ridge$finalModel, ridge$bestTune$lambda)
plot_var_imp <- tibble(features = row.names(df_varImp$importance),
                       coefs = coefs[-1],
                       importance = df_varImp$importance %>% unlist())

##########use this one for presencr
plot_var_imp$features <- factor(plot_var_imp$features,
                           levels = c('p53','RelA_2',
                                      'CREBZF', 'AR', 'ANAC013', 'EKLF',
                                      'WxxLF', 'stringent_9aa', 'moderate_9aa',
                                      'lenient_9aa', 'no_basic'))

plot_var_imp$colors <- c(rep('tan',10),'red')

#add auc to this graph
library(pROC)
predicted_prob<-predict.train(ridge, test_set,type = "prob") 
myroc<-roc(test_set$AD_set, predicted_prob$AD_negative)
myroc$auc

plot_var_imp$xaxis <- c(2:12)
plot_var_imp$xaxis <- factor(plot_var_imp$xaxis,
                             levels = c(2:12))

#graph feature importance
ggplot(plot_var_imp, aes(xaxis, coefs,fill = colors))+
  geom_col(show.legend = F,color='black')+
  geom_hline(yintercept = 0) + 
  scale_fill_manual(values = c('#ff4c4c','#dcbc67'))+
  theme_bw()+
  xlab('Motifs')+
  ylab('ML Feature Coefficient')+
  ggtitle('Gcn4', subtitle = 'ROC AUC: 0.732')
ggsave("Figures/new_Fig1c_Gcn4-Regression.tiff",
       height = 3.0, width = 2.6)







# AUC without 12



ridge3 <- train(
  AD_set ~., data = train_set[,-12], method = "glmnet",
  trControl = trainControl("cv", number = 5),
  tuneGrid = expand.grid(alpha = 0, lambda = lambda)
)
predicted_prob<-predict.train(ridge3, test_set,type = "prob") 
myroc<-roc(test_set$AD_set, predicted_prob$AD_negative)
myroc$auc

# AUC = 0.689











