random_library <- readRDS("Datasets/Datasets/random_library.rds")

#HSF ridge and LASSO
head(random_library)


HSF_all_aa<-random_library[,c('sequence', 'binary_stop')]
head(HSF_all_aa)



HSF_all_aa$R_count<-0
HSF_all_aa$H_count<-0
HSF_all_aa$K_count<-0
HSF_all_aa$D_count<-0  
HSF_all_aa$E_count<-0
HSF_all_aa$S_count<-0  
HSF_all_aa$T_count<-0
HSF_all_aa$N_count<-0  
HSF_all_aa$Q_count<-0
HSF_all_aa$C_count<-0  
HSF_all_aa$G_count<-0  
HSF_all_aa$P_count<-0  
HSF_all_aa$A_count<-0
HSF_all_aa$V_count<-0
HSF_all_aa$I_count<-0
HSF_all_aa$L_count<-0
HSF_all_aa$M_count<-0
HSF_all_aa$F_count<-0
HSF_all_aa$Y_count<-0
HSF_all_aa$W_count<-0


HSF_all_aa$R_count<-str_count(HSF_all_aa$sequence,"R")
HSF_all_aa$H_count<-str_count(HSF_all_aa$sequence,"H")
HSF_all_aa$K_count<-str_count(HSF_all_aa$sequence,"K")
HSF_all_aa$D_count<-str_count(HSF_all_aa$sequence,"D")  
HSF_all_aa$E_count<-str_count(HSF_all_aa$sequence,"E")
HSF_all_aa$S_count<-str_count(HSF_all_aa$sequence,"S")  
HSF_all_aa$T_count<-str_count(HSF_all_aa$sequence,"T")
HSF_all_aa$N_count<-str_count(HSF_all_aa$sequence,"N")  
HSF_all_aa$Q_count<-str_count(HSF_all_aa$sequence,"Q")
HSF_all_aa$C_count<-str_count(HSF_all_aa$sequence,"C")  
HSF_all_aa$G_count<-str_count(HSF_all_aa$sequence,"G")  
HSF_all_aa$P_count<-str_count(HSF_all_aa$sequence,"P")  
HSF_all_aa$A_count<-str_count(HSF_all_aa$sequence,"A")
HSF_all_aa$V_count<-str_count(HSF_all_aa$sequence,"V")
HSF_all_aa$I_count<-str_count(HSF_all_aa$sequence,"I")
HSF_all_aa$L_count<-str_count(HSF_all_aa$sequence,"L")
HSF_all_aa$M_count<-str_count(HSF_all_aa$sequence,"M")
HSF_all_aa$F_count<-str_count(HSF_all_aa$sequence,"F")
HSF_all_aa$Y_count<-str_count(HSF_all_aa$sequence,"Y")
HSF_all_aa$W_count<-str_count(HSF_all_aa$sequence,"W")


nrow(HSF_all_aa)
Y<-nrow(HSF_all_aa)

head(HSF_all_aa)


X<-which(HSF_all_aa$binary_stop=='live')
LX<-length(X)
live<-HSF_all_aa[X,]


head(live)
nrow(live)


Y<-which(HSF_all_aa$binary_stop=='die')
die<-HSF_all_aa[Y,]




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
  mutate(binary_stop = 'live')

die <- die %>% 
  as_tibble() %>% 
  select(-sequence) %>% 
  mutate(binary_stop = 'die')


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
train_set$binary_stop <- factor(train_set$binary_stop)

#tricky tricky can't have a column called NA with the following train() function
colnames(train_set)[which(colnames(train_set) == 'NA')] <- 'nA'
colnames(test_set)[which(colnames(test_set) == 'NA')] <- 'nA'

lambda <- 10^seq(-3, 3, length = 100)


# Build the model
ridge <- train(
  binary_stop ~., data = train_set, method = "glmnet",
  trControl = trainControl("cv", number = 10),
  tuneGrid = expand.grid(alpha = 0, lambda = lambda)
)


# Plot Variable importance----
df_varImp <- varImp(ridge)
plot_var_imp <- tibble(features = row.names(df_varImp$importance),
                       importance = df_varImp$importance %>% unlist())

plot_var_imp$coefficient <- coef(ridge$finalModel, ridge$bestTune$lambda)[2:21]

plot_var_imp$features <- factor(plot_var_imp$features,
                                levels = plot_var_imp$features[order(plot_var_imp$importance)])
levels(plot_var_imp$features) <- substr(levels(plot_var_imp$features), 1, 1)

predicted_prob<-predict.train(ridge, test_set,type = "prob") 
myroc<-roc(test_set$binary_stop, predicted_prob$die)

myroc$auc







#making the 'funnel' graph Dr. erkine wanted with the different AA having colored by class

AAcolors<-c('red', 'yellow', 'orange', 'blue', 'green', 'pink')



AAclass<-c(Y='Aromatic', W= 'Aromatic', F= 'Aromatic',
           I= 'Aliphatic', V='Aliphatic', L= 'Aliphatic', A= 'Aliphatic', M= 'Aliphatic',
           R= 'Basic', H= 'Basic', K= 'Basic',
           D= 'Acidic', E= 'Acidic',
           S= 'Polar', T= 'Polar', N= 'Polar', Q= 'Polar',
           C= 'Special', G = 'Special', P= 'Special')


plot_var_imp$AAclass=''


for(i in 1:length(AAclass)){
  plot_var_imp$AAclass[which(plot_var_imp$features== names(AAclass)[i])]<-AAclass[i]
}



#plot_var_imp$AAcolors<-AAcolors



plot_var_imp_ridge <- plot_var_imp


ggplot(plot_var_imp, aes(features, coefficient, fill = AAclass))+
  geom_col()+
  coord_flip()+
  scale_fill_manual(values = AAcolors)+ xlab("") +
#  ggtitle('Ridge Variable Importance HSF') +
  theme_bw() +
  theme(legend.title = element_blank())
ggsave("Figures/Fig2_Ridge_HSF.tiff", height = 3, width = 3,
       units = "in")



lasso <- train(
  binary_stop ~., data = train_set, method = "glmnet",
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
levels(plot_var_imp$features) <- substr(levels(plot_var_imp$features), 1, 1)



#making the variabel importance graph 
plot_var_imp$coefficient <- coef(lasso$finalModel, lasso$bestTune$lambda)[2:21]

plot_var_imp$features <- factor(plot_var_imp$features,
                                levels = plot_var_imp$features[order(plot_var_imp$importance)])


predicted_prob<-predict.train(lasso, test_set,type = "prob") 
myroc<-roc(test_set$binary_stop, predicted_prob$die)

myroc$auc



plot_var_imp$AAclass=''


for(i in 1:length(AAclass)){
  plot_var_imp$AAclass[which(plot_var_imp$features== names(AAclass)[i])]<-AAclass[i]
}



#plot_var_imp$AAcolors<-AAcolors






ggplot(plot_var_imp, aes(features, coefficient, fill = AAclass))+
  geom_col()+
  coord_flip()+
  scale_fill_manual(values = AAcolors)+ xlab("") +
#  ggtitle('Lasso Variable Importance HSF') +
  theme_bw() +
  theme(legend.title = element_blank())
ggsave("Figures/Fig2_Lasso_HSF.tiff", height = 3, width = 3,
       units = "in")


plot_var_imp_lasso <- plot_var_imp


save(plot_var_imp_ridge, plot_var_imp_lasso, file = "Output/HSF_regression.rData")






