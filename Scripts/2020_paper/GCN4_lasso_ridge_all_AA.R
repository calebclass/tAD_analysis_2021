hahn_tads <- readRDS("Datasets/Datasets/hahn_tads.rds")

#remaking the ridge and lasso but with all amino acids 
head(hahn_tads)

GCN4_all_aa<-hahn_tads[,c('sequence', 'AD_set')]
head(GCN4_all_aa)


GCN4_all_aa$R_count<-0
GCN4_all_aa$H_count<-0
GCN4_all_aa$K_count<-0
GCN4_all_aa$D_count<-0  
GCN4_all_aa$E_count<-0
GCN4_all_aa$S_count<-0  
GCN4_all_aa$T_count<-0
GCN4_all_aa$N_count<-0  
GCN4_all_aa$Q_count<-0
GCN4_all_aa$C_count<-0  
GCN4_all_aa$G_count<-0  
GCN4_all_aa$P_count<-0  
GCN4_all_aa$A_count<-0
GCN4_all_aa$V_count<-0
GCN4_all_aa$I_count<-0
GCN4_all_aa$L_count<-0
GCN4_all_aa$M_count<-0
GCN4_all_aa$F_count<-0
GCN4_all_aa$Y_count<-0
GCN4_all_aa$W_count<-0


GCN4_all_aa$R_count<-str_count(GCN4_all_aa$sequence,"R")
GCN4_all_aa$H_count<-str_count(GCN4_all_aa$sequence,"H")
GCN4_all_aa$K_count<-str_count(GCN4_all_aa$sequence,"K")
GCN4_all_aa$D_count<-str_count(GCN4_all_aa$sequence,"D")  
GCN4_all_aa$E_count<-str_count(GCN4_all_aa$sequence,"E")
GCN4_all_aa$S_count<-str_count(GCN4_all_aa$sequence,"S")  
GCN4_all_aa$T_count<-str_count(GCN4_all_aa$sequence,"T")
GCN4_all_aa$N_count<-str_count(GCN4_all_aa$sequence,"N")  
GCN4_all_aa$Q_count<-str_count(GCN4_all_aa$sequence,"Q")
GCN4_all_aa$C_count<-str_count(GCN4_all_aa$sequence,"C")  
GCN4_all_aa$G_count<-str_count(GCN4_all_aa$sequence,"G")  
GCN4_all_aa$P_count<-str_count(GCN4_all_aa$sequence,"P")  
GCN4_all_aa$A_count<-str_count(GCN4_all_aa$sequence,"A")
GCN4_all_aa$V_count<-str_count(GCN4_all_aa$sequence,"V")
GCN4_all_aa$I_count<-str_count(GCN4_all_aa$sequence,"I")
GCN4_all_aa$L_count<-str_count(GCN4_all_aa$sequence,"L")
GCN4_all_aa$M_count<-str_count(GCN4_all_aa$sequence,"M")
GCN4_all_aa$F_count<-str_count(GCN4_all_aa$sequence,"F")
GCN4_all_aa$Y_count<-str_count(GCN4_all_aa$sequence,"Y")
GCN4_all_aa$W_count<-str_count(GCN4_all_aa$sequence,"W")


nrow(GCN4_all_aa)
Y<-nrow(GCN4_all_aa)

head(GCN4_all_aa)


X<-which(GCN4_all_aa$AD_set=='AD_positive')
LX<-length(X)
live<-GCN4_all_aa[X,]


head(live)
nrow(live)


Y<-which(GCN4_all_aa$AD_set=='AD_negative')
die<-GCN4_all_aa[Y,]



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

plot_var_imp$coefficient <- coef(ridge$finalModel, ridge$bestTune$lambda)[2:21]

plot_var_imp$features <- factor(plot_var_imp$features,
                                levels = plot_var_imp$features[order(plot_var_imp$importance)])
levels(plot_var_imp$features) <- substr(levels(plot_var_imp$features), 1, 1)


#making the 'funnel' graph Dr. erkine wanted with the different AA having colored by class

colnames(GCN4_all_aa)


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
#  ggtitle('Ridge Variable Importance GCN4') +
  theme_bw() +
  theme(legend.title = element_blank())
ggsave("Figures/Fig2_Ridge_GCN4.tiff", height = 3, width = 3,
       units = "in")




rm(list = c('df_varImp', 'plot_var_imp'))

# Get ROC curve AUC results... graph if you dare. it takes a while ----
predicted_prob<-predict.train(ridge, test_set,type = "prob") 
myroc<-roc(test_set$AD_set, predicted_prob$AD_negative)

myroc$auc

#ROC AUC 0.9345

rm(list = c('predicted_prob', 'myroc'))
rm('mystats')

#now starting the LASSO



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
levels(plot_var_imp$features) <- substr(levels(plot_var_imp$features), 1, 1)



#making the variabel importance graph 
plot_var_imp$coefficient <- coef(lasso$finalModel, lasso$bestTune$lambda)[2:21]

plot_var_imp$features <- factor(plot_var_imp$features,
                                levels = plot_var_imp$features[order(plot_var_imp$importance)])

# Get ROC curve AUC results... graph if you dare. it takes a while ----
predicted_prob<-predict.train(lasso, test_set,type = "prob") 
myroc<-roc(test_set$AD_set, predicted_prob$AD_negative)

myroc$auc

#ROC: 0.9347


#AAcolors<-c('red', 'yellow', 'orange', 'blue', 'green', 'pink')

#AAclass<-c(Y_count='Aromatic', W_count= 'Aromatic', F_count= 'Aromatic',
#           I_count= 'Aliphatic', V_count='Aliphatic', L_count= 'Aliphatic', A_count= 'Aliphatic', M_count= 'Aliphatic',
#           R_count= 'Basic', H_count= 'Basic', K_count= 'Basic',
#           D_count= 'Acidic', E_count= 'Acidic',
#           S_count= 'Polar', T_count= 'Polar', N_count= 'Polar', Q_count= 'Polar',
#           C_count= 'Special', G_count = 'Special', P_count= 'Special')

plot_var_imp$AAclass=''


for(i in 1:length(AAclass)){
  plot_var_imp$AAclass[which(plot_var_imp$features== names(AAclass)[i])]<-AAclass[i]
}



#plot_var_imp$AAcolors<-AAcolors




ggplot(plot_var_imp, aes(features, coefficient, fill = AAclass))+
  geom_col()+
  coord_flip()+
  scale_fill_manual(values = AAcolors)+ xlab("") +
#  ggtitle('Lasso Variable Importance GCN4') +
  theme_bw() +
  theme(legend.title = element_blank())
ggsave("Figures/Fig2_Lasso_GCN4.tiff", height = 3, width = 3,
       units = "in")

plot_var_imp_lasso <- plot_var_imp

save(plot_var_imp_ridge, plot_var_imp_lasso, file = "Output/GCN4_regression.rData")






