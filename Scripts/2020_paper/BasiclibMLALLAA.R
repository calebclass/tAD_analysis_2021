#ridge and lasso on the different ML by AA
head(hahn_tads)

hahn_tads$acidic<-0
hahn_tads$basic<-0
hahn_tads$charge<-0

hahn_tads$acidic<-str_count(hahn_tads$sequence, '[DE]')
hahn_tads$basic<-str_count(hahn_tads$sequence, '[RK]')
hahn_tads$charge<-hahn_tads$basic - hahn_tads$acidic

x<-which(hahn_tads$charge > 0)

hahn_tads[x,]

Basiclib<-hahn_tads[x,c("sequence", "AD_set")]

head(Basiclib)


Basiclib$R_count<-0
Basiclib$H_count<-0
Basiclib$K_count<-0
Basiclib$D_count<-0  
Basiclib$E_count<-0
Basiclib$S_count<-0  
Basiclib$T_count<-0
Basiclib$N_count<-0  
Basiclib$Q_count<-0
Basiclib$C_count<-0  
Basiclib$G_count<-0  
Basiclib$P_count<-0  
Basiclib$A_count<-0
Basiclib$V_count<-0
Basiclib$I_count<-0
Basiclib$L_count<-0
Basiclib$M_count<-0
Basiclib$F_count<-0
Basiclib$Y_count<-0
Basiclib$W_count<-0


Basiclib$R_count<-str_count(Basiclib$sequence,"R")
Basiclib$H_count<-str_count(Basiclib$sequence,"H")
Basiclib$K_count<-str_count(Basiclib$sequence,"K")
Basiclib$D_count<-str_count(Basiclib$sequence,"D")  
Basiclib$E_count<-str_count(Basiclib$sequence,"E")
Basiclib$S_count<-str_count(Basiclib$sequence,"S")  
Basiclib$T_count<-str_count(Basiclib$sequence,"T")
Basiclib$N_count<-str_count(Basiclib$sequence,"N")  
Basiclib$Q_count<-str_count(Basiclib$sequence,"Q")
Basiclib$C_count<-str_count(Basiclib$sequence,"C")  
Basiclib$G_count<-str_count(Basiclib$sequence,"G")  
Basiclib$P_count<-str_count(Basiclib$sequence,"P")  
Basiclib$A_count<-str_count(Basiclib$sequence,"A")
Basiclib$V_count<-str_count(Basiclib$sequence,"V")
Basiclib$I_count<-str_count(Basiclib$sequence,"I")
Basiclib$L_count<-str_count(Basiclib$sequence,"L")
Basiclib$M_count<-str_count(Basiclib$sequence,"M")
Basiclib$F_count<-str_count(Basiclib$sequence,"F")
Basiclib$Y_count<-str_count(Basiclib$sequence,"Y")
Basiclib$W_count<-str_count(Basiclib$sequence,"W")


nrow(Basiclib)
Y<-nrow(Basiclib)

head(Basiclib)

bhold<-Basiclib

head(bhold)

X<-which(bhold$AD_set=='AD_positive')
LX<-length(X)
live<-bhold[X,]

LX/Y*100

head(live)
nrow(live)

Y<-which(bhold$AD_set=='AD_negative')
die<-bhold[Y,]



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


#make variable importance graph
plot_var_imp$coefficient <- coef(ridge$finalModel, ridge$bestTune$lambda)[2:21]

plot_var_imp$features <- factor(plot_var_imp$features,
                                levels = plot_var_imp$features[order(plot_var_imp$importance)])






#order features by importance
plot_var_imp$features <- factor(plot_var_imp$features,
                                levels = plot_var_imp$features[order(plot_var_imp$importance)])

#need to fix the colors 



AAcolors<-c('red', 'yellow', 'orange', 'blue', 'green', 'pink')








AAclass<-c(Y_count='Aromatic', W_count= 'Aromatic', F_count= 'Aromatic',
           I_count= 'Aliphatic', V_count='Aliphatic', L_count= 'Aliphatic', A_count= 'Aliphatic', M_count= 'Aliphatic',
           R_count= 'Basic', H_count= 'Basic', K_count= 'Basic',
           D_count= 'Acidic', E_count= 'Acidic',
           S_count= 'Polar', T_count= 'Polar', N_count= 'Polar', Q_count= 'Polar',
           C_count= 'Special', G_count = 'Special', P_count= 'Special')

plot_var_imp$AAclass=''


for(i in 1:length(AAclass)){
  plot_var_imp$AAclass[which(plot_var_imp$features== names(AAclass)[i])]<-AAclass[i]
}



plot_var_imp$AAcolors<-AAcolors






ggplot(plot_var_imp, aes(features, coefficient, fill = AAclass))+
  geom_col()+
  coord_flip()+
  scale_fill_manual(values = AAcolors)+
  ggtitle('Ridge BasicLib GCN4')




rm(list = c('df_varImp', 'plot_var_imp'))

# Get ROC curve AUC results... graph if you dare. it takes a while ----
predicted_prob<-predict.train(ridge, test_set,type = "prob") 
myroc<-roc(test_set$AD_set, predicted_prob$AD_negative)

myroc$auc



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




#making the variabel importance graph 
plot_var_imp$coefficient <- coef(ridge$finalModel, ridge$bestTune$lambda)[2:21]

plot_var_imp$features <- factor(plot_var_imp$features,
                                levels = plot_var_imp$features[order(plot_var_imp$importance)])






#order features by importance
plot_var_imp$features <- factor(plot_var_imp$features,
                                levels = plot_var_imp$features[order(plot_var_imp$importance)])




AAcolors<-c('red', 'yellow', 'orange', 'blue', 'green', 'pink')








AAclass<-c(Y_count='Aromatic', W_count= 'Aromatic', F_count= 'Aromatic',
           I_count= 'Aliphatic', V_count='Aliphatic', L_count= 'Aliphatic', A_count= 'Aliphatic', M_count= 'Aliphatic',
           R_count= 'Basic', H_count= 'Basic', K_count= 'Basic',
           D_count= 'Acidic', E_count= 'Acidic',
           S_count= 'Polar', T_count= 'Polar', N_count= 'Polar', Q_count= 'Polar',
           C_count= 'Special', G_count = 'Special', P_count= 'Special')

plot_var_imp$AAclass=''


for(i in 1:length(AAclass)){
  plot_var_imp$AAclass[which(plot_var_imp$features== names(AAclass)[i])]<-AAclass[i]
}



plot_var_imp$AAcolors<-AAcolors






ggplot(plot_var_imp, aes(features, coefficient, fill = AAclass))+
  geom_col()+
  coord_flip()+
  scale_fill_manual(values = AAcolors)+
  ggtitle('Lasso BasicLib GCN4')


#clean up environment
rm(list = c('df_varImp', 'plot_var_imp'))

# Get ROC curve AUC results... graph if you dare. it takes a while ----
predicted_prob<-predict.train(lasso, test_set,type = "prob") 
myroc<-roc(test_set$AD_set, predicted_prob$AD_negative)

myroc$auc









