hahn_tads <- readRDS("Datasets/Datasets/hahn_tads.rds")

#remaking the ridge and lasso but with all amino acids 
head(hahn_tads)

GCN4_all_aa<-hahn_tads[,c('sequence', 'AD_set')]
head(GCN4_all_aa)


motif_search <- c("R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W",
                  "[DE][DE]", "[DE].[DE]", "[DE]..[DE]", "[DE]...[DE]", "[DE]....[DE]",
                  "[WYF][WYF]", "[WYF].[WYF]", "[WYF]..[WYF]", "[WYF]...[WYF]", "[WYF]....[WYF]",
                  "[RHK][RHK]", "[RHK].[RHK]", "[RHK]..[RHK]", "[RHK]...[RHK]", "[RHK]....[RHK]",
                  "[DE][WYF]", "[DE].[WYF]", "[DE]..[WYF]", "[DE]...[WYF]", "[DE]....[WYF]",
                  "[WYF][DE]", "[WYF].[DE]", "[WYF]..[DE]", "[WYF]...[DE]", "[WYF]....[DE]")

names(motif_search) <- c(motif_search[1:20], 
                         "AA_0", "AA_1", "AA_2", "AA_3", "AA_4",
                         "RR_0", "RR_1", "RR_2", "RR_3", "RR_4",
                         "BB_0", "BB_1", "BB_2", "BB_3", "BB_4",
                         "AR_0", "AR_1", "AR_2", "AR_3", "AR_4",
                         "RA_0", "RA_1", "RA_2", "RA_3", "RA_4")


GCN4_features <- sapply(motif_search, 
                        function(i) stringr::str_count(GCN4_all_aa$sequence, i))

GCN4_all_aa <- cbind(GCN4_all_aa, GCN4_features)

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

plot_var_imp$coefficient <- coef(ridge$finalModel, ridge$bestTune$lambda)[-1]

plot_var_imp$features <- factor(plot_var_imp$features,
                                levels = plot_var_imp$features[order(plot_var_imp$importance)])
#levels(plot_var_imp$features) <- substr(levels(plot_var_imp$features), 1, 1)


# Save data for excel heatmap
write.table(plot_var_imp, "Output/GCN4_regression_aa+spacing.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

plot_var_imp <- plot_var_imp[grep("\\_", plot_var_imp$features),]

plot_var_imp$spacing <- gsub(".*\\_", "", plot_var_imp$features)
plot_var_imp$features <- gsub("\\_.*", "", plot_var_imp$features)

imp_matrix <- tidyr::spread(plot_var_imp[,-2], key = spacing, value = coefficient)
write.table(imp_matrix[c(3,1,5,2,4),], "Output/GCN4_regression_spacingMat.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)



#making the 'funnel' graph Dr. erkine wanted with the different AA having colored by class

colnames(GCN4_all_aa)


AAcolors<-c('red', 'yellow', 'orange', 'blue', 'green', 'pink',
            'dark')


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
#ggsave("Figures/Fig2_Ridge_GCN4.tiff", height = 3, width = 3,
#       units = "in")




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
#levels(plot_var_imp$features) <- substr(levels(plot_var_imp$features), 1, 1)



#making the variabel importance graph 
plot_var_imp$coefficient <- coef(lasso$finalModel, lasso$bestTune$lambda)[-1]

plot_var_imp$features <- factor(plot_var_imp$features,
                                levels = plot_var_imp$features[order(plot_var_imp$importance)])

write.table(plot_var_imp, "Output/GCN4_lassoRegression_aa+spacing.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

plot_var_imp <- plot_var_imp[grep("\\_", plot_var_imp$features),]

plot_var_imp$spacing <- gsub(".*\\_", "", plot_var_imp$features)
plot_var_imp$features <- gsub("\\_.*", "", plot_var_imp$features)

imp_matrix <- tidyr::spread(plot_var_imp[,-2], key = spacing, value = coefficient)
write.table(imp_matrix[c(3,1,5,2,4),], "Output/GCN4_lassoRegression_spacingMat.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


AAclass2 <- c(AAclass, c(AR = "Acidic-Aromatic", RA = "Aromatic-Acidic",
                         AA = "Acidic-Acidic", RR = "Aromatic-Aromatic", BB = "Basic-Basic"))

# Plot coefficients w/ bargraphs
var_imp <- read.table("Output/GCN4_lassoRegression_aa+spacing.txt",
                      header = TRUE)
var_imp$aas <- gsub("\\_.*", "", var_imp$features)
var_imp$spacing <- gsub(".*\\_", "", var_imp$features)
var_imp$spacing[nchar(var_imp$aas) == 1] <- "single"

var_imp$features <- factor(var_imp$features, 
                           levels = var_imp$features[order(abs(var_imp$coefficient), decreasing = FALSE)])
var_imp$aaclass <- AAclass2[var_imp$aas]

xmax <- max(var_imp$coefficient)
xmin <- min(var_imp$coefficient)

ggplot(var_imp[var_imp$spacing == "single",], aes(features, coefficient, fill = aaclass))+
  geom_col(colour = "black")+
  coord_flip()+
  scale_fill_manual(values = AAcolors)+ xlab("") +
  #  ggtitle('Lasso Variable Importance GCN4') +
  theme_bw() +
  ylim(xmin, xmax) +
  theme(legend.title = element_blank())
ggsave("Figure 5/Figures/Fig5c_Lasso_solos.tiff", height = 3, width = 3,
       units = "in")

var_imp2 <- var_imp[var_imp$spacing != "single",]
var_imp2$aaclass <- factor(var_imp2$aaclass, 
                           levels = c("Basic-Basic", "Aromatic-Aromatic", "Acidic-Acidic",
                                      "Aromatic-Acidic", "Acidic-Aromatic"))
var_imp2$lab <- ifelse(var_imp2$coefficient > 0, 
                       yes = -0.05, no = var_imp2$coefficient - 0.05)

ggplot(var_imp2 , aes(aaclass, coefficient, group = spacing, 
                      alpha = spacing, fill = aaclass))+
  geom_col(colour = "black", position = "dodge", width = 0.8)+
  geom_hline(yintercept = 0) +
  geom_text(aes(x = aaclass, y = lab, label = spacing), alpha = 1, size = 2,
            position = position_dodge(width = 0.8)) +
  coord_flip()+
  scale_fill_manual(values = c("blue", "orange", "red", "brown", "grey15"))+ xlab("") +
  #  ggtitle('Lasso Variable Importance GCN4') +
  theme_bw() +
  ylim(xmin, xmax) +
  theme(legend.position = "none")
ggsave("Figure 5/Figures/Fig5c_Lasso_motifs.tiff", height = 3, width = 2.8,
       units = "in")

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
ggsave("Figure 5/Figures/Fig2_Lasso_GCN4.tiff", height = 3, width = 3,
       units = "in")

plot_var_imp_lasso <- plot_var_imp

save(plot_var_imp_ridge, plot_var_imp_lasso, file = "Output/GCN4_regression.rData")




# Try a bargraph instead


