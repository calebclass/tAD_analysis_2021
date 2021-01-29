##ML on acidic aromatic and basic counts only
gcn4 <- readRDS('Datasets/Datasets/hahn_tads.rds') %>% as_tibble() %>% 
  select(sequence, AD_set)

#add acidic counts
gcn4$E <- str_count(gcn4$sequence,'E')
gcn4$D <- str_count(gcn4$sequence,'D')

#add basic counts
gcn4$R <- str_count(gcn4$sequence,'R')
gcn4$H <- str_count(gcn4$sequence,'H')
gcn4$K <- str_count(gcn4$sequence,'K')

#add aromatic counts
gcn4$W <- str_count(gcn4$sequence,'W')
gcn4$Y <- str_count(gcn4$sequence,'Y')
gcn4$F <- str_count(gcn4$sequence,'F')

#split into live and die
live <- gcn4 %>% 
  select(-sequence) %>% 
  filter(AD_set == 'AD_positive')

die <- gcn4 %>% 
  select(-sequence) %>% 
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


#grid of lamdas to search over
lambda <- 10^seq(-9, -3, length = 20)

# Build the model
ridge <- train(
  AD_set ~., data = train_set, method = "glmnet",
  trControl = trainControl("cv", number = 5),
  tuneGrid = expand.grid(alpha = 0, lambda = lambda)
)

#add auc to this graph

predicted_prob<-predict.train(ridge, test_set,type = "prob") 
myroc<-roc(test_set$AD_set, predicted_prob$AD_negative)
myroc$auc

#coef----
coef <- coef(ridge$finalModel, ridge$bestTune$lambda)
coef <- coef[,1]
coef <- coef[-1]

coef <- tibble(feature = names(coef),
               coef = coef)

coef$color <- c(rep('red',2),
                rep('blue',3),
                rep('orange',3))

coef$mag_coef <- sqrt(coef$coef * coef$coef)

#order by coeff
coef$feature <- factor(coef$feature,
                       levels = coef$feature[order(coef$mag_coef)])

ggplot(coef, aes(feature, coef,fill = color))+
  geom_col(show.legend = F,color='black')+
  scale_fill_manual(values = c('blue','orange','red'))+
  xlab('Motifs')+
  ylab('Feature Coefficient')+
  ggtitle('Gcn4', subtitle = 'Ridge ROC AUC: 0.921')+
  coord_flip()
