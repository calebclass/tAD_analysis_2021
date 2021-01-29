#making a lasso of the amino acid counts 
aminoacidlasso<- data.frame()

aminoacidlasso<- randomlibraryfeatures[, c("sequence", "binary_stop")]

#now make the string count for each sequence 



aminoacidlasso$R_count<-0
aminoacidlasso$H_count<-0
aminoacidlasso$K_count<-0
aminoacidlasso$D_count<-0  
aminoacidlasso$E_count<-0
aminoacidlasso$S_count<-0  
aminoacidlasso$T_count<-0
aminoacidlasso$N_count<-0  
aminoacidlasso$Q_count<-0
aminoacidlasso$C_count<-0  
aminoacidlasso$G_count<-0  
aminoacidlasso$P_count<-0  
aminoacidlasso$A_count<-0
aminoacidlasso$V_count<-0
aminoacidlasso$I_count<-0
aminoacidlasso$L_count<-0
aminoacidlasso$M_count<-0
aminoacidlasso$F_count<-0
aminoacidlasso$Y_count<-0
aminoacidlasso$W_count<-0



aminoacidlasso$R_count<-str_count(aminoacidlasso$sequence,"R")
aminoacidlasso$H_count<-str_count(aminoacidlasso$sequence,"H")
aminoacidlasso$K_count<-str_count(aminoacidlasso$sequence,"K")
aminoacidlasso$D_count<-str_count(aminoacidlasso$sequence,"D")  
aminoacidlasso$E_count<-str_count(aminoacidlasso$sequence,"E")
aminoacidlasso$S_count<-str_count(aminoacidlasso$sequence,"S")  
aminoacidlasso$T_count<-str_count(aminoacidlasso$sequence,"T")
aminoacidlasso$N_count<-str_count(aminoacidlasso$sequence,"N")  
aminoacidlasso$Q_count<-str_count(aminoacidlasso$sequence,"Q")
aminoacidlasso$C_count<-str_count(aminoacidlasso$sequence,"C")  
aminoacidlasso$G_count<-str_count(aminoacidlasso$sequence,"G")  
aminoacidlasso$P_count<-str_count(aminoacidlasso$sequence,"P")  
aminoacidlasso$A_count<-str_count(aminoacidlasso$sequence,"A")
aminoacidlasso$V_count<-str_count(aminoacidlasso$sequence,"V")
aminoacidlasso$I_count<-str_count(aminoacidlasso$sequence,"I")
aminoacidlasso$L_count<-str_count(aminoacidlasso$sequence,"L")
aminoacidlasso$M_count<-str_count(aminoacidlasso$sequence,"M")
aminoacidlasso$F_count<-str_count(aminoacidlasso$sequence,"F")
aminoacidlasso$Y_count<-str_count(aminoacidlasso$sequence,"Y")
aminoacidlasso$W_count<-str_count(aminoacidlasso$sequence,"W")

head(aminoacidlasso)

set.seed(6969)

aminoacidlasso<- subset(aminoacidlasso, select= -sequence)

aminoacidlasso$binary_stop<-factor(aminoacidlasso$binary_stop)

aminohold<-createDataPartition(aminoacidlasso$binary_stop, p=0.75, list = F)

lassotrain<-aminoacidlasso[aminohold,]
lassotest<-aminoacidlasso[-aminohold,]




myfeatures_plus<-paste(colnames(
  subset(
    lassotrain,select = -c(binary_stop)
  )),collapse = " + ")

myform<-as.formula(paste0("binary_stop ~",paste0(myfeatures_plus
)))

myhold<-findLinearCombos(
  model.matrix(myform,as.data.frame(lassotrain)))$remove
myhold<-colnames(lassotrain)[myhold]
lassotrain[,myhold]<-NULL
lassotest[,myhold]<-NULL






#subset the data, take out the factor, espcially because we are testing for it 
holdtrain<-subset(lassotrain, select = -binary_stop)

#downsampling 
lassotrain<- downSample(holdtrain, lassotrain$binary_stop)

cls.ctrl <- trainControl(method = "repeatedcv", #boot, cv, LOOCV, timeslice OR adaptive etc.
                         number = 10, #number of folds 
                         repeats = 5,
                         classProbs = TRUE, 
                         summaryFunction = twoClassSummary,
                         savePredictions = "final", 
                         allowParallel = TRUE)

set.seed(1997)
glm.fit <- train(Class ~ ., data = lassotrain, trControl = cls.ctrl,
                 method = "glm", #look up difference between GLM and GLMnet 
                 #this line changes the model, this is the package the machine learning is in. this will downlad and install the GLM package 
                 family = "binomial", metric = "ROC",
                 preProcess = c("nzv", "center", "scale"))

varImp(glm.fit)%>% plot()

Variableimport<- varImp(glm.fit)

Present<- ggplot(data=Variableimport, aes(Variableimport$importance))+
  ggtitle("Lasso Amino Acid Importance")


lassoraw<- predict.train(glm.fit, 
                         lassotest, type = "raw")

lassoprob <- predict.train(glm.fit, 
                           lassotest, type = "prob")

CrossTable(lassoraw, 
           lassotest$binary_stop)

myroc2<-roc(lassotest$binary_stop,
            lassoprob$die)

#AUC is 0.767
plot(myroc2)


mystats<-data.frame(myroc2$sensitivities,myroc2$specificities)
colnames(mystats)<-c('sensitivities','specificities')

method= 'Lasso'
graphhold<-mystats %>%
  ggplot(
    mapping = aes(x = 1 - specificities, y = sensitivities, color = method)
  ) +
  geom_rect(aes(xmin=0, xmax=0.25, ymin=0.75, ymax=1)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  geom_point(size = 0.1)+
  geom_line(size=1)+
  ggtitle("Amino Acid LASSO Roc")+
  annotate('rect',xmin = 0.625,xmax = 0.875,ymin = 0.1875,
           ymax = 0.3125,fill = 'gray')+
  annotate('text',x=0.75,y=0.25,label=paste('AUC:',round(myroc2$auc[1],4)),color = 'red')
myroc2$auc
ggplotly(graphhold,width = 800,height = 450)

#now I am running LASSO on a library that only has K, this is in the hopes of understanding the outliers 
####----

K<-which(aminoacidlasso$K_count > 1)

Rhold<- data.frame()

Rhold<-randomlibraryfeatures[K, c("sequence", "binary_stop","acidic_number", "tiny_number", "small_number", "aliphatic_number","aromatic_number","nonpolar_number","polar_number","charged_number","basic_number")]



set.seed(6969)

Rhold<- subset(Rhold, select= -sequence)

Rhold$binary_stop<-factor(Rhold$binary_stop)

Rholdhold<-createDataPartition(Rhold$binary_stop, p=0.75, list = F)

lassotrain<-Rhold[Rholdhold,]
lassotest<-Rhold[-Rholdhold,]




myfeatures_plus<-paste(colnames(
  subset(
    lassotrain,select = -c(binary_stop)
  )),collapse = " + ")

myform<-as.formula(paste0("binary_stop ~",paste0(myfeatures_plus
)))

myhold<-findLinearCombos(
  model.matrix(myform,as.data.frame(lassotrain)))$remove
myhold<-colnames(lassotrain)[myhold]
lassotrain[,myhold]<-NULL
lassotest[,myhold]<-NULL






#subset the data, take out the factor, espcially because we are testing for it 
holdtrain<-subset(lassotrain, select = -binary_stop)

#downsampling 
lassotrain<- downSample(holdtrain, lassotrain$binary_stop)

cls.ctrl <- trainControl(method = "repeatedcv", #boot, cv, LOOCV, timeslice OR adaptive etc.
                         number = 10, #number of folds 
                         repeats = 5,
                         classProbs = TRUE, 
                         summaryFunction = twoClassSummary,
                         savePredictions = "final", 
                         allowParallel = TRUE)

set.seed(1997)
glm.fit <- train(Class ~ ., data = lassotrain, trControl = cls.ctrl,
                 method = "glm", #look up difference between GLM and GLMnet 
                 #this line changes the model, this is the package the machine learning is in. this will downlad and install the GLM package 
                 family = "binomial", metric = "ROC",
                 preProcess = c("nzv", "center", "scale"))

varImp(glm.fit)%>% plot()

Variableimport<- varImp(glm.fit)

Present<- ggplot(data=Variableimport, aes(Variableimport$importance))+
  ggtitle("K library feature importance")


lassoraw<- predict.train(glm.fit, 
                         lassotest, type = "raw")

lassoprob <- predict.train(glm.fit, 
                           lassotest, type = "prob")

CrossTable(lassoraw, 
           lassotest$binary_stop)

myroc2<-roc(lassotest$binary_stop,
            lassoprob$die)

#AUC is 0.767
plot(myroc2)


mystats<-data.frame(myroc2$sensitivities,myroc2$specificities)
colnames(mystats)<-c('sensitivities','specificities')

method= 'Lasso'
graphhold<-mystats %>%
  ggplot(
    mapping = aes(x = 1 - specificities, y = sensitivities, color = method)
  ) +
  geom_rect(aes(xmin=0, xmax=0.25, ymin=0.75, ymax=1)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  geom_point(size = 0.1)+
  geom_line(size=1)+
  ggtitle("K LASSO ROC")+
  annotate('rect',xmin = 0.625,xmax = 0.875,ymin = 0.1875,
           ymax = 0.3125,fill = 'gray')+
  annotate('text',x=0.75,y=0.25,label=paste('AUC:',round(myroc2$auc[1],4)),color = 'red')
myroc2$auc
ggplotly(graphhold,width = 800,height = 450)


####--------
#now I am running LASSO on a library with both R and K  

basichold<-grep("[RHK]", randomlibraryfeatures$sequence)

Bhold<-randomlibraryfeatures[basichold, c("sequence", "binary_stop","acidic_number", "tiny_number", "small_number", "aliphatic_number","aromatic_number","nonpolar_number","polar_number","charged_number","basic_number")]


set.seed(6969)

Bhold<- subset(Bhold, select= -sequence)

Bhold$binary_stop<-factor(Bhold$binary_stop)

Bholdhold<-createDataPartition(Bhold$binary_stop, p=0.75, list = F)

lassotrain<-Bhold[Bholdhold,]
lassotest<-Bhold[-Bholdhold,]




myfeatures_plus<-paste(colnames(
  subset(
    lassotrain,select = -c(binary_stop)
  )),collapse = " + ")

myform<-as.formula(paste0("binary_stop ~",paste0(myfeatures_plus
)))

myhold<-findLinearCombos(
  model.matrix(myform,as.data.frame(lassotrain)))$remove
myhold<-colnames(lassotrain)[myhold]
lassotrain[,myhold]<-NULL
lassotest[,myhold]<-NULL






#subset the data, take out the factor, espcially because we are testing for it 
holdtrain<-subset(lassotrain, select = -binary_stop)

#downsampling 
lassotrain<- downSample(holdtrain, lassotrain$binary_stop)


cls.ctrl <- trainControl(method = "repeatedcv", #boot, cv, LOOCV, timeslice OR adaptive etc.
                         number = 10, #number of folds 
                         repeats = 5,
                         classProbs = TRUE, 
                         summaryFunction = twoClassSummary,
                         savePredictions = "final", 
                         allowParallel = TRUE)

set.seed(1997)
glm.fit <- train(Class ~ ., data = lassotrain, trControl = cls.ctrl,
                 method = "glm", #look up difference between GLM and GLMnet 
                 #this line changes the model, this is the package the machine learning is in. this will downlad and install the GLM package 
                 family = "binomial", metric = "ROC",
                 preProcess = c("nzv", "center", "scale"))

varImp(glm.fit)%>% plot()

Variableimport<- varImp(glm.fit)

Present<- ggplot(data=Variableimport, aes(Variableimport$importance))+
  ggtitle("Basic library feature importance")


lassoraw<- predict.train(glm.fit, 
                         lassotest, type = "raw")

lassoprob <- predict.train(glm.fit, 
                           lassotest, type = "prob")

CrossTable(lassoraw, 
           lassotest$binary_stop)

myroc2<-roc(lassotest$binary_stop,
            lassoprob$die)

#AUC is 0.767
plot(myroc2)


mystats<-data.frame(myroc2$sensitivities,myroc2$specificities)
colnames(mystats)<-c('sensitivities','specificities')

method= 'Lasso'
graphhold<-mystats %>%
  ggplot(
    mapping = aes(x = 1 - specificities, y = sensitivities, color = method)
  ) +
  geom_rect(aes(xmin=0, xmax=0.25, ymin=0.75, ymax=1)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  geom_point(size = 0.1)+
  geom_line(size=1)+
  ggtitle("Basic LASSO ROC")+
  annotate('rect',xmin = 0.625,xmax = 0.875,ymin = 0.1875,
           ymax = 0.3125,fill = 'gray')+
  annotate('text',x=0.75,y=0.25,label=paste('AUC:',round(myroc2$auc[1],4)),color = 'red')
myroc2$auc
ggplotly(graphhold,width = 800,height = 450)

####-----
#Running LASSO on the library without aromatic or acidic 
arohold<-grep("[DEWYF]", randomlibraryfeatures$sequence)

arohold<-randomlibraryfeatures[-arohold, c("peptide_length","sequence", "binary_stop","acidic_number", "tiny_number", "small_number", "aliphatic_number","aromatic_number","nonpolar_number","polar_number","charged_number","basic_number")]

set.seed(6969)

arohold<- subset(arohold, select= -sequence)

arohold$binary_stop<-factor(arohold$binary_stop)

aroholdhold<-createDataPartition(arohold$binary_stop, p=0.75, list = F)

lassotrain<-arohold[aroholdhold,]
lassotest<-arohold[-aroholdhold,]





myfeatures_plus<-paste(colnames(
  subset(
    lassotrain,select = -c(binary_stop)
  )),collapse = " + ")

myform<-as.formula(paste0("binary_stop ~",paste0(myfeatures_plus
)))

myhold<-findLinearCombos(
  model.matrix(myform,as.data.frame(lassotrain)))$remove
myhold<-colnames(lassotrain)[myhold]
lassotrain[,myhold]<-NULL
lassotest[,myhold]<-NULL






#subset the data, take out the factor, espcially because we are testing for it 
holdtrain<-subset(lassotrain, select = -binary_stop)

#downsampling 
lassotrain<- downSample(holdtrain, lassotrain$binary_stop)


cls.ctrl <- trainControl(method = "repeatedcv", #boot, cv, LOOCV, timeslice OR adaptive etc.
                         number = 10, #number of folds 
                         repeats = 5,
                         classProbs = TRUE, 
                         summaryFunction = twoClassSummary,
                         savePredictions = "final", 
                         allowParallel = TRUE)

set.seed(1997)
glm.fit <- train(Class ~ ., data = lassotrain, trControl = cls.ctrl,
                 method = "glm", #look up difference between GLM and GLMnet 
                 #this line changes the model, this is the package the machine learning is in. this will downlad and install the GLM package 
                 family = "binomial", metric = "ROC",
                 preProcess = c("nzv", "center", "scale"))

varImp(glm.fit)%>% plot()

Variableimport<- varImp(glm.fit)

Present<- ggplot(data=Variableimport, aes(Variableimport$importance))+
  ggtitle("[RHK] library feature importance")


lassoraw<- predict.train(glm.fit, 
                         lassotest, type = "raw")

lassoprob <- predict.train(glm.fit, 
                           lassotest, type = "prob")

CrossTable(lassoraw, 
           lassotest$binary_stop)

myroc2<-roc(lassotest$binary_stop,
            lassoprob$die)

#AUC is 0.767
plot(myroc2)


mystats<-data.frame(myroc2$sensitivities,myroc2$specificities)
colnames(mystats)<-c('sensitivities','specificities')

method= 'Lasso'
graphhold<-mystats %>%
  ggplot(
    mapping = aes(x = 1 - specificities, y = sensitivities, color = method)
  ) +
  geom_rect(aes(xmin=0, xmax=0.25, ymin=0.75, ymax=1)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  geom_point(size = 0.1)+
  geom_line(size=1)+
  ggtitle("W/O Aromatic LASSO ROC")+
  annotate('rect',xmin = 0.625,xmax = 0.875,ymin = 0.1875,
           ymax = 0.3125,fill = 'gray')+
  annotate('text',x=0.75,y=0.25,label=paste('AUC:',round(myroc2$auc[1],4)),color = 'red')
myroc2$auc
ggplotly(graphhold,width = 800,height = 450)













#looking into the W..LF motif, see if it lives in our library 

didntwork<-grep("W..LF", randomlibraryfeatures$sequence)

didnt<-data.frame()

didnt<-randomlibraryfeatures[didntwork,c("sequence", "binary_stop")]

summary(didnt$binary_stop)

















