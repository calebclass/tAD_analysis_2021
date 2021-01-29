#restarting after losing everything 

install.packages("tidyverse")
library(tidyverse)
install.packages('stringr')
library(stringr)
install.packages('fansi')
library(fansi)


head(hahn_tads)
hahn_hold<-hahn_tads[,c('sequence', 'AD_set', 'enrichmen_score')]

head(hahn_hold)

hahn_hold$Aromatic<-0
hahn_hold$Acidic<-0
hahn_hold$Aromatic<-str_count(hahn_hold$sequence, '[WYF]')
hahn_hold$Acidic<-str_count(hahn_tads$sequence, '[DE]')

head(hahn_hold)

hahn_hold$Aro_Aci_diff<-0
hahn_hold$Aro_Aci_diff<-hahn_hold$Aromatic - hahn_hold$Acidic

head(hahn_hold)

#now i will try to make the heat map
library(ggplot2)

graphhold<-ggplot(hahn_hold, aes(Aromatic, Acidic, fill = enrichmen_score))+
  geom_tile()+
  scale_fill_gradient(low= "white", high = "blue")+
  coord_flip()



#the heat map works, now I am going to make the histogram which follows the differece 
#plug in numbers ranging from -10 to 10 
#negative number <-more acidic than aromatic 
#postive number <- more aromatic than acidic 

diffhold<-subset(hahn_hold, hahn_hold$Aro_Aci_diff==10)
print(diffhold)
X<-nrow(diffhold)
print(X)
Y<-subset(diffhold, diffhold$AD_set=='AD_positive')
Y<-nrow(Y)
print(Y)
Y/X *100

#enter the output into excel for Dr.Erkine 
