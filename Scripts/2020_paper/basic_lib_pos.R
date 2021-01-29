#this script is investigating if the position of the amino acids by class have an impact in the basic library 
head(hahn_tads)
library(stringr)

hahn_tads$acidic<-0
hahn_tads$basic<-0
hahn_tads$charge<-0

hahn_tads$acidic<-str_count(hahn_tads$sequence, '[DE]')
hahn_tads$basic<-str_count(hahn_tads$sequence, '[RKH]')
hahn_tads$charge<-hahn_tads$basic - hahn_tads$acidic

x<-which(hahn_tads$charge > 0)

hahn_tads[x,]

Basiclib<-hahn_tads[x,c("sequence", "AD_set")]

nrow(Basiclib)

head(Basiclib)

x<-which(Basiclib$AD_set == 'AD_positive')

length(x)


#Basiclib$W_pos<-0

Basiclib$Aromatic_pos<-str_locate(Basiclib$sequence, '[WYF]')

which(Basiclib$Aromatic_pos == '30')

x<-which(Basiclib$Aromatic_pos == '30')

head(Basiclib)

Basiclib_new<-Basiclib[x, c("sequence", "AD_set")]

head(Basiclib_new)

nrow(Basiclib_new)

B<-nrow(Basiclib_new)

A<-which(Basiclib_new$AD_set == 'AD_positive')

length(A)

A<-length(A)

A/B*100




###

library(tidyverse)
install.packages('tibble')

library(stringr)


test<-Basiclib 

head(test)

test$Aro_pos<-0

test$position_TF<-F



test$Aro_pos<-str_locate_all(test$sequence, '[WYF]')

head(test)


which(test$Aro_pos == "1")

str_locate(test$Aro_pos, "1")


for ( i in 1:nrow(test)){
  test$position_TF[i]<-6 %in% test$Aro_pos[i][[1]][,1]
}



test$position_TF[1]<-6 %in% test$Aro_pos[1][[1]][,1]


 for ( i in 1:nrow(test)){
   hold1<-str_split(test$sequence[i], '')[[1]]
   
   test$position_TF[i]<- hold1[1] %in% c("W","Y","F")
 }


hold1<-str_split(test$sequence[1], '')[[1]]

test$position_TF[1]<- hold1[1] %in% c("W","Y","F")



#####
 
for (i in 1:30){
  test$hold<-''
  colnames(test)[which(colnames(test)== 'hold')]<-paste('pos', i , sep = '')
}
i = 1

test

for (i in 1:nrow(test)){
  hold1<-str_split(test$sequence[i], '')[[1]]
  test$pos1[i]<-hold1[1]
  test$pos2[i]<-hold1[2]
  test$pos3[i]<-hold1[3]
  test$pos4[i]<-hold1[4]
  test$pos5[i]<-hold1[5]
  test$pos6[i]<-hold1[6]
  test$pos7[i]<-hold1[7]
  test$pos8[i]<-hold1[8]
  test$pos9[i]<-hold1[9]
  test$pos10[i]<-hold1[10]
  test$pos11[i]<-hold1[11]
  test$pos12[i]<-hold1[12]
  test$pos13[i]<-hold1[13]
  test$pos14[i]<-hold1[14]
  test$pos15[i]<-hold1[15]
  test$pos16[i]<-hold1[16]
  test$pos17[i]<-hold1[17]
  test$pos18[i]<-hold1[18]
  test$pos19[i]<-hold1[19]
  test$pos20[i]<-hold1[20]
  test$pos21[i]<-hold1[21]
  test$pos22[i]<-hold1[22]
  test$pos23[i]<-hold1[23]
  test$pos24[i]<-hold1[24]
  test$pos25[i]<-hold1[25]
  test$pos26[i]<-hold1[26]
  test$pos27[i]<-hold1[27]
  test$pos28[i]<-hold1[28]
  test$pos29[i]<-hold1[29]
  test$pos30[i]<-hold1[30]
  
  if ((i%%100)==0){
    print(i)
  }
}








