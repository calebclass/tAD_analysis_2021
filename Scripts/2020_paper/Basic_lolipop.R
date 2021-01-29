#making a lollipop for the basic lib

Basiclib<-hahn_tads

Basiclib$Bacic<-str_count(Basiclib$sequence, '[RHK]')
Basiclib$Acidic<-str_count(Basiclib$sequence, '[DE]')
Basiclib$charge<-Basiclib$Bacic-Basiclib$Acidic

head(Basiclib)
B<- which(Basiclib$charge > 0)

Basiclib<-Basiclib[B,]

head(Basiclib)

Basiclib<-Basiclib[,c("sequence", "AD_set")]



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


head(Basiclib)

livehold<-which(Basiclib$AD_set=='AD_positive')
livehold<-Basiclib[livehold,]

diehold<-which(Basiclib$AD_set=='AD_negative')
diehold<-Basiclib[diehold,]


#adding up the counts for all of the amino acids 
sum(livehold$W_count)

livesum<-list()
colnames(livehold)

for(i in colnames(livehold)[-c(1,2)]){
  
  livesum[i] <-sum(livehold[,i]) 
  
}

diesum<- list()

for(i in colnames(diehold)[-c(1,2)]){
  
  diesum[i]<-sum(diehold[,i])
}

print(diesum)


liveVSdieDF<- data.frame(row.names = 1:40)


liveVSdieDF$AAcount<- names(livesum)

liveVSdieDF$livesum<-unlist(livesum)

liveVSdieDF$CheckAA<- names(diesum)



liveVSdieDF$diesum<-unlist(diesum)


liveVSdieDF$ratio<- liveVSdieDF$livesum/liveVSdieDF$diesum

head(liveVSdieDF)


#make a new col based off the percent live vs percent die 

#also make new col with the absolute value of die sum therefor we will not have to take the negative log 

colnames(liveVSdieDF)

SumLhold_weight<-sum(liveVSdieDF$livesum)

liveVSdieDF$DieSumAbs<- abs(liveVSdieDF$diesum)

#this is taking the absolute value of the #letter * the slope 
#(note it will always be negative because dying sequences will always have a negative slope)


liveVSdieDF$LivePercent<- liveVSdieDF$livesum/SumLhold_weight *100


SumDhold_weight<-sum(liveVSdieDF$DieSumAbs)
liveVSdieDF$DiePercent<- liveVSdieDF$DieSumAbs/SumDhold_weight *100  

liveVSdieDF$PercentRatio<- liveVSdieDF$LivePercent/liveVSdieDF$DiePercent

#making the graph in log2   
colnames(liveVSdieDF)



#initilizing the new log col
liveVSdieDF$Log2Ratio<-0
#log2 connects length and magnitude, 50% had the same magnitude as 200%
#to notmalize the ratios to the lives to the dies 


#liveVSdieDF$Log2Ratio<- log2(liveVSdieDF$Ratio)
#above gave me some weird numbers, now trying frequency 
liveVSdieDF$Log2Ratio<- log2(liveVSdieDF$PercentRatio)


#note: not sure if I factored this yet, so it might be really ugly 


Weightedgraph<- liveVSdieDF[c(21:40),]



Weightedgraph$CheckAA<-gsub('_count','', Weightedgraph$CheckAA)

Weightedgraph

Weightedgraph$CheckAA<- factor(Weightedgraph$CheckAA, levels=Weightedgraph[order(Weightedgraph$PercentRatio,decreasing=T),]$CheckAA)

#now start making the graph 



ggplot(data=Weightedgraph, mapping= aes(Log2Ratio, CheckAA)) +
  geom_segment(aes(x = 0, y = CheckAA, xend = Log2Ratio, yend = CheckAA), color = "black") +
  geom_vline(xintercept = 0)+
  geom_point(aes(colour= CheckAA), size=7, show.legend = FALSE) + 
  geom_point(shape = 1,size = 7,colour = "black")+
  geom_text(label=Weightedgraph$CheckAA)+
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf,
           alpha = .2) +
  labs(x= "log_2 Frequency Functional/ Frequency Non-Functional")+
  labs(y= "Amino Acid")+
  geom_text(x=-0.45, y= 4, label="Depleted")+
  geom_text(x=0.75, y= 16, label="Enriched")+
  ggtitle("Basiclib Lollipop")

























