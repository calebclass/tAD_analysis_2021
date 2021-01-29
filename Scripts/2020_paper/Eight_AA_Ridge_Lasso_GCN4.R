#Ridge and lasso with 8 amino acids (DE, WYF, KRH)
eighthold<-hahn_tads

eighthold$W_count<-0
eighthold$Y_count<-0
eighthold$F_count<-0
eighthold$D_count<-0
eighthold$E_count<-0
eighthold$K_count<-0
eighthold$H_count<-0
eighthold$R_count<-0

eighthold$W_count<-str_count(eighthold$sequence, 'W')
eighthold$Y_count<-str_count(eighthold$sequence, 'Y')
eighthold$F_count<-str_count(eighthold$sequence, 'F')
eighthold$D_count<-str_count(eighthold$sequence, "D")
eighthold$E_count<-str_count(eighthold$sequence, 'E')
eighthold$K_count<-str_count(eighthold$sequence, 'K')
eighthold$H_count<-str_count(eighthold$sequence, 'H')
eighthold$R_count<-str_count(eighthold$sequence, 'R')


X<-which(eighthold$AD_set=='AD_positive')
live<-eighthold[X,]

head(live)
nrow(live)

Y<-which(eighthold$AD_set=='AD_negative')
die<-eighthold[Y,]

head(die)
nrow(die)


#remove sequence column, add AD_set back to my datasets
live <- live %>% 
  as_tibble() %>% 
  select(-sequence) %>% 
  mutate(AD_set = 'AD_positive')

die <- die %>% 
  as_tibble() %>% 
  select(-sequence) %>% 
  mutate(AD_set = 'AD_negative')












