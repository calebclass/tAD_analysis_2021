hahn_tads <- readRDS("Datasets/Datasets/hahn_tads.rds")

head(hahn_tads)

anotherhold<-hahn_tads

anotherhold$Acidic_Count<-0
anotherhold$Aromatic_Count<-0
anotherhold$Aliphatic_Count<-0

anotherhold$Acidic_Count<-str_count(anotherhold$sequence, '[DE]')
anotherhold$Aromatic_Count<-str_count(anotherhold$sequence, '[WYF]')
anotherhold$Aliphatic_Count<-str_count(anotherhold$sequence, '[AVLI]')



anotherhold$Group_acidicANDaromatic= paste0(anotherhold$Acidic_Count, "_", anotherhold$Aromatic_Count)

head(anotherhold)


agg_library <- aggregate(anotherhold[,c("Acidic_Count", "Aromatic_Count", "enrichmen_score")], 
                         by = list(anotherhold$Group_acidicANDaromatic), FUN = mean)
head(agg_library)

ggplot(agg_library, aes(Acidic_Count, Aromatic_Count, fill = enrichmen_score))+
  geom_tile(colour = "grey50")+
  scale_fill_continuous(name = "Enrichment\nScore", low= 'white', high= 'blue')+
  coord_flip()+
  geom_abline(slope = 1, intercept = 0) +
  ylab("# Aromatic Residues") + xlab("# Acidic Residues") +
  guides(fill = guide_colourbar(barwidth = 0.6, barheight = 4)) +
  theme(legend.title = element_text(size = 5), 
        legend.text = element_text(size = 5),
        panel.background = element_rect(fill='grey70'))
ggsave("Figure 4/Figures/Fig4_GCN_AcidicAromaticHM.tiff", 
       height = 2.5, width = 3)


# bar graph for aromatic == acidic

df_eq <- data.frame(cts = 0:9,
                    seqs = NA,
                    num_func = NA)

for (i in 1:nrow(df_eq)) {
  df_eq$seqs[i] <- sum(anotherhold$Aromatic_Count == df_eq$cts[i] &
                         anotherhold$Acidic_Count == df_eq$cts[i])
  df_eq$num_func[i] <- sum(anotherhold$AD_set[anotherhold$Aromatic_Count == df_eq$cts[i] &
                                                  anotherhold$Acidic_Count == df_eq$cts[i]] == "AD_positive")
}
df_eq$perc_func <- df_eq$num_func / df_eq$seqs * 100

# Combine 8, 9
df_eq[9,] <- colSums(df_eq[9:10,])
df_eq <- df_eq[-10,]
df_eq$cts <- as.character(df_eq$cts)
df_eq$cts[9] <- "8+"
df_eq$perc_func[9] <- df_eq$num_func[9] / df_eq$seqs[9] * 100

ggplot(df_eq, aes(x = cts, y = perc_func)) +
  geom_col(fill = "darkgreen", colour = "black") +
  ylab("Percent Functional") + 
  xlab("# Acidic and Aromatic Residues") +
  theme_bw() +
#  scale_x_continuous(breaks = 0:9)  +
  theme(panel.grid.minor.x = element_blank())
ggsave("Figure 4/Figures/Fig4a_GCN_bargraph.tiff",
       height = 2, width = 3)


#filtering out all ad_negative sequences so that the heatmap shows true functionality of living sequneces, not brought down 
#by the non-lving sequences 
X<-which(anotherhold$AD_set=='AD_positive')

length(X)

anotherhold<-anotherhold[X,]

head(anotherhold)

nrow(anotherhold)

#making the heat map 

anotherhold$Group_acidicANDaromatic= paste0(anotherhold$Acidic_Count, "_", anotherhold$Aromatic_Count)

head(anotherhold)


agg_library <- aggregate(anotherhold[,c("Acidic_Count", "Aromatic_Count", "enrichmen_score")], 
                         by = list(anotherhold$Group_acidicANDaromatic), FUN = mean)
head(agg_library)





ggplot(agg_library, aes(Acidic_Count, Aromatic_Count, fill = enrichmen_score))+
  geom_tile()+
  scale_fill_continuous(low= 'dark blue', high= 'white')+
  coord_flip()+
  geom_abline()








#now i am making the same graph but with HSF
random_library <- readRDS("Datasets/Datasets/random_library.rds")

anotherhold<-random_library

anotherhold$Acidic_Count<-0
anotherhold$Aromatic_Count<-0
anotherhold$Aliphatic_Count<-0

anotherhold$Acidic_Count<-str_count(anotherhold$sequence, '[DE]')
anotherhold$Aromatic_Count<-str_count(anotherhold$sequence, '[WYF]')
anotherhold$Aliphatic_Count<-str_count(anotherhold$sequence, '[AVLI]')


#filtering out all ad_negative sequences so that the heatmap shows true functionality of living sequneces, not brought down 
#by the non-lving sequences 

# Caleb commenting out
#X<-which(anotherhold$binary_stop=='live')

#anotherhold<-anotherhold[X,]

head(anotherhold)



#making the heat map 

anotherhold$Group_acidicANDaromatic= paste0(anotherhold$Acidic_Count, "_", anotherhold$Aromatic_Count)

head(anotherhold)


#hist(anotherhold$estimate)
#boxplot(log2(anotherhold$estimate))

# Set outliers to maximum 10 (12 outliers)
anotherhold$estimate[anotherhold$estimate > 10] <- 10

agg_library <- aggregate(anotherhold[,c("Acidic_Count", "Aromatic_Count", "estimate")], 
                         by = list(anotherhold$Group_acidicANDaromatic), FUN = mean)
head(agg_library)


table(anotherhold$Group_acidicANDaromatic)


#ggplot(agg_library, aes(Acidic_Count, Aromatic_Count, fill = estimate))+
#  geom_tile()+
#  scale_fill_gradient(low= "white", high = "blue")+
#  coord_flip()+
#  geom_abline()

ggplot(agg_library, aes(Acidic_Count, Aromatic_Count, fill = estimate))+
  geom_tile(colour = "grey50")+
  scale_fill_continuous(name = "Growth\nEstimate", low= 'white', high= 'blue')+
  coord_flip()+
  geom_abline(slope = 1, intercept = 0) +
  ylab("# Aromatic Residues") + xlab("# Acidic Residues") +
  guides(fill = guide_colourbar(barwidth = 0.6, barheight = 4)) +
  theme(legend.title = element_text(size = 5), 
        legend.text = element_text(size = 5),
        panel.background = element_rect(fill='grey70')) +
  scale_x_continuous(breaks = seq(0, 8, by=2)) +
  scale_y_continuous(breaks = seq(0, 8, by=2))
ggsave("Figure 4/Figures/Fig4_HSF_AcidicAromaticHM.tiff", 
       height = 2.5, width = 3)




# bar graph for aromatic == acidic

df_eq <- data.frame(cts = 0:5,
                    seqs = NA,
                    num_func = NA)

for (i in 1:nrow(df_eq)) {
  df_eq$seqs[i] <- sum(anotherhold$Aromatic_Count == df_eq$cts[i] &
                         anotherhold$Acidic_Count == df_eq$cts[i])
  df_eq$num_func[i] <- sum(anotherhold$binary_stop[anotherhold$Aromatic_Count == df_eq$cts[i] &
                                                anotherhold$Acidic_Count == df_eq$cts[i]] == "live")
}

df_eq$perc_func <- df_eq$num_func / df_eq$seqs * 100
df_eq$cts <- as.character(df_eq$cts)
df_eq <- rbind(df_eq, c("4+", sum(df_eq$seq[5:6]), sum(df_eq$num_func[5:6]), NA))
df_eq$perc_func[7] <- as.numeric(df_eq$num_func[7]) / as.numeric(df_eq$seqs[7]) * 100
df_eq$perc_func <- as.numeric(df_eq$perc_func)

ggplot(df_eq[!(df_eq$cts %in% c("4", "5")),], aes(x = cts, y = perc_func)) +
  geom_col(fill = "darkgreen", colour = "black") +
  ylab("Percent Functional") + 
  xlab("# Acidic and Aromatic Residues") +
  theme_bw() +
#  scale_x_continuous(breaks = 0:9)  +
  theme(panel.grid.minor.x = element_blank())
ggsave("Figure 4/Figures/Fig4b_HSF_bargraph.tiff",
       height = 2, width = 3)















