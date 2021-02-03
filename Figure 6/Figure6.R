dat <- read.csv("Datasets/Datasets/Position_features.csv")


AAclass<-c(D= 'Acidic', E= 'Acidic',
           Y='Aromatic', W= 'Aromatic', F= 'Aromatic',
           I= 'Aliphatic', V='Aliphatic', L= 'Aliphatic', A= 'Aliphatic', M= 'Aliphatic',
           R= 'Basic', H= 'Basic', K= 'Basic',
           S= 'Polar', T= 'Polar', N= 'Polar', Q= 'Polar',
           C= 'Special', G = 'Special', P= 'Special')

dat$AAclass <- AAclass[dat$specific_composition]

library(ggplot2)

ggplot(data = dat, aes(x = Position, y = coefs, col = specific_composition)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ AAclass) + 
#  scale_colour_brewer(palette = "Set1") +
  theme_bw() +
  ylab("ML Feature Coefficient")
ggsave("Figures/FigS6_allCoefficients.tiff",
       height = 6, width = 8)

dat$specific_composition <- factor(dat$specific_composition,
                                   levels = names(AAclass))

ggplot(data = dat[dat$AAclass %in% c("Acidic", "Aromatic", "Basic"),], aes(x = Position, y = coefs, col = specific_composition)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ AAclass) + 
  scale_colour_brewer(palette = "Set1") +
  theme_bw() +
  theme(legend.title = element_blank()) +
  ylab("ML Feature Coefficient")
ggsave("Figures/Fig6_unsmoothed.tiff",
       height = 3, width = 8)


ggplot(data = dat[dat$AAclass %in% c("Acidic", "Aromatic", "Basic"),], aes(x = Position, y = coefs, col = specific_composition)) +
  geom_smooth(size = 1, method = "loess", se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ AAclass) + 
  scale_colour_brewer(palette = "Set1") +
  theme_bw() +
  theme(legend.title = element_blank()) +
  ylab("ML Feature Coefficient")
ggsave("Figures/Fig6_smoothed.tiff",
       height = 3, width = 8)


dat.sub <- dat[dat$AAclass %in% c("Acidic", "Aromatic", "Basic"),]
dat.sub$specific_composition <- factor(dat.sub$specific_composition, levels = c("W", "F", "Y", "D", "E", "R", "K", "H"))
dat.sub$AAclass <- factor(dat.sub$AAclass, levels = c("Aromatic", "Acidic", "Basic"))

ggplot(data = dat.sub, aes(x = Position, y = coefs, col = specific_composition)) +
  geom_smooth(size = 1, method = "loess", se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ AAclass, ncol = 1) + 
  scale_colour_manual(values = c("#ad7f00", "#ffbb00", "#fcda7c",
                                 "#ba0000", "#ff7070",
                                 "#13008f", "#2200ff", "#9686fc")) +
#  scale_colour_brewer(palette = "Set1") +
  theme_bw() +
  theme(legend.title = element_blank()) +
  ylab("ML Feature Coefficient")
ggsave("Figures/Fig6_smoothedVert-cols.tiff",
       height = 7, width = 4)


library(tidyverse)

dat.avg <- dat %>% group_by(AAclass, Position) %>%
  summarise(coef_mean = mean(coefs))

AAcolors<-c('red', 'yellow', 'orange', 'blue', 'green', 'pink')

ggplot(data = dat.avg, aes(x = Position, y = coef_mean, col = AAclass)) +
  geom_smooth(size = 1, method = "loess", se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_manual(values = AAcolors) +
  #scale_colour_brewer(palette = "Set1") +
  theme_bw() +
  theme(legend.title = element_blank()) +
  ylab("ML Feature Coefficient")
ggsave("Figures/Fig6_smoothedAgg-cols.tiff",
       height = 3, width = 5)


ggplot(data = dat, aes(x = Position, y = coefs, col = AAclass)) +
  geom_smooth(size = 1, method = "loess", se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_brewer(palette = "Set1") +
  theme_bw() +
  theme(legend.title = element_blank()) +
  ylab("ML Feature Coefficient")


# Live% for clusters
hahn_tads <- readRDS("Datasets/Datasets/hahn_tads.rds")
base_func <- mean(hahn_tads$AD_set == "AD_positive") * 100

hahn_tads$Basicstotal<- str_count(hahn_tads$sequence, "[HRK]")
hahn_tads$Acidstotal<- str_count(hahn_tads$sequence, "[ED]")
hahn_tads$Aromaticstotal<- str_count(hahn_tads$sequence, "[WYF]")
hahn_check <-hahn_tads %>% filter(Basicstotal >= 4 | Acidstotal >= 4 | Aromaticstotal >= 4)

clust_counts <- matrix(data = NA, nrow = 26, ncol = 3,
                      dimnames = list(1:26, c("Acidic", "Aromatic", "Basic")))
clust_check <- matrix(data = NA, nrow = 26, ncol = 3,
                      dimnames = list(1:26, c("Acidic", "Aromatic", "Basic")))

for (i in 1:nrow(clust_check)) {
  subseq <- substr(hahn_check$sequence, i, i+4)
  acidCt <- str_count(subseq, "[DE]")
  aromCt <- str_count(subseq, "[WYF]")
  baseCt <- str_count(subseq, "[HRK]")
  
  clust_counts[i,] <- c(sum(acidCt >= 4), sum(aromCt >= 4), sum(baseCt >= 4))
  clust_check[i,] <- c(mean(hahn_check$AD_set[acidCt >= 4] == "AD_positive"),
                       mean(hahn_check$AD_set[aromCt >= 4] == "AD_positive"),
                       mean(hahn_check$AD_set[baseCt >= 4] == "AD_positive"))
}

clust_df <- as.data.frame(clust_check)
clust_df$pos <- 1:nrow(clust_df) + 4
clust_plot <- gather(clust_df,
                     key = "AAclass", value = "func", -4)
clust_plot$funcPerc <- clust_plot$func * 100

ggplot(data = clust_plot, aes(x = pos, y = funcPerc, col = AAclass)) +
  geom_smooth(size = 1, method = "loess", se = FALSE) +
  geom_hline(yintercept = base_func, linetype = "dashed") +
  scale_colour_manual(values = c("red", "orange", "blue")) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  ylab("% Functional") +
  xlab("Cluster position")
ggsave("Figures/Fig6B_clustPercents.tiff",
       height = 3, width = 5)

