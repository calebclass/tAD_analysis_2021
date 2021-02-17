library(tidyverse)

gcn4 <- readRDS('Datasets/Datasets/hahn_tads.rds') %>% as_tibble()

gcn4$aro_count <- str_count(gcn4$sequence, '[WYF]')
gcn4$acid_count <- str_count(gcn4$sequence, '[DE]')

gcn4$total <- gcn4$aro_count + gcn4$acid_count

gcn4$diff <- gcn4$aro_count - gcn4$acid_count

plot_me <- expand.grid(-14:13,0:19) %>% as_tibble()
colnames(plot_me) <- c('diff','total')

plot_me <- plot_me %>% mutate(
  live_count = 0,
  die_count = 0
)

for(i in 1:nrow(plot_me)){
  lib <- gcn4 %>% filter(diff == plot_me$diff[i], total == plot_me$total[i])
  plot_me$live_count[i] <- lib %>% filter(AD_set == 'AD_positive') %>% nrow()
  plot_me$die_count[i] <- lib %>% filter(AD_set == 'AD_negative') %>% nrow()
}

plot_me$live_percent <- plot_me$live_count/(plot_me$live_count + plot_me$die_count)*100

#total needs to be a factor
plot_me$total <- factor(plot_me$total,
                        levels = 0:19)

plot_me %>% filter(die_count > 99) %>% 
  filter(total %in% c(6,8,10,12,14,16)) %>% 
  ggplot(., aes(diff, live_percent, color = total))+
  geom_vline(xintercept = 0, linetype = 'longdash')+
  geom_line(size = 1)+
  geom_point()+
  theme_bw()+
  xlab('#Aromatic - #Acidic')+
  ylab('Percent Functional')

#ggsave("Figure 5/Figures/new_aro_acidic_balance.tiff",
#       height = 2.5, width = 3.0)

