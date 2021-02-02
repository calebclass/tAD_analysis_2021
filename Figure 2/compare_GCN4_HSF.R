(load("Output/HSF_regression.rData"))

hsf.res <- list(ridge = plot_var_imp_ridge,
                lasso = plot_var_imp_lasso)

(load("Output/GCN4_regression.rData"))

gcn.res <- list(ridge = plot_var_imp_ridge,
                lasso = plot_var_imp_lasso)

all(gcn.res$ridge$features == hsf.res$ridge$features)

ridge.merged <- cbind(gcn.res$ridge, hsf.res$ridge[,2:3])
colnames(ridge.merged)[c(2,3,5,6)] <- c("imp.HSF", "HSF", "imp.GCN4", "GCN4")

AAcolors<-c('red', 'yellow', 'orange', 'blue', 'green', 'pink')

cor.test(ridge.merged$HSF, ridge.merged$GCN4, method = "spearman")

#library(ggplot2)
#setEPS("Figures/Fig2_Ridge_Compare.tiff", height = 5, width = 6.5)
postscript("Figures/Fig2_Ridge_Compare.eps", height = )
ggplot(ridge.merged, aes(x = GCN4, y = HSF)) +
  geom_point(aes(fill = AAclass), size = 8,shape = 21) +
  ggrepel::geom_text_repel(aes(label = features),
                           point.padding = 0, box.padding = 0.1, label.padding = 1) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  scale_fill_manual(values = AAcolors) +
  theme_bw() +
  theme(legend.title = element_blank())
#dev.off()
ggsave("Figure 2/Figures/Fig2_Ridge_Compare.tiff", height = 5, width = 6.5,
       units = "in")







# Replot
AAcolors<-c('red', 'yellow', 'orange', 'blue', 'green', 'pink')


AAclass<-c(Y='Aromatic', W= 'Aromatic', F= 'Aromatic',
           I= 'Aliphatic', V='Aliphatic', L= 'Aliphatic', A= 'Aliphatic', M= 'Aliphatic',
           R= 'Basic', H= 'Basic', K= 'Basic',
           D= 'Acidic', E= 'Acidic',
           S= 'Polar', T= 'Polar', N= 'Polar', Q= 'Polar',
           C= 'Special', G = 'Special', P= 'Special')



ggplot(gcn.res$lasso, aes(features, coefficient, fill = AAclass))+
  geom_col(col = "black")+
  coord_flip()+
  scale_fill_manual(values = AAcolors)+ xlab("") +
  #  ggtitle('Lasso Variable Importance GCN4') +
  theme_bw() +
  theme(legend.title = element_blank())
ggsave("Figure 2/Figures/Fig2_Lasso_GCN4.tiff", height = 3, width = 3,
       units = "in")


ggplot(gcn.res$ridge, aes(features, coefficient, fill = AAclass))+
  geom_col(col = "black")+
  coord_flip()+
  scale_fill_manual(values = AAcolors)+ xlab("") +
  #  ggtitle('Lasso Variable Importance GCN4') +
  theme_bw() +
  theme(legend.title = element_blank())
ggsave("Figure 2/Figures/Fig2_Ridge_GCN4.tiff", height = 3, width = 3,
       units = "in")

ggplot(hsf.res$lasso, aes(features, coefficient, fill = AAclass))+
  geom_col(col = "black")+
  coord_flip()+
  scale_fill_manual(values = AAcolors)+ xlab("") +
  #  ggtitle('Lasso Variable Importance GCN4') +
  theme_bw() +
  theme(legend.title = element_blank())
ggsave("Figure 2/Figures/Fig2_Lasso_HSF.tiff", height = 3, width = 3,
       units = "in")

ggplot(hsf.res$ridge, aes(features, coefficient, fill = AAclass))+
  geom_col(col = "black")+
  coord_flip()+
  scale_fill_manual(values = AAcolors)+ xlab("") +
  #  ggtitle('Lasso Variable Importance GCN4') +
  theme_bw() +
  theme(legend.title = element_blank())
ggsave("Figure 2/Figures/Fig2_Ridge_HSF.tiff", height = 3, width = 3,
       units = "in")