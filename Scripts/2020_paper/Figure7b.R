meths <- c("Lasso Regression", "Ridge Regression",
           "Fully Connected NN", "LSTM Network")
aucs <- data.frame(method = factor(meths, levels = rev(meths)),
                   AUC = c(0.935, 0.934, 0.975, 0.983))

library(ggplot2)

ggplot(aucs, aes(x = method, y = AUC)) +
  geom_col(fill = "darkgreen", colour = "black", width = 0.5) +
  scale_y_continuous(limits=c(0.5, 1),oob = scales::rescale_none) +
  coord_flip() +
  xlab("") +
  theme_bw()
ggsave("Figures/Fig7_aucBars.tiff",
       height = 2, width = 3)
