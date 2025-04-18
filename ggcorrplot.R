# Correlation plot with ggcorrplot
# E. Lamont
# 4/18/25

# http://www.sthda.com/english/wiki/ggcorrplot-visualization-of-a-correlation-matrix-using-ggplot2

# Coppola 2021 just uses Spearman's, so just do this one (is it because of the way they normalized?)
# Is it bad that there are lots of NAs? (Some of the sputum datasets have a lot fewer genes)

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=10),
        legend.title = element_text(size = 12),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=10), 
        axis.text.x = element_text(angle = 0, size=10, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10), 
        plot.subtitle = element_text(size=9), 
        plot.margin = margin(10, 10, 10, 20),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_blank()
  )

source("Most_Highly_Expressed_Genes.R") # To get AllSputum_RankExpression


###########################################################
################### SPEARMAN GGCORRPLOT ###################
# Only going to do Pearson because the number of genes is so high CLT applies and parametric tests can be used

# Need to make the gene names the rownames! 
AllSputum_RankExpression_2 <- AllSputum_RankExpression %>%
  column_to_rownames(var = "Gene")

# Make the correlation
corr <- cor(AllSputum_RankExpression_2, method = "spearman",
            use = "pairwise.complete.obs")
# Needs the use argument so it will ignore the NAs between each pairwise comparison

# Compute a matrix of correlation p-values
# p.mat <- cor_pmat(AllSputum_RankExpression_2,, method = "spearman")
# head(p.mat[, 1:4])

# min(corr) # -0.02468536
# my_min <- round(min(corr), 1)

# Plot
ggcorrplot_Spearman <- corr %>% 
  ggcorrplot(hc.order = F, 
             lab = TRUE, lab_size = 4,
             type = c("lower"),
             outline.col = "white") + 
  # scale_fill_gradient2(limit = c(my_min,1), low = "blue", high =  "red", mid = "white", midpoint = (((1-my_min)/2)+my_min)) + # Make sure to change based on the min!
  my_plot_themes + 
  scale_x_discrete(guide = guide_axis(angle = 45)) + 
  labs(title = "Spearman Correlation Literature sputum", 
       subtitle = NULL, 
       fill = "Correlation")
ggcorrplot_Spearman

# ggsave(ggcorrplot_SpearmanLog10,
#        file = "ggcorrplot_SpearmanLog10_v2.pdf",
#        path = "ggcorrplot_Figures",
#        width = 7, height = 6, units = "in")



###########################################################
#################### SPEARMAN COR.TEST ####################

# Using cor.test for significance because that is what Coppola 2021 does:
### "Since the datasets used in our analysis consist of a large sample size, all Spearmans rank correlation coefficients above 0.13 were significant (calculated using the cor.test function in R)."
## Don't really understand what this means....

cor.test(AllSputum_RankExpression_2$EllaW0, 
         AllSputum_RankExpression_2$Lai2021, 
         method = "spearman")

# Do them all at once
corr_pvalue_results <- list()
cols <- colnames(AllSputum_RankExpression_2)
for (i in 1:(length(cols)-1)) {
  for (j in (i+1):length(cols)) {
    x <- AllSputum_RankExpression_2[[cols[i]]]
    y <- AllSputum_RankExpression_2[[cols[j]]]
    complete <- complete.cases(x,y) # Check for missing values
    test <- cor.test(x[complete], y[complete], method = "spearman")
    corr_pvalue_results[[paste(cols[i], cols[j], sep = "_vs_")]] <- list(
      Corr_coeficient = test$estimate,
      p = test$p.value,
      n = sum(complete)
    )
  }
}

# Convert the results into a dataframe
corr_pvalue_results_df <- do.call(rbind, lapply(names(corr_pvalue_results), function(name) {
  data.frame(Comparison = name,
             Spearman_Corr_Coeffient = corr_pvalue_results[[name]]$Corr_coeficient,
             P_Value = corr_pvalue_results[[name]]$p,
             N = corr_pvalue_results[[name]]$n)
}))



###########################################################
########### SPEARMAN GGCORRPLOT WITH P-VALUES #############

# Do this a slightly different way for adding to ggcorrplot - THIS WAY IS CLEANER!!!!
n <- ncol(AllSputum_RankExpression_2)
corr <- matrix(NA, n, n)
p.mat <- matrix(NA, n, n)
colnames(corr) <- colnames(p.mat) <- colnames(AllSputum_RankExpression_2)
rownames(corr) <- rownames(p.mat) <- colnames(AllSputum_RankExpression_2)
for (i in 1:n) {
  for (j in 1:n) {
    test <- cor.test(AllSputum_RankExpression_2[[i]], AllSputum_RankExpression_2[[j]], method = "spearman",  use = "pairwise.complete.obs")
    corr[i, j] <- test$estimate
    p.mat[i, j] <- test$p.value
  }
}

ggcorrplot_Spearman_2 <- corr %>% 
  ggcorrplot(hc.order = F, 
             p.mat = p.mat,
             insig = "pch", 
             sig.level = 0.05,
             lab = TRUE, lab_size = 4,
             type = c("lower"),
             outline.col = "white") + 
  scale_fill_gradient2(limit = c(-0.03,1), low = "blue", high =  "red", mid = "white", midpoint = (((1-(-0.03))/2)+(-0.03))) + # Make sure to change based on the min!
  my_plot_themes + 
  scale_x_discrete(guide = guide_axis(angle = 45)) + 
  labs(title = "Spearman Correlation Literature sputum", 
       subtitle = "X's mean not significant", 
       fill = "Correlation")
ggcorrplot_Spearman_2

ggsave(ggcorrplot_Spearman_2,
       file = "ggcorrplot_AllSputum_Spearman_v2.pdf",
       path = "ggcorrplot_Figures",
       width = 7, height = 6, units = "in")


###########################################################
################### PEARSON CORRELATION ####################
# Just for fun - Basically the same as Spearman

n <- ncol(AllSputum_RankExpression_2)
corr <- matrix(NA, n, n)
p.mat <- matrix(NA, n, n)
colnames(corr) <- colnames(p.mat) <- colnames(AllSputum_RankExpression_2)
rownames(corr) <- rownames(p.mat) <- colnames(AllSputum_RankExpression_2)
for (i in 1:n) {
  for (j in 1:n) {
    test <- cor.test(AllSputum_RankExpression_2[[i]], AllSputum_RankExpression_2[[j]], method = "pearson",  use = "pairwise.complete.obs")
    corr[i, j] <- test$estimate
    p.mat[i, j] <- test$p.value
  }
}

ggcorrplot_Pearson <- corr %>% 
  ggcorrplot(hc.order = F, 
             p.mat = p.mat,
             insig = "pch", 
             sig.level = 0.05,
             lab = TRUE, lab_size = 4,
             type = c("lower"),
             outline.col = "white") + 
  # scale_fill_gradient2(limit = c(my_min,1), low = "blue", high =  "red", mid = "white", midpoint = (((1-my_min)/2)+my_min)) + # Make sure to change based on the min!
  my_plot_themes + 
  scale_x_discrete(guide = guide_axis(angle = 45)) + 
  labs(title = "Pearson Correlation Literature sputum", 
       subtitle = "X's mean not significant", 
       fill = "Correlation")
ggcorrplot_Pearson


