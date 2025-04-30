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

source("Most_Highly_Expressed_Genes.R") # To get AllSputum_RankExpression_W0RawReads

###########################################################
############### GGCORRPLOT NA ROWS REMOVED ################

# From Coppola2021 spearman correlation analysis: "The entire datasets-comparison consisted of the ranked expression of Mtb genes (n=1813) that were commonly investigated in all nine Mtb transcriptomes."
# So I think I have to keep just the genes that are in ALL that datasets (remove any row that contains NA) 
# I think this is the best way to do it, but not sure.....

# Remove any row that contains NA
AllSputum_RankExpression_W0RawReads_2_noNA <- AllSputum_RankExpression_W0RawReads_2 %>% na.omit()
# Now there are 1932 genes... still more than they have, but they had included non-sputum datasets so this may be where the descrepancy is....

# Make the correlation
corr <- cor(AllSputum_RankExpression_W0RawReads_2_noNA, method = "spearman")
# Needs the use argument so it will ignore the NAs between each pairwise comparison

# Compute a matrix of correlation p-values
p.mat <- cor_pmat(AllSputum_RankExpression_W0RawReads_2_noNA,, method = "spearman")
# head(p.mat[, 1:4])

# Plot
ggcorrplot_Spearman_noNA <- corr %>% 
  ggcorrplot(hc.order = F, 
             p.mat = p.mat,
             lab = TRUE, lab_size = 4,
             type = c("lower"),
             outline.col = "white") + 
  my_plot_themes + 
  scale_x_discrete(guide = guide_axis(angle = 45)) + 
  labs(title = "Spearman Correlation Literature sputum. All NA rows removed", 
       subtitle = "n=1932 genes, X's mean not significant", 
       fill = "Correlation")
ggcorrplot_Spearman_noNA

ggsave(ggcorrplot_Spearman_noNA,
       file = "ggcorrplot_Spearman_v2.pdf",
       path = "ggcorrplot_Figures",
       width = 7, height = 6, units = "in")


###########################################################
################### SPEARMAN GGCORRPLOT ###################
# Using all the data, even when there are NAs in some rows... Not sure the right use= for this in the cor() function! 

# Need to make the gene names the rownames! 
AllSputum_RankExpression_W0RawReads_2 <- AllSputum_RankExpression_W0RawReads %>%
  column_to_rownames(var = "Gene")

# Make the correlation
corr <- cor(AllSputum_RankExpression_W0RawReads_2, method = "spearman",
            use = "pairwise.complete.obs")
# Needs the use argument so it will ignore the NAs between each pairwise comparison

# Compute a matrix of correlation p-values
p.mat <- cor_pmat(AllSputum_RankExpression_W0RawReads_2,, method = "spearman")
head(p.mat[, 1:4])

# min(corr) # -0.02468536
# my_min <- round(min(corr), 1)

# Plot
ggcorrplot_Spearman <- corr %>% 
  ggcorrplot(hc.order = F, 
             p.mat = p.mat,
             lab = TRUE, lab_size = 4,
             type = c("lower"),
             outline.col = "white") + 
  # scale_fill_gradient2(limit = c(my_min,1), low = "blue", high =  "red", mid = "white", midpoint = (((1-my_min)/2)+my_min)) + # Make sure to change based on the min!
  my_plot_themes + 
  scale_x_discrete(guide = guide_axis(angle = 45)) + 
  labs(title = "Spearman Correlation Literature sputum", 
       subtitle = "X's mean not significant", 
       fill = "Correlation")
ggcorrplot_Spearman

ggsave(ggcorrplot_Spearman,
       file = "ggcorrplot_Spearman_v1.pdf",
       path = "ggcorrplot_Figures",
       width = 7, height = 6, units = "in")


###########################################################
########### SPEARMAN GGCORRPLOT WITH P-VALUES #############
# The same as above but just more complicated way to generate the correlations

# Using cor.test for significance because that is what Coppola 2021 does:
### "Since the datasets used in our analysis consist of a large sample size, all Spearmans rank correlation coefficients above 0.13 were significant (calculated using the cor.test function in R)."
## Don't really understand what this means....


# Do this a slightly different way for adding to ggcorrplot 
n <- ncol(AllSputum_RankExpression_W0RawReads_2)
corr <- matrix(NA, n, n)
p.mat <- matrix(NA, n, n)
colnames(corr) <- colnames(p.mat) <- colnames(AllSputum_RankExpression_W0RawReads_2)
rownames(corr) <- rownames(p.mat) <- colnames(AllSputum_RankExpression_W0RawReads_2)
for (i in 1:n) {
  for (j in 1:n) {
    test <- cor.test(AllSputum_RankExpression_W0RawReads_2[[i]], AllSputum_RankExpression_W0RawReads_2[[j]], method = "spearman",  use = "pairwise.complete.obs")
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
  # scale_fill_gradient2(limit = c(-0.03,1), low = "blue", high =  "red", mid = "white", midpoint = (((1-(-0.03))/2)+(-0.03))) + # Make sure to change based on the min!
  my_plot_themes + 
  scale_x_discrete(guide = guide_axis(angle = 45)) + 
  labs(title = "Spearman Correlation Literature sputum", 
       subtitle = "X's mean not significant", 
       fill = "Correlation")
ggcorrplot_Spearman_2

# ggsave(ggcorrplot_Spearman_2,
#        file = "ggcorrplot_AllSputum_Spearman_v2.pdf",
#        path = "ggcorrplot_Figures",
#        width = 7, height = 6, units = "in")


###########################################################
###### TOP 15% SPEARMAN GGCORRPLOT WITH P-VALUES ##########

# Just use the top 15% expressed genes (instead of all the genes)
# top15_AllSputum_RankExpression_W0RawReads

# Need to make the gene names the rownames! 
top15_AllSputum_RankExpression_W0RawReads_2 <- top15_AllSputum_RankExpression_W0RawReads %>%
  column_to_rownames(var = "Gene")


# Do this a slightly different way for adding to ggcorrplot - THIS WAY IS CLEANER!!!!
n <- ncol(top15_AllSputum_RankExpression_W0RawReads_2)
corr <- matrix(NA, n, n)
p.mat <- matrix(NA, n, n)
colnames(corr) <- colnames(p.mat) <- colnames(top15_AllSputum_RankExpression_W0RawReads_2)
rownames(corr) <- rownames(p.mat) <- colnames(top15_AllSputum_RankExpression_W0RawReads_2)
for (i in 1:n) {
  for (j in 1:n) {
    test <- cor.test(top15_AllSputum_RankExpression_W0RawReads_2[[i]], top15_AllSputum_RankExpression_W0RawReads_2[[j]], method = "spearman",  use = "pairwise.complete.obs")
    corr[i, j] <- test$estimate
    p.mat[i, j] <- test$p.value
  }
}

ggcorrplot_Spearman_top15 <- corr %>% 
  ggcorrplot(hc.order = F, 
             p.mat = p.mat,
             insig = "pch", 
             sig.level = 0.05,
             lab = TRUE, lab_size = 4,
             type = c("lower"),
             outline.col = "white") + 
  # scale_fill_gradient2(limit = c(-0.03,1), low = "blue", high =  "red", mid = "white", midpoint = (((1-(-0.03))/2)+(-0.03))) + # Make sure to change based on the min!
  my_plot_themes + 
  scale_x_discrete(guide = guide_axis(angle = 45)) + 
  labs(title = "Spearman Correlation Literature sputum; Top 15% expressed genes", 
       subtitle = "X's mean not significant", 
       fill = "Correlation")
ggcorrplot_Spearman_top15












