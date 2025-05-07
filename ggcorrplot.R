# Correlation plot with ggcorrplot - RANK EXPRESSION
# E. Lamont
# 4/18/25

# http://www.sthda.com/english/wiki/ggcorrplot-visualization-of-a-correlation-matrix-using-ggplot2

# Coppola 2021 just uses Spearman's, so just do this one (is it because of the way they normalized?)
# Is it bad that there are lots of NAs? (Some of the sputum datasets have a lot fewer genes)

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=16, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=16), 
        plot.subtitle = element_text(size=12), 
        # plot.margin = margin(10, 10, 10, 20),
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
AllSputum_RankExpression_W0RawReads_noNA <- AllSputum_RankExpression_W0RawReads %>% na.omit()
# Now there are 1932 genes... still more than they have, but they had included non-sputum datasets so this may be where the descrepancy is....

# See how many genes (not including NA's are in each column)
AllSputum_RankExpression_W0RawReads %>%
  summarise(across(everything(), ~sum(!is.na(.))))
#   Gene Walter2015 Garcia2016 Sharma2017 Lai2021 EllaW0
# 1 4597       2406       1970       3924    4111   4499

# Not sure why I had to process the data so much here
AllSputum_RankExpression_W0RawReads_noNA <- AllSputum_RankExpression_W0RawReads_noNA %>%
  as.data.frame()
rownames(AllSputum_RankExpression_W0RawReads_noNA) <- NULL
AllSputum_RankExpression_W0RawReads_noNA <- AllSputum_RankExpression_W0RawReads_noNA %>%
  column_to_rownames(var = "Gene")

# Change the name
AllSputum_RankExpression_W0RawReads_noNA <- AllSputum_RankExpression_W0RawReads_noNA %>% rename(PredictTB.W0 = EllaW0)

# Make the correlation
corr <- cor(AllSputum_RankExpression_W0RawReads_noNA, method = "spearman")
# Needs the use argument so it will ignore the NAs between each pairwise comparison

# Compute a matrix of correlation p-values
p.mat <- cor_pmat(AllSputum_RankExpression_W0RawReads_noNA,, method = "spearman")
# head(p.mat[, 1:4])

# Plot
ggcorrplot_Spearman_noNA <- corr %>% 
  ggcorrplot(hc.order = F, 
             p.mat = p.mat,
             lab = TRUE, lab_size = 6,
             type = c("lower"),
             # colors = c("#1B86B4", "white", "#AC204B"),
             legend.title = "Spearman\ncorrelation",
             outline.col = "white") + 
  my_plot_themes + 
  scale_x_discrete(guide = guide_axis(angle = 45)) + 
  labs(title = "Spearman Correlation Literature sputum. All NA rows removed", 
       subtitle = "n = 1932 genes, X = not significant", 
       x = NULL, y = NULL)
ggcorrplot_Spearman_noNA

ggsave(ggcorrplot_Spearman_noNA,
       file = "ggcorrplot_Spearman_v3.pdf",
       path = "ggcorrplot_Figures",
       width = 7, height = 6, units = "in")


###########################################################
################ GGCORRPLOT W0 vs LAI2021 #################

# What does the rank expression correlation look like when I can use all the data for these two datasets?

# Subset just the W0 and Lai2021
RankExpression_subset <- AllSputum_RankExpression_W0RawReads %>% select(all_of(c("Gene", "Lai2021", "EllaW0"))) # 4597 rows

# Remove any row that contains NA
RankExpression_subset_noNA <- RankExpression_subset %>% na.omit()
# Now there are 4026 genes

# Not sure why I had to process the data so much here
RankExpression_subset_noNA <- RankExpression_subset_noNA %>%
  as.data.frame()
rownames(RankExpression_subset_noNA) <- NULL
RankExpression_subset_noNA <- RankExpression_subset_noNA %>% 
  column_to_rownames(var = "Gene")

# Make the correlation
cor(RankExpression_subset_noNA, method = "spearman")
#           Lai2021    EllaW0
# Lai2021 1.0000000 0.7943109
# EllaW0  0.7943109 1.0000000
# Good, this is the same as above

# Compute a matrix of correlation p-values
p.mat <- cor_pmat(RankExpression_subset_noNA,, method = "spearman")
# head(p.mat[, 1:4])










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



