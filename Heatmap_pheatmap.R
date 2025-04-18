# Try a heatmap with pheatmap
# E. Lamont 
# 4/18/25

# https://rpubs.com/tgjohnst/heatmaps_testing_1

source("Most_Highly_Expressed_Genes.R") # To get AllSputum_RankExpression
# Using rank expression here, NOT TPM!

my_annotation_colors <- list(
  Week = c("0" = "#0072B2",  # Blue
           "2" = "#E66900",  # Orange
           "Broth" = "#999999")  # Grey
)


###########################################################
###################### PROCESS DATA #######################

# Need to make the gene names the rownames! 
AllSputum_RankExpression_2 <- AllSputum_RankExpression %>%
  column_to_rownames(var = "Gene")

# Just take the top 50 genes to make things easier when figuring it out
AllSputum_RankExpression_filtered20 <- AllSputum_RankExpression_2 %>%
  slice_head(n = 20)

###########################################################
######################## PHEATMAP #########################

pheatmap(AllSputum_RankExpression_filtered20, 
         # annotation_row = Gene_Category, 
         # annotation_col = my_pipeSummary["Week"],
         # annotation_colors = my_annotation_colors,
         scale = "none",
         display_numbers = T)












