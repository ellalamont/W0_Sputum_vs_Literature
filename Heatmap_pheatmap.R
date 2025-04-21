# Try a heatmap with pheatmap
# E. Lamont 
# 4/18/25

# https://rpubs.com/tgjohnst/heatmaps_testing_1

source("Most_Highly_Expressed_Genes.R") # To get AllSputum_RankExpression
# Using rank expression here, NOT TPM!


###########################################################
###################### PROCESS DATA #######################

# Need to make the gene names the rownames! 
AllSputum_RankExpression_2 <- AllSputum_RankExpression %>%
  column_to_rownames(var = "Gene")

# Just take the top 20 genes to make things easier when figuring it out
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



###########################################################
################ PHEATMAP WITH ALL GENES ##################

# Remove rows where there are at least 2 NAs, hope this will help with clustering
AllSputum_RankExpression_2_filtered <- AllSputum_RankExpression_2 %>% filter(rowSums(is.na(.)) < 2) # Trimmed down to 2383 rows

pheatmap(AllSputum_RankExpression_2_filtered, 
         # annotation_row = Gene_Category, 
         # annotation_col = my_pipeSummary["Week"],
         # annotation_colors = my_annotation_colors,
         scale = "none",
         display_numbers = F)

# 4/21/25: The Lai2021 is super different... many low numbers, not many high numbers... what is going on here??
## Because there were so many not detected or lowly detected genes in this RNAseq dataset... Mentioned in the Coppola2021 paper but not really explained how they dealt with it....


