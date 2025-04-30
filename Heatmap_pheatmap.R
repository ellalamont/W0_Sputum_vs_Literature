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


###########################################################
############### PHEATMAP WITH W0 AND LAI2021 ##############

# Combine the tpm and and metadata for the samples that I want to include (broth and W0 only)
# W0 samples: "S_250754", "S_355466", "S_503557"
# Broth samples: "H37Ra_Broth_4", "H37Ra_Broth_5", "H37Ra_Broth_6"
# Lai2021 W0 sputum: "SRR10125319", "SRR10125320", "SRR10125321", "SRR10125322", "SRR10125323", "SRR10125324"
# Lai2021 exponential broth: "SRR10125314", "SRR10125315"

Lai2021_SputumNames <- c("SRR10125319", "SRR10125320", "SRR10125321", "SRR10125322", "SRR10125323", "SRR10125324")
Lai2021_ExpNames <- c("SRR10125314", "SRR10125315")
W0_Sputum_SampleNames <- c("S_250754", "S_355466", "S_503557")
Ella_Broth_SampleNames <- c("H37Ra_Broth_4", "H37Ra_Broth_5", "H37Ra_Broth_6")















