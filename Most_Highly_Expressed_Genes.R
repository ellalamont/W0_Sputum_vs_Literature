# Look at most highly expressed genes between sputum transcriptomes
# 4/16/25

source("Import_data.R")
source("Import_literature_data.R")


###########################################################
################### MY BROTH TPM DATA #####################
# Average_tpm_broth
# # 4499 genes total, top 15% would be the top 674 genes

Broth_top_tpm_genes <- Average_tpm_broth %>%
  arrange(desc(Broth_tpm_Average)) %>%
  slice_head(n = 674) # %>%
# select(Gene)
write.csv(Broth_top_tpm_genes,
          file = "Top15PercentGeneLists/EllaBroth_top15Percent_genes.csv")



###########################################################
##################### MY W0 TPM DATA ######################
# Average_tpm_W0Sputum
# 4499 genes total, top 15% would be the top 674 genes

Ella_top_tpm_genes <- Average_tpm_W0Sputum %>%
  arrange(desc(Ella_W0_tpm_Average)) %>%
  slice_head(n = 674) # %>%
  # select(Gene)
write.csv(Ella_top_tpm_genes,
          file = "Top15PercentGeneLists/EllaW0Sputum_top15Percent_genes.csv")

###########################################################
########### WALTER2015 W0 MEDIAN EXPRESSION DATA ##########
# Walter2015_medianExpression
# 2412 genes total, top 15% would be the top 361

Walter2015_top_tpm_genes <- Walter2015_medianExpression %>%
  arrange(desc(Day0_MedianExpression)) %>%
  slice_head(n = 361)
write.csv(Walter2015_top_tpm_genes,
          file = "Top15PercentGeneLists/Walter2015_top15Percent_genes.csv")


###########################################################
############### GARCIA2016 W0 MEDIAN CT DATA ##############
# Garcia2016_medianCt
# 1970 genes total, top 15% would be 296

Garcia2016_top_tpm_genes <- Garcia2016_medianCt %>%
  arrange(Median.CT.Sputum) %>% # Needs to be ascending here because lower Ct is higher expression
  slice_head(n = 296)
write.csv(Garcia2016_top_tpm_genes,
          file = "Top15PercentGeneLists/Garcia2016_top_tpm_genes.csv")


###########################################################
######### SHARMA2017 W0 MEDIAN REALTIVE SCORE RANK ########
# Sharma2017_medianRelativeScoreRank
# 2924 genes total, top 15% would be 442

Sharma2017_top_tpm_genes <- Sharma2017_medianRelativeScoreRank %>%
  arrange(desc(Median_Relative_Score_Rank)) %>% # Needs to be ascending here because lower Ct is higher expression
  slice_head(n = 442)
write.csv(Sharma2017_top_tpm_genes,
          file = "Top15PercentGeneLists/Sharma2017_top_tpm_genes.csv")






