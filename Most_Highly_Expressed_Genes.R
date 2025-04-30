# Look at most highly expressed genes between sputum transcriptomes
# 4/16/25

source("Import_data.R")
source("Import_literature_data.R")
source("Quantile_Normalization.R")


###########################################################
################## MY W0 READS_M DATA #####################
# 4499 genes total, top 15% would be the top 674 genes

# Do this with the quantile normalization ranking think I did
# RANK_Average_RawReads_W0Sputum From Quantile Normalization! 
RANK_W0_RawReads_top15 <- RANK_Average_RawReads_W0Sputum %>%
  arrange(desc(RANK_Average_RawReads_W0Sputum)) %>%
  slice_head(n = round(nrow(RANK_Average_RawReads_W0Sputum)*0.15)) 
RANK_W0_RawReads_top15_GeneNames <- RANK_W0_RawReads_top15 %>%
  rownames()


###########################################################
############### USE DIRECT FROM COPPOLA2021 ###############
# I want the entire list here, not just the 15%, so remake everything and maybe change the colunm names while I'm at it 

# Coppola2021_Walter2015
Coppola2021_Walter2015_top15 <- Coppola2021_Walter2015 %>% 
  arrange(desc(Walter2015_MedianRelativeScoreRank)) %>% 
  slice_head(n = round(nrow(Coppola2021_Walter2015)*0.15))
Coppola2021_Walter2015_top15_GeneNames <- Coppola2021_Walter2015_top15 %>%
  pull(Walter2015_Gene)

# Coppola2021_Garcia2016
Coppola2021_Garcia2016_top15 <- Coppola2021_Garcia2016 %>% 
  arrange(desc(Garcia2016_PercentileRelativeScoreRank)) %>% 
  slice_head(n = round(nrow(Coppola2021_Garcia2016)*0.15)) 
Coppola2021_Garcia2016_top15_GeneNames <- Coppola2021_Garcia2016_top15 %>%
  pull(Garcia2016_Gene)

# Coppola2021_Sharma2017
Coppola2021_Sharma2017_top15 <- Coppola2021_Sharma2017 %>% 
  arrange(desc(Sharma2017_MedianRelativeScoreRank)) %>% 
  slice_head(n = round(nrow(Coppola2021_Sharma2017)*0.15)) 
Coppola2021_Sharma2017_top15_GeneNames <- Coppola2021_Sharma2017_top15 %>%
  pull(Sharma2017_Gene)

# Coppola2021_Lai2021
Coppola2021_Lai2021_top15 <- Coppola2021_Lai2021 %>% 
  arrange(desc(Lai2021_MedianRelativeScoreRank)) %>% 
  slice_head(n = round(nrow(Coppola2021_Lai2021)*0.15))
Coppola2021_Lai2021_top15_GeneNames <- Coppola2021_Lai2021_top15 %>%
  pull(Lai2021_Gene)

# Add my data


###########################################################
############# PUT EVERTHING IN ONE DATAFRAME ##############

combined_df <- full_join(Coppola2021_Walter2015 %>% rename(Gene = Walter2015_Gene),
                         Coppola2021_Garcia2016 %>% rename(Gene = Garcia2016_Gene),
                         by = "Gene")
combined_df <- full_join(combined_df,
                         Coppola2021_Sharma2017 %>% rename(Gene = Sharma2017_Gene),
                         by = "Gene")
combined_df <- full_join(combined_df,
                         Coppola2021_Lai2021 %>% rename(Gene = Lai2021_Gene),
                         by = "Gene")

AllSputum_RankExpression_W0RawReads <- full_join(combined_df,
                                     RANK_Average_RawReads_W0Sputum %>% rownames_to_column("Gene"),
                                    by = "Gene") %>%
  rename(EllaW0_RawReads_RankAverage = RANK_Average) %>%
  rename_with(~ str_replace(., "_.*", ""))



###########################################################
############## PUT TOP 15% IN ONE DATAFRAME ###############

top15_combined_df <- full_join(Coppola2021_Walter2015_top15 %>% rename(Gene = Walter2015_Gene),
                         Coppola2021_Garcia2016_top15 %>% rename(Gene = Garcia2016_Gene),
                         by = "Gene")
top15_combined_df <- full_join(top15_combined_df,
                         Coppola2021_Sharma2017_top15 %>% rename(Gene = Sharma2017_Gene),
                         by = "Gene")
top15_combined_df <- full_join(top15_combined_df,
                         Coppola2021_Lai2021_top15 %>% rename(Gene = Lai2021_Gene),
                         by = "Gene")
top15_AllSputum_RankExpression_W0RawReads <- full_join(top15_combined_df,
                                                       RANK_W0_RawReads_top15 %>% rownames_to_column("Gene"),
                                                 by = "Gene") %>%
  rename(EllaW0_RawReads_RankAverage = RANK_Average) %>%
  rename_with(~ str_replace(., "_.*", ""))





