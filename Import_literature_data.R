# Import Literature Data
# 4/16/25

# Lets not overthink this until we talk to Bob!!

# Want to import the literature data and process it in some way so I can compare....
# Maybe something like Coppola et al. (2021)? https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2021.763364/full#h13

# 4/25/25
# Added code to import the DEG data from the literature (compared to broth)

################################################
################# WALTER 2015 ##################
# Raw data is median expression at Day 0....
# I think everthing is ranked between 0 and 100...
# https://pmc.ncbi.nlm.nih.gov/articles/PMC4548467/#JIV149C20

# Walter2015_medianExpression <- read.csv("LiteratureData/Walter2015_EL.csv")

### DEG Compared to Broth ###
# This paper doesn't have this 

################################################
################# GARCIA 2016 ##################
# Data is Ct values, also there seem to be two different Ct values for the sputum, Coppola2021 is using the Ct values from the first tab I think
# Actually took the values from what Coppola2021 used

# Garcia2016_medianCt <- read.csv("LiteratureData/Garcia2016_EL.csv")

### DEG Compared to Broth ###
# The ratio column in Garcia2016 is 2^(-sputum-Aerobic) which is the fold change?
# I'm not sure if I need to convert this to log2fold change or keep as is.....
# ***<1 is up in sputum, >1 is up in broth****
Garcia2016_DEG <- read.csv("LiteratureData/DEG_ComparedToBroth/Garcia2016_DEG_EL.csv")

# Separate the up and down regulated genes
# ** These have been switched because <1 is up in sputum, so want UP to always be UP in sputum
Garcia2016_DEG_DOWN <- Garcia2016_DEG %>%
  filter(Significance == "TRUE") %>% # Need to remove these because they included the not significant
  filter(Ratio > 1) %>% # THIS IS NOT FOLD CHANGE!!!!
  select(Gene, Ratio)
Garcia2016_DEG_UP <- Garcia2016_DEG %>%
  filter(Significance == "TRUE") %>% # Need to remove these because they included the not significant
  filter(Ratio <= 1) %>% # THIS IS NOT FOLD CHANGE!!!!
  select(Gene, Ratio)

################################################
################# SHARMA 2017 ##################
# Actually took the values from what Coppola2021 used
# So I am just using the relative rank score, but I don't want to go back to the raw data until I really know what I am doing and talk to Bob!

# Sharma2017_medianRelativeScoreRank <- read.csv("LiteratureData/Sharma2017_EL.csv")

### DEG Compared to Broth ###
Sharma2017_DEG <- read.csv("LiteratureData/DEG_ComparedToBroth/Sharma2017_DEG_EL.csv")

# Separate the up and down regulated genes
Sharma2017_DEG_UP <- Sharma2017_DEG %>%
  filter(Fold.Change > 1) %>% #  Setting it to 1 because we will ignore the small fold changes
  select(Name, Fold.Change)
Sharma2017_DEG_DOWN <- Sharma2017_DEG %>%
  filter(Fold.Change < -1) %>% # Setting it to -1 because we will ignore the small fold changes
  mutate(Fold.Change = abs(Fold.Change)) %>% # So it's all positive to match other data
  select(Name, Fold.Change)



################################################
################### LAI 2021 ###################
# Actually took the values from what Coppola2021 used
# So I am just using the relative rank score, but I don't want to go back to the raw data until I really know what I am doing and talk to Bob!

Lai2021_medianRelativeScoreRank <- read.csv("LiteratureData/Lai2021_EL.csv")

### DEG Compared to Broth ###
Lai2021_DEG <- read.csv("LiteratureData/DEG_ComparedToBroth/Lai2021_DEG_EL.csv")

# Separate the up and down regulated genes
Lai2021_DEG_UP <- Lai2021_DEG %>%
  filter(log2FoldChange > 1) %>% # Setting it to 1 because we will ignore the small fold changes
  select(Gene, log2FoldChange)
Lai2021_DEG_DOWN <- Lai2021_DEG %>%
  filter(log2FoldChange < -1) %>% # Setting it to -1 because we will ignore the small fold changes
  mutate(log2FoldChange = abs(log2FoldChange)) %>% # So it's all positive to match other data
  select(Gene, log2FoldChange)


################################################
############### HONEYBORNE 2016 ################

### DEG Compared to Broth ###
Honeyborne_DEG <- read.csv("LiteratureData/DEG_ComparedToBroth/HoneyBorne2016_DEG_EL.csv")

# Separate the up and down regulated genes
Honeyborne_DEG_UP <- Honeyborne_DEG %>% 
  filter(Regulation == "up") %>%
  select(Gene, Fold.Change)
Honeyborne_DEG_DOWN <- Honeyborne_DEG %>% 
  filter(Regulation == "down") %>%
  select(Gene, Fold.Change)


################################################
############# ALL FROM COPPOLA 2021 ############
# Just take from the supplemental of Coppola2021 for all so it will be more consistent! 

Coppola2021_all <- read.csv("LiteratureData/Coppola2021_EL.csv")
Coppola2021_Walter2015 <- Coppola2021_all %>% 
  select(Walter2015_Gene, Walter2015_MedianRelativeScoreRank) %>% 
  drop_na()
Coppola2021_Garcia2016 <- Coppola2021_all %>% 
  select(Garcia2016_Gene, Garcia2016_PercentileRelativeScoreRank) %>% 
  drop_na()
Coppola2021_Sharma2017 <- Coppola2021_all %>% 
  select(Sharma2017_Gene, Sharma2017_MedianRelativeScoreRank) %>% 
  drop_na()
Coppola2021_Lai2021 <- Coppola2021_all %>% 
  select(Lai2021_Gene, Lai2021_MedianRelativeScoreRank) %>% 
  drop_na()



################################################
############### LOOK AT HISTOGRAMS #############
# hist(Coppola2021_Walter2015$Walter2015_MedianRelativeScoreRank)
# hist(Coppola2021_Garcia2016$Garcia2016_PercentileRelativeScoreRank)
# hist(Coppola2021_Sharma2017$Sharma2017_MedianRelativeScoreRank)
# hist(Coppola2021_Lai2021$Lai2021_MedianRelativeScoreRank)

# They all look good except for Lai! Why is this not normalized?
## Because there were so many not detected or lowly detected genes in this RNAseq dataset... Mentioned in the Coppola2021 paper but not really explained how they dealt with it....



