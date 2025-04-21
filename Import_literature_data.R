# Import Literature Data
# 4/16/25

# Lets not overthink this until we talk to Bob!!

# Want to import the literature data and process it in some way so I can compare....
# Maybe something like Coppola et al. (2021)? https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2021.763364/full#h13

################################################
################# WALTER 2015 ##################
# Raw data is median expression at Day 0....
# I think everthing is ranked between 0 and 100...
# https://pmc.ncbi.nlm.nih.gov/articles/PMC4548467/#JIV149C20

# Walter2015_medianExpression <- read.csv("LiteratureData/Walter2015_EL.csv")

################################################
################# GARCIA 2016 ##################
# Data is Ct values, also there seem to be two different Ct values for the sputum, Coppola2021 is using the Ct values from the first tab I think
# Actually took the values from what Coppola2021 used

# Garcia2016_medianCt <- read.csv("LiteratureData/Garcia2016_EL.csv")


################################################
################# SHARMA 2017 ##################
# Actually took the values from what Coppola2021 used
# So I am just using the relative rank score, but I don't want to go back to the raw data until I really know what I am doing and talk to Bob!

# Sharma2017_medianRelativeScoreRank <- read.csv("LiteratureData/Sharma2017_EL.csv")


################################################
################### LAI 2021 ###################
# Actually took the values from what Coppola2021 used
# So I am just using the relative rank score, but I don't want to go back to the raw data until I really know what I am doing and talk to Bob!

Lai2021_medianRelativeScoreRank <- read.csv("LiteratureData/Lai2021_EL.csv")


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



