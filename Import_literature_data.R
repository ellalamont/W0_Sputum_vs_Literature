# Import Literature Data
# 4/16/25

# Want to import the literature data and process it in some way so I can compare....
# Maybe something like Coppola et al. (2021)? https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2021.763364/full#h13

# 4/25/25
# Added code to import the DEG data from the literature (compared to broth)
# For DEG, set the log2fold thershold to 2.5, because that's what Garton2008 did. Honeyborne2016 is at 2. But I don't think I can do less than 2.5 because thats all I have for Garton2008

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
# The ratio column in Garcia2016 is 2^(-sputum-Aerobic)
# To get log2fold, take the log2 of the ratio column
# Negative means up in sputum
Garcia2016_DEG <- read.csv("LiteratureData/DEG_ComparedToBroth/Garcia2016_DEG_EL.csv")

# Add the log2fold change column
Garcia2016_DEG <- Garcia2016_DEG %>% mutate(Log2Fold = log2(Ratio))

# Separate the up and down regulated genes
# ** These have been switched because <1 is up in sputum, so want UP to always be UP in sputum
Garcia2016_DEG_DOWN <- Garcia2016_DEG %>%
  filter(Significance == "TRUE") %>% # Need to remove these because they included the not significant
  filter(Log2Fold > 2.5) %>% # Set the log2fold threshold to be 2.5
  select(Gene, Log2Fold)
Garcia2016_DEG_UP <- Garcia2016_DEG %>%
  filter(Significance == "TRUE") %>% # Need to remove these because they included the not significant
  filter(Log2Fold < -2.5) %>% # Set the log2fold threshold to be 2.5
  mutate(Log2Fold = abs(Log2Fold)) %>% # So it's all positive to match other data
  select(Gene, Log2Fold)

################################################
################# SHARMA 2017 ##################
# Actually took the values from what Coppola2021 used
# So I am just using the relative rank score, but I don't want to go back to the raw data until I really know what I am doing and talk to Bob!

# Sharma2017_medianRelativeScoreRank <- read.csv("LiteratureData/Sharma2017_EL.csv")

### DEG Compared to Broth ###
Sharma2017_DEG <- read.csv("LiteratureData/DEG_ComparedToBroth/Sharma2017_DEG_EL.csv")

# Already has been filtered to only be significant changes

# Separate the up and down regulated genes
Sharma2017_DEG_UP <- Sharma2017_DEG %>%
  filter(Fold.Change > 2.5) %>% # Set the log2fold threshold to be 2.5
  select(Name, Fold.Change)
Sharma2017_DEG_DOWN <- Sharma2017_DEG %>%
  filter(Fold.Change < -2.5) %>% # Set the log2fold threshold to be 2.5
  mutate(Fold.Change = abs(Fold.Change)) %>% # So it's all positive to match other data
  select(Name, Fold.Change)



################################################
################### LAI 2021 ###################
# Actually took the values from what Coppola2021 used
# Lai2021_medianRelativeScoreRank <- read.csv("LiteratureData/Lai2021_EL.csv")

### DEG Compared to Broth ###
Lai2021_DEG <- read.csv("LiteratureData/DEG_ComparedToBroth/Lai2021_DEG_EL.csv")
# Already only significant DEG

# Separate the up and down regulated genes
Lai2021_DEG_UP <- Lai2021_DEG %>%
  filter(log2FoldChange > 2.5) %>% # Set the log2fold threshold to be 2.5
  select(Gene, log2FoldChange)
Lai2021_DEG_DOWN <- Lai2021_DEG %>%
  filter(log2FoldChange < -2.5) %>% # Set the log2fold threshold to be 2.5 changes
  mutate(log2FoldChange = abs(log2FoldChange)) %>% # So it's all positive to match other data
  select(Gene, log2FoldChange)

# Not many genes left when using 2.5 log2fold change as the threshold...

# Lai2021 is RNAseq data so I can run it though Bob's pipeline and get the TPM values to compare...
Lai2021_tpm <- read.csv("LiteratureData/Lai2021_BobsPipeline/Mtb.Expression.Gene.Data.TPM_Lai2021.csv")
# Convert column to rownames
rownames(Lai2021_tpm) <- Lai2021_tpm[,1] # add the rownames

Lai2021_SputumNames <- c("SRR10125319", "SRR10125320", "SRR10125321", "SRR10125322", "SRR10125323", "SRR10125324")

# Import the metadata that I made in excel 
Lai2021_metadata <- read.csv("LiteratureData/Lai2021_BobsPipeline/Lai2021_metadata_EL.csv")

################################################
############### HONEYBORNE 2016 ################

### DEG Compared to Broth ###
Honeyborne_DEG <- read.csv("LiteratureData/DEG_ComparedToBroth/HoneyBorne2016_DEG_EL.csv")
# Already filtered to be fold change >2 (not sure if this is log2fold change...)
# Already only significant DEG

# Separate the up and down regulated genes
Honeyborne_DEG_UP <- Honeyborne_DEG %>% 
  filter(Regulation == "up") %>%
  filter(Fold.Change > 2.5) %>% # Set the log2fold threshold to be 2.5
  select(Gene, Fold.Change)
Honeyborne_DEG_DOWN <- Honeyborne_DEG %>% 
  filter(Regulation == "down") %>%
  filter(Fold.Change > 2.5) %>% # Set the log2fold threshold to be 2.5
  select(Gene, Fold.Change)


################################################
################# GARTON 2008 ##################

### DEG Compared to Broth ###
# Should not need to re-run the commented out ones, all for adding Rv numbers instead of gene names, re-load the csv below!
# Garton2008_DEG <- read.csv("LiteratureData/DEG_ComparedToBroth/Garton2008_DEG_EL.csv")

# Need to convert all the common names to Rv numbers
# gene_annot <- read.delim("H37Rv.txt")
# row.names(gene_annot) <- gene_annot$Locus.Tag

# Garton2008_DEG_2 <- Garton2008_DEG %>%
#   mutate(Gene = if_else(
#     grepl("^Rv", Common),
#     Common,
#     gene_annot$Locus.Tag[match(Common, gene_annot$Symbol)] # Grabs the Rv# associated with the gene name
#   ))

# Not all genes matched, will adjust the remainder manually:
# write.csv(Garton2008_DEG_2, "LiteratureData/DEG_ComparedToBroth/Garton2008_DEG_EL_2.csv")
# Went through and manually changed all the NAs in the Gene column to be the correct number
Garton2008_DEG <- read.csv("LiteratureData/DEG_ComparedToBroth/Garton2008_DEG_EL_2.csv")

# Only the significant 2.5 fold change has been kept in this datasheet (not log2?)
# Separate the up and down regulated genes
Garton2008_DEG_UP <- Garton2008_DEG %>% 
  filter(Direction.in.sputum == "UP") %>%
  select(Gene, Sputum)
Garton2008_DEG_DOWN <- Garton2008_DEG %>% 
  filter(Direction.in.sputum == "DOWN") %>%
  select(Gene, Sputum)








#############################################################
######## HIGHLY EXPRESSED GENES ALL FROM COPPOLA 2021 #######
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



