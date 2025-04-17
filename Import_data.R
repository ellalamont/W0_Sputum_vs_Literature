# Import Data
# 4/16/25

################################################
################ LOAD PACKAGES #################

library(ggplot2)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(knitr)
library(plotly)
library(ggprism) # for add_pvalue()
library(rstatix) # for adjust_pvalue
library(ggpmisc) # https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph
library(ggrepel)
library(pheatmap)
library(dendextend) # May need this for looking at pheatmap clustering
library(ggplotify) # To convert pheatmaps to ggplots
library(corrplot)
library(ggcorrplot)
library(ggfortify) # To make pca plots with plotly

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Biobase")

################################################
#################### COLORS ####################

cbPalette_1 <- c("#999999", "#E69F00") # Gold and Grey
cbPalette_1.5 <- c("#E69F00", "#999999") # Gold and Grey
cbPalette_2 <- c( "#0072B2", "#999999") # Blue and Grey
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette2 <-  c("#bfbfbf", "#56B4E9")
cbPalette3 <-  c("#bfbfbf", "#E69F00")
cbPalette4 <- c("#56B4E9", "#009E73", "#F0E442")
c25 <- c(
  "dodgerblue2", "#E31A1C", "green4",
  "#6A3D9A","#FF7F00","black", "gold1",
  "skyblue2", "#FB9A99","palegreen2","#CAB2D6",
  "#FDBF6F","gray70", "khaki2","maroon", "orchid1", "deeppink1", "blue1", "steelblue4","darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown"
)
c12 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "palegreen2", "gray70", "maroon", "orchid1", "darkturquoise", "darkorange4") 
c16 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black","gold1", "#FB9A99", "#CAB2D6", "palegreen2", "gray70", "maroon", "orchid1", "blue1", "darkturquoise", "darkorange4") 


# Stop scientific notation
options(scipen = 999) 
# options(scipen = 0) # To revert back to default


###########################################################
############ IMPORT AND PROCESS ALL TPM VALUES ############
# NOT scaled

# ProbeTest5_tpm <- read.csv("ProbeTest5_Mtb.Expression.Gene.Data.SCALED.TPM.csv")
ProbeTest5_tpm <- read.csv("EllaData/ProbeTest5_Mtb.Expression.Gene.Data.TPM_moreTrim.csv") # This has the 3' end trimmed 40bp to increase the number of reads aligning
ProbeTest4_tpm <- read.csv("EllaData/ProbeTest4_Mtb.Expression.Gene.Data.TPM.csv")
ProbeTest3_tpm <- read.csv("EllaData/ProbeTest3_Mtb.Expression.Gene.Data.TPM.csv")

# Need to remove the undetermined which all share names
ProbeTest5_tpm$Undetermined_S0 <- NULL
ProbeTest4_tpm$Undetermined_S0 <- NULL
ProbeTest3_tpm$Undetermined_S0 <- NULL

# Merge the 3 documents
All_tpm <- merge(ProbeTest5_tpm, ProbeTest4_tpm, all = T)
All_tpm <- merge(All_tpm, ProbeTest3_tpm, all = T)

# Adjust the names so they are slightly shorter
names(All_tpm) <- gsub(x = names(All_tpm), pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it (I think...)

rownames(All_tpm) <- All_tpm[,1] # add the rownames


###########################################################
#################### SUBSET TPM VALUES ####################
# Just grab the sputum and broth samples that I want

# Unique Sputum: 
# W0 samples: "S_250754", "S_355466", "S_503557" 
# W2 samples: "S_349941_Probe_3D_25", "S_503937", "S_575533_MtbrRNA", "S_577208"
# W4 samples: "S_351946_Probe_4A_100", "S_687338_Probe_4A_100"

# Unique Sputum above 1M reads
# W0 samples: "S_250754", "S_355466", "S_503557" 
# W2 samples: "S_503937", "S_575533_MtbrRNA", "S_577208"

# Uncaptured Ra broth samples
# "H37Ra_Broth_4", "H37Ra_Broth_5", "H37Ra_Broth_6"

# my_sample_names <- c("S_250754", "S_355466", "S_503557", "S_503937", "S_575533_MtbrRNA", "S_577208", "H37Ra_Broth_4", "H37Ra_Broth_5", "H37Ra_Broth_6")
# W0vsBroth_sample_names <- c("S_250754", "S_355466", "S_503557", "H37Ra_Broth_4", "H37Ra_Broth_5", "H37Ra_Broth_6")
W0_Sputum_SampleNames <- c("S_250754", "S_355466", "S_503557")

# my_tpm <- All_tpm %>% select(all_of(my_sample_names))
# my_tpm_W0vsBroth <- All_tpm %>% select(all_of(W0vsBroth_sample_names))
my_tpm_W0Sputum <- All_tpm %>% select(all_of(W0_Sputum_SampleNames))

# Get the average tpm for my 3 week 0s
Average_tpm_W0Sputum <- my_tpm_W0Sputum %>% 
  mutate(Ella_W0_tpm_Average = rowMeans(across(where(is.numeric)))) %>%
  select(Ella_W0_tpm_Average) # %>%
  # rownames_to_column(var = "Gene")

# Rank the averages from 0 to 100?

# Get the broth averages
Average_tpm_broth <- All_tpm %>%
  select(c("H37Ra_Broth_4", "H37Ra_Broth_5", "H37Ra_Broth_6")) %>%
  mutate(Broth_tpm_Average = rowMeans(across(where(is.numeric)))) %>%
  select(Broth_tpm_Average) # %>%



###########################################################
############### IMPORT PIPELINE SUMMARY DATA ##############
# Importing the ProbeTests 3 and 4 and 5 to get all the sputum samples I have done

# This has been edited to include more metadata!
ProbeTest5_pipeSummary <- read.csv("EllaData/ProbeTest5_Pipeline.Summary.Details_moreTrim.csv") # This has the 3' end trimmed 40bp to increase the number of reads aligning
ProbeTest4_pipeSummary <- read.csv("EllaData/ProbeTest4_Pipeline.Summary.Details.csv")
ProbeTest3_pipeSummary <- read.csv("EllaData/ProbeTest3_Pipeline.Summary.Details.csv")

# Merge the 3 documents
All_pipeSummary <- merge(ProbeTest5_pipeSummary, ProbeTest4_pipeSummary, all = T)
All_pipeSummary <- merge(All_pipeSummary, ProbeTest3_pipeSummary, all = T)

All_pipeSummary$X <- NULL
All_pipeSummary$X.1 <- NULL

All_pipeSummary$Hyb_Time <- as.character(All_pipeSummary$Hyb_Time)
ordered_Hyb_Time <- c("4", "16")
All_pipeSummary$Hyb_Time <- factor(All_pipeSummary$Hyb_Time, levels = ordered_Hyb_Time)

All_pipeSummary$Week <- as.character(All_pipeSummary$Week)
ordered_Week <- c("0", "2", "4")
All_pipeSummary$Week <- factor(All_pipeSummary$Week, levels = ordered_Week)

All_pipeSummary$EukrRNADep <- as.character(All_pipeSummary$EukrRNADep)
ordered_EukrRNADep <- c("MtbrRNA", "DualrRNA")
All_pipeSummary$EukrRNADep <- factor(All_pipeSummary$EukrRNADep, levels = ordered_EukrRNADep)

# Remove the undetermined
All_pipeSummary <- All_pipeSummary %>% filter(SampleID != "Undetermined_S0")

# Remove the marmoset and the high low THP1 samples
All_pipeSummary <- All_pipeSummary %>% filter(!Sample_Type %in% c("Marmoset", "High_Low_THP1"))

All_pipeSummary$SampleID <- gsub(x = All_pipeSummary$SampleID, pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it (I think...)

All_pipeSummary <- All_pipeSummary %>% mutate(Sputum_Number = str_extract(SampleID, "S_[0-9]+"))

# Just get the samples I am interested in
my_pipeSummary <- All_pipeSummary %>% filter(SampleID %in% (my_sample_names))
rownames(my_pipeSummary) <- my_pipeSummary$SampleID

# Change the NA to broth
my_pipeSummary$Week <- as.character(my_pipeSummary$Week)
my_pipeSummary$Week[is.na(my_pipeSummary$Week)] <- "Broth"
my_pipeSummary$Week <- as.factor(my_pipeSummary$Week)  # Convert back if needed

my_pipeSummary["Week"]










