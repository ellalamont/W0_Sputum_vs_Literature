# Batch Correction and MDS plot, following what Mark has done....
# E. Lamont
# 5/1/25

# Need Mark's help for understanding what he did for the MDS plot

source("Import_data.R")
source("Import_literature_data.R")

# BiocManager::install("sva")
library(sva) # For ComBat_seq batch correction
library(edgeR)


###########################################################
###################### PROCESS DATA #######################

# Mark starts with READS_M for this
# All_RawReads
# Lai2021_RawReads

Lai2021_SputumNames <- c("SRR10125319", "SRR10125320", "SRR10125321", "SRR10125322", "SRR10125323", "SRR10125324")
Lai2021_ExpNames <- c("SRR10125314", "SRR10125315")
W0_Sputum_SampleNames <- c("S_250754", "S_355466", "S_503557")
Ella_Broth_SampleNames <- c("H37Ra_Broth_4", "H37Ra_Broth_5", "H37Ra_Broth_6")

# Combine my metadata with Lai2021 metadata
# I Only need SampleID and Sample_Type
merged_metadata <- bind_rows(
  select(my_pipeSummary, all_of(c("SampleID", "Sample_Type"))),
  select(Lai2021_metadata, all_of(c("SampleID", "Sample_Type")))
)
rownames(merged_metadata) <- merged_metadata$SampleID

# Add Batch designation to the metadata:
merged_metadata <- merged_metadata %>%
  mutate(Batch = case_when(
    Sample_Type %in% c("Broth", "Sputum") ~ 1,
    Sample_Type %in% c("Rv_Exponential", "Rv_Stationary", "TBSputum_nonHIV") ~ 2
  ))

merged_metadata_subset <- merged_metadata %>% 
 filter(SampleID %in% c(W0_Sputum_SampleNames, Ella_Broth_SampleNames, Lai2021_SputumNames, Lai2021_ExpNames))

# Combine my tpm with Lai2021 tpm
merged_RawReads <- full_join(
  select(All_RawReads, all_of(c(W0_Sputum_SampleNames, Ella_Broth_SampleNames, "X"))),
  select(Lai2021_RawReads, all_of(c(Lai2021_SputumNames, Lai2021_ExpNames, "X"))),
  by = "X"
)
merged_RawReads_rn <- merged_RawReads %>%
  column_to_rownames(var = "X") %>%
  as.matrix()


###########################################################
#################### BATCH CORRECTION #####################
# This basically all taken from Mark and I'm not sure what's happening

# https://academic.oup.com/nargab/article/2/3/lqaa078/5909519?login=true


batch <- merged_metadata_subset$Batch
counts_corrected <- ComBat_seq(merged_RawReads_rn, batch = batch)


###########################################################
################### MDS PLOT FROM MARK ####################
# This basically all taken from Mark and I'm not sure what's happening
# Think there are a few normalization steps going on here....

y <- DGEList(counts_corrected)
z <- DGEList(merged_RawReads_rn)

# Adding in some metadata?
y$samples$group <- factor(merged_metadata_subset$Sample_Type)
y$samples$batch <- factor(merged_metadata_subset$Batch)

logcounts <- cpm(y,log=TRUE)
logcounts_notcorrected <- cpm(z, log=TRUE)

y <- calcNormFactors(y)
z <- calcNormFactors(z)
design <- model.matrix(~ 0 + group) # STUCK! What is the model.matrix? and group not found?
colnames(design) <- levels(group)


###########################################################
########################### PCA ###########################
# Just do the PCA the way I have been doing it

# Convert the batch corrected counts to counts per million (cpm)
merged_cpm <- cpm(counts_corrected)
# merged_cpm_log2 <- cpm(counts_corrected, log = T)

# transform the data 
merged_cpm_t <- as.data.frame(t(merged_cpm))
# merged_cpm_log2_t <- as.data.frame(t(merged_cpm_log2))

# Remove columns that are all zero so the scale works for prcomp
merged_cpm_t <- merged_cpm_t %>% select_if(colSums(.) != 0) # Down to 4494 genes 
# merged_cpm_log2_t <- merged_cpm_log2_t %>% select_if(colSums(.) != 0) # Stays at 4499 genes

# Make the actual PCA
# Need to choose here if using log2 or not
my_PCA <- prcomp(merged_cpm_t, scale = TRUE)
# my_PCA <- prcomp(merged_cpm_log2_t, scale = TRUE) #. This is throwing an error, not sure why

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 24.2% of variance
summary_PCA[2,1] # PC2 explains 17.0% of variance
summary_PCA[3,1] # PC3 explains 13.0% of variance

# Add the metadata
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, merged_metadata, by = "SampleID", )

my_PCA_df <- my_PCA_df %>%
  mutate(Sample_Type2 = case_when(Sample_Type == "Broth" ~ "PredictTB broth",
                                  Sample_Type == "Rv_Exponential"~ "Lai2021 broth",
                                  Sample_Type == "Sputum"~ "PredictTB sputum",
                                  Sample_Type == "TBSputum_nonHIV" ~ "Lai2021 sputum"))

# MAKE PCA PLOT with GGPLOT 
fig_PC1vsPC2 <- my_PCA_df %>%
  ggplot(aes(x = PC1, y = PC2, fill = Sample_Type2, shape = Sample_Type2)) + 
  geom_point(size = 6, alpha = 0.9, stroke = 0.8) +
  scale_fill_manual(values=c(`PredictTB broth` = "#999999", `Lai2021 broth` = "lightgrey", `PredictTB sputum`= "#0072B2", `Lai2021 sputum` = "#C6DDF0")) +  
  scale_shape_manual(values=c(`PredictTB sputum` = 21, `Lai2021 sputum` = 21, `PredictTB broth`= 23, `Lai2021 broth` = 23)) + 
  # geom_text_repel(aes(label = EukrRNADep), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(title = "Ella W0 and H37Ra broth, with Lai2021 sputum and exponential Rv broth",
       subtitle = "Batch corrected raw reads -> CPM",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%"),
       fill = NULL, shape = NULL) +
  my_plot_themes
fig_PC1vsPC2
ggsave(fig_PC1vsPC2,
       file = "PCA_EllavsLai2021_BatchCorrect.CPM_v2.pdf",
       path = "PCA_Figures",
       width = 8, height = 5, units = "in")
ggsave(fig_PC1vsPC2,
       file = "PCA_EllavsLai2021_BatchCorrect.CPM_v2.png",
       path = "PCA_Figures",
       dpi = 300,
       width = 8, height = 5, units = "in")


###########################################################
####################### MDS MY WAY ########################

# Convert the batch corrected counts to counts per million (cpm)
merged_cpm <- cpm(counts_corrected)

# transform the data 
merged_cpm_t <- as.data.frame(t(merged_cpm))

# Cmpute MDS
mds <- cmdscale(dist(merged_cpm_t))
colnames(mds) <- c("Dim.1", "Dim.2")

# Add the metadata, just hoping that the rows are in the same order....??
# Add the metadata
my_MDS_df <- as.data.frame(mds) %>%
  rownames_to_column(var = "SampleID")
my_MDS_df <- merge(my_MDS_df, merged_metadata_subset, by = "SampleID")


fig_MDS <- my_MDS_df %>%
  ggplot(aes(x = Dim.1, y = Dim.2, fill = Sample_Type, shape = Sample_Type)) + 
  geom_point(size = 6, alpha = 0.9, stroke = 0.8) +
  scale_fill_manual(values=c(`Broth` = "#999999", `Rv_Exponential` = "lightgrey", `Sputum`= "#0072B2", `TBSputum_nonHIV` = "#C6DDF0")) +  
  scale_shape_manual(values=c(`Sputum` = 21, `TBSputum_nonHIV` = 22, `Broth`= 24, `Rv_Exponential` = 25)) + 
  # geom_text_repel(aes(label = EukrRNADep), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(title = "MDS: Ella W0 and H37Ra broth, with Lai2021 sputum and exponential Rv broth",
       subtitle = "Batch corrected raw reads -> CPM",
       x = "Dimension 1",
       y = "Dimension 2") +
  my_plot_themes
# my_plot_themes_thumbnail
# fig_MDS
# ggsave(fig_MDS,
#        file = "MDS_EllavsLai2021_BatchCorrect.CPM_v1.pdf",
#        # file = "PCA_EllavsLai2021_BatchCorrect.LOG2.CPM_v1.pdf",
#        path = "MDS_Figures",
#        width = 8, height = 5, units = "in")






