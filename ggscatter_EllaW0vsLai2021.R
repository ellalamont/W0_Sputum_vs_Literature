# Correlation scatterplot between Ella W0 and Lai2021
# E. Lamont 
# 5/2/25

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=9), 
        # plot.margin = margin(10, 10, 10, 20),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_blank()
  )

# Stop scientific notation
# options(scipen = 999) 
options(scipen = 0) # To revert back to default


###########################################################
################# PROCESS AND MERGE TPM ###################

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

# Combine my tpm with Lai2021 tpm
merged_tpm <- full_join(
  select(All_tpm, all_of(c(W0_Sputum_SampleNames, Ella_Broth_SampleNames, "X"))),
  select(Lai2021_tpm, all_of(c(Lai2021_SputumNames, Lai2021_ExpNames, "X"))),
  by = "X"
)
merged_tpm_2 <- merged_tpm %>% rename(Gene = X)
rownames(merged_tpm_2) <- merged_tpm_2$Gene

merged_tpm_3 <- merged_tpm_2 %>% select(all_of(c(W0_Sputum_SampleNames, Lai2021_SputumNames, "Gene")))

# Log10 transform the data
merged_tpm_Log10 <- merged_tpm_3 %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values

# Add columns for averages
merged_tpm_Log10 <- merged_tpm_Log10 %>% 
  mutate(
    AVERAGE_EllaW0 = rowMeans(select(., all_of(W0_Sputum_SampleNames)), na.rm = TRUE),
    AVERAGE_Lai2021 = rowMeans(select(., all_of(Lai2021_SputumNames)), na.rm = TRUE),
  )

# Columns for averages not log10 transformed
merged_tpm_4 <- merged_tpm_3 %>% 
  mutate(
    AVERAGE_EllaW0 = rowMeans(select(., all_of(W0_Sputum_SampleNames)), na.rm = TRUE),
    AVERAGE_Lai2021 = rowMeans(select(., all_of(Lai2021_SputumNames)), na.rm = TRUE),
  )


###########################################################
########### TPM SCATTERPLOT CORR WITH LAI2021 #############

Sample1 <- "AVERAGE_EllaW0" # THP1 spiked Captured
Sample2 <- "AVERAGE_Lai2021" # Broth Not Captured
ScatterCorr <- merged_tpm_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.6, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "Pearson correlation; TPM",
       x = paste0("Log10(TPM+1) EllaW0 (n=3)"), y = paste0("Log10(TPM+1) Lai2021 (n=6)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
# ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = paste0("EllaW0_vs_Lai2021_TPM_v1.pdf"),
       path = "CorrelationScatter_Figures",
       width = 7, height = 5, units = "in")


###########################################################
############# PROCESS BATCH CORRECTION CPM ################

source("BatchCorrection_MDS_PCA.R") # for merged_cpm

# Add gene column
merged_cpm_2 <- merged_cpm %>% 
  as.data.frame(merged_cpm) %>% 
  mutate(Gene = rownames(merged_cpm))
rownames(merged_cpm_2) <- merged_cpm_2$Gene

merged_cpm_3 <- merged_cpm_2 %>% select(all_of(c(W0_Sputum_SampleNames, Lai2021_SputumNames, "Gene")))

# Log10 transform the data
merged_cpm_Log10 <- merged_cpm_3 %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values

# Add columns for averages
merged_cpm_Log10 <- merged_cpm_Log10 %>% 
  mutate(
    AVERAGE_EllaW0 = rowMeans(select(., all_of(W0_Sputum_SampleNames)), na.rm = TRUE),
    AVERAGE_Lai2021 = rowMeans(select(., all_of(Lai2021_SputumNames)), na.rm = TRUE),
  )

# Columns for averages not log10 transformed
merged_cpm_4 <- merged_cpm_3 %>% 
  mutate(
    AVERAGE_EllaW0 = rowMeans(select(., all_of(W0_Sputum_SampleNames)), na.rm = TRUE),
    AVERAGE_Lai2021 = rowMeans(select(., all_of(Lai2021_SputumNames)), na.rm = TRUE),
  )


################################################################
###### BATCH CORRECTED CPM SCATTERPLOT CORR WITH LAI2021 #######

Sample1 <- "AVERAGE_EllaW0" # THP1 spiked Captured
Sample2 <- "AVERAGE_Lai2021" # Broth Not Captured
ScatterCorr <- merged_cpm_Log10 %>% 
  ggplot(aes(x = .data[[Sample1]], y = .data[[Sample2]])) + 
  geom_point(aes(text = Gene), alpha = 0.6, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue") + 
  labs(title = paste0(Sample1, " vs ", Sample2),
       subtitle = "Pearson correlation; Batch corrected -> CPM",
       x = paste0("Log10(CPM+1) PredictTB sputum (n=3)"), y = paste0("Log10(CPM+1) Lai2021 sputum (n=6)")) + 
  stat_cor(method="pearson") + # add a correlation to the plot
  my_plot_themes
ScatterCorr
# ggplotly(ScatterCorr)
ggsave(ScatterCorr,
       file = paste0("EllaW0_vs_Lai2021_BC.CPM_v1.jpeg"),
       path = "CorrelationScatter_Figures",
       width = 7, height = 5, units = "in")






