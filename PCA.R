# PCA plot with my W0 and broth and Lai2021 W0 and broth
# E. Lamont
# 4/30/25

# Look into ggbiplot for more PCA stuff??
# https://cran.r-project.org/web/packages/ggbiplot/readme/README.html

# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

# Two options in base R, prcomp() and princomp()
# prcomp() is preferred according to the website above

source("Import_data.R")
source("Import_literature_data.R")

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
        plot.margin = margin(10, 10, 10, 20),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_blank()
  )


Lai2021_SputumNames <- c("SRR10125319", "SRR10125320", "SRR10125321", "SRR10125322", "SRR10125323", "SRR10125324")
Lai2021_ExpNames <- c("SRR10125314", "SRR10125315")
W0_Sputum_SampleNames <- c("S_250754", "S_355466", "S_503557")
Ella_Broth_SampleNames <- c("H37Ra_Broth_4", "H37Ra_Broth_5", "H37Ra_Broth_6")

###########################################################
################# ALL W0, BROTH, LAI2021 ##################

# Combine my metadata with Lai2021 metadata
# I Only need SampleID and Sample_Type
merged_metadata <- bind_rows(
  select(my_pipeSummary, all_of(c("SampleID", "Sample_Type"))),
  select(Lai2021_metadata, all_of(c("SampleID", "Sample_Type")))
)

# Combine my tpm with Lai2021 tpm
merged_tpm <- full_join(
  select(All_tpm, all_of(c(W0_Sputum_SampleNames, Ella_Broth_SampleNames, "X"))),
  select(Lai2021_tpm, all_of(c(Lai2021_SputumNames, Lai2021_ExpNames, "X"))),
  by = "X"
  )

# add rownames to the tpm 
merged_tpm <- merged_tpm %>% column_to_rownames(var = "X")

# transform the data 
merged_tpm_t <- as.data.frame(t(merged_tpm))

# Remove columns that are all zero so the scale works for prcomp
merged_tpm_t <- merged_tpm_t %>% select_if(colSums(.) != 0) # Down to 4494 genes 

# Make the actual PCA
my_PCA <- prcomp(merged_tpm_t, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 32.2% of variance
summary_PCA[2,1] # PC2 explains 16.0% of variance
summary_PCA[3,1] # PC3 explains 11.8% of variance

# Add the metadata
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, merged_metadata, by = "SampleID", )

# MAKE PCA PLOT with GGPLOT 
fig_PC1vsPC2 <- my_PCA_df %>%
  ggplot(aes(x = PC1, y = PC2, fill = Sample_Type, shape = Sample_Type)) + 
  geom_point(size = 6, alpha = 0.9, stroke = 0.8) +
  scale_fill_manual(values=c(`Broth` = "#999999", `Rv_Exponential` = "lightgrey", `Sputum`= "#0072B2", `TBSputum_nonHIV` = "#C6DDF0")) +  
  scale_shape_manual(values=c(`Sputum` = 21, `TBSputum_nonHIV` = 22, `Broth`= 24, `Rv_Exponential` = 25)) + 
  # geom_text_repel(aes(label = EukrRNADep), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(title = "Ella W0 and H37Ra broth, with Lai2021 sputum and exponential Rv broth",
       subtitle = "Using TPM",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
# my_plot_themes_thumbnail
fig_PC1vsPC2
# ggplotly(fig_PC1vsPC2)
ggsave(fig_PC1vsPC2,
       file = "PCA_EllavsLai2021_TPM_v1.pdf",
       path = "PCA_Figures",
       width = 8, height = 5, units = "in")


# 3D plot
# https://plotly.com/r/pca-visualization/
PCA_3D <- plot_ly(my_PCA_df, x = ~PC1, y = ~PC2, z = ~PC3,
                  type = "scatter3d", mode = "markers",
                  color = ~Sample_Type# , 
                  # colors = c12,
                  # text = ~Replicate
)
# PCA_3D
# htmlwidgets::saveWidget(as_widget(PCA_3D), "PCA_3D.html")

###########################################################
######## LOG TRANSFORMED ALL W0, BROTH, LAI2021 ###########

# Combine my metadata with Lai2021 metadata
# I Only need SampleID and Sample_Type
merged_metadata <- bind_rows(
  select(my_pipeSummary, all_of(c("SampleID", "Sample_Type"))),
  select(Lai2021_metadata, all_of(c("SampleID", "Sample_Type")))
)

# Combine my tpm with Lai2021 tpm
merged_tpm <- full_join(
  select(All_tpm, all_of(c(W0_Sputum_SampleNames, Ella_Broth_SampleNames, "X"))),
  select(Lai2021_tpm, all_of(c(Lai2021_SputumNames, Lai2021_ExpNames, "X"))),
  by = "X"
)

# add rownames to the tpm 
merged_tpm <- merged_tpm %>% column_to_rownames(var = "X")

# Log10 transform the data
merged_tpm_log10 <- merged_tpm %>% 
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values

# transform the data 
merged_tpm_t <- as.data.frame(t(merged_tpm_log10))

# Remove columns that are all zero so the scale works for prcomp
merged_tpm_t <- merged_tpm_t %>% select_if(colSums(.) != 0) # Down to 4494 genes 

# Make the actual PCA
my_PCA <- prcomp(merged_tpm_t, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 25.5% of variance
summary_PCA[2,1] # PC2 explains 19.7% of variance
summary_PCA[3,1] # PC3 explains 16.0% of variance

# Add the metadata
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, merged_metadata, by = "SampleID", )

# MAKE PCA PLOT with GGPLOT 
fig_PC1vsPC2 <- my_PCA_df %>%
  ggplot(aes(x = PC1, y = PC2, fill = Sample_Type, shape = Sample_Type)) + 
  geom_point(size = 6, alpha = 0.9, stroke = 0.8) +
  scale_fill_manual(values=c(`Broth` = "#999999", `Rv_Exponential` = "lightgrey", `Sputum`= "#0072B2", `TBSputum_nonHIV` = "#C6DDF0")) +  
  scale_shape_manual(values=c(`Sputum` = 21, `TBSputum_nonHIV` = 22, `Broth`= 24, `Rv_Exponential` = 25)) + 
  # geom_text_repel(aes(label = EukrRNADep), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  labs(title = "Ella W0 and H37Ra broth, with Lai2021 sputum and exponential Rv broth",
       subtitle = "Using log10(TPM+1)",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
# my_plot_themes_thumbnail
fig_PC1vsPC2
ggsave(fig_PC1vsPC2,
       file = "PCA_EllavsLai2021_LOG10.TPM_v1.pdf",
       path = "PCA_Figures",
       width = 8, height = 5, units = "in")


