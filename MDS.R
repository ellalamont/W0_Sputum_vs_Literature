# MDS plot with my W0 and broth and Lai2021 W0 and broth
# E. Lamont
# 5/1/25

# https://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/122-multidimensional-scaling-essentials-algorithms-and-r-code/
# https://www.geeksforgeeks.org/multidimensional-scaling-using-r/



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


###########################################################
###################### PROCESS DATA #######################

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
merged_metadata_subset <- merged_metadata %>% 
  filter(SampleID %in% c(W0_Sputum_SampleNames, Ella_Broth_SampleNames, Lai2021_SputumNames, Lai2021_ExpNames))


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

###########################################################
##################### CLASSICAL MDS #######################

# Cmpute MDS
mds <- cmdscale(dist(merged_tpm_t))
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
       subtitle = "Using TPM",
       x = "Dimension 1",
       y = "Dimension 2") +
  my_plot_themes
# my_plot_themes_thumbnail
fig_MDS
ggsave(fig_MDS,
       file = "MDS_EllavsLai2021_TPM_v1.pdf",
       # file = "PCA_EllavsLai2021_BatchCorrect.LOG2.CPM_v1.pdf",
       path = "MDS_Figures",
       width = 8, height = 5, units = "in")



