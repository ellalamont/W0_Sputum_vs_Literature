# Making a VennDiagram
# E. Lamont
# 4/18/25

# https://gaospecial.github.io/ggVennDiagram/articles/using-ggVennDiagram.html


source("Import_data.R")
source("Import_literature_data.R")
source("Most_Highly_Expressed_Genes.R") # top15_AllSputum_RankExpression_W0RawReads



# Plot basics
my_plot_themes <- theme_grey() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=10),
        legend.title = element_text(size = 12),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=10), 
        axis.text.x = element_text(angle = 0, size=10, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10), 
        plot.subtitle = element_text(size=9), 
        plot.margin = margin(10, 10, 10, 20),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_blank()
  )



###########################################################
######### MAKE THE VENN DIAGRAM - FROM READS_M ############

# These are the ranked genes pulled directly from Coppola2021, then my W0 averages that I tried to rank the READS_M on based on how Coppola did it.

gene_list_top15 <- list("Ella_W0\nn=675" = RANK_W0_RawReads_top15_GeneNames,
                  "Walter2015\nn=361" = Coppola2021_Walter2015_top15_GeneNames, 
                  "Garcia2016\nn=296" = Coppola2021_Garcia2016_top15_GeneNames, 
                  "Sharma2017\nn=589" = Coppola2021_Sharma2017_top15_GeneNames,
                  "Lai2021\nn=617" = Coppola2021_Lai2021_top15_GeneNames)

Venn_FromRawReads <- ggVennDiagram(gene_list_top15,
                    show_intersect = F,
                    label = "both",
                    set_size = 4,
                    label_size = 3,
                    label_alpha = 0) + # Text background
  scale_fill_distiller(palette = "RdBu", limits = c(0, 420),) +
  # scale_fill_distiller(palette = "Reds", direction = 1) +
  scale_x_continuous(expand = expansion(mult = .1)) +
  labs(title = "Top 15% expressed genes in literature sputum",
       subtitle = "W0 data started with READS_M, literature data from Coppola2021",
       fill = "Number of genes") 
Venn_FromRawReads

# Make it plotly! 
ggVennDiagram(gene_list2,
              show_intersect = T)

# ggsave(Venn_FromRawReads,
#        file = "W0vsLitSputum_Top15PercentGenes_wREADS_M_v1.pdf",
#        path = "VennDiagrams_Figures",
#        width = 9, height = 6, units = "in")


###########################################################
######### VENN - MY RAW READS WITH LAI2021 ONLY ###########

gene_list3 <- list("Ella_W0\nn=674" = RANK_W0_RawReads_top15_GeneNames,
                   "Lai2021\nn=617" = Coppola2021_Lai2021_top15_GeneNames)

Venn_EllaRawReads_Lai2021 <- ggVennDiagram(gene_list3,
                                   show_intersect = F,
                                   label = "both",
                                   set_size = 4,
                                   label_size = 3,
                                   label_alpha = 0) + # Text background
  scale_fill_distiller(palette = "RdBu", limits = c(0, 440),) +
  # scale_fill_distiller(palette = "Reds", direction = 1) +
  scale_x_continuous(expand = expansion(mult = .1)) +
  labs(title = "Top 15% expressed genes Ella W0 vs Lai2021",
       subtitle = "W0 data started with READS_M, literature data from Coppola2021",
       fill = "Number of genes") 
Venn_EllaRawReads_Lai2021

ggsave(Venn_EllaRawReads_Lai2021,
       file = "W0vsLai2021_Top15PercentGenes_wREADS_M_v1.pdf",
       path = "VennDiagrams_Figures",
       width = 9, height = 6, units = "in")



###########################################################
########## VENN - MY RAW READS W/OUT SHARMA2017 ###########

gene_list4 <- list("Ella_W0\nn=674" = RANK_W0_RawReads_top15_GeneNames,
                   "Walter2015\nn=361" = Coppola2021_Walter2015_top15_GeneNames, 
                   "Garcia2016\nn=296" = Coppola2021_Garcia2016_top15_GeneNames, 
                   # "Sharma2017\nn=589" = Coppola2021_Sharma2017_top15_GeneNames,
                   "Lai2021\nn=617" = Coppola2021_Lai2021_top15_GeneNames)

Venn_FromRawReads_NoSharma <- ggVennDiagram(gene_list_top15[],
                                   show_intersect = F,
                                   label = "both",
                                   set_size = 4,
                                   label_size = 3,
                                   label_alpha = 0) + # Text background
  scale_fill_distiller(palette = "RdBu", limits = c(0, 310),) +
  # scale_fill_distiller(palette = "Reds", direction = 1) +
  scale_x_continuous(expand = expansion(mult = .1)) +
  labs(title = "Top 15% expressed genes in literature sputum",
       subtitle = "W0 data started with READS_M, literature data from Coppola2021",
       fill = "Number of genes") 
Venn_FromRawReads_NoSharma



###########################################################
########## INDIVIDUAL VENN DIAGRAM COMPARISONS ############

Venns_plot <- ggVennDiagram(gene_list_top15[c(1,5)],
                             show_intersect = F,
                             label = "both",
                             set_size = 4,
                             label_size = 3,
                             label_alpha = 0) + # Text background
  scale_fill_distiller(palette = "RdBu", limits = c(0, 525),) +
  # scale_fill_distiller(palette = "Reds", direction = 1) +
  scale_x_continuous(expand = expansion(mult = .1)) +
  labs(title = "Upregulated genes in literature sputum compared to broth",
       subtitle = "DEGs compared to broth from literature, P<0.05, Log2fold>2.5",
       fill = "Number of genes") 
Venns_plot


