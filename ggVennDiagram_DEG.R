# Making a VennDiagram of the DEG genes
# E. Lamont
# 4/25/25

# Pulled the DEG of sputum compared to broth for all the literature ("Import_literature_data.R") and want to do separate venn diagrams of the UP and DOWN significantly regulated genes. Needs to be separate in case they are super different
# Note for DEG, the genes that are always high in all samples will not show up, this method will show differences between samples more than similarities(?)

# https://gaospecial.github.io/ggVennDiagram/articles/using-ggVennDiagram.html


source("Import_data.R")
source("Import_literature_data.R")
# source("Most_Highly_Expressed_Genes.R")

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
#################### PROCESS THE DATA #####################

# Getting a list of all the DEG gene names for all the data
genelist_DEG_UP <- list("Ella_W0\nn=472" = EllaW0_DEG_UP %>% pull(GENE_ID),
                        "Honeyborne2016\nn=528" = Honeyborne_DEG_UP %>% pull(Gene),
                        "Garcia2016\nn=302" = Garcia2016_DEG_UP %>% pull(Gene),
                        "Sharma2017\nn=333" = Sharma2017_DEG_UP %>% pull(Name),
                        "Lai2021\nn=114" = Lai2021_DEG_UP %>% pull(Gene)
                        )

genelist_DEG_DOWN <- list("Ella_W0\nn=339" = EllaW0_DEG_DOWN %>% pull(GENE_ID),
                          "Honeyborne2016\nn=555" = Honeyborne_DEG_DOWN %>% pull(Gene),
                          "Garcia2016\nn=511" = Garcia2016_DEG_DOWN %>% pull(Gene),
                          "Sharma2017\nn=469" = Sharma2017_DEG_DOWN %>% pull(Name),
                          "Lai2021\nn=77" = Lai2021_DEG_DOWN %>% pull(Gene)
                          )


###########################################################
########### VENN DIAGRAM - UP REGULATED GENES #############

Venn_DEG_UP <- ggVennDiagram(genelist_DEG_UP,
                                   show_intersect = F,
                                   label = "both",
                                   set_size = 4,
                                   label_size = 3,
                                   label_alpha = 0) + # Text background
  scale_fill_distiller(palette = "RdBu", limits = c(0, 400),) +
  # scale_fill_distiller(palette = "Reds", direction = 1) +
  scale_x_continuous(expand = expansion(mult = .1)) +
  labs(title = "Upregulated genes in literature sputum compared to broth",
       subtitle = "DEGs compared to broth from literature, P<0.05, Log2fold>1 (Garcia stil as ratio)",
       fill = "Number of genes") 
Venn_DEG_UP
ggsave(Venn_DEG_UP,
       file = "W0vsLitSputum_DEG_UP_v1.pdf",
       path = "VennDiagrams_Figures",
       width = 9, height = 6, units = "in")

# Make it plotly! 
ggVennDiagram(genelist_DEG_UP,
              show_intersect = T)


###########################################################
########### VENN DIAGRAM - UP REGULATED GENES #############

Venn_DEG_DOWN <- ggVennDiagram(genelist_DEG_DOWN,
                             show_intersect = F,
                             label = "both",
                             set_size = 4,
                             label_size = 3,
                             label_alpha = 0) + # Text background
  scale_fill_distiller(palette = "RdBu", limits = c(0, 340),) +
  # scale_fill_distiller(palette = "Reds", direction = 1) +
  scale_x_continuous(expand = expansion(mult = .1)) +
  labs(title = "Downregulated genes in literature sputum compared to broth",
       subtitle = "DEGs compared to broth from literature, P<0.05, Log2fold<-1 (Garcia stil as ratio)",
       fill = "Number of genes") 
Venn_DEG_DOWN
ggsave(Venn_DEG_DOWN,
       file = "W0vsLitSputum_DEG_DOWN_v1.pdf",
       path = "VennDiagrams_Figures",
       width = 9, height = 6, units = "in")

# Make it plotly! 
ggVennDiagram(genelist_DEG_DOWN,
              show_intersect = T)
