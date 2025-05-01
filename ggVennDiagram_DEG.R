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

# These are all with p<0.05 and log2fold (or just fold, for some, it's unclear) of at least 2.5 

cat(paste0("Ella_W0\nn=", nrow(EllaW0_DEG_UP)))
paste0("Ella_W0\nn=", nrow(EllaW0_DEG_UP))

# Getting a list of all the DEG gene names for all the data
genelist_DEG_UP <- list("Ella_W0" = EllaW0_DEG_UP %>% pull(GENE_ID),
                        "Honeyborne2016" = Honeyborne_DEG_UP %>% pull(Gene),
                        "Garcia2016" = Garcia2016_DEG_UP %>% pull(Gene),
                        "Sharma2017" = Sharma2017_DEG_UP %>% pull(Name),
                        "Garton2008" = Garton2008_DEG_UP %>% pull(Gene),
                        "Lai2021" = Lai2021_DEG_UP %>% pull(Gene)
                        )
# Set the names to include the number of genes
names(genelist_DEG_UP) <- c(
  paste0("Ella_W0\nn=", nrow(EllaW0_DEG_UP)),
  paste0("Honeyborne2016\nn=", nrow(Honeyborne_DEG_UP)),
  paste0("Garcia2016\nn=", nrow(Garcia2016_DEG_UP)),
  paste0("Sharma2017\nn=", nrow(Sharma2017_DEG_UP)),
  paste0("Garton2008\nn=", nrow(Garton2008_DEG_UP)),
  paste0("Lai2021\nn=", nrow(Lai2021_DEG_UP))
)

genelist_DEG_DOWN <- list("Ella_W0" = EllaW0_DEG_DOWN %>% pull(GENE_ID),
                          "Honeyborne2016" = Honeyborne_DEG_DOWN %>% pull(Gene),
                          "Garcia2016" = Garcia2016_DEG_DOWN %>% pull(Gene),
                          "Sharma2017" = Sharma2017_DEG_DOWN %>% pull(Name),
                          "Garton2008" = Garton2008_DEG_DOWN %>% pull(Gene),
                          "Lai2021" = Lai2021_DEG_DOWN %>% pull(Gene)
                          )
# Set the names to include the number of genes
names(genelist_DEG_DOWN) <- c(
  paste0("Ella_W0\nn=", nrow(EllaW0_DEG_DOWN)),
  paste0("Honeyborne2016\nn=", nrow(Honeyborne_DEG_DOWN)),
  paste0("Garcia2016\nn=", nrow(Garcia2016_DEG_DOWN)),
  paste0("Sharma2017\nn=", nrow(Sharma2017_DEG_DOWN)),
  paste0("Garton2008\nn=", nrow(Garton2008_DEG_DOWN)),
  paste0("Lai2021\nn=", nrow(Lai2021_DEG_DOWN))
)

###########################################################
########### VENN DIAGRAM - UP REGULATED GENES #############

Venn_DEG_UP <- ggVennDiagram(genelist_DEG_UP,
                                   show_intersect = F,
                                   label = "both",
                                   set_size = 4,
                                   label_size = 3,
                                   label_alpha = 0) + # Text background
  scale_fill_distiller(palette = "RdBu", limits = c(0, 320),) +
  # scale_fill_distiller(palette = "Reds", direction = 1) +
  scale_x_continuous(expand = expansion(mult = .1)) +
  labs(title = "Upregulated genes in literature sputum compared to broth",
       subtitle = "DEGs compared to broth from literature, P<0.05, Log2fold>2.5",
       fill = "Number of genes") 
Venn_DEG_UP
ggsave(Venn_DEG_UP,
       file = "W0vsLitSputum_DEG_UP_v2.pdf",
       path = "VennDiagrams_Figures",
       width = 9, height = 6, units = "in")

# Make it plotly! 
# ggVennDiagram(genelist_DEG_UP,
#               show_intersect = T)


###########################################################
########## VENN DIAGRAM - DOWN REGULATED GENES ############

Venn_DEG_DOWN <- ggVennDiagram(genelist_DEG_DOWN,
                             show_intersect = F,
                             label = "both",
                             set_size = 4,
                             label_size = 3,
                             label_alpha = 0) + # Text background
  scale_fill_distiller(palette = "RdBu", limits = c(0, 210),) +
  # scale_fill_distiller(palette = "Reds", direction = 1) +
  scale_x_continuous(expand = expansion(mult = .1)) +
  labs(title = "Downregulated genes in literature sputum compared to broth",
       subtitle = "DEGs compared to broth from literature, P<0.05, Log2fold<-2.5",
       fill = "Number of genes") 
Venn_DEG_DOWN
ggsave(Venn_DEG_DOWN,
       file = "W0vsLitSputum_DEG_DOWN_v2.pdf",
       path = "VennDiagrams_Figures",
       width = 9, height = 6, units = "in")

# Make it plotly! 
# ggVennDiagram(genelist_DEG_DOWN,
#               show_intersect = T)


###########################################################
########## INDIVIDUAL VENN DIAGRAM COMPARISONS ############

Venn_DEG_UP <- ggVennDiagram(genelist_DEG_UP[c(4,5)],
                             show_intersect = F,
                             label = "both",
                             set_size = 4,
                             label_size = 3,
                             label_alpha = 0) + # Text background
  scale_fill_distiller(palette = "RdBu", limits = c(0, 400),) +
  # scale_fill_distiller(palette = "Reds", direction = 1) +
  scale_x_continuous(expand = expansion(mult = .1)) +
  labs(title = "Upregulated genes in literature sputum compared to broth",
       subtitle = "DEGs compared to broth from literature, P<0.05, Log2fold>2.5",
       fill = "Number of genes") 
Venn_DEG_UP

