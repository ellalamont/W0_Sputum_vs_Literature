library(shiny)

# dev.off() # Sometimes need this to get the shiny to show the pheatmap?

source("Import_data.R")
source("Import_literature_data.R")

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none",legend.text=element_text(size=12),
        legend.title = element_text(size = 12),
        plot.title = element_text(size=12), 
        axis.title.x = element_text(size=12), 
        axis.text.x = element_text(angle = 0, size=12, vjust=1, hjust=0.5),
        axis.title.y = element_text(size=12),
        axis.text.y = element_text(size=12), 
        plot.subtitle = element_text(size=12), 
        plot.margin = margin(10, 10, 10, 20))

my_annotation_colors <- list(
  Sample_Type = c("Sputum" = "#0072B2",  # Blue
           "Broth" = "#999999",
           "Rv_Exponential" = "lightgrey",
           "TBSputum_nonHIV" = "#C6DDF0")
)

# Need to merge the tpm and metadata because It's not done in the import data files
merged_metadata <- bind_rows(
  select(my_pipeSummary, all_of(c("SampleID", "Sample_Type"))),
  select(Lai2021_metadata, all_of(c("SampleID", "Sample_Type")))
)
rownames(merged_metadata) <- merged_metadata$SampleID

merged_tpm <- full_join(
  select(All_tpm, all_of(c(W0_Sputum_SampleNames, Ella_Broth_SampleNames, "X"))),
  Lai2021_tpm,
  by = "X"
)
merged_tpm <- merged_tpm %>% column_to_rownames(var = "X")


# Define UI ----
ui <- fluidPage(
  titlePanel("Ella W0 and H37Ra broth vs Lai2021 W0 and Exponential broth"),
  
  fluidRow(
    
    column(width = 4,
           
           # Which samples to plot
           checkboxGroupInput("timepoints",
                              label = "Select Samples",
                              choices = c("Broth", "Sputum", "Rv_Exponential", "TBSputum_nonHIV"),
                              selected = c("Broth", "Sputum", "Rv_Exponential", "TBSputum_nonHIV")),
           
           # Dropdown for selecting which rda file (gene set source)
           selectInput("my_GeneSetSource",
                       label = "Gene Set Source",
                       choices = names(allGeneSetList)),
           # Dropdown for selecting the gene set within the chosen rda file.
           selectInput("my_GeneSet",
                       label = "Gene Set",
                       choices = NULL),
           
           # Select genes manually
           textInput("manual_genes",
                     label = "Enter Rv# (comma separated)",
                     placeholder = "Rv1473A, Rv2011c, Rv0494"),
           # Add row clustering options
           numericInput("cutree_rows", 
                        label = "Number of Row Clusters", 
                        value = 1, min = 1, step = 1),
           # Add column clustering options
           numericInput("cutree_cols", 
                        label = "Number of Column Clusters", 
                        value = 1, min = 1, step = 1),
           selectInput("my_scaling",
                       label = "How to scale",
                       choices = (c("row", "column", "none"))),
           # Add checkbox to toggle display_numbers in heatmap
           checkboxInput("show_numbers", label = "Show Values", value = FALSE),
           textInput("my_GeneID", 
                     label = "Pick a gene to link to mycobrowser",
                     placeholder = "Rv..."),
           uiOutput("gene_link")  # New UI output for the link
    ),
    
    column(width = 8, # Max is 12...
           uiOutput("dynamic_pheatmap")
           # plotOutput("pheatmap", width = "200%", height = "600")
           # plotOutput("pheatmap", width = "100%", height = "600px")
    )
    
  )
  
)

# Define server logic ----
server <- function(input, output, session) {
  
  # Gene Link
  output$gene_link <- renderUI({
    req(input$my_GeneID)  # Ensure there's a valid input
    url <- paste0("https://mycobrowser.epfl.ch/genes/", input$my_GeneID)
    tags$a(href = url, target = "_blank", paste0("View Details of ", input$my_GeneID, " on Mycobrowser"))
  })
  
  # When a new gene set source is selected, update the gene set dropdown
  observeEvent(input$my_GeneSetSource, {
    updateSelectInput(session, "my_GeneSet",
                      choices = names(allGeneSetList[[input$my_GeneSetSource]]),
                      selected = NULL)
  })
  
  # Function to get the selected genes
  get_selected_genes <- reactive({
    if (input$manual_genes != "") {
      # Process manual input: remove extra spaces, split by commas or spaces
      genes <- unlist(strsplit(input$manual_genes, "[,\\s]+"))
      genes <- trimws(genes)  # Trim whitespace
      genes <- genes[genes != ""]  # Remove empty entries
    } else if (!is.null(input$my_GeneSet) && input$my_GeneSet %in% names(allGeneSetList[[input$my_GeneSetSource]])) {
      genes <- allGeneSetList[[input$my_GeneSetSource]][[input$my_GeneSet]]
    } else {
      genes <- character(0)  # Empty vector if nothing is selected
    }
    
    return(genes)
  })
  
  # Dynamic UI for heatmap height
  output$dynamic_pheatmap <- renderUI({
    # req(input$my_GeneSetSource, input$my_GeneSet)
    
    # Get the list of genes from either dropdown selection or manual input
    selected_genes <- get_selected_genes()
    
    # Count the number of genes in the selected set
    num_genes <- sum(rownames(merged_tpm) %in% allGeneSetList[[input$my_GeneSetSource]][[input$my_GeneSet]])
    
    # Dynamically set plot height (base height + extra space per gene)
    plot_height <- max(400, min(2500, num_genes * 40))  # Adjust as needed
    # plot_width <- max(2000, min(2000)) # Still need to figure out how to change this! 
    
    plotOutput("pheatmap", height = paste0(plot_height, "px"))
  })
  
  
  # Render the pheatmap
  output$pheatmap <- renderPlot({
    
    # Make a list of the columns for each timepoint
    timepoint_columns <- list(
      "Broth" = c("H37Ra_Broth_4", "H37Ra_Broth_5", "H37Ra_Broth_6"),
      "Sputum" = c("S_250754", "S_355466", "S_503557"),
      "Rv_Exponential" = c("SRR10125314", "SRR10125315"),
      "TBSputum_nonHIV" = c("SRR10125319", "SRR10125320", "SRR10125321", "SRR10125322", "SRR10125323", "SRR10125324")
    )
    
    selected_genes <- get_selected_genes()
    
    # Filter data based on selected genes
    my_data <- merged_tpm[rownames(merged_tpm) %in% selected_genes, , drop = FALSE]
    
    # Now filter columns based on which timepoints are checked
    if (!is.null(input$timepoints)) {
      selected_cols <- unlist(timepoint_columns[input$timepoints])
      selected_cols <- intersect(selected_cols, colnames(my_data))  # only keep existing columns
      my_data <- my_data[, selected_cols, drop = FALSE]
    }
    
    # If no columns left, show error
    if (ncol(my_data) == 0) {
      showNotification("No columns match the selected timepoints.", type = "error")
      return(NULL)
    }
    
    # Check if we have at least 2 genes
    if (nrow(my_data) < 2) {
      showNotification("At least two valid genes are required for clustering.", type = "error")
      return(NULL)
    }
    
    
    p <- pheatmap(my_data, 
                  annotation_col = merged_metadata["Sample_Type"], 
                  annotation_colors = my_annotation_colors,
                  # col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                  scale = input$my_scaling, 
                  display_numbers = input$show_numbers,
                  fontsize_number = 8,
                  cutree_rows = input$cutree_rows,
                  cutree_cols = input$cutree_cols,
                  fontsize = 18)
    p
    
    # pheatmap returns a complex grid object; use grid.draw() to render it in Shiny.
    grid::grid.newpage()
    grid::grid.draw(p$gtable)
  })
  
}

# Run the app ----
shinyApp(ui = ui, server = server)# , options = list(launch.browser = TRUE))