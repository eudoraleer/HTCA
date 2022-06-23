library(shiny)
library(Seurat)
library(shinyjs)

colorschemes <- c("Default","Colorful","Tableau1","Tableau2","Tenx","Mark","Bright","Cold","Warm")

shinyUI(fluidPage(
    useShinyjs(),
    conditionalPanel(
        "false",
        downloadButton("dl", label = "Download Post-integration Seurat Object"),
    ),
    br(),
    titlePanel("scRNA-Seq: Data Integration"),
    tabsetPanel(
        tabPanel("Upload File",
                 titlePanel("Upload Post-filtering Seurat .RDS file"),
                 sidebarLayout(
                     sidebarPanel(
                         fileInput('file1', 'Upload .RDS File',
                                   accept=c('.rds', 'application/RDS','.RDS')),
                         selectInput('integrationmethod', 'Select an Integration Method', choices = c("Seurat V3","Harmony"), selected = "Seurat V3"),
                         downloadButton("downloadRDS", label = "Example Data"),
                         tags$a(href="http://dx.doi.org/10.1038/nbt.4042", "(Data Source)", target="_blank"),
                         h6("Alternatively, right click on the button and save link as."),
                     ),
                     mainPanel(
                         tableOutput('contents')
                     )
                 ),
        ),
        tabPanel("Plot",
                 pageWithSidebar(
                     headerPanel('Visualization'),
                     sidebarPanel(
                       sliderInput("pointsize", 
                                   label = "Choose Point Size",
                                   min = 0.1, max = 5, value = 0.3, step = 0.1),
                       # sliderInput("colornum", 
                       #             label = "Choose Color Scheme Index",
                       #             min = 1, max = 100, value = 1, step = 1),
                         sliderInput("width", 
                                     label = "Plot Width",
                                     min = 400, max = 8000, value = 880, step = 10),
                         sliderInput("height", 
                                     label = "Plot Height",
                                     min = 400, max = 8000, value = 400, step = 10),
                         actionButton("go", "Download post-integration file")
                    ),
                     # actionButton("DimReduction", "Screenshot"),
                     mainPanel(
                         plotOutput('umap', width = "100%"),
                         plotOutput('tsne', width = "100%")
                     )
                 )
        )
        
    )
)
)
