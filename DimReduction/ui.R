library(shiny)
library(Seurat)
library(shinyjs)

colorschemes <- c("Default","Colorful","Tableau1","Tableau2","Tenx","Mark","Bright","Cold","Warm")

shinyUI(fluidPage(
    useShinyjs(),
    conditionalPanel(
        "false",
        downloadButton("dl", label = "Download Dimension Reduction Coordinates"),
        downloadButton("dl2", label = "Download Seurat Object"),
    ),
    br(),
    titlePanel("scRNA-Seq: Dimension Reduction"),
    tabsetPanel(
        tabPanel("Upload File",
                 titlePanel("Upload"),
                 h5("Please upload the post-filtering (for one-sample project) or post-integration (for multiple-samples project) file in .RDS format"),
                 sidebarLayout(
                     sidebarPanel(
                         fileInput('file1', 'Upload RDS File',
                                   accept=c('.rds', 
                                            '.RDS', 
                                            'application/rds', 'application/RDS')),
                         downloadButton("downloadRDS", label = "Example Data"),
                         tags$a(href="http://dx.doi.org/10.1038/nbt.4042", "(Data Source)", target="_blank"),
                         h6("Alternatively, right click on the button and save link as."),
                         checkboxInput('runPCA', 'Run PCA', TRUE),
                         checkboxInput('runUMAP', 'Run UMAP', FALSE),
                         checkboxInput('runTSNE', 'Run tSNE (longer run time)', FALSE),
                         actionButton("run", "Start Dimension Reduction"),
                     ),
                     mainPanel(
                         tableOutput('contents')
                     )
                 ),
        ),
        tabPanel("Plot",
                 pageWithSidebar(
                     headerPanel('Plot'),
                     sidebarPanel(
                       selectInput('dimreduction', 'Dimension Reduction', "", selected = ""),
                       sliderInput("colorindex", 
                                   label = "Color Scheme Index",
                                   min = 1, max = 100, value = 1, step = 1),
                       sliderInput("pointsize", 
                                   label = "Point Size",
                                   min = 0.1, max = 10, value = 0.5, step = 0.1),
                       sliderInput("width", 
                                     label = "Plot Width",
                                     min = 400, max = 8000, value = 840, step = 10),
                         sliderInput("height", 
                                     label = "Plot Height",
                                     min = 400, max = 8000, value = 460, step = 10),
                       actionButton("go", "Download post-dimensional reduction files")),
                     mainPanel(
                         plotOutput('ScatteredPlot', width = "100%")
                     )
                 )
        )
        
    )
)
)
