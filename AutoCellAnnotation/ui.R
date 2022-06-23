library(shiny)
library(Seurat)
library(shinyjs)

shinyUI(fluidPage(
    useShinyjs(),
    conditionalPanel(
        "false", # always hide the download button
        downloadButton("dl", label = "Download Cell Type Coordinates"),
        downloadButton("dl2", label = "Download Seurat Object"),
    ),
    br(),
    titlePanel("scRNA-Seq: Auto-Cell Annotation"),
    tabsetPanel(
        tabPanel("Upload File",
                 h4("Upload Post-Clustering Seurat Object in RDS Format"),
                 sidebarLayout(
                     sidebarPanel(
                         fileInput('file1', 'Upload RDS File',
                                   accept=c('application/rds', 'application/RDS')),
                         downloadButton("downloadRDS", label = "Download Example"),
                         tags$a(href="http://dx.doi.org/10.1038/nbt.4042", "(Data Source)", target="_blank"),
                         h6("Alternatively, right click on the button and save link as."),
                         actionButton("exampledata", "View an Example"),
                         br(),
                         h5("Default Annotation Reference: Human Primary Cell Atlas"),
                         downloadButton("downloadHPCA", label = "Download Default Annotation Reference Provided by celldex"),
                         tags$a(href="http://dx.doi.org/10.1038/s41590-018-0276-y", "(Reference Index)", target="_blank"),
                         h5("Alternatively, please provide an alternative reference file similar to the default reference file format:"),
                         fileInput('altref', 'Upload Alternative Reference File (.RDS)',
                                   accept=c('application/rds', 'application/RDS')),
                         tags$br(),
                     ),
                     mainPanel(
                         tableOutput('contents')
                     )
                 ),
        ),
        tabPanel("Plot",
                 pageWithSidebar(
                     headerPanel('Auto-Cell Annotation'),
                     sidebarPanel(
                         selectInput('dimreduction', 'Dimension Reduction', "", selected = ""),
                         selectInput('group', 'Separate by Group', "", selected = ""),
                         actionButton("go", "Run DE Analysis for Cell Types"),
                         checkboxInput('label', 'Label Cell Types', TRUE),
                         # selectInput('annotate', 'Label Type', choices = annotate, selected = annotate[2]),
                         sliderInput("pointsize", 
                                     label = "Size of Points",
                                     min = 0, max = 10, value = 0.5, step = 0.1),
                         sliderInput("labelsize", 
                                     label = "Size of Labels",
                                     min = 0, max = 10, value = 4, step = 0.1),
                         sliderInput("colindex",
                                     label = "Choose A Color Scheme Index",
                                     min = 1, max = 100, value = 1, step = 1),
                         sliderInput("ncols",
                                     label = "Number of Columns (for multiple groups)",
                                     min = 1, max = 10, value = 1, step = 1),
                         sliderInput("width", 
                                     label = "Plot Width",
                                     min = 400, max = 8000, value = 780, step = 10),
                         sliderInput("height", 
                                     label = "Plot Height",
                                     min = 400, max = 8000, value = 480, step = 10)),
                     mainPanel(
                         plotOutput('ScatteredPlot', width = "100%")
                     )
                 )
        )
        
    )
)
)
