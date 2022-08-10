library(shiny)
library(Seurat)
library(shinyjs)
library(shinyscreenshot)
library(DT)
library(plotly)

shinyUI(fluidPage(
    useShinyjs(),
    br(),
    titlePanel("scRNA-Seq: Manual Annotation"),
    tabsetPanel(
        tabPanel("Upload File",
                 h4("Upload Post-Auto-Cell Annotation or Auto-Clustering Seurat Object in RDS Format"),
                 sidebarLayout(
                     sidebarPanel(
                         fileInput('file1', 'Upload RDS File',
                                   accept=c('application/rds', 'application/RDS')),
                         downloadButton("downloadRDS", label = "Download Example"),
                         tags$a(href="http://dx.doi.org/10.1038/nbt.4042", "(Data Source)", target="_blank"),
                         h6("Alternatively, right click on the button and save link as."),
                         actionButton("exampledata", "View an Example"),
                         br()
                     ),
                     mainPanel(
                         tableOutput('contents')
                     )
                 ),
        ),
        tabPanel("Plot",
                 pageWithSidebar(
                     headerPanel('Plotting'),
                     sidebarPanel(
                       selectInput('dimreduction', 'Dimension Reduction', "", selected = ""),
                       selectInput('clusters', label = 'Choose an Annotation Group to View', choices = "", selected = ""),
                       sliderInput('res', label = 'Resolution Level for Auto-Clustering', min = 0, max = 10, value = 0.8, step = 0.05),
                       actionButton("go", "Run DE Analysis for the Chosen Annotation Group"),
                       br(),
                       br(),
                       selectInput('group', 'Separate by Group', choices = "", selected = ""),
                       sliderInput('pointsize', label = 'Size of Points', min = 0.1, max = 10, value = 0.5, step = 0.1),
                       checkboxInput('label', label = 'Label Cluster/Cell Type (For Static Plot)', value = TRUE),
                       sliderInput('labelsize', label = 'Size of Labels', min = 0.1, max = 20, value = 10, step = 0.1),
                       sliderInput('axislabelsize', label = 'Size of Axis Labels', min = 10, max = 50, value = 20, step = 0.1),
                       sliderInput("colindex", label = "Choose A Color Scheme Index", min = 1, max = 100, value = 1, step = 1),
                       sliderInput("ncols", label = "Number of Columns (for multiple groups)", min = 1, max = 10, value = 2, step = 1),
                       sliderInput('padj', label = 'Keep Genes with P-Value Smaller Than', min = 0, max = 1, value = 0.05, step = 0.01),
                       sliderInput('log2fc', label = 'Keep Genes with Log2FC Greater Than', min = -20, max = 20, value = 0.25, step = 0.1),
                       sliderInput('width', label = 'Plot Width', min = 400, max = 8000, value = 900, step = 10),
                       sliderInput('height', label = 'Plot Height', min = 200, max = 8000, value = 550, step = 10),
                         downloadButton("dlpng", "Download PNG"),
                         downloadButton('dlpdf', 'Download PDF'),
                         br(),
                         br(),
                       downloadButton('dl1', 'Download Current Seurat Object'),
                       h5("(If no manual annotation is provided, post-auto-clustering Seurat object will be downloaded)"),
                       downloadButton('dl2', 'Download Chosen Annotation Group DE Table*'),
                       downloadButton('dl3', 'Download Manual-Annotated DE Table*'),
                       h5("*Full Table and No Filter")),
                     mainPanel(
                       # h4("Please carry out manual annotation via the 'Lasso Select' icon or via clusters:"),
                       plotlyOutput('ScatteredPlotly', width = "100%", height = "100%", inline = TRUE),
                       tabsetPanel(
                         tabPanel("DE Result", DTOutput("table")),
                         tabPanel("Manual Annotation", DTOutput("form"),                     actionButton("updateManual", "Update Manual Annotation to Plot"),
       actionButton("runManualDE", "Run Current Annotation DE Analysis")),
       tabPanel("DE Result (Manual Annotation)", DTOutput("manualtable")),
       tabPanel("Static Plot", plotOutput("ScatteredPlot", width = "100%"))
                       )
                     )
                 )
        )
        
    )
)
)
