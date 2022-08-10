library(shiny)
library(Seurat)
library(shinyjs)
library(shinyscreenshot)
library(DT)
library(plotly)
library(shinyalert)

shinyUI(fluidPage(
    useShinyjs(),
    br(),
    titlePanel("scRNA-Seq: Data Splicing"),
    tabsetPanel(
        tabPanel("Upload File",
                 h4("Upload Seurat Object in RDS Format"),
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
                       selectInput('group', 'Separate by Group', choices = "", selected = ""),
                       sliderInput('pointsize', label = 'Size of Points', min = 0.1, max = 10, value = 0.5, step = 0.1),
                       checkboxInput('label', label = 'Label Cluster/Cell Type (For Static Plot)', value = TRUE),
                       sliderInput('labelsize', label = 'Size of Labels', min = 0.1, max = 20, value = 10, step = 0.1),
                       sliderInput('axislabelsize', label = 'Size of Axis Labels', min = 10, max = 50, value = 20, step = 0.1),
                       sliderInput("colindex", label = "Choose A Color Scheme Index", min = 1, max = 100, value = 1, step = 1),
                       sliderInput("ncols", label = "Number of Columns (for multiple groups)", min = 1, max = 10, value = 2, step = 1),
                       sliderInput('width', label = 'Plot Width', min = 400, max = 8000, value = 900, step = 10),
                       sliderInput('height', label = 'Plot Height', min = 200, max = 8000, value = 550, step = 10),
                         downloadButton("dlpng", "Download PNG"),
                         downloadButton('dlpdf', 'Download PDF'),
                         br(),
                         br(),
                       downloadButton('dl1', 'Download Selected Cells (Seurat Object)')),
                     mainPanel(
                       h4("Data splicing can be done via the 'Lasso Select' icon at the top right menu bar of the plot"),
                       plotlyOutput('ScatteredPlotly', width = "100%", height = "100%", inline = TRUE),
                       verbatimTextOutput("selecting"),
                       verbatimTextOutput("selected")
                     )
                 )
        )
        
    )
)
)
