library(shiny)
library(Seurat)
library(shinyjs)
library(shinyscreenshot)

shinyUI(fluidPage(
    useShinyjs(),
    conditionalPanel(
        "false",
        downloadButton("dl", label = "Download Clustering Coordinates"),
        downloadButton("dl2", label = "Download Seurat Object"),
    ),
    br(),
    titlePanel("Imputation Procedure for scRNA-Seq Data"),
    tabsetPanel(
        tabPanel("Upload File",
                 h4("Upload Post-Quality Control Seurat Object in RDS Format"),
                 tags$a(href="https://www.nature.com/articles/s41467-021-27729-z", "(Imputation method: ALRA)", target="_blank"),
                 br(),
                 br(),
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
                         selectInput('sample', label = 'Choose a Sample to View', choices = "", selected = ""),
                         sliderInput('pointsize', label = 'Size of Points', min = 0.1, max = 10, value = 1.5, step = 0.1),
                         sliderInput('labelsize', label = 'Size of Labels', min = 10, max = 50, value = 20, step = 0.1),
                         sliderInput('width', label = 'Plot Width', min = 400, max = 8000, value = 900, step = 10),
                         sliderInput('height', label = 'Plot Height', min = 200, max = 8000, value = 400, step = 10),
                         downloadButton("dlpng", "Download PNG"),
                         downloadButton('dlpdf', 'Download PDF'),
                         br(),
                         br(),
                         downloadButton('dlRDS', 'Download Post-Imputation Seurat Object')),
                     mainPanel(
                         plotOutput('ScatteredPlot', width = "100%")
                     )
                 )
        )
        
    )
)
)
