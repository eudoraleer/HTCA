library(shiny)
library(Seurat)
library(shinyjs)
library(shinyscreenshot)
library(liana)
library(DT)

shinyUI(fluidPage(
  useShinyjs(),
  conditionalPanel(
    "false",
    downloadButton("dl", label = "Download Clustering Coordinates"),
    downloadButton("dl2", label = "Download Seurat Object")
  ),
  br(),
  titlePanel("scRNA-Seq: Cell-Cell Communication Analysis"),
  tabsetPanel(
    tabPanel("Upload File",
             h4("Upload Post-Annotation or Post-Clustering Seurat Object in RDS Format"),
             sidebarLayout(
               sidebarPanel(
                 fileInput('file1', 'Upload RDS File',
                           accept=c('application/rds', 'application/RDS')),
                 downloadButton("downloadRDS", label = "Download Example"),
                 tags$a(href="http://dx.doi.org/10.1038/nbt.4042", "(Data Source)", target="_blank"),
                 h6("Alternatively, right click on the button and save link as."),
                 actionButton("exampledata", "View an Example"),
                 br(),
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
                 checkboxGroupInput('methods', 'Select method(s)', selected = ""),
                 checkboxGroupInput('dbs', 'Select knowledge resource(s)', selected = ""),
                 selectInput('celllabels', label = 'Select the cell labels to perform the analysis', choices = "", selected = ""),
                 selectInput('groups', label = 'Select a group to perform the analysis separately (only for data with multiple groups, for example, diseased and control groups. Ignored for single-group data)', choices = "", selected = NULL),
                 actionButton("go", "Run LIANA Cell-Cell Communication Analysis"),
                 br(),
                 tags$a(href="https://doi.org/10.1038/s41467-022-30755-0", "Reference to LIANA", target="_blank"),
                 br(),
                 br(),
                 checkboxInput('set3d', label = '3D Vertex', value = TRUE),
                 sliderInput('labelsize', label = 'Size of Labels', min = 0.1, max = 20, value = 2, step = 0.1),
                 sliderInput('labeldist', label = 'Distance of Labels', min = 0.1, max = 10, value = 3, step = 0.1),
                 sliderInput('arrowthickness', label = 'Thickness of Arrows', min = 0.1, max = 10, value = 0.5, step = 0.1),
                 sliderInput('vertexsize', label = 'Size of Vertex', min = 10, max = 50, value = 14, step = 1),
                 sliderInput('curvity', label = 'Curvity of Arrows', min = 0.1, max = 2, value = 0.8, step = 0.01),
                 sliderInput("colindex",
                             label = "Choose A Color Scheme Index",
                             min = 1, max = 100, value = 1, step = 1),
                 sliderInput("width", 
                             label = "Plot Width",
                             min = 400, max = 8000, value = 900, step = 10),
                 sliderInput("height", 
                             label = "Plot Height",
                             min = 400, max = 8000, value = 760, step = 10),
                 downloadButton("dlpng", "Download PNG"),
                 downloadButton('dlpdf', 'Download PDF'),
                 br(),
                 br(),
                 downloadButton("dlresult", "Download Result Tables")),
               mainPanel(
                 selectInput('viewmethod', label = 'Method(s)', choices = "", selected = ""),
                 selectInput('viewdb', label = 'Resource(s)', choices = "", selected = ""),
                 plotOutput('Graph', width = "100%", inline = TRUE),
                 fluidRow(column(12, DTOutput('table')))
               )
             )
    )
    
  )
)
)
