library(shiny)
library(Seurat)
library(shinyjs)

shinyUI(fluidPage(
    useShinyjs(),
    conditionalPanel(
        "false", # always hide the download button
        downloadButton("dl", label = "Download Clustering Coordinates"),
        downloadButton("dl2", label = "Download Seurat Object"),
        # downloadButton("SeuratRDS", label = "Download Seurat Output"),
        # downloadButton("DimCoord", label = "Download Coordinates")
    ),
    br(),
    titlePanel("scRNA-Seq: Auto-Clustering and DE Analysis"),
    tabsetPanel(
        tabPanel("Upload File",
                 h4("Upload Post-Dimensional Reduction Seurat Object in RDS Format"),
                 # h5("Acceptable Formats: Post-dimension reduction Seurat object file saved in RDS format"),
                 # h5("CSV/TXT file: Gene-Cell (Row-Column) Matrix CSV or TXT File (Un-Normalized Counts)"),
                 
                 sidebarLayout(
                     sidebarPanel(
                         fileInput('file1', 'Upload RDS File',
                                   accept=c('application/rds', 'application/RDS')),
                         
                         # added interface for uploading data from
                         # http://shiny.rstudio.com/gallery/file-upload.html
                         # tags$br(),
                         # checkboxInput('header', 'Header (Option for CSV/TXT file)', TRUE),
                         # checkboxInput('rownames', 'Row Names', TRUE),
                         # radioButtons('sep', 'Separator (Option for CSV/TXT file)',
                         #              c(Comma=',',
                         #                Semicolon=';',
                         #                Tab='\t'),
                         #              ','),
                         # radioButtons('quote', 'Quote (Option for CSV/TXT file)',
                         #              c(None='',
                         #                'Double Quote'='"',
                         #                'Single Quote'="'"),
                         #              '"'),
                         downloadButton("downloadRDS", label = "Download Example"),
                         tags$a(href="http://dx.doi.org/10.1038/nbt.4042", "(Data Source)", target="_blank"),
                         # downloadButton("downloadH5", label = "10X H5 Example"),
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
                         # "Empty inputs" - they will be updated after the data is uploaded
                         selectInput('dimreduction', 'Dimension Reduction', "", selected = ""),
                         # selectInput('ycol', 'Dimension 2', "", selected = ""),
                         # selectInput('group', 'Group', "", selected = ""),
                         sliderInput("res", 
                                     label = "Resolution Level for Cluster Numbers",
                                     min = 0, max = 10, value = 0.8, step = 0.05),
                         # sliderInput("nfeatureshigher", 
                         #             label = "Filter Genes Detected/Cell: Upper Bound",
                         #             min = 0, max = 25000, value = 2500, step = 10),
                         # sliderInput("mito", 
                         #             label = "Remove Mitochondria/Cell by Percent",
                         #             min = 0, max = 100, value = 5, step = 1),
                         # checkboxInput('runClusters', 'Run Auto-Clustering', TRUE),
                         # checkboxInput('runDE', 'DE Analysis', FALSE),
                         # checkboxInput('runTSNE', 'Run tSNE (longer run time)', FALSE),
                         actionButton("go", "Run DE Analysis"),
                         checkboxInput('label', 'Label Clusters', TRUE),
                         # selectInput('annotate', 'Label Type', choices = annotate, selected = annotate[2]),
                         sliderInput("pointsize", 
                                     label = "Size of Points",
                                     min = 0, max = 10, value = 0.5, step = 0.1),
                         sliderInput("labelsize", 
                                     label = "Size of Labels",
                                     min = 0, max = 10, value = 5, step = 0.1),
                         sliderInput("colindex",
                                     label = "Choose A Color Scheme Index",
                                     min = 1, max = 100, value = 1, step = 1),
                         sliderInput("width", 
                                     label = "Plot Width",
                                     min = 400, max = 8000, value = 780, step = 10),
                         sliderInput("height", 
                                     label = "Plot Height",
                                     min = 400, max = 8000, value = 480, step = 10)),
                         # selectInput('annotate', 'Label Type', choices = annotate, selected = annotate[2]),
                         # selectInput("color",
                         #             label = "Choose A Color Scheme",
                         #             choices = colorschemes,
                         #             selected = colorschemes[5]),
                         # checkboxInput('iflog', label = 'Log Scaling Violin Plot', value = TRUE),
                         # colourpicker::colourInput("col", "Choose Node Color", "#B995C2"),
                         # actionButton("DimReduction", "Screenshot"),
                     mainPanel(
                         plotOutput('ScatteredPlot', width = "100%")
                     )
                 )
        )
        
    )
)
)
