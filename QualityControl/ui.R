library(shiny)
library(Seurat)
library(colourpicker)
library(shinyjs)

shinyUI(fluidPage(
    useShinyjs(),
    conditionalPanel(
        "false",
        downloadButton("dl", label = "Download Seurat Object"),
    ),
    br(),
    titlePanel("scRNA-Seq: Quality Control"),
    tabsetPanel(
        tabPanel("Upload File",
                 titlePanel("Upload 10X H5/CSV/TXT File"),
                 h5("Acceptable Formats:"),
                 h5("H5 file: 10X Filtered Matrix H5 File"),
                 h5("CSV/TXT file: Gene-Cell (Row-Column) Matrix CSV or TXT File (Un-Normalized Counts)"),
                 h5("RDS file: Seurat object"),
                 sidebarLayout(
                     sidebarPanel(
                         fileInput('file1', '1.Upload H5/CSV/TXT/RDS Files',
                                   multiple = TRUE,
                                   accept=c('text/csv', '.h5',
                                            'text/comma-separated-values,text/plain', 'application/h5','application/rds','application/RDS')),
                         checkboxInput('header', 'Header (Option for CSV/TXT file)', TRUE),
                         checkboxInput('rownames', 'Row Names (Option for CSV/TXT file)', TRUE),
                         radioButtons('sep', 'Separator (Option for CSV/TXT file)',
                                      c(Comma=',',
                                        Semicolon=';',
                                        Tab='\t'),
                                      ','),
                         radioButtons('quote', 'Quote (Option for CSV/TXT file)',
                                      c(None='',
                                        'Double Quote'='"',
                                        'Single Quote'="'"),
                                      '"'),
                         tags$a(href="http://dx.doi.org/10.1038/nbt.4042", "Downsampling the data from this source as an example:", target="_blank"),
                         downloadButton("downloadCSV", label = "Example Data (CSV)"),
                         downloadButton("downloadH5", label = "Example Data (H5)"),
                         h6("Alternatively, right click on the button and save link as."),
          br(),
                         fileInput('metafile', '2.Upload a CSV/TXT meta file (Optional. Only for multiple-samples project)',
                                   accept=c('text/csv', 
                                            'text/comma-separated-values,text/plain', 
                                            '.csv')),
          checkboxInput('nometa', 'Click here if no metafile is to be provided', value = FALSE),
          radioButtons('metasep', 'Separator (Option for meta file)',
                       c(Comma=',',
                         Semicolon=';',
                         Tab='\t'),
                       ','),
                         downloadButton("downloadMeta", label = "Example Meta File"),
                     ),
                     mainPanel(
                         tableOutput('contents')
                     )
                 ),
        ),
        tabPanel("Plot",
                 pageWithSidebar(
                     headerPanel('Calculations'),
                     sidebarPanel(
                         sliderInput("nfeatureslower", 
                                     label = "Filter Genes Detected/Cell: Lower Bound",
                                     min = 0, max = 1000, value = 200, step = 10),
                         sliderInput("nfeatureshigher", 
                                     label = "Filter Genes Detected/Cell: Upper Bound",
                                     min = 0, max = 25000, value = 2500, step = 10),
                         sliderInput("mito", 
                                     label = "Remove Mitochondria/Cell by Percent",min = 0, max = 100, value = 5, step = 1),
                         actionButton("go", "Create filtered dataset"),
                         colourpicker::colourInput("col", "Choose Plot Color Theme", "#B995C2"),
                         sliderInput("width", 
                                     label = "Plot Width",
                                     min = 400, max = 8000, value = 780, step = 10),
                         sliderInput("height", 
                                     label = "Plot Height",
                                     min = 400, max = 8000, value = 620, step = 10)),
                     mainPanel(
                         plotOutput('VlnPlot', width = "100%")
                     )
                 )
        )
        
    )
)
)
