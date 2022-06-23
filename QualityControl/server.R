library(BiocManager)
options(repos = BiocManager::repositories())

library("shiny")
library("ggplot2")
library("Seurat")
library("colourpicker")
library("shinyjs")
library("gt")
library("hdf5r")

options(shiny.maxRequestSize=8000*1024^2)
exampleCSV <- "DB/QC_Example_Dataset.csv"
exampleH5 <- "DB/QC_Example_Dataset.h5"
examplemeta <- "DB/QC_Metadata_Example.csv"
sys_time <- format(Sys.time(), "%Y%m%d%H%M%S")

shinyServer(function(input, output, session) {
    output$downloadCSV <- downloadHandler(
        filename = "QC_Example_Dataset.csv",
        content = function(file) {
            file.copy(exampleCSV, file)
        },
        contentType = "text/csv"
    )
    
    output$downloadH5 <- downloadHandler(
        filename = "HTCA_QC_Example_Dataset.h5",
        content <- function(file) {
            file.copy(exampleH5, file)
        },
        contentType = "application/h5"
    )
    
    output$downloadMeta <- downloadHandler(
      filename = "HTCA_QC_Metadata_Example.csv",
      content <- function(file) {
        file.copy(examplemeta, file)
      },
      contentType = "text/csv"
    )
    
    output$dl <- downloadHandler(
        filename <- paste("HTCA_Post_QC_Seurat_Object_",sys_time,".RDS", sep = ""),
        content <- function(file) {
            x <- qualityControl()$x
            saveRDS(x,file)
        },
        contentType = "application/RDS"
    )
    
    meta <- reactive({
      if(!is.null(input$metafile)){
        m <- read.table(input$metafile$datapath, sep = input$metasep, header = T)
        colnames(m) <- c("File_Name","Sample_ID","Batch","Group")
      }else{
        m <- NULL
      }
      return(m)
    })
    
    data <- reactive({
        req(input$file1)
        req((input$nometa == FALSE & !is.null(input$metafile)) |
            (input$nometa == TRUE & is.null(input$metafile)))
        showModal(modalDialog("Running filtering..please proceed to plots after filtering is done", footer=NULL))
        inFile <- input$file1
        current <- NULL
          for(i in 1:nrow(inFile)){
            cfile <- inFile$datapath[i]
            df <- NULL
        if(length(grep("\\.csv$", cfile, ignore.case = T)) > 0){
        df <- read.csv(cfile, header = input$header, sep = input$sep, quote = input$quote, row.names = if(input$rownames == TRUE) newobj <- 1 else newobj <- NULL)
        df <- CreateSeuratObject(counts = df, min.cells = 3, min.features = 200)
        }else if(length(grep("\\.h5$", cfile, ignore.case = T)) > 0){
        df <- Read10X_h5(cfile)
        if((length(df) == 1) | (length(df) > 10)){
            df <- CreateSeuratObject(counts = df, min.cells = 3, min.features = 200)
        }else{
            df <- CreateSeuratObject(counts = df$`Gene Expression`, min.cells = 3, min.features = 200)
        }
        }
        crow <- NULL
        if(!is.null(meta())){
              crow <- which(gsub("(.*)\\s+$","\\1",meta()$File_Name) == inFile$name[i])
        df$orig.ident <- meta()[crow,"Sample_ID"]
        df$Batch <- meta()[crow,"Batch"]
        df$Group <- meta()[crow,"Group"]
        }else{
          df$orig.ident <- "Sample"
        }
        df$Percent_Mito <- PercentageFeatureSet(df, pattern = "^MT-")
        current[[i]] <- df
          }
        if(length(current) > 1){
          current <- merge(current[[1]], current[c(2:length(current))])
        }else{
          current <- current[[1]]
        }
        
        removeModal()
        
        updateSliderInput(session, inputId = 'nfeatureslower', label = 'Filter Genes Detected/Cell (Lower Bound)', min = 0, max = 1000, value = 200, step = 10)
        updateSliderInput(session, inputId = 'nfeatureshigher', label = 'Filter Genes Detected/Cell (Upper Bound)', min = 0, max = 25000, value = 2500, step = 10)
        updateSliderInput(session, inputId = 'mito', label = 'Remove Mitochondria/Cell by Percent', min = 0, max = 100, value = 5, step = 1)
        updateSliderInput(session, inputId = 'width', label = 'Plot Width', min = 400, max = 8000, value = 780, step = 10)
        updateSliderInput(session, inputId = 'height', label = 'Plot Height', min = 400, max = 8000, value = 620, step = 10)
        return(current)
    })
    
    output$contents <- renderPrint({
            data() 
    })
    
    w <- reactive({
        req(input$width)
        input$width
    })
    
    h <- reactive({
        req(input$height)
        input$height
    })
    
    qualityControl <- reactiveVal()
    observeEvent(input$go, {
        showModal(modalDialog("Creating filtered dataset..", footer=NULL))
        x <- data()
        x <- subset(x, subset = nFeature_RNA > input$nfeatureslower & nFeature_RNA < input$nfeatureshigher & Percent_Mito < input$mito)
        seuratObj <- NULL
        seuratObj$x <- x
        print("Completed the creation of filtered object!")
        removeModal()
        qualityControl(seuratObj)
        runjs("$('#dl')[0].click();")
    })
    
    output$VlnPlot <- renderPlot({
        x <- data()
        x <- subset(x, subset = nFeature_RNA > input$nfeatureslower & nFeature_RNA < input$nfeatureshigher & Percent_Mito < input$mito)
        csamples <- length(unique(x$orig.ident))
        if(csamples > 1){
          ccols <- adjust_luminance(input$col, steps = 1.5)
          ccols <- colorRampPalette(c(ccols,input$col))(csamples)
        }else{
          ccols <- input$col
        }
        print(VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", "Percent_Mito"), ncol = 3,pt.size = 0.01,cols = ccols, group.by = "orig.ident"))
    }, width = w, height = h)
})










