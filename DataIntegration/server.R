library(BiocManager)
options(repos = BiocManager::repositories())

library("shiny")
library("ggplot2")
library("randomcoloR")
library("Seurat")
library("shinyjs")
library("harmony")

options(shiny.maxRequestSize=8000*1024^2)

exampleRDS <- "DB/DataIntegration_Example_Dataset.RDS"
sys_time <- format(Sys.time(), "%Y%m%d%H%M%S")

ScatteredPlot <- function(plotx, x, y, group = "Batch", point_size = 1, cols, facet = T, facet_group = "DATA_TYPE"){
  
    p <- ggplot(plotx, aes_string(x = x, y = y, color = group)) +
      geom_point(size = point_size)
  
  p <- p + 
    scale_color_manual(values = cols)+
    theme_classic(base_size = 25) +
    xlab(x) + ylab(y)
  
  if(facet == T){
    p <- p+facet_wrap(facet_group)
  }
  
  return(p)
  
}

shinyServer(function(input, output, session) {
    
    output$downloadRDS <- downloadHandler(
        filename = "HTCA_DataIntegration_Example_Dataset.RDS",
        content = function(file) {
            file.copy(exampleRDS, file)
        },
        contentType = "application/RDS"
    )
    
    output$dl <- downloadHandler(
      filename <- paste("HTCA_Post_Integration_Seurat_Object_",sys_time,".RDS", sep = ""),
      content <- function(file) {
        x <- postIntegrate()$x
        saveRDS(x,file)
      },
      contentType = "application/RDS"
    )
    
    data <- reactive({
        req(input$file1)
        showModal(modalDialog("Uploading..", footer=NULL))
        inFile <- input$file1
        removeModal()

        showModal(modalDialog("Running integration..", footer=NULL))
        x <- readRDS(inFile$datapath)
        x <- NormalizeData(x, verbose = F)
        x <- FindVariableFeatures(x, verbose = F)
        x <- ScaleData(x, verbose = F)
        m <- ifelse(ncol(x) < 50, ncol(x) - 2, 50)
        x <- RunPCA(x, npcs = m)
        x <- RunUMAP(x, reduction = "pca", dims = 1:m)
        x <- RunTSNE(x, reduction = "pca", dims = 1:m, check_duplicates = FALSE)
        plotx <- data.frame(DATA_TYPE = "Before Integration",
                            UMAP_1 = x@reductions$umap@cell.embeddings[,"UMAP_1"],
                            UMAP_2 = x@reductions$umap@cell.embeddings[,"UMAP_2"],
                            tSNE_1 = x@reductions$tsne@cell.embeddings[,"tSNE_1"],
                            tSNE_2 = x@reductions$tsne@cell.embeddings[,"tSNE_2"],
                            Batch = x$Batch)
        
        if(input$integrationmethod == "Seurat V3"){
          x <- SplitObject(x, split.by = "Batch")
          x <- lapply(x, function(y) {
            y <- NormalizeData(y)
            y <- FindVariableFeatures(y)
          })
          
          canchor <- FindIntegrationAnchors(x, anchor.features = SelectIntegrationFeatures(x))
          x <- IntegrateData(canchor)
          DefaultAssay(x) <- "integrated"
          x <- ScaleData(x)
          x <- RunPCA(x, npcs = m)
          x <- RunUMAP(x, reduction = "pca", dims = 1:m)
          x <- RunTSNE(x, reduction = "pca", dims = 1:m, check_duplicates = FALSE)
        }else{
          x <- RunHarmony(x, group.by.vars = "Batch")
          x <- RunUMAP(x, reduction = "harmony", dims = 1:m)
          x <- RunTSNE(x, reduction = "harmony", dims = 1:m, check_duplicates = FALSE)
        }
        
        plotx <- rbind(plotx,data.frame(DATA_TYPE = "After Integration",
                            UMAP_1 = x@reductions$umap@cell.embeddings[,"UMAP_1"],
                            UMAP_2 = x@reductions$umap@cell.embeddings[,"UMAP_2"],
                            tSNE_1 = x@reductions$tsne@cell.embeddings[,"tSNE_1"],
                            tSNE_2 = x@reductions$tsne@cell.embeddings[,"tSNE_2"],
                            Batch = x$Batch))
        plotx$DATA_TYPE <- factor(plotx$DATA_TYPE, levels = c("Before Integration",
                                                              "After Integration"))
        removeModal()
        
        updateSliderInput(session, inputId = 'pointsize', label = 'Choose Point Size', min = 0.1, max = 5, value = 0.3, step = 0.1)
        # updateSliderInput(session, inputId = 'colornum', label = "Choose Color Scheme Index", min = 1, max = 100, value = 1, step = 1)
        updateSliderInput(session, inputId = 'height', label = 'Plot Height', min = 400, max = 8000, value = 400, step = 10)
        updateSliderInput(session, inputId = 'height', label = 'Plot Height', min = 400, max = 8000, value = 400, step = 10)
        
        data <- NULL
        data$x <- x
        data$plotx <- plotx
        return(data)
    })
    
    output$contents <- renderPrint({
            data()$x
    })
    
    w <- reactive({
        req(input$width)
        input$width
    })
    
    h <- reactive({
        req(input$height)
        input$height
    })
    
    postIntegrate <- reactiveVal()
    observeEvent(input$go, {
      seuratObj <- NULL
      seuratObj$x <- data()$x
        postIntegrate(seuratObj)
        runjs("$('#dl')[0].click();")
    })
    
    output$umap <- renderPlot({
        plotx <- data()$plotx
        set.seed(9)
        ccols <- distinctColorPalette(k = length(unique(plotx$Batch)))
        print(ScatteredPlot(plotx, x = "UMAP_1", y = "UMAP_2",
                                group = "Batch", facet = T, facet_group = "DATA_TYPE",
                                point_size = input$pointsize, cols = ccols))
    }, width = w, height = h)
    
    output$tsne <- renderPlot({
      plotx <- data()$plotx
      set.seed(9)
      ccols <- distinctColorPalette(k = length(unique(plotx$Batch)))
      print(ScatteredPlot(plotx, x = "tSNE_1", y = "tSNE_2",
                          group = "Batch", facet = T, facet_group = "DATA_TYPE",
                          point_size = input$pointsize, cols = ccols))
    }, width = w, height = h)
})










