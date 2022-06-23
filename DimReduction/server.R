library(BiocManager)
options(repos = BiocManager::repositories())

library("shiny")
library("ggplot2")
library("Seurat")
library("shinyjs")
library("randomcoloR")

options(shiny.maxRequestSize=8000*1024^2)

exampleRDS <- "DB/DimReduction_Example_10X_Dataset.RDS"
sys_time <- format(Sys.time(), "%Y%m%d%H%M%S")

ScatteredPlot <- function(plotx, x, y, group = "orig.ident", point_size = 1, colindex = 1){
  set.seed(colindex)
  ccols <- randomColor(length(unique(plotx[,group])))
    p <- ggplot(plotx, aes_string(x = x, y = y, color = group)) +
      geom_point(size = point_size)
  p <- p + 
    scale_color_manual(values = ccols)+
    theme_classic(base_size = 25) +
    xlab(x) + ylab(y)
  return(p)
}

shinyServer(function(input, output, session) {
    
    output$downloadRDS <- downloadHandler(
        filename = "HTCA_DimReduction_Example_10X_Dataset.RDS",
        content = function(file) {
            file.copy(exampleRDS, file)
        },
        contentType = "application/RDS"
    )
    
    output$dl2 <- downloadHandler(
        filename <- paste("HTCA_DimReduction_Seurat_Object_",sys_time,".RDS", sep = ""),
        content <- function(file) {
            x <- dimReduct()$x
            saveRDS(x,file)
        },
        contentType = "application/RDS"
    )
    
    output$dl <- downloadHandler(
        filename <- paste("HTCA_DimReduction_Coordinates_",sys_time,".csv", sep = ""),
        content <- function(file) {
            x <- dimReduct()$plotx
            write.csv(x,file, row.names = T)
        },
        contentType = "text/csv"
    )
    
    data <- reactive({
        req(input$file1)
        req(input$run)
      showModal(modalDialog("Computing..(if tSNE is selected, calculations may take much longer time)..Download for dimension reduction coordinates CSV output and Seurat object RDS file will start automatically once the calculations are done.", footer=NULL))
      inFile <- input$file1
        x <- readRDS(inFile$datapath)
        if(length(grep("integrated", Assays(x), ignore.case = T)) > 0){
          DefaultAssay(x) <- "integrated"
        }
        cpca <- "pca"
        if(!is.null(x@reductions$harmony)){
          cpca <- "harmony"
        }
        x <- ScaleData(x)
        m <- ifelse(ncol(x) < 50, ncol(x) - 2, 50)
        x <- RunPCA(x, npcs = m)
        plotx <- data.frame(PC_1 = x@reductions$pca@cell.embeddings[,"PC_1"],
                            PC_2 = x@reductions$pca@cell.embeddings[,"PC_2"])
        methods <- "PCA"
        if(cpca == "harmony"){
          plotx <- cbind(plotx,
                         HARMONY_1 = x@reductions$harmony@cell.embeddings[,"harmony_1"],
                         HARMONY_2 = x@reductions$harmony@cell.embeddings[,"harmony_2"])
          methods <- c(methods, "HARMONY")
        }
        if(input$runUMAP == TRUE){
          x <- RunUMAP(x, reduction = cpca, dims = 1:m)
          plotx <- cbind(plotx,
                         UMAP_1 = x@reductions$umap@cell.embeddings[,"UMAP_1"],
                         UMAP_2 = x@reductions$umap@cell.embeddings[,"UMAP_2"])
          methods <- c(methods, "UMAP")
        }
        if(input$runTSNE == TRUE){
          x <- RunTSNE(x, reduction = cpca, dims = 1:m, check_duplicates = FALSE)
          plotx <- cbind(plotx,
                         tSNE_1 = x@reductions$tsne@cell.embeddings[,"tSNE_1"],
                         tSNE_2 = x@reductions$tsne@cell.embeddings[,"tSNE_2"])
          methods <- c(methods, "tSNE")
        }
        
        plotx <- cbind(plotx, x@meta.data)
        removeModal()
        
        updateSelectInput(session, inputId = 'dimreduction', label = 'Dimension Reduction', choices = methods, selected = methods[length(methods)])
        updateSliderInput(session, inputId = "colorindex", label = "Color Scheme Index", min = 1, max = 100, value = 1, step = 1)
        updateSliderInput(session, inputId = 'pointsize', label = 'Point Size', min = 0.1, max = 10, value = 0.5, step = 0.1)
        updateSliderInput(session, inputId = 'width', label = 'Plot Width', min = 400, max = 8000, value = 840, step = 10)
        updateSliderInput(session, inputId = 'height', label = 'Plot Height', min = 400, max = 8000, value = 460, step = 10)
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
    
    dimReduct <- reactiveVal()
    observeEvent(input$go, {
        print("Preparing to download post-dimensional reduction data and coordinates CSV file..")
        seuratObj <- NULL
        seuratObj$x <- data()$x
        seuratObj$plotx <- data()$plotx
        removeModal()
        dimReduct(seuratObj)
        runjs("$('#dl')[0].click();")
        runjs("$('#dl2')[0].click();")
    })
    
    output$ScatteredPlot <- renderPlot({
      plotx <- data()$plotx
      plotx <- plotx[,c(grep("^PC_1$|^PC_2$|^HARMONY_1$|^HARMONY_2$|^tSNE_1$|^tSNE_2$|^UMAP_1$|^UMAP_2$", colnames(plotx), ignore.case = T), grep("^PC_1$|^PC_2$|^HARMONY_1$|^HARMONY_2$|^tSNE_1$|^tSNE_2$|^UMAP_1$|^UMAP_2$", colnames(plotx), ignore.case = T, invert = T))]
      cmethod <- gsub("PCA","PC",input$dimreduction, ignore.case = T)
      print(ScatteredPlot(plotx, x = paste(cmethod,"_1",sep = ""), y = paste(cmethod,"_2",sep = ""), point_size = input$pointsize, colindex = input$colorindex))
    }, width = w, height = h)
    
})










