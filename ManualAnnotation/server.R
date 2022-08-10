library(BiocManager)
options(repos = BiocManager::repositories())

library("shiny")
library("shinyscreenshot")
library("shinyalert")
library("ggplot2")
library("ggrepel")
library("Seurat")
library("shinyjs")
library("randomcoloR")
library("DT")
library("plotly")

options(shiny.maxRequestSize=3000*1024^2)

example <- "DB/Manual_Annotation_Example_10X_Dataset.RDS"

sys_time <- format(Sys.time(), "%Y%m%d%H%M%S")

create_centroids <- function(plotx, x, y, group){
  centroids <- split(plotx, plotx[,group])
  centroids <- lapply(centroids, function(cent){
    cent <- data.frame(Cluster = ifelse(nrow(cent) > 0,as.character(unique(cent[,group])), NA),
                       Centroid_X = ifelse(nrow(cent) > 0, median(cent[,x]), NA),
                       Centroid_Y = ifelse(nrow(cent) > 0, median(cent[,y]), NA))})
  centroids <- do.call(rbind.data.frame, centroids)
  colnames(centroids) <- c(group,x,y)
  centroids <- centroids[!is.na(centroids[,group]),]
  return(centroids)
}

adjust_theme <- function(p, xangle = 0,legend = "right", title_size = 20, xsize=20){
  p <- p+ theme_classic()+
    theme(axis.text.x = element_text(size = xsize-5, angle = xangle, hjust=1,vjust = 1),
          axis.text.y = element_text(size = xsize-5), legend.position = legend,
          axis.title.x = element_text(size = xsize, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = xsize, margin=margin(0,10,0,0)),
          legend.title = element_text(size =xsize, face = "bold"),
          legend.text = element_text(size = xsize-5),
          strip.text.x = element_text(size = xsize),
          strip.background = element_blank(),
          plot.title = element_text(size =title_size, face = "bold", hjust = 0.5))
  return(p)
}

ScatteredPlot <- function(plotx, x, y, group = "CELL_TYPE", facet = "None", ncol = 1, label = TRUE, labelsize = 3, axislabelsize = 20, annotate = "Non-Numeric Annotation", point_size = 1, colindex = 1){
  
  set.seed(colindex)
  ccols <- distinctColorPalette(length(unique(plotx[,group])))
  if(label == TRUE){
    if(annotate == "Numeric Annotation"){
      plotx[,group] <- factor(plotx[,group], levels = sort(unique(as.numeric(as.character(plotx[,group])))))
    }else{
      plotx[,group] <- factor(plotx[,group], levels = sort(unique(as.character(plotx[,group]))))
    }
    centroids <- create_centroids(plotx, x, y, group)
    centroids$Col <- ccols
    p <- ggplot(plotx, aes_string(x = x, y = y, color = group)) +
      geom_point(size = point_size)
    p <- p + geom_text_repel(data = centroids, aes_string(x = x, y = y, label = group), box.padding = 2, max.overlaps = Inf, size = labelsize, show.legend = F)
  }else{
    p <- ggplot(plotx, aes_string(x = x, y = y, color = group)) +
      geom_point(size = point_size)
  }
  
  p <- p + 
    scale_color_manual(values = ccols)+
    theme_classic() +
    xlab(x) + ylab(y)
  
  p <- adjust_theme(p, title_size = axislabelsize, xsize = axislabelsize)
  if(facet != "None"){
    p <- p+facet_wrap(facet, ncol = ncol)
  }
  
  return(p)
  
}

shinyServer(function(input, output, session) {
    
    inFile <- reactive({
      if(input$exampledata){
        x <- example
      }else{
        x <- input$file1$datapath
      }
    })
    
    data <- reactive({
        req(inFile())
        showModal(modalDialog("Initializing..Proceed to tab 'Plot' after this process.", footer=NULL))
        x <- readRDS(inFile())
        DefaultAssay(x) <- "RNA"
        clusters <- colnames(x@meta.data)[grep("CELL_TYPE|seurat_clusters|CELL_TYPE_MANUAL", colnames(x@meta.data), ignore.case = T)]
        groups <- colnames(x@meta.data)[grep("Group|orig.ident|Batch|seurat_clusters|CELL_TYPE", colnames(x@meta.data), ignore.case = T)]
        Idents(x) <- "seurat_clusters"
        if(is.null(x@assays$integrated)){
          DefaultAssay(x) <- "integrated"
        }else{
          DefaultAssay(x) <- "RNA"
        }
        
        methods <- NULL
        if(!is.null(x@reductions$pca)){methods <- c(methods, "PCA")}
        if(!is.null(x@reductions$harmony)){methods <- c(methods, "HARMONY")}
        if(!is.null(x@reductions$tsne)){methods <- c(methods, "tSNE")}
        if(!is.null(x@reductions$umap)){methods <- c(methods, "UMAP")}
        removeModal()
        
        updateSelectInput(session, inputId = 'dimreduction', label = 'View by Dimension Reduction Type', choices = methods, selected = methods[length(methods)])
        updateSelectInput(session, inputId = 'clusters', label = 'Choose an Annotation Group to View', choices = clusters, selected = clusters[1])
        updateSliderInput(session, inputId = 'res', label = 'Resolution Level for Auto-Clustering', min = 0, max = 10, value = 0.8, step = 0.05)
        updateSelectInput(session, inputId = 'group', 'Separate by Group', choices = c("None",groups), selected = "None")
        updateSliderInput(session, inputId = 'padj', label = 'Keep Genes with Adjusted P-Value Smaller Than', min = 0, max = 1, value = 0.05, step = 0.01)
        updateSliderInput(session, inputId = 'log2fc', label = 'Keep Genes with Log2FC Greater Than', min = -20, max = 20, value = 0.25, step = 0.1)
        updateSliderInput(session, inputId = 'pointsize', label = 'Size of Points', min = 0.1, max = 10, value = 0.5, step = 0.1)
        updateCheckboxInput(session, inputId = 'label', label = 'Label Cluster/Cell Type (For Static Plot)', value = TRUE)
        updateSliderInput(session, inputId = 'labelsize', label = 'Size of Labels', min = 0.1, max = 20, value = 10, step = 0.1)
        updateSliderInput(session, inputId = 'axislabelsize', label = 'Size of Axis Labels', min = 10, max = 50, value = 20, step = 0.1)
        updateSliderInput(session, inputId = 'colindex', label = 'Choose A Color Scheme Index', min = 1, max = 100, value = 1, step = 1)
        updateSliderInput(session, inputId = 'ncols', label = 'Number of Columns (for multiple groups)', min = 1, max = 10, value = 2, step = 1)
        updateSliderInput(session, inputId = 'width', label = 'Plot Width', min = 400, max = 8000, value = 900, step = 10)
        updateSliderInput(session, inputId = 'height', label = 'Plot Height', min = 200, max = 8000, value = 550, step = 10)
        return(x)
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
    
    DEAnalysis <- reactiveValues(DEresults = NULL)
    observeEvent(input$go, {
      showModal(modalDialog("Differential expression analysis of the selection annotation group is in process (May take some time depending on sample size)..", footer=NULL))
      x <- data()
      DefaultAssay(x) <- "RNA"
      Idents(x) <- input$clusters
      deout <- FindAllMarkers(x, min.pct = 0.25, logfc.threshold = 0.25)
      deout[,c(grep("cluster|gene",colnames(deout), ignore.case = T),grep("cluster|gene",colnames(deout), ignore.case = T, invert = T))]
      deout <- deout[order(deout$avg_log2FC, decreasing = F),]
      print("Completed Auto-Cell Annotation and DE Analysis!")
      removeModal()
      
      DEAnalysis$DEresults$x <- x
      DEAnalysis$DEresults$deout <- deout
      
    })
    
    ManualDE <- reactiveValues(DEresults = NULL)
    observeEvent(input$runManualDE, {
      showModal(modalDialog("Differential expression analysis of the manually annotated group is in process (May take some time depending on sample size)..", footer=NULL))
      if(input$go){
        x <- DEAnalysis$DEresults$x
      }else{
        x <- data()
      }
      DefaultAssay(x) <- "RNA"
      cannot <- ManualAnnot$annot
      colnames(cannot) <- c("cluster","manual")
      cannot[is.na(cannot)] <- "NA"
      x$CELL_TYPE_MANUAL <- cannot[match(x$seurat_clusters, cannot$cluster),"manual"]
      Idents(x) <- "CELL_TYPE_MANUAL"
      deout <- FindAllMarkers(x, min.pct = 0.25, logfc.threshold = 0.25)
      deout[,c(grep("cluster|gene",colnames(deout), ignore.case = T),grep("cluster|gene",colnames(deout), ignore.case = T, invert = T))]
      deout <- deout[order(deout$avg_log2FC, decreasing = F),]
      print("Completed Auto-Cell Annotation and DE Analysis!")
      removeModal()
      
      ManualDE$DEresults$x <- x
      ManualDE$DEresults$deout <- deout
    })
    
    ManualUpdate <- reactiveValues(x = NULL)
    observeEvent(input$updateManual, {
      showModal(modalDialog("Updating annotation..", footer=NULL))
      cannot <- ManualAnnot$annot
      if(!all(cannot[,"Manual Annotation (Double Click to Fill in Accordingly)"] == "")){
        if(input$go){
          x <- DEAnalysis$DEresults$x
        }else{
          x <- data()
        }
        
        DefaultAssay(x) <- "RNA"
        cannot <- ManualAnnot$annot
        colnames(cannot) <- c("cluster","manual")
        cannot[is.na(cannot)] <- "NA"
        x$CELL_TYPE_MANUAL <- cannot[match(x$seurat_clusters, cannot$cluster),"manual"]
        Idents(x) <- "CELL_TYPE_MANUAL"
        print("Completed Updating Annotation!")
        ManualUpdate$x <- x
      }else{
        shinyalert("No Mannual Annotation Provided", "Something went wrong.", type = "error")
      }
      removeModal()
      })
    
    ScatteredPlotObject <- reactive({
      cclust <- input$clusters
      if(!is.null(ManualUpdate$x)){
        x <- ManualUpdate$x
        cclust <- "CELL_TYPE_MANUAL"
      }else if(input$go){
        x <- DEAnalysis$DEresults$x
      }else{
        x <- data()
      }
      
      if(!is.null(x@assays$integrated)){
        DefaultAssay(x) <- "integrated"
      }
      Idents(x) <- "orig.ident"
      x <- FindClusters(x, resolution = input$res)
      plotx <- x@meta.data
      if(length(grep("^PC_1$|^PC_2$|^HARMONY_1$|^HARMONY_2$|^tSNE_1$|^tSNE_2$|^UMAP_1$|^UMAP_2$", colnames(plotx), ignore.case = T)) == 0){
        if(!is.null(x@reductions$pca)){
          plotx <- cbind(plotx,data.frame(PC_1 = x@reductions$pca@cell.embeddings[,"PC_1"],
                                          PC_2 = x@reductions$pca@cell.embeddings[,"PC_2"]))
        }
        
        if(!is.null(x@reductions$tsne)){
          plotx <- cbind(plotx,data.frame(tSNE_1 = x@reductions$tsne@cell.embeddings[,"tSNE_1"],
                                          tSNE_2 = x@reductions$tsne@cell.embeddings[,"tSNE_2"]))
        }
        
        if(!is.null(x@reductions$umap)){
          plotx <- cbind(plotx,data.frame(UMAP_1 = x@reductions$umap@cell.embeddings[,"UMAP_1"],
                                          UMAP_2 = x@reductions$umap@cell.embeddings[,"UMAP_2"]))
        }
      }
      
      cmethod <- gsub("PCA","PC",input$dimreduction, ignore.case = T)
      if(cclust != "seurat_clusters"){
        cannot <- "Non-Numeric Annotation"
      }else{
        cannot <- "Numeric Annotation"
      }
      p <- ScatteredPlot(plotx, x = paste(cmethod,"_1",sep = ""), y = paste(cmethod,"_2",sep = ""),
                         label = input$label, labelsize = input$labelsize,
                         facet = input$group, ncol = input$ncols, group = cclust,
                         annotate = cannot, axislabelsize = input$axislabelsize,
                         point_size = input$pointsize, colindex = input$colindex)
    })
    
    output$ScatteredPlot <- renderPlot({
      print(ScatteredPlotObject())
    }, width = w, height = h)
    
    output$ScatteredPlotly <- renderPlotly({
      return(ggplotly(ScatteredPlotObject(), width = w(), height = h()))
    })
    
    output$dlpng <- downloadHandler(
      filename = function() { paste("HTCA_Manual_Annotation_",input$dimreduction,"_",input$clusters,"_",sys_time,".png", sep = "")},
      content <- function(file) {
        png(file, width=w()*5, height=h()*5, res = 300)
        print(ScatteredPlotObject())
        dev.off()
      }
    )
    
    output$dlpdf <- downloadHandler(
      filename = function() { paste("HTCA_Manual_Annotation_",input$dimreduction,"_",input$clusters,"_",sys_time,".pdf", sep = "")},
      content <- function(file) {
        pdf(file, width=w()/72, height=h()/72)
        print(ScatteredPlotObject())
        dev.off()
      }
    )
    
    output$downloadRDS <- downloadHandler(
      filename = "Manual_Annotation_Example_10X_Dataset.RDS",
      content = function(file) {
        file.copy(example, file)
      },
      contentType = "application/RDS"
    )
    
    output$dl1 <- downloadHandler(
      filename <- paste("HTCA_Post_Manual_Annotation_Seurat_Object_",sys_time,".RDS", sep = ""),
      content <- function(file) {
        showModal(modalDialog("Create download object..", footer=NULL))
        if(!is.null(ManualUpdate$x)){
          x <- ManualUpdate$x
          cclust <- "CELL_TYPE_MANUAL"
        }else if(input$go){
          x <- DEAnalysis$DEresults$x
        }else{
          x <- data()
        }
        saveRDS(x,file)
        removeModal()
      },
      contentType = "application/RDS"
    )
    
    output$dl2 <- downloadHandler(
      filename <- paste("HTCA_Post_Manual_Annotation_DE_Analysis_By_Cluster_Resolution_",input$res,"_",sys_time,".csv", sep = ""),
      content <- function(file) {
        if(input$go){
          x <- DEAnalysis$DEresults$deout
          write.csv(x,file, row.names = F, quote = F)
        }else{
          shinyalert("DE Analysis Not Run", "Something went wrong.", type = "error")
        }
      },
      contentType = "text/csv"
    )
    
    output$dl3 <- downloadHandler(
      filename <- paste("HTCA_Post_Manual_Annotation_DE_Analysis_By_Manual_Annotation_",sys_time,".csv", sep = ""),
      content <- function(file) {
        if(input$runManualDE){
          x <- ManualDE$DEresults$deout
          write.csv(x,file, row.names = T, quote = F)
        }else{
          shinyalert("DE Analysis for Manual Annotation Not Run", "Something went wrong.", type = "error")
        }
      },
      contentType = "text/csv"
    )

    deout <- reactive({
      x <- DEAnalysis$DEresults$deout
      x <- x[which(x$p_val_adj < input$padj & x$avg_log2FC > input$log2fc),]
      return(x)
    })
    
    manualdeout <- reactive({
      x <- ManualDE$DEresults$deout
      x <- x[which(x$p_val_adj < input$padj & x$avg_log2FC > input$log2fc),]
      return(x)
    })
    
    output$manualtable <- renderDT(manualdeout(),
                             filter = "top",
                             style="bootstrap",
                             rownames = F,
                             options = list(
                               pageLength = 5)
    )
    
    output$table <- renderDT(deout(),
                             filter = "top",
                             style="bootstrap",
                             rownames = F,
                             options = list(
                               pageLength = 5)
    )
    
    output$form <- renderDT({
      x <- deform()
      datatable(x, editable = TRUE, rownames = F, style = "bootstrap4")
    })
    
    deform <- reactive({
      if(input$go){
        x <- DEAnalysis$DEresults$x
      }else{
        x <- data()
        if(!is.null(x@assays$integrated)){
          DefaultAssay(x) <- "integrated"
        }
        Idents(x) <- "orig.ident"
        x <- FindClusters(x, resolution = input$res)
      }
      x <- data.frame(Cluster = sort(unique(x@meta.data[,input$clusters])), Manual_Annotation = "")
      colnames(x)[which(colnames(x) == "Manual_Annotation")] <- "Manual Annotation (Double Click to Fill in Accordingly)"
      return(x)
    })
    
    ManualAnnot <- reactiveValues(annot = NULL)
    observeEvent(input$form_cell_edit, {
      if(is.null(ManualAnnot$annot)){
        ManualAnnot$annot <- deform()
      }
      row  <- input$form_cell_edit$row
      ManualAnnot$annot[row, "Manual Annotation (Double Click to Fill in Accordingly)"] <- input$form_cell_edit$value
      print(ManualAnnot$annot)
    })

})



