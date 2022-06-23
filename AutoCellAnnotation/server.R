library("BiocManager")
options(repos = BiocManager::repositories())
library("shiny")
library("ggplot2")
library("ggrepel")
library("Seurat")
library("shinyjs")
library("SingleR")
library("SingleCellExperiment")
library("randomcoloR")
options(shiny.maxRequestSize=8000*1024^2)

create_centroids <- function(plotx, x, y, group){
    
    centroids <- split(plotx, plotx[,group])
    centroids <- lapply(centroids, function(cent){
        cent <- data.frame(Cluster = ifelse(nrow(cent) > 0,as.character(unique(cent[,group])), NA),
                           Centroid_X = ifelse(nrow(cent) > 0, median(cent[,x]), NA),
                           Centroid_Y = ifelse(nrow(cent) > 0, median(cent[,y]), NA))})
    centroids <- do.call(rbind.data.frame, centroids)
    centroids <- centroids[!is.na(centroids[,"Cluster"]),]
    return(centroids)
    
}

adjust_theme <- function(p, xangle = 0,legend = "right", title_size = 20, xsize=20){
    p <- p+ theme_classic()+
        theme(axis.text.x = element_text(size = xsize, angle = xangle, hjust=1,vjust = 1),
              axis.text.y = element_text(size = 20), legend.position = legend,
              axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
              axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)),
              legend.title = element_text(size =20, face = "bold"),
              legend.text = element_text(size = 15),
              strip.text.x = element_text(size = 20),
              strip.background = element_blank(),
              plot.title = element_text(size =title_size, face = "bold", hjust = 0.5))
    return(p)
}

ScatteredPlot <- function(plotx, x, y, group = "CELL_TYPE", facet = "None", ncol = 1, label = TRUE, labelsize = 3, annotate = "Non-Numeric Annotation", point_size = 1, colindex = 1){
    
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
        
        plotx$Label <- ""
        for(i in 1:nrow(centroids)){
            plotx[nrow(plotx)+1,x] <- centroids[i,"Centroid_X"]
            plotx[nrow(plotx),y] <- centroids[i,"Centroid_Y"]
            plotx[nrow(plotx),group] <- centroids[i,"Cluster"]
            plotx[nrow(plotx),"Label"] <- centroids[i,"Cluster"]
        }
        p <- ggplot(plotx, aes_string(x = x, y = y, color = group, label = "Label")) +
            geom_point(alpha = ifelse(plotx$Label != "", 0, 1), size = point_size)
    }else{
        p <- ggplot(plotx, aes_string(x = x, y = y, color = group)) +
            geom_point(size = point_size)
    }
    
    p <- p + 
        scale_color_manual(values = ccols)+
        theme_classic() +
        xlab(x) + ylab(y)
    
    if(label == TRUE){
        if(annotate == "Non-Numeric Annotation"){
            p <- p + geom_text_repel(box.padding = 2, max.overlaps = Inf, size = labelsize, show.legend = F) # colour = "black"
        }else{
            p <- p + annotate("text", x=centroids$Centroid_X, y=centroids$Centroid_Y, label= centroids$Cluster, size = labelsize, hjust = 0, fontface =1)
        }
    }
    
    p <- adjust_theme(p)
    if(facet != "None"){
      p <- p+facet_wrap(facet, ncol = ncol)
    }
    
    return(p)
    
}

example <- "DB/AutoCellAnnotation_Example.RDS"
ref <- "DB/HumanPrimaryCellAtlas_Reference.RDS"
hpca.se <- readRDS("DB/HumanPrimaryCellAtlas_Reference.RDS")

sys_time <- format(Sys.time(), "%Y%m%d%H%M%S")

shinyServer(function(input, output, session) {
    
    output$downloadRDS <- downloadHandler(
        filename = "HTCA_AutoCellAnnotation_Example.RDS",
        content = function(file) {
            file.copy(example, file)
        },
        contentType = "application/RDS"
    )
    
    output$downloadHPCA <- downloadHandler(
      filename = "HPCA_HumanPrimaryCellAtlas_Reference.RDS",
      content = function(file) {
        file.copy(ref, file)
      },
      contentType = "application/RDS"
    )
    
    
    output$dl2 <- downloadHandler(
        filename <- paste("HPCA_Cell_Type_DE_Analysis_",sys_time,".csv", sep = ""),
        content <- function(file) {
            x <- autoAnnot()$deout
            write.csv(x,file, row.names = T)
        },
        contentType = "text/csv"
    )
    
    output$dl <- downloadHandler(
        filename <- paste("HPCA_AutoCellAnnotation_Coordinates_",sys_time,".csv", sep = ""),
        content <- function(file) {
            x <- autoAnnot()$plotx
            write.csv(x,file, row.names = T)
        },
        contentType = "text/csv"
    )
    
    inFile <- reactive({
      if(input$exampledata){
        x <- example
      }else{
        x <- input$file1$datapath
      }
    })
    
    data <- reactive({
        req(inFile())
        showModal(modalDialog("Annotating..", footer=NULL))
        df <- readRDS(inFile())
        Idents(df) <- "seurat_clusters"
        DefaultAssay(df) <- "RNA"
        if(!is.null(input$altref)){
          cref <- readRDS(input$altref$datapath)
        }else{
          cref <- hpca.se
        }
        clu_ann <- SingleR(test = as.SingleCellExperiment(DietSeurat(df)),
                           clusters =  df$seurat_clusters,
                           ref = cref, assay.type.test=1,
                           labels = cref$label.main)
        df$CELL_TYPE <- clu_ann$labels[match(df$seurat_clusters,row.names(clu_ann))]
        df@meta.data[which(is.na(df$CELL_TYPE)),"CELL_TYPE"] <- "Unidentifiable"
        
        if(!is.null(df@reductions$pca)){
            df@meta.data <- cbind(df@meta.data,data.frame(PC_1 = df@reductions$pca@cell.embeddings[,"PC_1"],
                                            PC_2 = df@reductions$pca@cell.embeddings[,"PC_2"]))
        }
        
        if(!is.null(df@reductions$harmony)){
          df@meta.data <- cbind(df@meta.data,data.frame(HARMONY_1 = df@reductions$harmony@cell.embeddings[,"harmony_1"],
                                                        HARMONY_2 = df@reductions$harmony@cell.embeddings[,"harmony_2"]))
        }
        if(!is.null(df@reductions$tsne)){
            df@meta.data <- cbind(df@meta.data,data.frame(tSNE_1 = df@reductions$tsne@cell.embeddings[,"tSNE_1"],
                                            tSNE_2 = df@reductions$tsne@cell.embeddings[,"tSNE_2"]))
        }
        
        if(!is.null(df@reductions$umap)){
            df@meta.data <- cbind(df@meta.data,data.frame(UMAP_1 = df@reductions$umap@cell.embeddings[,"UMAP_1"],
                                            UMAP_2 = df@reductions$umap@cell.embeddings[,"UMAP_2"]))
        }
        removeModal()
        
        methods <- NULL
        if(!is.null(df@reductions$pca)){
            methods <- c(methods, "PCA")
        }
        
        if(!is.null(df@reductions$harmony)){
          methods <- c(methods, "HARMONY")
        }
        
        if(!is.null(df@reductions$tsne)){
            methods <- c(methods, "tSNE")
        }
        
        if(!is.null(df@reductions$umap)){
            methods <- c(methods, "UMAP")
        }
        
        groups <- c("None","Group","orig.ident","Batch","seurat_clusters")
        
        updateSelectInput(session, inputId = 'dimreduction', label = 'Dimension Reduction', choices = methods, selected = methods[length(methods)])
        updateSelectInput(session, inputId = 'group', 'Separate by Group', choices = groups, selected = groups[1])
        updateSliderInput(session, inputId = 'pointsize', label = 'Size of Points', min = 0, max = 10, value = 2, step = 0.1)
        updateSliderInput(session, inputId = 'labelsize', label = 'Size of Labels', min = 0, max = 20, value = 8, step = 0.1)
        updateSliderInput(session, inputId = 'colindex', label = 'Choose A Color Scheme Index', min = 1, max = 100, value = 1, step = 1)
        updateSliderInput(session, inputId = 'ncols', label = 'Number of Columns (for multiple groups)', min = 1, max = 10, value = 1, step = 1)
        updateCheckboxInput(session, inputId = 'label', label = 'Label Cell Types', value = TRUE)
        updateSliderInput(session, inputId = 'width', label = 'Plot Width', min = 400, max = 8000, value = 780, step = 10)
        updateSliderInput(session, inputId = 'height', label = 'Plot Height', min = 400, max = 8000, value = 480, step = 10)
        
        return(df)
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
    
    autoAnnot <- reactiveVal()
    observeEvent(input$go, {
        showModal(modalDialog("Differential expression analysis is in process to identify top markers specific to each cell types (May take some time depending on sample size)..Download for post-cell annotation cell type coordinates CSV file and DE genes CSV file will start automatically once the calculations are done.", footer=NULL))
        x <- data()
        plotx <- x@meta.data
        plotx <- plotx[,c(grep("^PC_1$|^PC_2$|^HARMONY_1$|^HARMONY_2$|^tSNE_1$|^tSNE_2$|^UMAP_1$|^UMAP_2$", colnames(plotx), ignore.case = T), grep("^PC_1$|^PC_2$|^HARMONY_1$|^HARMONY_2$|^tSNE_1$|^tSNE_2$|^UMAP_1$|^UMAP_2$", colnames(plotx), ignore.case = T, invert = T))]
        
        DefaultAssay(x) <- "RNA"
        Idents(x) <- "CELL_TYPE"
        deout <- FindAllMarkers(x, min.pct = 0.25, logfc.threshold = 0.25)
        deout <- deout[order(deout$p_val_adj, decreasing = F),]
        
        seuratObj <- NULL
        seuratObj$plotx <- plotx
        seuratObj$deout <- deout
        print("Completed Auto-Cell Annotation and DE Analysis!")
        removeModal()
        # return(seuratObj)
        # autoClust(plotx)
        autoAnnot(seuratObj)
        
        runjs("$('#dl')[0].click();")
        runjs("$('#dl2')[0].click();")
    })
    
    output$ScatteredPlot <- renderPlot({
        x <- data()
        plotx <- x@meta.data
        cmethod <- gsub("PCA","PC",input$dimreduction, ignore.case = T)
        print(ScatteredPlot(plotx, x = paste(cmethod,"_1",sep = ""), y = paste(cmethod,"_2",sep = ""),
                                label = input$label, labelsize = input$labelsize,
                                facet = input$group, ncol = input$ncols,
                                point_size = input$pointsize, colindex = input$colindex))
    }, width = w, height = h)
})










