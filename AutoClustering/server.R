library(BiocManager)
options(repos = BiocManager::repositories())

library("shiny")
library("ggplot2")
library("ggrepel")
library("Seurat")
library("shinyjs")
library("randomcoloR")

options(shiny.maxRequestSize=3000*1024^2)

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
        # ggtitle(title) +
        # scale_fill_manual(values = col) +
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

sca_ScatteredPlot <- function(plotx, x, y, group = "CLUSTERS", label = TRUE, labelsize = 3, annotate = "Numeric Annotation", point_size = 2, colindex = 1){
    
    if(label == TRUE){
        if(annotate == "Numeric Annotation"){
            plotx[,group] <- factor(plotx[,group], levels = sort(unique(as.numeric(as.character(plotx[,group])))))
        }else{
            plotx[,group] <- factor(plotx[,group], levels = sort(unique(as.character(plotx[,group]))))
        }
        centroids <- create_centroids(plotx, x, y, group)
        set.seed(colindex)
        ccols <- distinctColorPalette(length(unique(plotx[,group])))
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
            p <- p + geom_text_repel(max.overlaps = Inf, size = labelsize) # colour = "black"
        }else{
            p <- p + annotate("text", x=centroids$Centroid_X, y=centroids$Centroid_Y, label= centroids$Cluster, size = labelsize, hjust = 0, fontface =1)
        }
    }
    
    p <- adjust_theme(p)
    
    return(p)
    
}

example <- "DB/AutoClustering_Example_10X_Dataset.RDS"

sys_time <- format(Sys.time(), "%Y%m%d%H%M%S")

shinyServer(function(input, output, session) {
    
    output$downloadRDS <- downloadHandler(
        filename = "AutoClustering_Example_10X_Dataset.RDS",
        content = function(file) {
            file.copy(example, file)
        },
        contentType = "application/RDS"
    )
    
    output$dl2 <- downloadHandler(
        filename <- paste("HTCA_AutoClustering_Seurat_Object_",sys_time,".RDS", sep = ""),
        content <- function(file) {
            x <- autoClust()$x
            saveRDS(x,file)
        },
        contentType = "application/RDS"
    )
    
    output$dl <- downloadHandler(
        filename <- paste("HTCA_Top_Markers_Post_Clustering_",sys_time,".csv", sep = ""),
        content <- function(file) {
            x <- autoClust()$deout
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
        showModal(modalDialog("Processing..", footer=NULL))
        df <- readRDS(inFile())
        cpca <- "pca"
        if(!is.null(df@reductions$harmony)){
          cpca <- "harmony"
        }
        m <- ncol(df@reductions[[cpca]]@cell.embeddings)
        df <- FindNeighbors(df, reduction = cpca, dims = 1:m)
        removeModal()
        
        methods <- NULL
        if(!is.null(df@reductions$pca)){
            methods <- c(methods, "PCA")
        }
        
        if(!is.null(df@reductions$tsne)){
            methods <- c(methods, "tSNE")
        }
        
        if(!is.null(df@reductions$umap)){
            methods <- c(methods, "UMAP")
        }
        
        updateSelectInput(session, inputId = 'dimreduction', label = 'Dimension Reduction', choices = methods, selected = methods[length(methods)])
        updateSliderInput(session, inputId = 'res', label = 'Resolution Level for Cluster Numbers', min = 0, max = 10, value = 0.8, step = 0.05)
        # updateCheckboxInput(session, inputId = 'runClusters', label = 'Run Auto-Clustering', value = TRUE)
        # updateCheckboxInput(session, inputId = 'runDE', label = 'DE Analysis', value = FALSE)
        # updateCheckboxInput(session, inputId = 'runTSNE', label = 'Run tSNE (longer run time)', value = FALSE)
        updateSliderInput(session, inputId = 'pointsize', label = 'Size of Points', min = 0, max = 10, value = 2, step = 0.1)
        updateSliderInput(session, inputId = 'labelsize', label = 'Size of Labels', min = 0, max = 20, value = 8, step = 0.1)
        updateCheckboxInput(session, inputId = 'label', label = 'Label Clusters', value = TRUE)
        # updateCheckboxInput(session, inputId = 'iflog', label = 'Log Scaling Violin Plot', value = TRUE)
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
    
    autoClust <- reactiveVal()
    observeEvent(input$go, {
        showModal(modalDialog("Differential expression analysis is in process to identify top markers specific to each clusters(May take some time depending on the number of clusters and sample size)..Download for post-clustering Seurat object RDS file and DE genes CSV file will start automatically once the calculations are done.", footer=NULL))
        x <- data()
        if(!is.null(x@assays$integrated)){
          DefaultAssay(x) <- "integrated"
        }
        Idents(x) <- "orig.ident"
        x <- FindClusters(x, resolution = input$res)
        DefaultAssay(x) <- "RNA"
        Idents(x) <- "seurat_clusters"
        deout <- FindAllMarkers(x, min.pct = 0.25, logfc.threshold = 0.25)
        deout <- deout[order(deout$p_val_adj, decreasing = F),]
        
        seuratObj <- NULL
        seuratObj$x <- x
        seuratObj$deout <- deout
        print("Completed Auto-clustering and DE Analysis!")
        removeModal()
        # return(seuratObj)
        # autoClust(plotx)
        autoClust(seuratObj)
        
        runjs("$('#dl')[0].click();")
        runjs("$('#dl2')[0].click();")
    })
    
    output$ScatteredPlot <- renderPlot({
        x <- data()
        if(!is.null(x@assays$integrated)){
          DefaultAssay(x) <- "integrated"
        }
        Idents(x) <- "orig.ident"
        x <- FindClusters(x, resolution = input$res)
        plotx <- x@meta.data
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
        
        plotx <- plotx[,c(grep("^PC_1$|^PC_2$|^HARMONY_1$|^HARMONY_2$|^tSNE_1$|^tSNE_2$|^UMAP_1$|^UMAP_2$", colnames(plotx), ignore.case = T), grep("^PC_1$|^PC_2$|^HARMONY_1$|^HARMONY_2$|^tSNE_1$|^tSNE_2$|^UMAP_1$|^UMAP_2$", colnames(plotx), ignore.case = T, invert = T))]
        plotx$CLUSTERS <- plotx$seurat_clusters
        cmethod <- gsub("PCA","PC",input$dimreduction, ignore.case = T)
        print(sca_ScatteredPlot(plotx, x = paste(cmethod,"_1",sep = ""), y = paste(cmethod,"_2",sep = ""),
                                label = input$label, labelsize = input$labelsize,
                                point_size = input$pointsize, colindex = input$colindex))
    }, width = w, height = h)
})










