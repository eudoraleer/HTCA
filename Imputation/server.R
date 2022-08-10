library(BiocManager)
options(repos = BiocManager::repositories())

library("shiny")
library("shinyscreenshot")
library("ggplot2")
library("ggrepel")
library("Seurat")
library("patchwork")
library("shinyjs")
library("randomcoloR")
# library("fastRPCA")

source("DB/alra.R")

options(shiny.maxRequestSize=3000*1024^2)

example <- "DB/Imputation_Example_10X_Dataset.RDS"

sys_time <- format(Sys.time(), "%Y%m%d%H%M%S")

shinyServer(function(input, output, session) {
    
    output$downloadRDS <- downloadHandler(
        filename = "Imputation_Example_10X_Dataset.RDS",
        content = function(file) {
            file.copy(example, file)
        },
        contentType = "application/RDS"
    )
    
    output$dl <- downloadHandler(
        filename <- paste("HTCA_Post_Imputation_Seurat_Object_",sys_time,".RDS", sep = ""),
        content <- function(file) {
            x <- Imputation()$x
            saveRDS(x,file)
        },
        contentType = "application/RDS"
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
        showModal(modalDialog("Imputation in process. For multiple samples data, imputation is carried out separately for each sample (may take some time depending on the number of samples)..Proceed to tab 'Plot' after the imputation.", footer=NULL))
        x <- readRDS(inFile())
        DefaultAssay(x) <- "RNA"
        toremove <- colSums(x@assays$RNA@counts)
        toremove <- toremove[which(toremove == 0)]
        if(length(toremove) > 0){
          x <- subset(x, cells = row.names(x@meta.data)[which(!row.names(x@meta.data) %in% toremove)])
        }
        
        cchoices <- NULL
        x <- SplitObject(x, split.by = "orig.ident")
        for(i in 1:length(x)){
          current <- x[[i]]
          current <- NormalizeData(current)
          current_norm <- t(as.matrix(current@assays$RNA@data))
          k_choice <- choose_k(current_norm)
          cchoices[[length(cchoices)+1]] <- k_choice
          names(cchoices)[length(cchoices)] <- unique(current$orig.ident)
          current_norm <- alra(current_norm,k=k_choice$k)[[3]]
          current_norm <- t(current_norm)
          colnames(current_norm) <- colnames(current)
          current@assays$RNA@data <- current_norm
          x[[i]] <- current
        }
        
        if(length(x) > 1){
          x <- merge(x[[1]], x[[2:length(x)]])
        }else{
          x <- x[[1]]
        }
        out <- NULL
        out$x <- x
        out$cchoices <- cchoices
        removeModal()
        
        updateSelectInput(session, inputId = 'sample', label = 'Choose a Sample to View', choices = names(cchoices), selected = names(cchoices)[1])
        updateSliderInput(session, inputId = 'pointsize', label = 'Size of Points', min = 0.1, max = 10, value = 1.5, step = 0.1)
        updateSliderInput(session, inputId = 'labelsize', label = 'Size of Labels', min = 10, max = 50, value = 20, step = 0.1)
        updateSliderInput(session, inputId = 'width', label = 'Plot Width', min = 400, max = 8000, value = 900, step = 10)
        updateSliderInput(session, inputId = 'height', label = 'Plot Height', min = 200, max = 8000, value = 400, step = 10)
        
        return(out)
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
    
    ScatteredPlotObject <- reactive({
        k_choice <- data()$cchoices[[input$sample]]
        df <- data.frame(x=1:100,y=k_choice$d)
        p1<-ggplot(df,aes(x=x,y=y),) + geom_point(size=input$pointsize) + geom_point(data = df[which(df$x == k_choice$k),], aes(x,y), color = "darkred", size = input$pointsize+1) + geom_line(size=0.5)+ geom_vline(xintercept=k_choice$k)  +theme_classic(base_size = input$labelsize) + theme( axis.title.x=element_blank() ) + scale_x_continuous(breaks=seq(10,100,10)) + ylab('s_i') + ggtitle("Singular values")
        df <- data.frame(x=2:100,y=diff(k_choice$d))[3:99,]
        p2<-ggplot(df,aes(x=x,y=y),) + geom_point(size=input$pointsize)  + geom_point(data = df[which(df$x == k_choice$k),], aes(x,y), color = "darkred", size = input$pointsize+1) + geom_line(size=0.5)+ geom_vline(xintercept=k_choice$k+1) + theme_classic(base_size = input$labelsize) + theme(axis.title.x=element_blank() ) + scale_x_continuous(breaks=seq(10,100,10)) + ylab('s_{i} - s_{i-1}') + ggtitle('Singular value spacings')
        p <- p1+p2+plot_annotation(title = paste('Sample:',input$sample,', Choosen K=',k_choice$k, sep =""), theme = theme(plot.title = element_text(size = input$labelsize, face = "bold", hjust = 0.5)))
    })
    
    output$ScatteredPlot <- renderPlot({
      print(ScatteredPlotObject())
    }, width = w, height = h)
    
    output$dlpng <- downloadHandler(
      filename = function() { paste("HTCA_Imputation_",sys_time,".png", sep = "")},
      content <- function(file) {
        png(file, width=w()*5, height=h()*5, res = 300)
        print(ScatteredPlotObject())
        dev.off()
      }
    )
    
    output$dlpdf <- downloadHandler(
      filename = function() { paste("HTCA_Imputation_",sys_time,".pdf", sep = "")},
      content <- function(file) {
        pdf(file, width=w()/72, height=h()/72)
        print(ScatteredPlotObject())
        dev.off()
      }
    )
    
    output$dlRDS <- downloadHandler(
      filename <- paste("HTCA_Post_Imputation_Seurat_Object_",sys_time,".RDS", sep = ""),
      content <- function(file) {
        x <- data()$x
        saveRDS(x,file)
      },
      contentType = "application/RDS"
    )
})



