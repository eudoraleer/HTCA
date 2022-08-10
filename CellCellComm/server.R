library(BiocManager)
options(repos = BiocManager::repositories())

library("shiny")
library("shinyscreenshot")
library("ggplot2")
library("ggrepel")
library("Seurat")
library("shinyjs")
library("randomcoloR")
library("Connectome")
library("entropy")
library("igraph")
library("qgraph")
library("dplyr")
library("DT")
library("liana")

methods <- readRDS("DB/CellCellComm_Methods.RDS")
resources <- readRDS("DB/CellCellComm_Resources.RDS")

options(shiny.maxRequestSize=3000*1024^2)

example <- "DB/CellCellComm_Example_10X_Dataset.RDS"
sys_time <- format(Sys.time(), "%Y%m%d%H%M%S")

shinyServer(function(input, output, session) {
    
    output$downloadRDS <- downloadHandler(
        filename = "CellCellComm_Example_10X_Dataset.RDS",
        content = function(file) {
            file.copy(example, file)
        },
        contentType = "application/RDS"
    )
    
    output$dlresult <- downloadHandler(
        filename <- paste("HTCA_CellCellComm_Result_",sys_time,".RDS", sep = ""),
        content <- function(file) {
            x <- cellcell$allresults
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
        showModal(modalDialog("Uploading data..Proceed to tab 'Plot' once the upload is done.", footer=NULL))
        x <- readRDS(inFile())
        DefaultAssay(x) <- "RNA"
        groups <- colnames(x@meta.data)[grep("Group",colnames(x@meta.data), ignore.case = T)]
        if(length(groups) == 0){
          groups <- "None"
        }else{
          groups <- sort(unique(x@meta.data[,groups]))
          if(length(groups) == 1){
            groups <- "None"
          }else{
            groups <- c("None", groups)
          }
        }
        celllabels <- colnames(x@meta.data)[grep("seurat_clusters|cell_type",colnames(x@meta.data), ignore.case = T)]
        cviewmethods <- input$methods
        cviewdbs <- input$dbs
        removeModal()
        updateCheckboxGroupInput(session, inputId = 'methods', label = 'Select method(s)', choices = methods, selected = methods[1])
        updateCheckboxGroupInput(session, inputId = 'dbs', label = 'Select knowledge resource(s)', choices = resources, selected = resources[1])
        updateSelectInput(session, inputId = 'celllabels', label = 'Select the cell labels to perform the analysis', choices = celllabels, selected = celllabels[length(celllabels)])
        updateSelectInput(session, inputId = 'groups', label = 'Select a group to perform the analysis separately (for data with multiple groups, for example, diseased and control groups. Ignored for single-group data)', choices = groups, selected = groups[1])
        updateCheckboxInput(session, inputId = 'set3d', label = '3D Vertex',value = TRUE)
        updateSliderInput(session, inputId = 'labelsize', label = 'Size of Labels', min = 0.1, max = 20, value = 2, step = 0.1)
        updateSliderInput(session, inputId = 'labeldist', label = 'Distance of Labels', min = 0.1, max = 10, value = 3, step = 0.1)
        updateSliderInput(session, inputId = 'arrowthickness', label = 'Thickness of Arrows', min = 0.1, max = 10, value = 0.5, step = 0.1)
        updateSliderInput(session, inputId = 'vertexsize', label = 'Size of Vertex', min = 10, max = 50, value = 14, step = 1)
        updateSliderInput(session, inputId = 'curvity', label = 'Curvity of Arrows', min = 0.1, max = 2, value = 0.8, step = 0.01)
        updateSliderInput(session, inputId = "colindex", label = "Choose A Color Scheme Index", min = 1, max = 100, value = 1, step = 1)
        updateSliderInput(session, inputId = 'width', label = 'Plot Width', min = 400, max = 8000, value = 900, step = 10)
        updateSliderInput(session, inputId = 'height', label = 'Plot Height', min = 400, max = 8000, value = 760, step = 10)
        updateSelectInput(session, inputId = 'viewmethod', label = 'Method(s)', choices = cviewmethods, selected = cviewmethods[1])
        updateSelectInput(session, inputId = 'viewdb', label = 'Resource(s)', choices = cviewdbs, selected = cviewdbs[1])
        
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
    
    cellcell <- reactiveValues(allresults = NULL)
    observeEvent(input$go, {
        showModal(modalDialog("Cell-cell communication analysis is in process (may take some time depending on the sample size and the number of methods/resources/celltypes to run)...", footer=NULL, easyClose = F))
        x <- data()
        if(input$groups != "None"){
          cx <- subset(x, subset = "Group" == input$groups)
          showModal(modalDialog("Cell-cell communication analysis is in process (may take some time depending on the sample size and the number of methods/resources/celltypes to run)...", footer=NULL, easyClose = F))
          cellresult <- liana_wrap(cx, method = input$methods, resource = input$dbs, idents_col = input$celllabels)
          removeModal()
        }else{
          showModal(modalDialog("Cell-cell communication analysis is in process (may take some time depending on the sample size and the number of methods/resources/celltypes to run)...", footer=NULL, easyClose = F))
          cellresult <- liana_wrap(x, method = input$methods, resource = input$dbs, idents_col = input$celllabels)
          removeModal()
          print(paste("Done with liana_wrap"))
        }

        cellcell$allresults <- cellresult
        print("Completed Cell-Cell Communication Analysis!")
        removeModal()
    })
    
    GraphObject <- function(){
      cellresult <- cellcell$allresults
      if(class(cellresult)[1] == "list" & length(cellresult) != 1){
        if(length(cellresult[[input$viewdb]]) > 0){
          current <- cellresult[[input$viewdb]]
          if(length(current[[input$viewmethod]]) > 0){
            cresult <- current[[input$viewmethod]]
          }else{
            cresult <- current
          }
        }else{
          current <- cellresult[[input$viewmethod]]
          if(length(current[[input$viewdb]]) > 0){
            cresult <- current[[input$viewdb]]
          }else{
            cresult <- current
          }
        }
      }else{
        cresult <- cellresult
      }
      mynet <- data.frame(table(cresult[,c("source","target")]))
      mynet <- mynet[which(mynet$Freq > 0),]
      net<- graph_from_data_frame(mynet)
      deg <- igraph::degree(net, mode="all")
      
      set.seed(input$colindex)
      ccols <- distinctColorPalette(length(unique(c(mynet$source, mynet$target))))
      names(ccols) <- unique(c(mynet$source, mynet$target))
      karate_groups <- cluster_optimal(net)
      coords <- layout_in_circle(net, order = order(membership(karate_groups)))
      E(net)$width  <- E(net)$Freq*input$arrowthickness
      for (i in 1: length(unique(mynet$source))){
        E(net)[purrr::map(unique(mynet$source),function(x) {
          get.edge.ids(net,vp = c(unique(mynet$source)[i],x))
        })%>% unlist()]$color <- ccols[unique(mynet$source)[i]]
      }
      plotx <- NULL
      plotx$net <- net
      plotx$ccols <- ccols
      plotx$coords <- coords
      plotx$cresult <- cresult[,which(!colnames(cresult) %in% c("ligand.complex","receptor.complex"))]
      return(plotx)
      removeModal()
    }
    
    output$Graph <- renderPlot({
      req(input$go)
      plot(GraphObject()$net, edge.arrow.size=.4, 
           vertex.size=input$vertexsize,
           edge.curved=input$curvity,
           vertex.color=GraphObject()$ccols,
           vertex.shape=ifelse(input$set3d == TRUE, "sphere", "circle"),
           vertex.frame.color="#555555",
           vertex.label.color=GraphObject()$ccols,
           layout = GraphObject()$coords,
           vertex.label.family="Helvetica",
           vertex.label.degree=45,
           vertex.label.dist=input$labeldist,
           vertex.label.cex=input$labelsize)
      title(paste(input$viewmethod, input$viewdb, sep = " - "),cex.main=input$labelsize,col.main="black")
    }, width = w, height = h)
    
    output$dlpng <- downloadHandler(
      filename = function() { paste("HTCA_CellCellCommunication_",sys_time,".png", sep = "")},
      content <- function(file) {
        png(file, width=w()*5, height=h()*5, res = 300)
        plot(GraphObject()$net, edge.arrow.size=.4, 
             vertex.size=input$vertexsize,
             edge.curved=input$curvity,
             vertex.color=GraphObject()$ccols,
             vertex.shape="sphere",
             vertex.frame.color="#555555",
             vertex.label.color=GraphObject()$ccols,
             layout = GraphObject()$coords,
             vertex.label.family="Helvetica",
             vertex.label.degree=45,
             vertex.label.dist=3,
             vertex.label.cex=input$labelsize)
        title(paste(input$viewmethod, input$viewdb, sep = " - "),cex.main=input$labelsize,col.main="black")
        dev.off()
      }
    )
    
    output$dlpdf <- downloadHandler(
      filename = function() { paste("HTCA_CellCellCommunication_",sys_time,".pdf", sep = "")},
      content <- function(file) {
        pdf(file, width=w()/72, height=h()/72)
        plot(GraphObject()$net, edge.arrow.size=.4, 
             vertex.size=input$vertexsize,
             edge.curved=input$curvity,
             vertex.color=GraphObject()$ccols,
             vertex.shape="sphere",
             vertex.frame.color="#555555",
             vertex.label.color=GraphObject()$ccols,
             layout = GraphObject()$coords,
             vertex.label.family="Helvetica",
             vertex.label.degree=45,
             vertex.label.dist=3,
             vertex.label.cex=input$labelsize)
        title(paste(input$viewmethod, input$viewdb, sep = " - "),cex.main=input$labelsize,col.main="black")
        dev.off()
      }
    )
    
    selectTable <- reactive({
      req(input$go)
      return(GraphObject()$cresult)
    })
    
    output$table <- renderDT(selectTable(),
                             filter = "top",
                             style="bootstrap",
                             rownames = F,
                             options = list(
                               pageLength = 10)
    )
})



