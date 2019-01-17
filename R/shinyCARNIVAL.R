#'\code{shinyCARNIVAL}
#'
#' This function takes results from ILP optimisation and display the results on an R-Shiny application
#'
#' @param dir_name The directory name that stores CARNIVAL results
#' @param nrIntActInc Maximum number of included interactions in the graph
#' @param thresNodeAll Minimal node activity to be shown on the graph
#' @param thresEdgeAll Minimal edge presence percentage to be shown on the graph
#' @param graphLayout Defined igraph-style layout: 1=tree, 2=circle, 3=grid, 4=sphere, other=default[nicely]
#' @param UP2GS Convert UniprotID to Gene Symbols to facilitate the interpretation
#' @param inverseCR Set parameter to TRUE to remove the perturbation node from the network upon display
#'
#' @import shiny
#' @import visNetwork
#' @import tidyverse
#'
#' @return Launched Rshiny app with plotted results
#' @export

shinyCARNIVAL <- function(dir_name, nrIntActInc=1000, thresNodeAll=NULL, thresEdgeAll=NULL, graphLayout=0, UP2GS=F, inverseCR=F){
  
  library(shiny)
  # library(rcytoscapejs)
  library(visNetwork)
  # library(threejs)
  library(tidyverse)
  
  if (graphLayout==1) {graphLayoutName <- "layout_as_tree"
  } else if (graphLayout==2) {graphLayoutName <- "layout_in_circle"
  } else if (graphLayout==3) {graphLayoutName <- "layout_on_grid"
  } else if (graphLayout==4) {graphLayoutName <- "layout_on_sphere"
  } else {graphLayoutName <- "layout_nicely"}
  
  ui = shinyUI(
    fluidPage(    
      titlePanel("CARNIVAL results on R-Shiny"),
      mainPanel(
        # rcytoscapejsOutput("g3plot") 
        visNetworkOutput("g3plot") 
      )
    )
  )
  
  server = function(input, output) {
    
    nodeDataReactive <- reactiveValues()
    
    observe({ 
      net = read.table(file = paste0(dir_name,"/weightedModel_1.txt"),header = T,sep = "\t",stringsAsFactors = F)
      nodeAttr = read.table(file = paste0(dir_name,"/nodesAttributes_1.txt"),header = T,sep = "\t",stringsAsFactors = F)
      
      if (inverseCR) {
        net <- net[-which(net[,1]=="Perturbation"),]
      }
      
      if (UP2GS) {
        # Select mapping file
        IDmap <- read.table(file = system.file("HUMAN_9606_idmapping_onlyGeneName.dat",package="CARNIVAL"),header = F,sep = "\t",stringsAsFactors = F)
        
        print("Mapping Uniprot IDs to Gene Symbols....")
        
        # Map each element in 'res' -> always take the first ID if there are many hit entries and collect unmapped nodes
        Unmapped <- NULL
        
        # Common SIF
        for (counter in 1:length(net[,1])) {
          if (length(IDmap[which(IDmap[,1] == net[counter,1]),3])>0) {
            net[counter,1] <- IDmap[which(IDmap[,1] == net[counter,1]),3][1]
          } else {
            Unmapped <- c(Unmapped,net[counter,1])
          }
          if (length(IDmap[which(IDmap[,1] == net[counter,3]),3])>0) {
            net[counter,3] <- IDmap[which(IDmap[,1] == net[counter,3]),3][1]
          } else {
            Unmapped <- c(Unmapped,net[counter,3])
          }
        }
        
        # Common node activity
        for (counter in 1:length(nodeAttr[,1])) {
          if (length(IDmap[which(IDmap[,1] == nodeAttr[counter,1]),3])>0) {
            nodeAttr[counter,1] <- IDmap[which(IDmap[,1] == nodeAttr[counter,1]),3][1]
          } else {
            Unmapped <- c(Unmapped,nodeAttr[counter,1])
          }
        }
        
        print("- Done -")
        
      }
      
      # Remove edges based on edge threshold (if specified)
      if (!is.null(thresEdgeAll)) {
        if (length(which(as.numeric(net[,4])<thresEdgeAll))!=0) {
          net <- net[-which(as.numeric(net[,4])<thresEdgeAll),]
        }
      }
      
      # Remove nodes (and corresponding edges) based on node threshold (if specified)
      if (!is.null(thresNodeAll)) {
        nodesBelowThres <- nodeAttr[which(abs(as.numeric(nodeAttr[,5]))<thresNodeAll & abs(as.numeric(nodeAttr[,5]))!=0),1]
        if (length(nodesBelowThres)>0) {
          IdxEdge2rm <- NULL
          for (counter in 1:length(nodesBelowThres)) {
            if (nodesBelowThres[counter] %in% c(net[,1],net[,3])) {
              IdxEdge2rm <- c(IdxEdge2rm,which(nodesBelowThres[counter] == net[,1]),which(nodesBelowThres[counter] == net[,3]))
            }
          }
          if (!is.null(IdxEdge2rm)) {
            net <- net[-unique(IdxEdge2rm),]
          }
        }
      }
      
      # Identify input TFs
      inputTFs <- nodeAttr[which(nodeAttr[,6]=="P"),1]
      
      # ======================= #
      # "rcytoscapejs" pipeline #
      # ======================= #
      
      # # Extract and map network information
      # netInputs <- setdiff(net[,1],net[,3])
      # 
      # network <- as.data.frame(cbind(net[,1],net[,3],net[,2],net[,4]),stringsAsFactors=F); colnames(network) <- c("Source","Target","Type","Weight")
      # 
      # allNodes <- unique(c(network[,1],network[,2]))
      # nodeAct <- rep(NA,length(allNodes))
      # for (counter in 1:length(nodeAct)) {
      #   nodeAct[counter] <- nodeAttr$AvgAct[which(allNodes[counter]==nodeAttr$Node)]
      # }
      # 
      # nodeData <- data.frame(id=unique(c(network[,1],network[,2])),name=unique(c(network[,1],network[,2])),stringsAsFactors = F)
      # nodeColorDefault <- "#DBDBDB"
      # nodeData$color <- rep(nodeColorDefault,nrow(nodeData))
      # # Blue set
      # nodeData$color[which(nodeAct>90)] <- "#8488FF"
      # nodeData$color[which(nodeAct>70 & nodeAct<90)] <- "#A1A4FF"
      # nodeData$color[which(nodeAct>50 & nodeAct<70)] <- "#BDBFFF"
      # nodeData$color[which(nodeAct>30 & nodeAct<50)] <- "#DADBFF"
      # nodeData$color[which(nodeAct>10 & nodeAct<30)] <- "#F6F6FF"
      # # # Green set
      # # nodeData$color[which(nodeAct>90)] <- "#4EE968"
      # # nodeData$color[which(nodeAct>70 & nodeAct<90)] <- "#68E97D"
      # # nodeData$color[which(nodeAct>50 & nodeAct<70)] <- "#82E993"
      # # nodeData$color[which(nodeAct>30 & nodeAct<50)] <- "#9BE9A8"
      # # nodeData$color[which(nodeAct>10 & nodeAct<30)] <- "#CFE9D4"
      # # Red set
      # nodeData$color[which(nodeAct< -90)] <- "#FF7070"
      # nodeData$color[which(nodeAct< -70 & nodeAct> -90)] <- "#FF8E8E"
      # nodeData$color[which(nodeAct< -50 & nodeAct> -70)] <- "#FFA5A5"
      # nodeData$color[which(nodeAct< -30 & nodeAct> -50)] <- "#FFBCBC"
      # nodeData$color[which(nodeAct< -10 & nodeAct> -30)] <- "#FFECEC"
      # 
      # nodeShapeDefault <- "ellipse"
      # nodeData$shape <- rep(nodeShapeDefault, nrow(nodeData))
      # nodeData$shape[which(nodeData[,1] %in% inputTFs)] <- "triangle"
      # nodeData$shape[which(nodeData[,1] %in% netInputs)] <- "diamond"
      # 
      # edgeData <- data.frame(source=network[,1], target=network[,2], stringsAsFactors=FALSE)
      # edgeColorDefault <- "#DBDBDB"
      # edgeData$color <- rep(edgeColorDefault,nrow(edgeData))
      # 
      # network[,4] <- as.numeric(network[,4])
      # 
      # # Red set    
      # edgeData$color[which(network[,3]==-1 & network[,4]>90)] <- "#FF5050"
      # edgeData$color[which(network[,3]==-1 & network[,4]>70 & network[,4]<90)] <- "#FF8E8E"
      # edgeData$color[which(network[,3]==-1 & network[,4]>50 & network[,4]<70)] <- "#FFA5A5"
      # edgeData$color[which(network[,3]==-1 & network[,4]>30 & network[,4]<50)] <- "#FFBCBC"
      # edgeData$color[which(network[,3]==-1 & network[,4]>10 & network[,4]<30)] <- "#FFECEC"
      # 
      # # Blue Set
      # edgeData$color[which(network[,3]==1 & network[,4]>90)] <- "#5650FF"
      # edgeData$color[which(network[,3]==1 & network[,4]>70 & network[,4]<90)] <- "#A1A4FF"
      # edgeData$color[which(network[,3]==1 & network[,4]>50 & network[,4]<70)] <- "#BDBFFF"
      # edgeData$color[which(network[,3]==1 & network[,4]>30 & network[,4]<50)] <- "#DADBFF"
      # edgeData$color[which(network[,3]==1 & network[,4]>10 & network[,4]<30)] <- "#F6F6FF"
      # 
      # nodeDataReactive$nodeData <- nodeData
      # nodeDataReactive$edgeData <- edgeData
      # nodeDataReactive$network <- network
      
      # =================== #
      # visNetwork pipeline #
      # =================== #
      
      # Initialise matrices
      AllNodes <- sort(unique(c(net[,1],net[,3])))
      nodesVis <- as.data.frame(matrix(NA,length(AllNodes),4))
      edgesVis <- as.data.frame(matrix(NA,nrow(net),4))
      colnames(nodesVis) <- c("id","label","level","value")
      colnames(edgesVis) <- c("from","to","sign","value")
      
      # Mapping
      # Node
      nodesVis[,1] <- as.numeric(1:length(AllNodes))
      nodesVis[,2] <- AllNodes
      
      nodesVis[,3] <- 2
      netInputs <- setdiff(net[,1],net[,3])
      for (counter in 1:length(netInputs)) {
        nodesVis[which(netInputs[counter]==AllNodes),3] <- 1
      }
      inputTFs <- nodeAttr[which(nodeAttr[,6]=="P"),1]
      for (counter in 1:length(inputTFs)) {
        nodesVis[which(inputTFs[counter]==AllNodes),3] <- 3
      }
      
      for (counter in 1:length(AllNodes)) {
        nodesVis[counter,4] <- nodeAttr$AvgAct[which(nodesVis[counter,2]==nodeAttr$Node)]
      }
      
      nodesVis$color.background <- rep(NA,nrow(nodesVis))
      for (counter in 1:nrow(nodesVis)) {
        nodesVis$color.background[counter] <- 
          if (nodesVis$val[counter]>0) {"lightgreen"} 
        else if (nodesVis$val[counter]<0) {"tomato"}
        else if (nodesVis$val[counter]==0) {"grey"}
      }
      nodesVis$color.highlight.background <- rep(NA,nrow(nodesVis))
      
      nodesVis$shape <- rep(NA,nrow(nodesVis))
      for (counter in 1:nrow(nodesVis)) {
        nodesVis$shape[counter] <- 
          if (nodesVis$level[counter]==1) {"diamond"} 
        else if (nodesVis$level[counter]==2) {"circle"}
        else if (nodesVis$level[counter]==3) {"triangle"}
      }
      
      nodesVis$value <- abs(nodesVis$value*10)
      nodesVis <- select(nodesVis,-level)
      
      # Edge
      for (counter in 1:nrow(net)) {
        edgesVis[counter,1] <- nodesVis$id[which(net$Node1[counter]==nodesVis$label)]
        edgesVis[counter,2] <- nodesVis$id[which(net$Node2[counter]==nodesVis$label)]
        edgesVis[counter,3] <- net$Sign[counter]
        edgesVis[counter,4] <- net$Weight[counter]
      }
      
      edgesVis$color <- rep(NA,nrow(edgesVis))
      for (counter in 1:length(edgesVis$color)) {
        edgesVis$color[counter] <- ifelse(test = edgesVis$sign[counter]==1,yes = "green",no = "red")
      }
      
      edgesVis$value <- edgesVis$value/10
      

      nodeDataReactive$nodeDataVis <- nodesVis
      nodeDataReactive$edgeDataVis <- edgesVis
      
    })
    
    # output$g3plot = renderRcytoscapejs({
    #   if (nrow(nodeDataReactive$nodeData) > 0){
    #     networkJS <- createCytoscapeJsNetwork(nodeData = nodeDataReactive$nodeData, edgeData = nodeDataReactive$edgeData,nodeLabelColor = "#000000")
    #     rcytoscapejs(nodeEntries = networkJS$nodes,edgeEntries = networkJS$edges, layout = graphLayout)
    #   } else {
    #     stop("No network to be shown")
    #     }
    # })
    
    output$g3plot = renderVisNetwork({
      if (nrow(nodeDataReactive$nodeDataVis) > 0){
        if (graphLayout==1) {
          visNetwork(nodes = nodeDataReactive$nodeDataVis, edges = nodeDataReactive$edgeDataVis, width="100%") %>% 
          visEdges(arrows = 'to') %>% 
          visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
          visNodes(color = list(border="grey")) %>%
          visIgraphLayout(layout = graphLayoutName,flip.y=F) 
        } else {
          visNetwork(nodes = nodeDataReactive$nodeDataVis, edges = nodeDataReactive$edgeDataVis, width="100%") %>% 
            visEdges(arrows = 'to') %>% 
            visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
            visNodes(color = list(border="grey")) %>%
            visIgraphLayout(layout = graphLayoutName) 
        }
      } else {
        stop("No network to be shown")
      }
    })
    
  }
  
  ShinyInput <- list(ui=ui,server=server)
  
  runApp(ShinyInput)
  
}
