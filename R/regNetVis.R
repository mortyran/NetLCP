#' @title Regulatory Network Visual Function
#' @description Regulatory network visualization
#'
#' @param regData standard out of the function binaryRegulation or multieleRegulation.
#' @param filterDegree an integer for filtering the nodes.
#' @param selectCREs a vector of elements or CREs.
#' @param netLayout the layout of network
#'
#' @examples
#' regNetVis(regData = NULL, filterDegree = 20, selectCREs = "9978", netLayout = "layout_in_circle")
#' @export
regNetVis = function(regData = NULL, filterDegree = 40, selectCREs = NULL, netLayout = "layout_nicely"){

  if(is.null(regData)){
    return("Please input the valid regData......")
  }
  nodeInfo = data.frame(table(c(regData$node1, regData$node2)))
  colnames(nodeInfo) = c("label", "value")
  nodeInfo$id = 1:nrow(nodeInfo)
  nodeInfo$group = sapply(nodeInfo$label, function(x){
    if(grepl("MIMAT", x)){
      return("miRNA")
    }else if(grepl("ENSG", x)){
      return("lncRNA")
    }else if(grepl("hsa_circ", x)){
      return("circRNA")
    }else if(grepl("R-HSA", x)){
      return("Reactome")
    }else if(grepl("WP", x)){
      return("Wikipathway")
    }else if(grepl("hsa", x)){
      return("KEGG")
    }else{
      return("mRNA")
    }
  })
  nodeInfo$title = nodeInfo$group
  edgeInfo = regData
  colnames(edgeInfo) = c("from", "to", "title", "group")
  edgeInfo$from = sapply(edgeInfo$from, function(x){
    return(nodeInfo$id[nodeInfo$label == x])
  })
  edgeInfo$to = sapply(edgeInfo$to, function(x){
    return(nodeInfo$id[nodeInfo$label == x])
  })

  if(is.null(selectCREs)){
    nodeFilter = nodeInfo$id[nodeInfo$value >= filterDegree]
    edgeInfo = edgeInfo[edgeInfo$from %in% nodeFilter & edgeInfo$to %in% nodeFilter,]
    nodeInfo = nodeInfo[nodeInfo$id %in% unique(c(edgeInfo$from, edgeInfo$to)),]
    if(nrow(nodeInfo) > 0){
      return(visNetwork::visNetwork(nodes = nodeInfo, edges = edgeInfo) %>% visNetwork::visIgraphLayout(layout = netLayout) %>% visNetwork::visOptions(manipulation = TRUE, highlightNearest = TRUE,selectedBy = list(variable = "group", highlight = TRUE), nodesIdSelection = TRUE))
    }

  }else{

    if(any(grepl("_", selectCREs))){
      elements_all = c()
      for(i in selectCREs){
        elements_all = c(unlist(strsplit(i, "_")), elements_all)
      }
      selectCREs = unique(elements_all)
    }

    if(length(selectCREs) == 1){
      nodeFilter = nodeInfo$id[nodeInfo$label %in% selectCREs]
      edgeInfo = edgeInfo[edgeInfo$from %in% nodeFilter | edgeInfo$to %in% nodeFilter,]
      nodeInfo = nodeInfo[nodeInfo$id %in% unique(c(edgeInfo$from, edgeInfo$to)),]
    }else{
      nodeFilter = nodeInfo$id[nodeInfo$label %in% selectCREs]
      edgeInfo = edgeInfo[edgeInfo$from %in% nodeFilter & edgeInfo$to %in% nodeFilter,]
      nodeInfo = nodeInfo[nodeInfo$id %in% unique(c(edgeInfo$from, edgeInfo$to)),]
    }
    if(nrow(nodeInfo) > 0){
      return(visNetwork::visNetwork(nodes = nodeInfo, edges = edgeInfo) %>% visNetwork::visIgraphLayout(layout = netLayout) %>% visNetwork::visOptions(manipulation = TRUE, highlightNearest = TRUE,selectedBy = list(variable = "group", highlight = TRUE), nodesIdSelection = TRUE))
    }
  }
}
