#' @title eQTLs Regulatory Network Visual Function
#' @description eQTLs Regulatory network visualization
#'
#' @param regData the standard out of the function binaryRegulationExtract or multieleRegulation.
#' @param eQTLsData the standard out of the function eQTLsDetection.
#' @param filterDegree an integer for filtering the nodes.
#' @param selectCREs a vector of elements or CREs.
#' @param netLayout layout of network, including "layout_in_circle" or "layout_nicely".
#'
#' @examples
#' eQTLsNetVis(regData = NULL, eQTLsData = NULL, filterDegree = 100,selectCREs = "MIMAT0000423",netLayout = "layout_in_circle")
#' @export
eQTLsNetVis = function(regData = NULL, eQTLsData = NULL, filterDegree = 40, selectCREs = NULL, netLayout = "layout_in_circle"){

  if(is.null(regData) | is.null(eQTLsData)){
    return("Please input the valid regData and eQTLsData......")
  }

  colnames(regData) = c("node1", "node2", "title", "group")
  if(is.null(selectCREs)){
    nodeInfo_temp = data.frame(table(c(regData$node1, regData$node2)))
    nodeFilter = nodeInfo_temp$Var1[nodeInfo_temp$Freq >= filterDegree]
    regData = regData[regData$node1 %in% nodeFilter & regData$node2 %in% nodeFilter,]
  }else{

    if(any(grepl("_", selectCREs))){
      elements_all = c()
      for(i in selectCREs){
        elements_all = c(unlist(strsplit(i, "_")), elements_all)
      }
      selectCREs = unique(elements_all)
    }

    if(length(selectCREs) == 1){
      regData = regData[regData$node1 %in% selectCREs | regData$node2 %in% selectCREs, ]
    }else{
      regData = regData[regData$node1 %in% selectCREs & regData$node2 %in% selectCREs, ]
    }
  }

  if(nrow(regData) == 0){
    return()
  }

  eqtls_edgeInfo = eQTLsData[,c(1,2,3,5)]
  colnames(eqtls_edgeInfo) = c("node1", "node2", "title", "group")
  eqtls_edgeInfo = eqtls_edgeInfo[eqtls_edgeInfo$node1 %in% unique(c(regData$node1, regData$node2)) | eqtls_edgeInfo$node2 %in% unique(c(regData$node1, regData$node2)), ]
  edgeInfo = rbind(regData, eqtls_edgeInfo) %>% distinct()

  nodeInfo = data.frame(table(c(edgeInfo$node1, edgeInfo$node2)))
  colnames(nodeInfo) = c("label", "value")
  eqtl_label = eQTLsData[,c(2, 5)] %>% distinct()
  colnames(eqtl_label) = c("label", "group")
  nodeInfo = nodeInfo %>% left_join(eqtl_label, by = "label")
  nodeInfo$group[is.na(nodeInfo$group)] = sapply(nodeInfo$label[is.na(nodeInfo$group)], function(x){
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
  nodeInfo = nodeInfo %>% dplyr::distinct()
  nodeInfo$id = 1:nrow(nodeInfo)
  nodeInfo$title = nodeInfo$group

  nodeInfo_1 = nodeInfo[,c(1, 4)] %>% dplyr::rename(node1=label,from=id)
  nodeInfo_2 = nodeInfo[,c(1, 4)] %>% dplyr::rename(node2=label,to=id)

  edgeInfo = left_join(edgeInfo, nodeInfo_1, by = "node1")
  edgeInfo = left_join(edgeInfo, nodeInfo_2, by = "node2")
  edgeInfo = edgeInfo[,-c(1,2)]

  return(visNetwork::visNetwork(nodes = nodeInfo, edges = edgeInfo) %>% visNetwork::visIgraphLayout(layout = netLayout) %>% visNetwork::visOptions(manipulation = TRUE, highlightNearest = TRUE,selectedBy = list(variable = "group", highlight = TRUE), nodesIdSelection = TRUE))
}
