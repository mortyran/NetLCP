#' @title biological elements of CREs statistics
#' @description biological elements of CREs statistics
#'
#' @param regData the standard out of the function binaryRegulation or multieleRegulation.
#' @param filterDegree an integer for filtering the nodes.
#' @param selectCREs a vector of elements or CREs.
#'
#' @examples
#' regStat(regData = NULL, filterDegree = 50)
#' @export
regStat = function(regData = NULL, filterDegree = 40, selectCREs = NULL){

  if(is.null(regData)){
    return("Please input the valid regData......")
  }

  if(is.null(selectCREs)){
    nodeInfo = data.frame(table(c(regData$node1, regData$node2)))
    nodeFilter = nodeInfo$Var1[nodeInfo$Freq >= filterDegree]
    regData = regData[regData$node1 %in% nodeFilter | regData$node2 %in% nodeFilter,]
  }else{

    if(any(grepl("_", selectCREs))){
      elements_all = c()
      for(i in selectCREs){
        elements_all = c(unlist(strsplit(i, "_")), elements_all)
      }
      selectCREs = unique(elements_all)
    }

    if(length(selectCREs) == 1){
      regData = regData[regData$node1 %in% selectCREs | regData$node2 %in% selectCREs,]
    }else{
      regData = regData[regData$node1 %in% selectCREs & regData$node2 %in% selectCREs,]
    }
  }

  A = data.frame(table(c(regData$node1, regData$node2)))
  A$type = sapply(A$Var1, function(x){
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
  A$color = sapply(A$Var1, function(x){
    if(grepl("MIMAT", x)){
      return("#de425b")
    }else if(grepl("ENSG", x)){
      return("#eb8c66")
    }else if(grepl("hsa_circ", x)){
      return("#f0c897")
    }else if(grepl("R-HSA", x)){
      return("#cdd59e")
    }else if(grepl("WP", x)){
      return("#91b164")
    }else if(grepl("hsa", x)){
      return("#488f31")
    }else{
      return("#80a856")
    }
  }
  )
  return(plotly::plot_ly(x = A$Var1, y = A$Freq, marker = list(color = A$color), type = "bar") %>% plotly::layout(title = "Element Degree", xaxis = list(title = "Element"), yaxis = list(title = "Degree")))
}
