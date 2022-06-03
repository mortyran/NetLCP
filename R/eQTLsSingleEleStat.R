#' @title single element eQTLs statistics
#' @description single element eQTLs statistics
#'
#' @param regData the standard out of the function binaryRegulationExtract or multieleRegulation.
#' @param eQTLsData the standard out of the function eQTLsDetection.
#' @param filterDegree an integer for filtering the nodes.
#' @param selectCREs select a or a set of nodes to perform statistics.
#'
#' @examples
#' eQTLsSingleEleStat(regData = NULL, eQTLsData = NULL, filterDegree = 100, selectCREs = NULL)
#' @export
eQTLsSingleEleStat = function(regData = NULL, eQTLsData = NULL, filterDegree = 40, selectCREs = NULL){

  if(is.null(regData) | is.null(eQTLsData)){
    return("Please input the valid regData and eQTLsData......")
  }
  nodeInfo = data.frame(table(c(regData$node1, regData$node2)))

  if(is.null(selectCREs)){
    nodeFilter = nodeInfo$Var1[nodeInfo$Freq >= filterDegree]
    eqtls = eQTLsData[eQTLsData$Elements %in% nodeFilter,]

    if(nrow(eqtls) > 0){
      eqtls_stat = table(eqtls$Elements, eqtls$RegType)
      eqtls_stat = data.frame(Element=rownames(eqtls_stat), cis_eQTLs=eqtls_stat[,1], trans_eQTLs=eqtls_stat[,2])

      return(plotly::plot_ly(eqtls_stat, x = ~Element, y = ~cis_eQTLs, type = 'bar', name = 'cis-eQTLs') %>% plotly::add_trace(y = ~trans_eQTLs, name = 'trans-eQTLs') %>% plotly::layout(yaxis = list(title = 'eQTLs count'), barmode = 'stack'))
    }

  }else{

    if(any(grepl("_", selectCREs))){
      elements_all = c()
      for(i in selectCREs){
        elements_all = c(unlist(strsplit(i, "_")), elements_all)
      }
      selectCREs = unique(elements_all)
    }

    eqtls = eQTLsData[eQTLsData$Elements %in% selectCREs,]
    if(nrow(eqtls) > 0){
      eqtls_stat = table(eqtls$Elements, eqtls$RegType)
      if(ncol(eqtls_stat) == 1){
        return(plotly::plot_ly(x = ~rownames(eqtls_stat), y = ~eqtls_stat[,1], type = 'bar', name = colnames(eqtls_stat)[1]) %>% plotly::layout(xaxis = list(title = 'Element'), yaxis = list(title = 'eQTLs count'), barmode = 'stack'))
      }else{
        return(plotly::plot_ly(x = ~rownames(eqtls_stat), y = ~eqtls_stat[,1], type = 'bar', name = colnames(eqtls_stat)[1]) %>% plotly::add_trace(y = ~eqtls_stat[,2], name = colnames(eqtls_stat)[2]) %>% plotly::layout(xaxis = list(title = 'Element'), yaxis = list(title = 'eQTLs count'), barmode = 'stack'))
      }
    }
  }
}
