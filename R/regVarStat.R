#' @title Variant 'switches' Statistics
#' @description Variant 'switches' Statistics
#'
#' @param regVar the standard out of the function regVarDetection.
#' @param regulationType a character representing the CREs type, it should be "miRNA-mRNA", "miRNA-mRNA-pathway", "lncRNA-miRNA-mRNA", "circRNA-miRNA-mRNA", "lncRNA-miRNA-mRNA-pathway" or "circRNA-miRNA-mRNA-pathway".
#' @param selectCREs a vector of elements or CREs.
#'
#' @examples
#' regVarStat(regVar = NULL, regulationType = "lncRNA-miRNA-mRNA", selectCREs = sampleData)
#' @export
regVarStat = function(regVar = NULL, regulationType = NULL, selectCREs = NULL){

  if(is.null(regVar) | is.character(regVar)){
    return(NULL)
  }
  if(is.null(regulationType) | !regulationType %in% c("miRNA-mRNA", "miRNA-mRNA-pathway", "lncRNA-miRNA-mRNA", "circRNA-miRNA-mRNA", "lncRNA-miRNA-mRNA-pathway", "circRNA-miRNA-mRNA-pathway")){
    return("Please choose regulation type from miRNA-mRNA, miRNA-mRNA-pathway, lncRNA-miRNA-mRNA, circRNA-miRNA-mRNA, lncRNA-miRNA-mRNA-pathway or circRNA-miRNA-mRNA-pathway")
  }

  if(is.null(selectCREs)){
    return("selectCREs is empty!")
  }

  if(any(grepl("_", selectCREs))){
    elements_all = c()
    for(i in selectCREs){
      elements_all = c(unlist(strsplit(i, "_")), elements_all)
    }
    selectCREs = unique(elements_all)
  }

  miRNAID = c()
  lnc_cirRNAID = c()
  pathwayID = c()
  mRNAID = c()
  for(i in selectCREs){
    if(grepl("MIMAT", i)){
      miRNAID = append(miRNAID, i)
    }else if(grepl("ENSG", i) | grepl("hsa_circ", i)){
      lnc_cirRNAID = append(lnc_cirRNAID, i)
    }else if(grepl("R-HSA", i) | grepl("WP", i) | grepl("hsa", i)){
      pathwayID = append(pathwayID, i)
    }else{
      mRNAID = append(mRNAID, i)
    }
  }

  if(regulationType == "lncRNA-miRNA-mRNA" | regulationType == "circRNA-miRNA-mRNA"){
    regVar_filter = regVar[regVar$ceRNA %in% lnc_cirRNAID & regVar$miRNA %in% miRNAID & regVar$mRNA %in% mRNAID, ]
    regVar_filter$reg = paste(regVar_filter$ceRNA, regVar_filter$miRNA, regVar_filter$mRNA, sep = "_")
  }else if(regulationType == "lncRNA-miRNA-mRNA-pathway" | regulationType == "circRNA-miRNA-mRNA-pathway"){
    regVar_filter = regVar[regVar$ceRNA %in% lnc_cirRNAID & regVar$miRNA %in% miRNAID & regVar$mRNA %in% mRNAID & regVar$pathway %in% pathwayID, ]
    regVar_filter$reg = paste(regVar_filter$ceRNA, regVar_filter$miRNA, regVar_filter$mRNA, regVar_filter$pathway, sep = "_")
  }else if(regulationType == "miRNA-mRNA-pathway"){
    regVar_filter = regVar[regVar$miRNAID %in% miRNAID & regVar$geneID %in% mRNAID & regVar$pathway %in% pathwayID, ]
    regVar_filter$reg = paste(regVar_filter$miRNAID, regVar_filter$geneID, regVar_filter$pathway, sep = "_")
  }else if(regulationType == "miRNA-mRNA"){
    regVar_filter = regVar[regVar$miRNAID %in% miRNAID & regVar$geneID %in% mRNAID, ]
    regVar_filter$reg = paste(regVar_filter$miRNAID, regVar_filter$geneID, sep = "_")
  }

  if(nrow(regVar_filter) == 0){
    return(NULL)
  }
  if(regulationType == "lncRNA-miRNA-mRNA-pathway" | regulationType == "circRNA-miRNA-mRNA-pathway"){
    regVar_filter_lncvar = regVar_filter[,c(4, 6, 17)] %>% distinct()
    regVar_filter_mimvar = regVar_filter[,c(11, 12, 17)] %>% distinct()
    vartable = cbind(table(regVar_filter_lncvar$reg, regVar_filter_lncvar$lncMutType), table(regVar_filter_mimvar$reg, regVar_filter_mimvar$mimMutType))
  }else if(regulationType == "lncRNA-miRNA-mRNA" | regulationType == "circRNA-miRNA-mRNA"){
    regVar_filter_lncvar = regVar_filter[,c(4, 6, 14)] %>% distinct()
    regVar_filter_mimvar = regVar_filter[,c(11, 12, 14)] %>% distinct()
    vartable = cbind(table(regVar_filter_lncvar$reg, regVar_filter_lncvar$lncMutType), table(regVar_filter_mimvar$reg, regVar_filter_mimvar$mimMutType))
  }else{
    vartable = table(regVar_filter$reg, regVar_filter$MutType)
  }
  if(ncol(vartable) == 6){
    return(plotly::plot_ly(x = ~rownames(vartable), y = ~vartable[,1], type = 'bar', name = colnames(vartable)[1]) %>% plotly::add_trace(y = ~vartable[,2], name = colnames(vartable)[2]) %>% add_trace(y = ~vartable[,3], name = colnames(vartable)[3]) %>% plotly::add_trace(y = ~vartable[,4], name = colnames(vartable)[4])%>% plotly::add_trace(y = ~vartable[,5], name = colnames(vartable)[5]) %>% plotly::add_trace(y = ~vartable[,6], name = colnames(vartable)[6]) %>% plotly::layout(xaxis = list(title = 'CREs'), yaxis = list(title = 'Variant switches count'), barmode = 'stack'))
  }else if(ncol(vartable) == 5){
    return(plotly::plot_ly(x = ~rownames(vartable), y = ~vartable[,1], type = 'bar', name = colnames(vartable)[1]) %>% plotly::add_trace(y = ~vartable[,2], name = colnames(vartable)[2]) %>% add_trace(y = ~vartable[,3], name = colnames(vartable)[3]) %>% plotly::add_trace(y = ~vartable[,4], name = colnames(vartable)[4])%>% plotly::add_trace(y = ~vartable[,5], name = colnames(vartable)[5]) %>% plotly::layout(xaxis = list(title = 'CREs'), yaxis = list(title = 'Variant switches count'), barmode = 'stack'))
  }else if(ncol(vartable) == 4){
    return(plotly::plot_ly(x = ~rownames(vartable), y = ~vartable[,1], type = 'bar', name = colnames(vartable)[1]) %>% plotly::add_trace(y = ~vartable[,2], name = colnames(vartable)[2]) %>% add_trace(y = ~vartable[,3], name = colnames(vartable)[3]) %>% plotly::add_trace(y = ~vartable[,4], name = colnames(vartable)[4])%>% plotly::layout(xaxis = list(title = 'CREs'), yaxis = list(title = 'Variant switches count'), barmode = 'stack'))
  }else if(ncol(vartable) == 3){
    return(plotly::plot_ly(x = ~rownames(vartable), y = ~vartable[,1], type = 'bar', name = colnames(vartable)[1]) %>% plotly::add_trace(y = ~vartable[,2], name = colnames(vartable)[2]) %>% add_trace(y = ~vartable[,3], name = colnames(vartable)[3]) %>% plotly::layout(xaxis = list(title = 'CREs'), yaxis = list(title = 'Variant switches count'), barmode = 'stack'))
  }else if(ncol(vartable) == 2){
    return(plotly::plot_ly(x = ~rownames(vartable), y = ~vartable[,1], type = 'bar', name = colnames(vartable)[1]) %>% plotly::add_trace(y = ~vartable[,2], name = colnames(vartable)[2]) %>% layout(xaxis = list(title = 'CREs'), yaxis = list(title = 'Variant switches count'), barmode = 'stack'))
  }else if(ncol(vartable) == 1){
    return(plotly::plot_ly(x = ~rownames(vartable), y = ~vartable[,1], type = 'bar', name = colnames(vartable)[1]) %>% plotly::layout(xaxis = list(title = 'CREs'), yaxis = list(title = 'Variant switches count'), barmode = 'stack'))
  }
}
