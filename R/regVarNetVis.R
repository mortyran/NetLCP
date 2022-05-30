#' @title Variant 'switches' Network Visualization
#' @description Variant 'switches' network visualization
#'
#' @param regVar the standard out of the function regVarDetection.
#' @param regulationType a character representing the CREs type, it should be "miRNA-mRNA", "miRNA-mRNA-pathway", "lncRNA-miRNA-mRNA", "circRNA-miRNA-mRNA", "lncRNA-miRNA-mRNA-pathway" or "circRNA-miRNA-mRNA-pathway".
#' @param selectNode select a set of nodes to visualize in the network.
#'
#' @examples
#' regVarNetVis(regVar = regVar, regulationType = "miRNA-mRNA-pathway", selectNode = c("MIMAT0000717","6774", "hsa04550"))
#' @export
regVarNetVis = function(regVar = NULL, regulationType = c("miRNA-mRNA", "miRNA-mRNA-pathway", "lncRNA-miRNA-mRNA", "circRNA-miRNA-mRNA", "lncRNA-miRNA-mRNA-pathway", "circRNA-miRNA-mRNA-pathway"), selectNode = NULL){

  if(is.null(regVar) | nrow(regVar) == 0){
    return(NULL)
  }
  if(!regulationType %in% c("miRNA-mRNA", "miRNA-mRNA-pathway", "lncRNA-miRNA-mRNA", "circRNA-miRNA-mRNA", "lncRNA-miRNA-mRNA-pathway", "circRNA-miRNA-mRNA-pathway")){
    return("Please choose regulation type from miRNA-mRNA, miRNA-mRNA-pathway, lncRNA-miRNA-mRNA, circRNA-miRNA-mRNA, lncRNA-miRNA-mRNA-pathway or circRNA-miRNA-mRNA-pathway")
  }
  if(is.null(selectNode)){
    return("selectNode is empty!")
  }

  miRNAID = c()
  lnc_cirRNAID = c()
  pathwayID = c()
  mRNAID = c()
  for(i in selectNode){
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
    if(length(lnc_cirRNAID) == 0 | length(miRNAID) == 0 | length(mRNAID) == 0){
      return(NULL)
    }
    regVar_filter = regVar[regVar$ceRNA %in% lnc_cirRNAID & regVar$miRNA %in% miRNAID & regVar$mRNA %in% mRNAID, ]
    if(nrow(regVar_filter) == 0){
      return(NULL)
    }
    lnc_mim = regVar_filter[,c(1,2)] %>% dplyr::rename(node1=ceRNA, node2=miRNA)
    lnc_mim$group = "lncRNA-miRNA"
    mim_m = regVar_filter[,c(2,3)] %>% dplyr::rename(node1=miRNA, node2=mRNA)
    mim_m$group = "miRNA-mRNA"
    lnc_var = regVar_filter[,c(1,4,6)] %>% dplyr::rename(node1=ceRNA, node2=Location, group=lncMutType)
    mi_var = regVar_filter[,c(2,11,12)] %>% dplyr::rename(node1=miRNA, node2=SNPID, group=mimMutType) %>% filter(group == "Seed_SNP")
    m_var = regVar_filter[,c(3,11,12)] %>% dplyr::rename(node1=mRNA, node2=SNPID, group=mimMutType) %>% filter(group == "UTR3_SNP")

    regVar_edgeInfo = rbind(lnc_mim, mim_m, lnc_var, mi_var, m_var) %>% distinct()
    regVar_edgeInfo$title = regVar_edgeInfo$group
    colnames(regVar_edgeInfo) = c("node1", "node2", "group", "title")
  }else if(regulationType == "lncRNA-miRNA-mRNA-pathway" | regulationType == "circRNA-miRNA-mRNA-pathway"){
    if(length(lnc_cirRNAID) == 0 | length(miRNAID) == 0 | length(mRNAID) == 0 | length(pathwayID) == 0){
      return(NULL)
    }
    regVar_filter = regVar[regVar$ceRNA %in% lnc_cirRNAID & regVar$miRNA %in% miRNAID & regVar$mRNA %in% mRNAID & regVar$pathway %in% pathwayID, ]
    if(nrow(regVar_filter) == 0){
      return(NULL)
    }
    lnc_mim = regVar_filter[,c(1,2)] %>% dplyr::rename(node1=ceRNA, node2=miRNA)
    lnc_mim$group = "lncRNA-miRNA"
    mim_m = regVar_filter[,c(2,3)] %>% dplyr::rename(node1=miRNA, node2=mRNA)
    mim_m$group = "miRNA-mRNA"
    m_path = regVar_filter[,c(3,14)] %>% dplyr::rename(node1=mRNA, node2=pathway)
    m_path$group = "mRNA-pathway"
    lnc_var = regVar_filter[,c(1,4,6)] %>% dplyr::rename(node1=ceRNA, node2=Location, group=lncMutType)
    mi_var = regVar_filter[,c(2,11,12)] %>% dplyr::rename(node1=miRNA, node2=SNPID, group=mimMutType) %>% filter(group == "Seed_SNP")
    m_var = regVar_filter[,c(3,11,12)] %>% dplyr::rename(node1=mRNA, node2=SNPID, group=mimMutType) %>% filter(group == "UTR3_SNP")

    regVar_edgeInfo = rbind(lnc_mim, mim_m, m_path, lnc_var, mi_var, m_var) %>% distinct()
    regVar_edgeInfo$title = regVar_edgeInfo$group
    colnames(regVar_edgeInfo) = c("node1", "node2", "group", "title")
  }else if(regulationType == "miRNA-mRNA-pathway"){
    if(length(miRNAID) == 0 | length(mRNAID) == 0 | length(pathwayID) == 0){
      return(NULL)
    }
    regVar_filter = regVar[regVar$miRNAID %in% miRNAID & regVar$geneID %in% mRNAID & regVar$pathway %in% pathwayID, ]
    if(nrow(regVar_filter) == 0){
      return(NULL)
    }
    mim_m = regVar_filter[,c(1,2)] %>% dplyr::rename(node1=miRNAID, node2=geneID)
    mim_m$group = "miRNA-mRNA"
    m_path = regVar_filter[,c(2,6)] %>% dplyr::rename(node1=geneID, node2=pathway)
    m_path$group = "mRNA-pathway"
    mi_var = regVar_filter[,c(1,3,4)] %>% dplyr::rename(node1=miRNAID, node2=SNPID, group=MutType) %>% filter(group == "Seed_SNP")
    m_var = regVar_filter[,c(2,3,4)] %>% dplyr::rename(node1=geneID, node2=SNPID, group=MutType) %>% filter(group == "UTR3_SNP")

    regVar_edgeInfo = rbind(mim_m, m_path, mi_var, m_var) %>% distinct()
    regVar_edgeInfo$title = regVar_edgeInfo$group
    colnames(regVar_edgeInfo) = c("node1", "node2", "group", "title")
  }else if(regulationType == "miRNA-mRNA"){
    if(length(miRNAID) == 0 | length(mRNAID) == 0){
      return(NULL)
    }
    regVar_filter = regVar[regVar$miRNAID %in% miRNAID & regVar$geneID %in% mRNAID, ]
    if(nrow(regVar_filter) == 0){
      return(NULL)
    }
    mim_m = regVar_filter[,c(1,2)] %>% dplyr::rename(node1=miRNAID, node2=geneID)
    mim_m$group = "miRNA-mRNA"
    mi_var = regVar_filter[,c(1,3,4)] %>% dplyr::rename(node1=miRNAID, node2=SNPID, group=MutType) %>% filter(group == "Seed_SNP")
    m_var = regVar_filter[,c(2,3,4)] %>% dplyr::rename(node1=geneID, node2=SNPID, group=MutType) %>% filter(group == "UTR3_SNP")
    regVar_edgeInfo = rbind(mim_m, mi_var, m_var) %>% distinct()
    regVar_edgeInfo$title = regVar_edgeInfo$group
    colnames(regVar_edgeInfo) = c("node1", "node2", "group", "title")
  }

  edgeInfo = regVar_edgeInfo %>% distinct()

  nodeInfo = data.frame(table(c(edgeInfo$node1, edgeInfo$node2)))
  colnames( nodeInfo) = c("label", "value")

  regVar_label = regVar_edgeInfo[,c(2, 3)] %>% distinct()
  colnames(regVar_label) = c("label", "group")
   nodeInfo =  nodeInfo %>% dplyr::left_join(regVar_label, by = "label")
   nodeInfo$group[is.na( nodeInfo$group)] = sapply( nodeInfo$label[is.na( nodeInfo$group)], function(x){
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
  nodeInfo =  nodeInfo %>% dplyr::distinct()
  nodeInfo$id = 1:nrow( nodeInfo)
  nodeInfo$title =  nodeInfo$group

  nodeInfo_1 = nodeInfo[,c(1, 4)] %>% dplyr::rename(node1=label,from=id)
  nodeInfo_2 = nodeInfo[,c(1, 4)] %>% dplyr::rename(node2=label,to=id)

  edgeInfo = dplyr::left_join(edgeInfo, nodeInfo_1, by = "node1")
  edgeInfo = dplyr::left_join(edgeInfo, nodeInfo_2, by = "node2")
  edgeInfo = edgeInfo[,-c(1,2)]
  return(visNetwork::visNetwork(nodes = nodeInfo, edges = edgeInfo) %>% visNetwork::visIgraphLayout(layout = "layout_nicely") %>% visNetwork::visOptions(manipulation = TRUE, highlightNearest = TRUE,selectedBy = list(variable = "group", highlight = TRUE), nodesIdSelection = TRUE))
}

