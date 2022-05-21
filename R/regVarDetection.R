#' @title Variants on regulations Detection
#' @description Detect the Variants on regulations, which are integrated from LnCeVar and miRNASNPv3.
#'
#' @param regData the dataframe of output from binaryRegulation and multieleRegulation functions,
#'                accepting miRNA-mRNA, miRNA-mRNA-pathway, lncRNA-miRNA-mRNA, circRNA-miRNA-mRNA, lncRNA-miRNA-mRNA-pathway, circRNA-miRNA-mRNA-pathway.
#' @param regulationType a character representing the extracted regulation data type, it should be "miRNA-mRNA", "miRNA-mRNA-pathway", "lncRNA-miRNA-mRNA", "circRNA-miRNA-mRNA", "lncRNA-miRNA-mRNA-pathway" or "circRNA-miRNA-mRNA-pathway".
#' @examples
#' regVarDetection(regData = NULL)
#' @export
regVarDetection = function(regData = NULL, regulationType = c("miRNA-mRNA", "miRNA-mRNA-pathway", "lncRNA-miRNA-mRNA", "circRNA-miRNA-mRNA", "lncRNA-miRNA-mRNA-pathway", "circRNA-miRNA-mRNA-pathway")){
  if(is.null(regulationType) | !(regulationType %in% c("miRNA-mRNA", "miRNA-mRNA-pathway", "lncRNA-miRNA-mRNA", "circRNA-miRNA-mRNA", "lncRNA-miRNA-mRNA-pathway", "circRNA-miRNA-mRNA-pathway"))){
    return("Please choose a regulation type from miRNA-mRNA, miRNA-mRNA-pathway, lncRNA-miRNA-mRNA, circRNA-miRNA-mRNA, lncRNA-miRNA-mRNA-pathway, circRNA-miRNA-mRNA-pathway")
  }

  if(is.null(regData)){
    return("Please enter at least one element......")
  }else{
    print("Variants on regulations extracting begins......")
  }
  load(paste0(system.file("extdata", package = "NetLCP"), "/extraction.rda"))

  if(regulationType == "lncRNA-miRNA-mRNA" | regulationType == "circRNA-miRNA-mRNA"){
    all_node = unique(c(regData$node1, regData$node2))
    ce_var = CEREGVAR[(CEREGVAR$ceRNA %in% all_node) & (CEREGVAR$miRNA %in% all_node) & (CEREGVAR$mRNA %in% all_node), ]
    ce_var$mim = paste(ce_var$miRNA, ce_var$mRNA, sep = "-")
    mim_var = MIMSNP[(MIMSNP$miRNAID %in% all_node) & (MIMSNP$geneID %in% all_node), ]
    mim_var$mim = paste(mim_var$miRNAID, mim_var$geneID, sep = "-")
    ce_var = ce_var %>% dplyr::left_join(mim_var, by = "mim") %>% distinct() %>% na.omit()
    colnames(ce_var)[c(6,7,12,13)] = c("lncMutType", "Source1", "mimMutType", "Source2")
    if(nrow(ce_var) > 0){
      return(ce_var)
    }else{
      return("Sorry, search result is empty")
    }

  }else if(regulationType == "lncRNA-miRNA-mRNA-pathway" | regulationType == "circRNA-miRNA-mRNA-pathway"){
    all_node = unique(c(regData$node1, regData$node2))
    ce_var = CEREGVAR[(CEREGVAR$ceRNA %in% all_node) & (CEREGVAR$miRNA %in% all_node) & (CEREGVAR$mRNA %in% all_node), ]
    ce_var$mim = paste(ce_var$miRNA, ce_var$mRNA, sep = "-")

    mim_var = MIMSNP[(MIMSNP$miRNAID %in% all_node) & (MIMSNP$geneID %in% all_node), ]
    mim_var$mim = paste(mim_var$miRNAID, mim_var$geneID, sep = "-")

    ce_var = ce_var %>% dplyr::left_join(mim_var, by = "mim") %>% dplyr::distinct() %>% na.omit()
    colnames(ce_var)[c(6,7,12,13)] = c("lncMutType", "Source1", "mimMutType", "Source2")

    m_path = regData %>% dplyr::filter(regType == "mRNA-pathway")
    colnames(m_path) = c("pathway", "mRNA", "pathwaySource", "regType")
    ce_var$mRNA = as.character(ce_var$mRNA)
    ce_var = dplyr::left_join(ce_var, m_path, by = "mRNA") %>% dplyr::distinct() %>% na.omit()

    if(nrow(ce_var) > 0){
      return(ce_var)
    }else{
      return("Sorry, search result is empty")
    }
  }else if(regulationType == "miRNA-mRNA-pathway"){
    all_node = unique(c(regData$node1, regData$node2))
    mim_var = MIMSNP[(MIMSNP$miRNAID %in% all_node) & (MIMSNP$geneID %in% all_node), ]
    m_path = regData %>% dplyr::filter(regType == "mRNA-pathway")
    colnames(m_path) = c("pathway", "geneID", "pathwaySource", "regType")
    mim_var$geneID = as.character(mim_var$geneID)
    mim_var = dplyr::left_join(mim_var, m_path, by = "geneID") %>% dplyr::distinct() %>% na.omit()

    if(nrow(mim_var) > 0){
      return(mim_var)
    }else{
      return("Sorry, search result is empty")
    }

  }else if(regulationType == "miRNA-mRNA"){
    all_node = unique(c(regData$node1, regData$node2))
    mim_var = MIMSNP[(MIMSNP$miRNAID %in% all_node) & (MIMSNP$geneID %in% all_node), ]
    if(nrow(mim_var) > 0){
      return(mim_var)
    }else{
      return("Sorry, search result is empty")
    }
  }
}
