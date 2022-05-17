#' @title Single Biological Elements eQTLs Detection
#' @description Detect the eQTLs of single biological elements, which are integrated from ncRNA-eQTLs and eQTLGen.
#'
#' @param regData the dataframe of output from binaryRegulation and multieleRegulation functions,
#'
#' @examples
#' eQTLsDetection(regData = NULL)
#' @export
eQTLsDetection = function(regData = NULL){

  if(is.null(regData)){
    return("Please enter at least one element......")
  }else{
    print("Single biological elements eQTLs extracting extracting begins......")
  }
  load(paste0(system.file("extdata", package = "NetLCP"), "/extraction.rda"))
  SinEleVar = ELEVAR %>% dplyr::filter(Elements %in% unique(c(regData$node1, regData$node2)))
  if(nrow(SinEleVar) == 0){
    return("Sorry, search result is empty, please try to check the ID types of input transcriptome or the regulationtype......")
  }else{
    return(SinEleVar)
  }
}
