#' @title Binary CREs in Local Heterogeneous Network
#' @description Locate the binary CREs in regional Heterogeneous network which are experimentally verified.
#'
#' @param elementList a vector consist of elements including circRNA, lncRNA, miRNA ,mRNA or KEGG/Reactome/Wikipathway pathway IDs.
#' @param regulationType character string naming the type of regulation which can be located including "circRNA-miRNA", "lncRNA-miRNA", "lncRNA-mRNA", "miRNA-mRNA", "miRNA-pathway", "mRNA-pathway".
#' @param allRegulation a logical value indicates exploring all related CREs from the storage (TRUE) or internal CREs between input elements (FALSE).
#'
#' @examples
#' binaryRegulation(elementList = c("2309", "3838", "MIMAT0000255"), regulationType = "miRNA-mRNA",  allRegulation = FALSE)
#' @export
binaryRegulation = function(elementList = NULL, regulationType = NULL, allRegulation = FALSE) {

  if(is.null(elementList)){
    return("Please enter elements......")
  }

  if(is.null(regulationType) | !(regulationType %in% c("circRNA-miRNA", "lncRNA-miRNA", "lncRNA-mRNA", "miRNA-mRNA", "miRNA-pathway", "mRNA-pathway"))){
    return("Please choose a regulation type from circRNA-miRNA, lncRNA-miRNA, lncRNA-mRNA, miRNA-mRNA, miRNA-pathway, mRNA-pathway")
  }
  load(system.file("extdata", "netelements.rda", package = "NetLCP"))
  if(!all(elementList %in% NETELEMENTS)){
    print(paste0("Filtering the missing input elements in input network......"))
    print(paste0(paste0(elementList[!(elementList %in% NETELEMENTS)], collapse = "/"), " have been filtered....."))
    elementList = elementList[elementList %in% NETELEMENTS]
    print(paste0("Now remain ", length(elementList)))
  }

  print("Binary Regulation extraction begins, please wait while we do something......")
  load(system.file("extdata", "extraction.rda", package = "NetLCP"))

  REG =  REG %>% dplyr::filter(regType == regulationType)

  if(allRegulation == TRUE){
    ExtractedData = REG[REG$node1 %in% elementList | REG$node2 %in% elementList,]
    if(nrow(ExtractedData) == 0){
      return("Sorry, search result is empty, please try to check the ID types of input elements or the regulationtype......")
    }else{
      return(ExtractedData)
    }
  }else{
    ExtractedData = REG[REG$node1 %in% elementList & REG$node2 %in% elementList,]
    if(nrow(ExtractedData) == 0){
      return("Sorry, search result is empty, please try to check the ID types of input elements or the regulationtype......")
    }else{
      return(ExtractedData)
    }
  }
}
