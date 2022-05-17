#' @title Binary Regulatory in Local Heterogeneous network
#' @description Locate the experimental binary regulatory data among biological elements in regional Heterogeneous network which is integrated from a wide range of bioinformatics databases.
#'
#' @param transcriptomeList a vector consist of transcriptome including circRNA, lncRNA, miRNA ,mRNA or KEGG/Reactome/Wikipathway pathway IDs.
#' @param regulationType character string naming the type of regulation which can be located including "circRNA-miRNA", "lncRNA-miRNA", "lncRNA-mRNA", "miRNA-mRNA", "miRNA-pathway", "mRNA-pathway".
#' @param allRegulation logical value indicates locating all related regulations from the storage or internal regulations between input elements.
#'
#' @examples
#' binaryRegulation(transcriptomeList = c("2309", "3838", "MIMAT0000255"), regulationType = "miRNA-mRNA",  allRegulation = FALSE)
#' @export
binaryRegulation = function(transcriptomeList = NULL, regulationType = c("circRNA-miRNA", "lncRNA-miRNA", "lncRNA-mRNA", "miRNA-mRNA", "miRNA-pathway", "mRNA-pathway"), allRegulation = FALSE) {

  if(is.null(transcriptomeList)){
    return("Please enter at least one element......")
  }

  if(!(regulationType %in% c("circRNA-miRNA", "lncRNA-miRNA", "lncRNA-mRNA", "miRNA-mRNA", "miRNA-pathway", "mRNA-pathway"))){
    return("Please choose a regulation type from circRNA-miRNA, lncRNA-miRNA, lncRNA-mRNA, miRNA-mRNA, miRNA-pathway, mRNA-pathway")
  }
  load(system.file("extdata", "netelements.rda", package = "NetLCP"))
  if(!all(transcriptomeList %in% NETELEMENTS)){
    print(paste0("Filtering the missing input transcriptome in input network......"))
    print(paste0(paste0(transcriptomeList[!(transcriptomeList %in% NETELEMENTS)], collapse = "/"), " have been filtered....."))
    transcriptomeList = transcriptomeList[transcriptomeList %in% NETELEMENTS]
    print(paste0("Now remain ", length(transcriptomeList)))
  }

  print("Binary Regulation extraction begins, please wait while we do something......")
  load(system.file("extdata", "extraction.rda", package = "NetLCP"))

  REG =  REG %>% dplyr::filter(regType == regulationType)

  if(allRegulation == TRUE){
    ExtractedData = REG[REG$node1 %in% transcriptomeList | REG$node2 %in% transcriptomeList,]
    if(nrow(ExtractedData) == 0){
      return("Sorry, search result is empty, please try to check the ID types of input transcriptome or the regulationtype......")
    }else{
      return(ExtractedData)
    }
  }else{
    ExtractedData = REG[REG$node1 %in% transcriptomeList & REG$node2 %in% transcriptomeList,]
    if(nrow(ExtractedData) == 0){
      return("Sorry, search result is empty, please try to check the ID types of input transcriptome or the regulationtype......")
    }else{
      return(ExtractedData)
    }
  }
}
