#' @title Multielements CREs in Local Heterogeneous Network
#' @description Locate the multielements CREs in regional heterogeneous network which are experimentally verified.
#'
#' @param elementList a vector consist of elements including circRNA, lncRNA, miRNA, mRNA IDs.
#' @param regulationType character string naming the type of CREs including "circRNA-miRNA-mRNA", "lncRNA-miRNA-mRNA", "lncRNA-miRNA-pathway" and "circRNA-miRNA-pathway".
#' @param allRegulation a logical value indicates exploring all related CREs from the storage or internal CREs between input elements.
#'
#' @examples
#' multieleRegulation(elementList = c("2309", "3838", "ENSG00000259366", "MIMAT0000255"), regulationType = "lncRNA-miRNA-mRNA",  allRegulation = FALSE)
#' @export
multieleRegulation = function(elementList = NULL, regulationType = NULL, allRegulation = FALSE) {

  if(is.null(elementList)){
    return("Please enter at least one element......")
  }

  if(is.null(regulationType) | !(regulationType %in% c("circRNA-miRNA-mRNA", "lncRNA-miRNA-mRNA", "lncRNA-miRNA-mRNA-pathway", "circRNA-miRNA-mRNA-pathway", "miRNA-mRNA-pathway"))){
    return("Please choose a regulation type from circRNA-miRNA-mRNA, lncRNA-miRNA-mRNA, miRNA-mRNA-pathway, lncRNA-miRNA-mRNA-pathway, circRNA-miRNA-mRNA-pathway")
  }
  load(system.file("extdata", "netelements.rda", package = "NetLCP"))
  if(!all(elementList %in% NETELEMENTS)){
    print(paste0("Filtering the missing input elements in input network......"))
    print(paste0(paste0(elementList[!(elementList %in% NETELEMENTS)], collapse = "/"), " have been filtered....."))
    elementList = elementList[elementList %in% NETELEMENTS]
    print(paste0("Now remain ", length(elementList)))
  }

  print("Multielement regulation extraction begins, please wait while we do something......")
  # setwd(as.character(utils::installed.packages()[,"LibPath"][names(utils::installed.packages()[,"LibPath"]) == "NetLCP"]))
  load(system.file("extdata", "extraction.rda", package = "NetLCP"))

  if(regulationType == "circRNA-miRNA-mRNA"){
    REG =  REG %>% filter(regType == "circRNA-miRNA" | regType == "miRNA-mRNA")
    if(allRegulation == TRUE){
      ExtractedData = REG[REG$node1 %in% elementList | REG$node2 %in% elementList,]
      if("circRNA-miRNA" %in% ExtractedData$regType & "miRNA-mRNA" %in% ExtractedData$regType){
        ExtractedData_1 = ExtractedData %>% dplyr::filter(regType == "circRNA-miRNA")
        ExtractedData_2 = ExtractedData %>% dplyr::filter(regType == "miRNA-mRNA")
        intersect_miRNA = intersect(ExtractedData_1$node1, ExtractedData_2$node1)
        if(length(intersect_miRNA) > 0){
          ExtractedData = ExtractedData[ExtractedData$node1 %in% intersect_miRNA, ]
          return(ExtractedData)
        }else{
          return("Sorry, search result is empty, please try to check the ID types of input elements or the regulationtype......")
        }
      }else{
        return("Sorry, search result is empty, please try to check the ID types of input elements or the regulationtype......")
      }
    }else{
      ExtractedData = REG[REG$node1 %in% elementList & REG$node2 %in% elementList,]
      if("circRNA-miRNA" %in% ExtractedData$regType & "miRNA-mRNA" %in% ExtractedData$regType){
        return(ExtractedData)
      }else{
        return("Sorry, search result is empty, please try to check the ID types of input elements or the regulationtype......")
      }
    }
  }else if(regulationType == "lncRNA-miRNA-mRNA"){
    REG =  REG %>% dplyr::filter(regType == "lncRNA-miRNA" | regType == "miRNA-mRNA")
    if(allRegulation == TRUE){
      ExtractedData = REG[REG$node1 %in% elementList | REG$node2 %in% elementList,]
      if("lncRNA-miRNA" %in% ExtractedData$regType & "miRNA-mRNA" %in% ExtractedData$regType){
        ExtractedData_1 = ExtractedData %>% dplyr::filter(regType == "lncRNA-miRNA")
        ExtractedData_2 = ExtractedData %>% dplyr::filter(regType == "miRNA-mRNA")
        intersect_miRNA = intersect(ExtractedData_1$node1, ExtractedData_2$node1)
        if(length(intersect_miRNA) > 0){
          ExtractedData = ExtractedData[ExtractedData$node1 %in% intersect_miRNA, ]
          return(ExtractedData)
        }else{
          return("Sorry, search result is empty, please try to check the ID types of input elements or the regulationtype......")
        }
      }else{
        return("Sorry, search result is empty, please try to check the ID types of input elements or the regulationtype......")
      }
    }else{
      ExtractedData = REG[REG$node1 %in% elementList & REG$node2 %in% elementList,]
      if("lncRNA-miRNA" %in% ExtractedData$regType & "miRNA-mRNA" %in% ExtractedData$regType){
        return(ExtractedData)
      }else{
        return("Sorry, search result is empty, please try to check the ID types of input elements or the regulationtype......")
      }
    }
  }else if(regulationType == "lncRNA-miRNA-mRNA-pathway"){
    REG =  REG %>% filter(regType == "lncRNA-miRNA" | regType == "miRNA-mRNA" | regType == "mRNA-pathway")
    if(allRegulation == TRUE){
      # print("HERE1......")
      ExtractedData = REG[REG$node1 %in% elementList | REG$node2 %in% elementList,]
      if("lncRNA-miRNA" %in% ExtractedData$regType & "miRNA-mRNA" %in% ExtractedData$regType & "mRNA-pathway" %in% ExtractedData$regType){
        # print("HERE2......")
        ExtractedData_1 = ExtractedData %>% dplyr::filter(regType == "lncRNA-miRNA")
        ExtractedData_2 = ExtractedData %>% dplyr::filter(regType == "miRNA-mRNA")
        ExtractedData_3 =  ExtractedData %>% dplyr::filter(regType == "mRNA-pathway")
        intersect_miRNA = intersect(ExtractedData_1$node1, ExtractedData_2$node1)
        if(length(intersect_miRNA) > 0){
          # print("HERE3......")
          ExtractedData_4 = ExtractedData[ExtractedData$node1 %in% intersect_miRNA, ]
          intersect_mRNA = intersect(ExtractedData_4$node2, ExtractedData_3$node2)
          if(length(intersect_mRNA) > 0){
            # print("HERE4......")
            ExtractedData = rbind(ExtractedData_4[ExtractedData_4$node2 %in% intersect_mRNA,], ExtractedData_3[ExtractedData_3$node2 %in% intersect_mRNA,])
            return(ExtractedData)
          }else{
            return("Sorry, search result is empty, please try to check the ID types of input elements or the regulationtype......")
          }
        }else{
          return("Sorry, search result is empty, please try to check the ID types of input elements or the regulationtype......")
        }
      }else{
        return("Sorry, search result is empty, please try to check the ID types of input elements or the regulationtype......")
      }
    }else{
      ExtractedData = REG[REG$node1 %in% elementList & REG$node2 %in% elementList,]
      if("lncRNA-miRNA" %in% ExtractedData$regType & "miRNA-mRNA" %in% ExtractedData$regType & "mRNA-pathway" %in% ExtractedData$regType){
        ExtractedData_1 = ExtractedData %>% dplyr::filter(regType == "lncRNA-miRNA")
        ExtractedData_2 = ExtractedData %>% dplyr::filter(regType == "miRNA-mRNA")
        ExtractedData_3 =  ExtractedData %>% dplyr::filter(regType == "mRNA-pathway")
        intersect_mRNA = intersect(ExtractedData_2$node2, ExtractedData_3$node2)
        if(length(intersect_mRNA) > 0){
          ExtractedData = rbind(ExtractedData_1, ExtractedData[ExtractedData$node2 %in% intersect_mRNA, ])
          return(ExtractedData)
        }else{
          return("Sorry, search result is empty, please try to check the ID types of input elements or the regulationtype......")
        }
      }else{
        return("Sorry, search result is empty, please try to check the ID types of input elements or the regulationtype......")
      }
    }
  }else if(regulationType == "circRNA-miRNA-mRNA-pathway"){
    REG =  REG %>% filter(regType == "circRNA-miRNA" | regType == "miRNA-mRNA" | regType == "mRNA-pathway")
    if(allRegulation == TRUE){
      ExtractedData = REG[REG$node1 %in% elementList | REG$node2 %in% elementList,]
      if("circRNA-miRNA" %in% ExtractedData$regType & "miRNA-mRNA" %in% ExtractedData$regType & "mRNA-pathway" %in% ExtractedData$regType){
        ExtractedData_1 = ExtractedData %>% dplyr::filter(regType == "circRNA-miRNA")
        ExtractedData_2 = ExtractedData %>% dplyr::filter(regType == "miRNA-mRNA")
        ExtractedData_3 =  ExtractedData %>% dplyr::filter(regType == "mRNA-pathway")
        intersect_miRNA = intersect(ExtractedData_1$node1, ExtractedData_2$node1)
        if(length(intersect_miRNA) > 0){
          ExtractedData_4 = ExtractedData[ExtractedData$node1 %in% intersect_miRNA, ]
          intersect_mRNA = intersect(ExtractedData_4$node2, ExtractedData_3$node2)
          if(length(intersect_mRNA) > 0){
            ExtractedData = rbind(ExtractedData_4[ExtractedData_4$node2 %in% intersect_mRNA,], ExtractedData_3[ExtractedData_3$node2 %in% intersect_mRNA,])
            return(ExtractedData)
          }else{
            return("Sorry, search result is empty, please try to check the ID types of input elements or the regulationtype......")
          }
        }else{
          return("Sorry, search result is empty, please try to check the ID types of input elements or the regulationtype......")
        }
      }else{
        return("Sorry, search result is empty, please try to check the ID types of input elements or the regulationtype......")
      }
    }else{
      ExtractedData = REG[REG$node1 %in% elementList & REG$node2 %in% elementList,]
      if("circRNA-miRNA" %in% ExtractedData$regType & "miRNA-mRNA" %in% ExtractedData$regType & "mRNA-pathway" %in% ExtractedData$regType){
        ExtractedData_1 = ExtractedData %>% dplyr::filter(regType == "circRNA-miRNA")
        ExtractedData_2 = ExtractedData %>% dplyr::filter(regType == "miRNA-mRNA")
        ExtractedData_3 =  ExtractedData %>% dplyr::filter(regType == "mRNA-pathway")
        intersect_mRNA = intersect(ExtractedData_2$node2, ExtractedData_3$node2)
        if(length(intersect_mRNA) > 0){
          ExtractedData = rbind(ExtractedData_1, ExtractedData[ExtractedData$node2 %in% intersect_mRNA, ])
          return(ExtractedData)
        }else{
          return("Sorry, search result is empty, please try to check the ID types of input elements or the regulationtype......")
        }
      }else{
        return("Sorry, search result is empty, please try to check the ID types of input elements or the regulationtype......")
      }
    }
  }else if(regulationType == "miRNA-mRNA-pathway"){
    REG =  REG %>% filter(regType == "miRNA-mRNA" | regType == "mRNA-pathway")
    if(allRegulation == TRUE){
      ExtractedData = REG[REG$node1 %in% elementList | REG$node2 %in% elementList,]
      if("miRNA-mRNA" %in% ExtractedData$regType & "mRNA-pathway" %in% ExtractedData$regType){
        ExtractedData_1 = ExtractedData %>% dplyr::filter(regType == "miRNA-mRNA")
        ExtractedData_2 = ExtractedData %>% dplyr::filter(regType == "mRNA-pathway")
        intersect_mRNA = intersect(ExtractedData_1$node2, ExtractedData_2$node1)
        if(length(intersect_mRNA) > 0){
          ExtractedData = ExtractedData[(ExtractedData$node1 %in% intersect_mRNA) | (ExtractedData$node2 %in% intersect_mRNA), ]
          return(ExtractedData)
        }else{
          return("Sorry, search result is empty, please try to check the ID types of input elements or the regulationtype......")
        }
      }else{
        return("Sorry, search result is empty, please try to check the ID types of input elements or the regulationtype......")
      }
    }else{
      ExtractedData = REG[REG$node1 %in% elementList & REG$node2 %in% elementList,]
      if("miRNA-mRNA" %in% ExtractedData$regType & "mRNA-pathway" %in% ExtractedData$regType){
        return(ExtractedData)
      }else{
        return("Sorry, search result is empty, please try to check the ID types of input elements or the regulationtype......")
      }
    }
  }
}
