#' @title eQTLs of CREs statistics
#' @description eQTLs of CREs statistics
#'
#' @param regData the standard out of the function binaryRegulationExtract or multieleRegulationExtract.
#' @param eQTLsData the standard out of the function eQTLsDetection.
#' @param regulationType a character representing the CREs type, it should be one of "circRNA-miRNA", "lncRNA-miRNA", "lncRNA-mRNA", "miRNA-mRNA", "miRNA-pathway", "mRNA-pathway", "circRNA-miRNA-mRNA", "lncRNA-miRNA-mRNA", "miRNA-mRNA-pathway", "lncRNA-miRNA-mRNA-pathway", "circRNA-miRNA-mRNA-pathway".
#' @param topCREs an integer representing the number of the top prioritized CREs to exhibit.
#' @param filterDegree an integer for filtering the nodes.
#' @param selectCREs a vector of elements or CREs.
#'
#' @examples
#' eQTLsRegStat(regData = NULL, eQTLsData = NULL, regulationType = "miRNA-mRNA", filterDegree = 40, topCREs = 10, selectCREs = NULL)
#' @export
eQTLsRegStat = function(regData = NULL, eQTLsData = NULL, regulationType = NULL, topCREs = NULL, filterDegree = 40, selectCREs = NULL){

  if(is.null(regData) | is.null(eQTLsData) | is.null(regulationType)){
    return("Please input the valid parameter......")
  }

  if(is.null(regulationType) | !(regulationType %in% c("circRNA-miRNA", "lncRNA-miRNA", "lncRNA-mRNA", "miRNA-mRNA", "miRNA-pathway", "mRNA-pathway", "circRNA-miRNA-mRNA", "lncRNA-miRNA-mRNA", "miRNA-mRNA-pathway", "lncRNA-miRNA-mRNA-pathway", "circRNA-miRNA-mRNA-pathway"))){
    return("Please choose a regulation type from circRNA-miRNA, lncRNA-miRNA, lncRNA-mRNA, miRNA-mRNA, miRNA-pathway, mRNA-pathway, circRNA-miRNA-mRNA, lncRNA-miRNA-mRNA, miRNA-mRNA-pathway, lncRNA-miRNA-mRNA-pathway, circRNA-miRNA-mRNA-pathway")
  }

  if(is.null(selectCREs)){
    nodeInfo = data.frame(table(c(regData$node1, regData$node2)))
    nodeFilter = nodeInfo$Var1[nodeInfo$Freq >= filterDegree]
    regData = regData[regData$node1 %in% nodeFilter & regData$node2 %in% nodeFilter,]
    eqtls = eQTLsData[eQTLsData$Elements %in% nodeFilter,]
  }else{

    if(any(grepl("_", selectCREs))){
      elements_all = c()
      for(i in selectCREs){
        elements_all = c(unlist(strsplit(i, "_")), elements_all)
      }
      selectCREs = unique(elements_all)
    }

    regData = regData[regData$node1 %in% selectCREs & regData$node2 %in% selectCREs,]
    eqtls = eQTLsData[eQTLsData$Elements %in% selectCREs,]
  }

  if(regulationType == "lncRNA-miRNA-mRNA-pathway" | regulationType == "circRNA-miRNA-mRNA-pathway"){
    all_regtypte = unique(regData$regType)
    if(length(all_regtypte) == 3){
      A1 = regData %>% filter(regType == all_regtypte[1]) %>% dplyr::select(node1, node2) %>% dplyr::rename(node1=node2, node2=node1)
      A2 = regData %>% filter(regType == all_regtypte[2]) %>% dplyr::select(node1, node2) %>% dplyr::rename(node2=node1, node3=node2)
      A3 = regData %>% filter(regType == all_regtypte[3]) %>% dplyr::select(node1, node2) %>% dplyr::rename(node4=node1, node3=node2)
      regulation = A1[,c(2,1)] %>% dplyr::left_join(A2, by = "node2") %>% dplyr::left_join(A3, by = "node3") %>% dplyr::distinct() %>% na.omit()
      regulation$reg = paste(regulation$node1, regulation$node2, regulation$node3, regulation$node4, sep = "_")

      eqtls_stat = table(eqtls$Elements)
      regulation$eQTLs = apply(regulation[,c(1,2,3,4)], 1, function(x){
        sum(eqtls_stat[rownames(eqtls_stat) %in% x])
      })

      regulation = regulation %>% dplyr::arrange(desc(eQTLs))


      if(is.null(topCREs)){
        top_num = nrow(regulation)
      }else{
        if(topCREs <= nrow(regulation) & topCREs > 0){
          top_num = topCREs
        }else{
          top_num = nrow(regulation)
        }
      }

      regulation = regulation[1:top_num,]

      print("saving prioritized CREs eQTLs......")
      write.csv(regulation %>% select(reg, eQTLs), file = "CREs_eQTLs.csv", quote = F, row.names = F)
      print("Prioritized CREs eQTLs has been output !")

      return(regulation %>% plotly::plot_ly(x = ~reg, y = ~eQTLs, type = 'bar', name = 'eQTLs') %>% plotly::layout(xaxis = list(title = "CREs", categoryorder = "trace"), yaxis = list(title = 'eQTLs Count'), barmode = 'stack'))

      # if(ncol(eqtls_stat) == 2){
      #   regulation$eQTLs = apply(regulation[,c(1,2,3,4)], 1, function(x){
      #     sum(eqtls_stat[rownames(eqtls_stat) %in% x])
      #   })
      #   regulation$trans_eQTLs = apply(regulation[,c(1,2,3,4)], 1, function(x){
      #     sum(eqtls_stat[,2][rownames(eqtls_stat) %in% x])
      #   })
      #   return(regulation %>% plotly::plot_ly(regulation, x = ~reg, y = ~cis_eQTLs, type = 'bar', name = 'cis-eQTLs') %>% plotly::add_trace(y = ~trans_eQTLs, name = 'trans-eQTLs') %>% plotly::layout(xaxis = list(title = 'Regulation'), yaxis = list(title = 'eQTLs'), barmode = 'stack'))
      # }else{
      #   regulation$cis_eQTLs = apply(regulation[,c(1,2,3,4)], 1, function(x){
      #     sum(eqtls_stat[,1][rownames(eqtls_stat) %in% x])
      #   })
      #   return(plotly::plot_ly(regulation, x = ~reg, y = ~cis_eQTLs, type = 'bar', name = 'cis-eQTLs') %>% plotly::layout(xaxis = list(title = 'Regulation'), yaxis = list(title = 'eQTLs'), barmode = 'stack'))
      # }
    }

  }else if(regulationType == "lncRNA-miRNA-mRNA" | regulationType == "circRNA-miRNA-mRNA"){
    all_regtypte = unique(regData$regType)
    if(length(all_regtypte) == 2){
      A1 = regData %>% filter(regType == all_regtypte[1]) %>% select(node1, node2) %>% dplyr::rename(node1=node2, node2=node1)
      A2 = regData %>% filter(regType == all_regtypte[2]) %>% select(node1, node2) %>% dplyr::rename(node2=node1,node3=node2)
      regulation = A1[,c(2,1)] %>% left_join(A2, by = "node2") %>% dplyr::distinct() %>% na.omit()
      regulation$reg = paste(regulation$node1, regulation$node2, regulation$node3, sep = "_")

      eqtls_stat = table(eqtls$Elements)
      regulation$eQTLs = apply(regulation[,c(1,2,3,4)], 1, function(x){
        sum(eqtls_stat[rownames(eqtls_stat) %in% x])
      })
      regulation = regulation %>% dplyr::arrange(desc(eQTLs))



      if(is.null(topCREs)){
        top_num = nrow(regulation)
      }else{
        if(topCREs <= nrow(regulation) & topCREs > 0){
          top_num = topCREs
        }else{
          top_num = nrow(regulation)
        }
      }

      regulation = regulation[1:top_num,]

      print("saving prioritized CREs eQTLs......")
      write.csv(regulation %>% select(reg, eQTLs), file = "CREs_eQTLs.csv", quote = F, row.names = F)
      print("Prioritized CREs eQTLs has been output !")

      return(regulation %>% plotly::plot_ly(x = ~reg, y = ~eQTLs, type = 'bar', name = 'eQTLs') %>% plotly::layout(xaxis = list(title = "CREs", categoryorder = "trace"), yaxis = list(title = 'eQTLs Count'), barmode = 'stack'))
      #
      # eqtls_stat = table(eqtls$Elements, eqtls$RegType)
      # if(ncol(eqtls_stat) == 2){
      #   regulation$cis_eQTLs = apply(regulation[,c(1,2,3)], 1, function(x){
      #     sum(eqtls_stat[,1][rownames(eqtls_stat) %in% x])
      #   })
      #   regulation$trans_eQTLs = apply(regulation[,c(1,2,3)], 1, function(x){
      #     sum(eqtls_stat[,2][rownames(eqtls_stat) %in% x])
      #   })
      #   return(plotly::plot_ly(regulation, x = ~reg, y = ~cis_eQTLs, type = 'bar', name = 'cis-eQTLs') %>% plotly::add_trace(y = ~trans_eQTLs, name = 'trans-eQTLs') %>% plotly::layout(xaxis = list(title = 'Regulation'), yaxis = list(title = 'eQTLs'), barmode = 'stack'))
      # }else{
      #   regulation$cis_eQTLs = apply(regulation[,c(1,2,3)], 1, function(x){
      #     sum(eqtls_stat[,1][rownames(eqtls_stat) %in% x])
      #   })
      #   return(plotly::plot_ly(regulation, x = ~reg, y = ~cis_eQTLs, type = 'bar', name = 'cis-eQTLs') %>% plotly::layout(xaxis = list(title = 'Regulation'), yaxis = list(title = 'eQTLs'), barmode = 'stack'))
      # }

    }

  }else if(regulationType == "miRNA-mRNA-pathway"){
    all_regtypte = unique(regData$regType)
    if(length(all_regtypte) == 2){
      A1 = regData %>% dplyr::filter(regType == all_regtypte[1]) %>% dplyr::select(node1, node2)
      A2 = regData %>% dplyr::filter(regType == all_regtypte[2]) %>% dplyr::select(node1, node2) %>% dplyr::rename(node3=node1)
      regulation = A1 %>% dplyr::left_join(A2, by = "node2") %>% dplyr::distinct() %>% na.omit()
      regulation$reg = paste(regulation$node1, regulation$node2, regulation$node3, sep = "_")

      eqtls_stat = table(eqtls$Elements)
      regulation$eQTLs = apply(regulation[,c(1,2,3,4)], 1, function(x){
        sum(eqtls_stat[rownames(eqtls_stat) %in% x])
      })
      regulation = regulation %>% dplyr::arrange(desc(eQTLs))

      if(is.null(topCREs)){
        top_num = nrow(regulation)
      }else{
        if(topCREs <= nrow(regulation) & topCREs > 0){
          top_num = topCREs
        }else{
          top_num = nrow(regulation)
        }
      }

      regulation = regulation[1:top_num,]

      print("saving prioritized CREs eQTLs......")
      write.csv(regulation %>% select(reg, eQTLs), file = "CREs_eQTLs.csv", quote = F, row.names = F)
      print("Prioritized CREs eQTLs has been output !")


      return(regulation %>% plotly::plot_ly(x = ~reg, y = ~eQTLs, type = 'bar', name = 'eQTLs') %>% plotly::layout(xaxis = list(title = "CREs", categoryorder = "trace"), yaxis = list(title = 'eQTLs Count'), barmode = 'stack'))
      #
      #       eqtls_stat = table(eqtls$Elements, eqtls$RegType)
      #       if(ncol(eqtls_stat) == 2){
      #         regulation$cis_eQTLs = apply(regulation[,c(1,2,3)], 1, function(x){
      #           sum(eqtls_stat[,1][rownames(eqtls_stat) %in% x])
      #         })
      #         regulation$trans_eQTLs = apply(regulation[,c(1,2,3)], 1, function(x){
      #           sum(eqtls_stat[,2][rownames(eqtls_stat) %in% x])
      #         })
      #         return(plotly::plot_ly(regulation, x = ~reg, y = ~cis_eQTLs, type = 'bar', name = 'cis-eQTLs') %>% plotly::add_trace(y = ~trans_eQTLs, name = 'trans-eQTLs') %>% plotly::layout(xaxis = list(title = 'Regulation'), yaxis = list(title = 'eQTLs'), barmode = 'stack'))
      #       }else{
      #         regulation$cis_eQTLs = apply(regulation[,c(1,2,3)], 1, function(x){
      #           sum(eqtls_stat[,1][rownames(eqtls_stat) %in% x])
      #         })
      #         return(plotly::plot_ly(regulation, x = ~reg, y = ~cis_eQTLs, type = 'bar', name = 'cis-eQTLs') %>% plotly::layout(xaxis = list(title = 'Regulation'), yaxis = list(title = 'eQTLs'), barmode = 'stack'))
      #       }
    }
  }else{
    regulation = regData[,c(1,2)] %>% dplyr::distinct() %>% na.omit()
    if(nrow(regulation) > 0){
      regulation$reg = paste(regulation$node1, regulation$node2, sep = "_")

      eqtls_stat = table(eqtls$Elements)
      regulation$eQTLs = apply(regulation[,c(1,2)], 1, function(x){
        sum(eqtls_stat[rownames(eqtls_stat) %in% x])
      })
      regulation = regulation %>% dplyr::arrange(desc(eQTLs))

      if(is.null(topCREs)){
        top_num = nrow(regulation)
      }else{
        if(topCREs <= nrow(regulation) & topCREs > 0){
          top_num = topCREs
        }else{
          top_num = nrow(regulation)
        }
      }

      regulation = regulation[1:top_num,]


      print("saving prioritized CREs eQTLs......")
      write.csv(regulation %>% select(reg, eQTLs), file = "CREs_eQTLs.csv", quote = F, row.names = F)
      print("Prioritized CREs eQTLs has been output !")


      return(regulation %>% plotly::plot_ly(x = ~reg, y = ~eQTLs, type = 'bar', name = 'eQTLs') %>% plotly::layout(xaxis = list(title = "CREs", categoryorder = "trace"), yaxis = list(title = 'eQTLs Count'), barmode = 'stack'))
      # eqtls_stat = table(eqtls$Elements, eqtls$RegType)
      # if(ncol(eqtls_stat) == 2){
      #   regulation$cis_eQTLs = apply(regulation[,c(1,2)], 1, function(x){
      #     sum(eqtls_stat[,1][rownames(eqtls_stat) %in% x])
      #   })
      #   regulation$trans_eQTLs = apply(regulation[,c(1,2)], 1, function(x){
      #     sum(eqtls_stat[,2][rownames(eqtls_stat) %in% x])
      #   })
      #   regulation$reg = paste(regulation$node1, regulation$node2, sep = "_")
      #   return(plotly::plot_ly(regulation, x = ~reg, y = ~cis_eQTLs, type = 'bar', name = 'cis-eQTLs') %>% plotly::add_trace(y = ~trans_eQTLs, name = 'trans-eQTLs') %>% plotly::layout(xaxis = list(title = 'Regulation'), yaxis = list(title = 'eQTLs'), barmode = 'stack'))
      # }else{
      #   regulation$cis_eQTLs = apply(regulation[,c(1,2)], 1, function(x){
      #     sum(eqtls_stat[,1][rownames(eqtls_stat) %in% x])
      #   })
      #   return(plotly::plot_ly(regulation, x = ~reg, y = ~cis_eQTLs, type = 'bar', name = 'cis-eQTLs') %>% plotly::layout(xaxis = list(title = 'Regulation'), yaxis = list(title = 'eQTLs'), barmode = 'stack'))
      # }
    }
  }
}


