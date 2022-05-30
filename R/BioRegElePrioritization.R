#' @title Regulatory Elements Prioritization.
#' @description Prioritize the regulatory elements through decoding the heterogeneous network according to input transcriptome.
#'
#' @param transcriptomeList a vector consist transcriptome including circRNA, lncRNA, miRNA or mRNA.
#' @param prioType character string naming the prioritizing element types including circRNA, lncRNA, KEGG, Reactome, Wikipathway.
#' @param empiricalPvalue logical value enabling or not enabling empirical pvalue computing.
#'
#' @examples
#' BioRegElePrioritization(transcriptomeList = c("2309", "3838", "MIMAT0000255"), prioType = "KEGG", empiricalPvalue = FALSE)
#' @export
BioRegElePrioritization = function(transcriptomeList = NULL, prioType = c("circRNA", "lncRNA", "KEGG", "Reactome", "Wikipathway"), empiricalPvalue = FALSE){

  warnings(FALSE)

  if(!(prioType %in% c("circRNA", "lncRNA", "KEGG", "Reactome", "Wikipathway"))){
    return("Please choose prioritization type from circRNA, lncRNA, KEGG, Reactome, Wikipathway")
  }

  load(system.file("extdata", "inputelements.rda", package = "NetLCP"))
  if(!all(transcriptomeList %in% INPUTELEMENTS)){
    print(paste0("Filtering the missing elements of transcriptomeList in the input network......"))
    print(paste0("Element ", paste0(transcriptomeList[!(transcriptomeList %in% INPUTELEMENTS)], collapse = "/"), " have been filtered....."))
    transcriptomeList = transcriptomeList[transcriptomeList %in% INPUTELEMENTS]
    print(paste0("Now remain ", length(transcriptomeList)))
  }

  if(length(unique(transcriptomeList)) < 50){
    return("Please enter at least 50 unique elements......")
  }else{
    print("Prioritization begins, please wait while we do something......")
  }


  utils::memory.limit(400000)
  HeterogeneousNetworkNormalizing = function(SimilarityNetwork=NULL, AssociationNetworkA=NULL, AssociationNetworkB=NULL, theta=NULL){

    if(!is.matrix(SimilarityNetwork)){

      AssociationCombination = sum(AssociationNetworkA) + sum(AssociationNetworkB)

    }else{

      AssociationCombination = Matrix::rowSums(as.matrix(AssociationNetworkA)) + Matrix::rowSums(as.matrix(AssociationNetworkB))

    }

    if(!is.matrix(SimilarityNetwork)){

      SimilarityNetwork_Nor = 0
      CoefficientMat = ifelse(AssociationCombination > 0, 1 / AssociationCombination, 0)
      AssociationNetworkA_Nor = AssociationNetworkA * as.vector(CoefficientMat)
      AssociationNetworkB_Nor = AssociationNetworkB * as.vector(CoefficientMat)

      return(c(SimilarityNetwork_Nor, AssociationNetworkA_Nor, AssociationNetworkB_Nor))
    }else{

      SimilarityNetwork_ColSum = colSums(SimilarityNetwork)

      CoefficientMat = ifelse(AssociationCombination != 0 & SimilarityNetwork_ColSum != 0, (1 - theta) /  SimilarityNetwork_ColSum, 0)
      CoefficientMat1 = ifelse(AssociationCombination == 0 & SimilarityNetwork_ColSum != 0, 1 / SimilarityNetwork_ColSum, 0)
      CoefficientMat = CoefficientMat + CoefficientMat1
      SimilarityNetwork_Nor = SimilarityNetwork * as.vector(CoefficientMat)

      CoefficientMat = ifelse(AssociationCombination != 0 & SimilarityNetwork_ColSum != 0, theta / AssociationCombination, 0)
      CoefficientMat1 = ifelse(AssociationCombination != 0 & SimilarityNetwork_ColSum == 0, 1 / AssociationCombination, 0)
      CoefficientMat = CoefficientMat + CoefficientMat1
      AssociationNetworkA_Nor = AssociationNetworkA * as.vector(CoefficientMat)
      AssociationNetworkB_Nor = AssociationNetworkB * as.vector(CoefficientMat)

      return(cbind(SimilarityNetwork_Nor, AssociationNetworkA_Nor, AssociationNetworkB_Nor))
    }

  }
  NormalizeHeterogeneousNetworkGenerator = function(HeterogeneousNetwork = NULL, Transcriptome = NULL, NetworkJumpProbability = NULL){

    HNet = HeterogeneousNetwork %>% as.data.frame()

    NodeTypes = data.frame(NodeName = colnames(HNet)[-1], NodeType = HNet$NodeType)

    # generally filtering
    pred_tar = NodeTypes$NodeName[NodeTypes$NodeType == "predicting"]
    EfficHNetNodeInd = NodeTypes$NodeName %in% c(pred_tar, Transcriptome)
    NodeTypes = NodeTypes[EfficHNetNodeInd,]

    HNet = HNet %>% select(-NodeType) %>% as.matrix()
    HNet = HNet[EfficHNetNodeInd, EfficHNetNodeInd]

    #
    pred_tar = NodeTypes$NodeName[NodeTypes$NodeType == "predicting"]
    trans_in = NodeTypes$NodeName[NodeTypes$NodeType != "predicting"]
    EfficHNetNodeInd = c(rep(TRUE, length(pred_tar)), Matrix::rowSums(HNet[(length(pred_tar) + 1) : nrow(HNet), (length(pred_tar) + 1) : nrow(HNet)]) != 0)
    HNet = HNet[EfficHNetNodeInd, EfficHNetNodeInd]
    EfficHNetNodeInfo = NodeTypes[EfficHNetNodeInd,]
    NodeTypes = EfficHNetNodeInfo

    trans_in = NodeTypes$NodeName[NodeTypes$NodeType != "predicting"]
    EfficHNetNodeInd = c(colSums(HNet[(length(pred_tar) + 1) : nrow(HNet), 1:length(pred_tar)]) != 0, rep(TRUE, length(trans_in)))
    HNet = HNet[EfficHNetNodeInd, EfficHNetNodeInd]
    EfficHNetNodeInfo = NodeTypes[EfficHNetNodeInd,]
    NodeTypes = EfficHNetNodeInfo

    if(all(NodeTypes$NodeType %in% c("miRNA", "mRNA"))){
      return("Sorry, please enlarge your input transcriptome data......")
    }

    #
    #
    # NETWORK JUMP PROBABILITY
    theta = as.numeric(NetworkJumpProbability)

    # extract subnetwork from HNet
    NodeType1_index = which(NodeTypes$NodeType == unique(NodeTypes$NodeType)[1])
    NodeType2_index = which(NodeTypes$NodeType == unique(NodeTypes$NodeType)[2])
    NodeType3_index = which(NodeTypes$NodeType == unique(NodeTypes$NodeType)[3])
    (length(NodeType1_index) + length(NodeType2_index) + length(NodeType3_index)) == dim(HNet)[1]

    HNet_AA_pri = HNet[NodeType1_index, NodeType1_index]
    HNet_AB_pri = HNet[NodeType1_index, NodeType2_index]
    HNet_AC_pri = HNet[NodeType1_index, NodeType3_index]
    HNet_A_Nor = HeterogeneousNetworkNormalizing(SimilarityNetwork = HNet_AA_pri, AssociationNetworkA = HNet_AB_pri, AssociationNetworkB = HNet_AC_pri, theta = theta)

    HNet_BA_pri = HNet[NodeType2_index, NodeType1_index]
    HNet_BB_pri = HNet[NodeType2_index, NodeType2_index]
    HNet_BC_pri = HNet[NodeType2_index, NodeType3_index]
    HNet_B_Nor = HeterogeneousNetworkNormalizing(SimilarityNetwork = HNet_BB_pri, AssociationNetworkA = HNet_BA_pri, AssociationNetworkB = HNet_BC_pri, theta = theta)

    HNet_CA_pri = HNet[NodeType3_index, NodeType1_index]
    HNet_CB_pri = HNet[NodeType3_index, NodeType2_index]
    HNet_CC_pri = HNet[NodeType3_index, NodeType3_index]
    HNet_C_Nor = HeterogeneousNetworkNormalizing(SimilarityNetwork = HNet_CC_pri, AssociationNetworkA = HNet_CA_pri, AssociationNetworkB = HNet_CB_pri, theta = theta)

    HNet_Nor = rbind(HNet_A_Nor, HNet_B_Nor, HNet_C_Nor) %>% t()# colSum = 1

    gc()
    return(list(HNet_Nor, NodeTypes))
  }
  Random_Walk_Restart = function(Network_Matrix, r, seed_score_matrix){

    Threeshold = 1e-12
    NetworkSize = ncol(Network_Matrix)

    residue = 1
    iter = 1
    restart_vector =  seed_score_matrix

    while(residue >= Threeshold){

      if(iter > 50){
        break
      }

      old_seed_score_matrix = seed_score_matrix
      seed_score_matrix = (1 - r) * MatrixMultiplication(Network_Matrix, seed_score_matrix) + r * restart_vector
      residue = sqrt(sum((seed_score_matrix - old_seed_score_matrix)^2))

      if(is.na(as.logical(residue))){
        return("Sorry, please enlarge your input transcriptome")
      }

      iter = iter + 1
      gc()
    }
    return(seed_score_matrix)
  }

  RandomNormalizeHeterogeneousNetworkGenerator = function(HeterogeneousNetwork = NULL, Transcriptome = NULL, NetworkJumpProbability = NULL){

    HNet = HeterogeneousNetwork %>% as.data.frame()
    NodeTypes = data.frame(NodeName = colnames(HNet)[-1], NodeType = HNet$NodeType)

    # generally filtering
    pred_tar = NodeTypes$NodeName[NodeTypes$NodeType == "predicting"]
    EfficHNetNodeInd = NodeTypes$NodeName %in% c(pred_tar, Transcriptome)
    NodeTypes = NodeTypes[EfficHNetNodeInd,]

    HNet = HNet %>% dplyr::select(-NodeType) %>% as.matrix()
    HNet = HNet[EfficHNetNodeInd, EfficHNetNodeInd]

    pred_tar = NodeTypes$NodeName[NodeTypes$NodeType == "predicting"]
    trans_in = NodeTypes$NodeName[NodeTypes$NodeType != "predicting"]
    EfficHNetNodeInd = c(rep(TRUE, length(pred_tar)), Matrix::rowSums(HNet[(length(pred_tar) + 1) : nrow(HNet), (length(pred_tar) + 1) : nrow(HNet)]) != 0)
    HNet = HNet[EfficHNetNodeInd, EfficHNetNodeInd]
    EfficHNetNodeInfo = NodeTypes[EfficHNetNodeInd,]
    NodeTypes = EfficHNetNodeInfo

    trans_in = NodeTypes$NodeName[NodeTypes$NodeType != "predicting"]
    EfficHNetNodeInd = c(colSums(HNet[(length(pred_tar) + 1) : nrow(HNet), 1:length(pred_tar)]) != 0, rep(TRUE, length(trans_in)))
    HNet = HNet[EfficHNetNodeInd, EfficHNetNodeInd]
    EfficHNetNodeInfo = NodeTypes[EfficHNetNodeInd,]
    NodeTypes = EfficHNetNodeInfo
    HNet = apply(HNet, 1, function(x){ind = sample(1:length(x), size = length(x));return(x[ind])})

    theta = as.numeric(NetworkJumpProbability)

    # extract subnetwork from HNet
    NodeType1_index = which(NodeTypes$NodeType == unique(NodeTypes$NodeType)[1])
    NodeType2_index = which(NodeTypes$NodeType == unique(NodeTypes$NodeType)[2])
    NodeType3_index = which(NodeTypes$NodeType == unique(NodeTypes$NodeType)[3])
    (length(NodeType1_index) + length(NodeType2_index) + length(NodeType3_index)) == dim(HNet)[1]

    HNet_AA_pri = HNet[NodeType1_index, NodeType1_index]
    HNet_AB_pri = HNet[NodeType1_index, NodeType2_index]
    HNet_AC_pri = HNet[NodeType1_index, NodeType3_index]

    HNet_A_Nor = HeterogeneousNetworkNormalizing(SimilarityNetwork = HNet_AA_pri, AssociationNetworkA = HNet_AB_pri, AssociationNetworkB = HNet_AC_pri, theta = theta)

    HNet_BA_pri = HNet[NodeType2_index, NodeType1_index]
    HNet_BB_pri = HNet[NodeType2_index, NodeType2_index]
    HNet_BC_pri = HNet[NodeType2_index, NodeType3_index]

    HNet_B_Nor = HeterogeneousNetworkNormalizing(SimilarityNetwork = HNet_BB_pri, AssociationNetworkA = HNet_BA_pri, AssociationNetworkB = HNet_BC_pri, theta = theta)

    HNet_CA_pri = HNet[NodeType3_index, NodeType1_index]
    HNet_CB_pri = HNet[NodeType3_index, NodeType2_index]
    HNet_CC_pri = HNet[NodeType3_index, NodeType3_index]

    HNet_C_Nor = HeterogeneousNetworkNormalizing(SimilarityNetwork = HNet_CC_pri, AssociationNetworkA = HNet_CA_pri, AssociationNetworkB = HNet_CB_pri, theta = theta)

    HNet_Nor = rbind(HNet_A_Nor, HNet_B_Nor, HNet_C_Nor) %>% t()# colSum = 1

    gc()
    return(list(HNet_Nor, NodeTypes))
  }
  transcriptomeList = transcriptomeList
  Rcpp::sourceCpp(paste0(system.file("extdata", package = "NetLCP"), "/MatrixFunctions.cpp"))
  load(system.file("extdata", "prioritization.rda", package = "NetLCP"))
  if(prioType == "lncRNA"){
    HNet_index = LMIM
  }else if(prioType == "circRNA"){
    HNet_index = CMIM
  }else if(prioType == "KEGG"){
    HNet_index = KPMIM
  }else if(prioType == "Reactome"){
    HNet_index = RPMIM
  }else if(prioType == "Wikipathway"){
    HNet_index = WPMIM
  }

  HNet = NormalizeHeterogeneousNetworkGenerator(HeterogeneousNetwork = HNet_index, Transcriptome = transcriptomeList, NetworkJumpProbability = 0.9)
  if(is.character(HNet)){
    return("Sorry, please enlarge your input transcriptome data......")
  }

  seed_score_matrix = matrix(0, nrow = nrow(HNet[[2]]), ncol = nrow(HNet[[2]]))
  diag(seed_score_matrix) = 1

  gc()

  HNet_RWR_Mat = Random_Walk_Restart(HNet[[1]], 0.7, seed_score_matrix)
  HNet_RWR_Mat_svd = svd(HNet_RWR_Mat)

  svd_reduction_num = round(0.1 * nrow(HNet_RWR_Mat))

  singular_value_redution = HNet_RWR_Mat_svd$d[1:svd_reduction_num]
  singular_value_matrix = matrix(0, svd_reduction_num, svd_reduction_num)
  diag(singular_value_matrix) = singular_value_redution

  HNet_new_RWR_matrix = MatrixMultiplication(singular_value_matrix, t(HNet_RWR_Mat_svd$v[ ,1:svd_reduction_num]))
  colnames(HNet_new_RWR_matrix) = HNet[[2]]$NodeName

  Pred_Risk_Score = c()
  pred_tar_vector = HNet_new_RWR_matrix[ ,1:sum(HNet[[2]]$NodeType == "predicting")]

  for(i in 1 : sum(HNet[[2]]$NodeType == "predicting")){
    sum_pred_tar = 0
    for(j in (sum(HNet[[2]]$NodeType == "predicting") + 1) : ncol(HNet_new_RWR_matrix)){
      sum_pred_tar = sum_pred_tar + sum(pred_tar_vector[ ,i] * HNet_new_RWR_matrix[ ,j])  / (sqrt(sum(pred_tar_vector[ ,i] ^ 2)) * sqrt(sum(HNet_new_RWR_matrix[, j] ^ 2)))
    }
    Pred_Risk_Score = append(Pred_Risk_Score, sum_pred_tar)
  }

  PredTransVecInfo = HNet[[2]][HNet[[2]]$NodeType == "predicting", ]
  PredTransVecInfo$RiskScore = Pred_Risk_Score
  PredTransVecInfo = PredTransVecInfo %>% as.data.frame()
  PredTransVecInfo = PredTransVecInfo[order(PredTransVecInfo$RiskScore, decreasing = T), ]
  PredTransVecInfo$Ranking = rank(sort(PredTransVecInfo$RiskScore))
  OfficialName = OFFINAME

  PredTransVecInfo = dplyr::left_join(PredTransVecInfo, OfficialName, by = "NodeName")
  result = PredTransVecInfo
  result$NorRiskScore = sapply(result$RiskScore, function(x){(x - mean(result$RiskScore)) / sd(result$RiskScore)})
  result = result %>% filter(NorRiskScore >= 1)
  result = result[,-c(2,3)]
  gc()

  print("Prioritization finished......")

  if(empiricalPvalue == FALSE){
    return(result)
  }else if(empiricalPvalue == TRUE){
    set.seed(156001)
    print("Computint empirical pvalue, it could be take a while......")
    pb = txtProgressBar(style=3)

    star_time = Sys.time()

    random_result = data.frame(NodeName=NA, NodeType=NA, RiskScore=NA, Ranking=NA)
    miRNA_in_trans = sum(grepl("MIMAT", transcriptomeList))
    mRNA_in_trans = length(transcriptomeList) - miRNA_in_trans

    allmRNA = MRNALIST
    allmiRNA = MIRNALIST
    allmRNA = allmRNA[!allmRNA$mRNAID %in% transcriptomeList]
    allmiRNA = allmiRNA[!allmiRNA$miRNAID %in% transcriptomeList]

    for (k in 1:200){
      transcriptomeList = c(sample(allmRNA$mRNAID, mRNA_in_trans), sample(allmiRNA$miRNAID, miRNA_in_trans))
      HNet = RandomNormalizeHeterogeneousNetworkGenerator(HeterogeneousNetwork = HNet_index, Transcriptome = transcriptomeList, NetworkJumpProbability = 0.9)

      all(round(colSums(HNet[[1]]), 0) == 1)
      gc()

      seed_score_matrix = matrix(0, nrow = nrow(HNet[[2]]), ncol = nrow(HNet[[2]]))
      diag(seed_score_matrix) = 1
      dim(seed_score_matrix)

      gc()

      HNet_RWR_Mat = Random_Walk_Restart(HNet[[1]], 0.7, seed_score_matrix)

      HNet_RWR_Mat_svd = svd(HNet_RWR_Mat)

      svd_reduction_num = round(0.1 * nrow(HNet_RWR_Mat))

      singular_value_redution = HNet_RWR_Mat_svd$d[1:svd_reduction_num]
      singular_value_matrix = matrix(0, svd_reduction_num, svd_reduction_num)
      diag(singular_value_matrix) = singular_value_redution

      HNet_new_RWR_matrix = MatrixMultiplication(singular_value_matrix, t(HNet_RWR_Mat_svd$v[ ,1:svd_reduction_num]))
      colnames(HNet_new_RWR_matrix) = HNet[[2]]$NodeName

      Pred_Risk_Score = c()
      pred_tar_vector = HNet_new_RWR_matrix[ ,1:sum(HNet[[2]]$NodeType == "predicting")]

      for(i in 1 : sum(HNet[[2]]$NodeType == "predicting")){
        sum_pred_tar = 0
        for(j in (sum(HNet[[2]]$NodeType == "predicting") + 1) : ncol(HNet_new_RWR_matrix)){
          sum_pred_tar = sum_pred_tar + sum(pred_tar_vector[ ,i] * HNet_new_RWR_matrix[ ,j])  / (sqrt(sum(pred_tar_vector[ ,i] ^ 2)) * sqrt(sum(HNet_new_RWR_matrix[, j] ^ 2)))
        }
        Pred_Risk_Score = append(Pred_Risk_Score, sum_pred_tar)
      }

      PredTransVecInfo = HNet[[2]][HNet[[2]]$NodeType == "predicting", ]
      PredTransVecInfo$RiskScore = Pred_Risk_Score
      PredTransVecInfo = PredTransVecInfo %>% as.data.frame()
      PredTransVecInfo = PredTransVecInfo[order(PredTransVecInfo$RiskScore, decreasing = T), ]
      PredTransVecInfo$Ranking = rank(sort(PredTransVecInfo$RiskScore))


      random_result = rbind(random_result, PredTransVecInfo)
      gc()

      setTxtProgressBar(pb, k/200)

    }
    end_time = Sys.time()
    close(pb)

    run_time = end_time - star_time
    print(paste0("Empirical pvalue has been computed! running time: ", run_time))

    random_result = random_result[-1,]

    score_nor = sapply(result$NodeName, function(x){
      random_score = random_result$RiskScore[which(random_result$NodeName == x)]
      random_score =  random_score[random_score < quantile(random_score)[4] * 4]
      score_diff = pnorm((result$RiskScore[which(result$NodeName == x)] - mean(random_score) + log2(result$RiskScore[which(result$NodeName == x) + 1])) / sd(random_score), lower.tail = F)
      return(score_diff)
    })

    result$P.value = score_nor
    return(result)
  }
}
