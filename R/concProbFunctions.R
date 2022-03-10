countInversions <- function(y, w = NULL, nu = 0, verb)
{
  result <- countInversionsCpp(y, w = w, nu = nu, verb)
  return(result);
}

concProb_bin_bruteForce <- function(obs, pred, w = NULL, tiesPred = FALSE, tiesObs = FALSE){
  if(is.null(w)) w <- rep(1,length(obs))
  if(tiesObs) tiesPred <- T
  
  if(!((length(w) == length(obs)) || (length(w) == length(pred)))) stop('obs, pred and w (if defined) should have the same length.')
  if(any(w < 0)) stop('all weights in w should be positive.')
    
  a <- Sys.time()
  ind1 <- obs == 1
  ind0 <- obs == 0
  pred1 <- pred[ind1]
  pred0 <- pred[ind0]
  w1 <- w[ind1]
  w0 <- w[ind0]
  
  conc <- 0
  nTiesObs <- 0
  nTiesObsPred <- 0
  tot<- 0
  nTiesPred <- 0
  
  if(length(pred0) < length(pred1)){
    for(i in 1:length(pred0)){
      concInd <- which(pred1 > pred0[i])
      conc <- conc + sum(w1[concInd]) * w0[i]
      
      totInd <- which(pred0[i] != pred1)
      tot <- tot + sum(w1[totInd]) * w0[i]
      if(tiesPred){
        totIndTie <- (1:length(pred1))[!((1:length(pred1)) %in% totInd)] 
        nTiesPred <- nTiesPred + sum(w1[totIndTie]) * w0[i]
      } 
    }
  }else{
    for(i in 1:length(pred1)){
      concInd <- which(pred0 < pred1[i])
      conc <- conc + sum(w0[concInd]) * w1[i]
      tot.Ind <- which(pred0 != pred1[i])
      tot <- tot + sum(w0[tot.Ind]) * w1[i]
      if(tiesPred){
        totIndTie <- (1:length(pred0))[! (1:length(pred0)) %in% tot.Ind ]
        nTiesPred <- nTiesPred + sum(w0[totIndTie]) * w1[i]
      } 
    }
  }
  
  if(tiesObs){
    combMat0 <- combn(w0, 2)
    combMat1 <- combn(w1, 2)
    nTiesObs <- sum(combMat0[1,] * combMat0[2,]) + sum(combMat1[1,] * combMat1[2,])
    # #Code below if one wants a difference between ties in observations and ties in observations & predictions
    # group <- vec_group_loc(pred0)
    # numbUniqTies <- which(lengths(group$loc) > 1) #We know already the number of ties
    # if(length(numbUniqTies) > 0){
    #   for(i in 1:length(numbUniqTies)){
    #     indTie <- group$loc[[numbUniqTies[i]]] #already ordered
    #     combMat0Tie <- combn(w0[indTie], 2)
    #     nTiesObsPred <- nTiesObsPred + sum(combMat0Tie[1,] * combMat0Tie[2,])
    #   }
    # }
    # 
    # group <- vec_group_loc(pred1)
    # numbUniqTies <- which(lengths(group$loc)>1) #We know already the number of ties
    # if(length(numbUniqTies)>0){
    #   for(i in 1:length(numbUniqTies)){
    #     indTie <- group$loc[[numbUniqTies[i]]] #already ordered
    #     combMat1Tie <- combn(w1[indTie], 2)
    #     nTiesObsPred <- nTiesObsPred + sum(combMat1Tie[1,] * combMat1Tie[2,])
    #   }
    # }
    # 
    # nTiesObs <- nTiesObs - nTiesObsPred
  }
  b <- Sys.time()
  return(list(concProb = (2 * nb.conc + nb.ties.pred + nb.ties.obs + nb.ties.pred.obs) / (2 * (nb.conc + nb.disc + nb.ties.obs + nb.ties.pred + nb.ties.pred.obs)),
              time = b - a, nb.conc = nb.conc, nb.disc = nb.disc, nb.ties.pred = nb.ties.pred, nb.ties.obs = nb.ties.obs, nb.ties.pred.obs = nb.ties.pred.obs))
}

concProb_bin_fast <- function(obs, pred, w = NULL, tiesPred = FALSE, tiesObs = FALSE){
  if(is.null(w)) w <- rep(1,length(obs))
  if(tiesObs) tiesPred <- T
  
  if(!((length(w) == length(obs)) || (length(w) == length(pred)))) stop('obs, pred and w (if defined) should have the same length.')
  if(any(w < 0)) stop('all weights in w should be positive.')
  
  a <- Sys.time()
  # First order the observations according to the predictions
  order.o <- order(pred)
  obs.o <- obs[order.o] # order the observed (0,1) based on predicted
  obs.o1 <- obs.o == 1
  obs.o0 <- obs.o == 0
  nb.disc <- nb.conc <- 0
  # #each value of 1 in obs.o causes a number of non-concordant pairs
  # #equal to the number of zeroes to the right of that 1:
  # nb.disc <- as.numeric(sum(rev(cumsum(rev(obs.o == 0)))[which(obs.o == 1)]))
 
  #with weights, this reduces to:
  w.sel <- cumsum(rev(w[obs.o0]))
  nb.zero <- rev(cumsum(rev(obs.o0)))[obs.o1]
  nb.zero <- nb.zero[nb.zero > 0]
  w.sel <- w.sel[nb.zero]
  if(length(w.sel) > 0 ) nb.disc <- w[obs.o1][1:length(w.sel)]%*%w.sel
  
  # #conversely, each value of 0 in obs.o causes a number of concordant
  # #pairs equal to the number of ones to the right of that 0
  # nb.conc  <- as.numeric(sum(rev(cumsum(rev(obs.o == 1)))[which(obs.o == 0)]))
  
  #with weights, this reduces to:
  w.sel <- cumsum(rev(w[obs.o1]))
  nb.one <- rev(cumsum(rev(obs.o1)))[obs.o0]
  nb.one <- nb.one[nb.one > 0]
  w.sel <- w.sel[nb.one]
  if(length(w.sel) > 0 ) nb.conc <- w[obs.o0][1:length(w.sel)]%*%w.sel
  
  # check for duplicates:
  no.nb.conc <- as.numeric(0); no.nb.disc <- as.numeric(0)
  nb.ties.pred <- as.numeric(0);nb.ties.obs <- as.numeric(0); nb.ties.pred.obs <- as.numeric(0)
  group <- vec_group_loc(pred)
  numbUniqTies <- which(lengths(group$loc) > 1) #We know already the number of ties
  if(length(numbUniqTies) > 0){
    for(i in 1:length(numbUniqTies)){
      # determine number of concordant pairs and number of discordant pairs that are not valid
      indTie <- group$loc[[numbUniqTies[i]]]
      obsTie <- obs[indTie]
      
      # # no weights
      # no.nb.disc <- no.nb.disc + sum(rev(cumsum(rev(obsTie == 0)))[which(obsTie == 1)])
      # no.nb.conc  <- no.nb.conc + sum(rev(cumsum(rev(obsTie == 1)))[which(obsTie == 0)])
      # if(tiesObs){
      #   nb.ties.pred.obs <- nb.ties.pred.obs + choose(sum(obsTie == 0), 2) + choose(sum(obsTie == 1), 2)
      # }
      
      obsTie0 <- obsTie == 0
      obsTie1 <- obsTie == 1
      wTie <- w[indTie]
      
      w.selTie <- cumsum(rev(wTie[obsTie0]))
      nb.zero <- rev(cumsum(rev(obsTie0)))[obsTie1]
      nb.zero <- nb.zero[nb.zero > 0]
      w.selTie <- w.selTie[nb.zero]
      if(length(w.selTie)>0) no.nb.disc <- no.nb.disc + wTie[obsTie1][1:length(w.selTie)]%*%w.selTie
      
      w.selTie <- cumsum(rev(wTie[obsTie1]))
      nb.one <- rev(cumsum(rev(obsTie1)))[obsTie0]
      nb.one <- nb.one[nb.one > 0]
      w.selTie <- w.selTie[nb.one]
      if(length(w.selTie)>0) no.nb.conc <- no.nb.conc + wTie[obsTie0][1:length(w.selTie)]%*%w.selTie
    }
    nb.conc <- nb.conc - no.nb.conc
    nb.disc <- nb.disc - no.nb.disc
  }
  
  if(tiesPred){
    nb.ties.pred <- no.nb.disc + no.nb.conc
  }
  
  if(tiesObs){
    # #no weights
    # nb.ties.obs <- choose(sum(obs == 1), 2) + choose(sum(obs == 0), 2) - nb.ties.pred.obs
    
    #with weights
    combMat0 <- combn(w[obs.o0], 2)
    combMat1 <- combn(w[obs.o1], 2)
    nb.ties.obs <- sum(combMat0[1,] * combMat0[2,]) + sum(combMat1[1,] * combMat1[2,])
    
    # #Code below if one wants a difference between ties in observations and ties in observations & predictions
    # group <- vec_group_loc(pred0)
    # numbUniqTies <- which(lengths(group$loc) > 1) #We know already the number of ties
    # if(length(numbUniqTies) > 0){
    #   for(i in 1:length(numbUniqTies)){
    #     indTie <- group$loc[[numbUniqTies[i]]] #already ordered
    #     combMat0Tie <- combn(w0[indTie], 2)
    #     nb.ties.pred.obs <- nb.ties.pred.obs + sum(combMat0Tie[1,] * combMat0Tie[2,])
    #   }
    # }
    # 
    # group <- vec_group_loc(pred1)
    # numbUniqTies <- which(lengths(group$loc)>1) #We know already the number of ties
    # if(length(numbUniqTies)>0){
    #   for(i in 1:length(numbUniqTies)){
    #     indTie <- group$loc[[numbUniqTies[i]]] #already ordered
    #     combMat1Tie <- combn(w1[indTie], 2)
    #     nb.ties.pred.obs <- nb.ties.pred.obs + sum(combMat1Tie[1,] * combMat1Tie[2,])
    #   }
    # }
    # 
    # nb.ties.obs <- nb.ties.obs - nb.ties.pred.obs
  }
  
  b <- Sys.time()
  # Using the formula of the C-index
  return(list(concProb = (2 * nb.conc + nb.ties.pred + nb.ties.obs + nb.ties.pred.obs) / (2 * (nb.conc + nb.disc + nb.ties.obs + nb.ties.pred + nb.ties.pred.obs)),
              time = b - a, nb.conc = nb.conc, nb.disc = nb.disc, nb.ties.pred = nb.ties.pred, nb.ties.obs = nb.ties.obs, nb.ties.pred.obs = nb.ties.pred.obs))
}

concProb_cont_fast <- function(obs, pred, w = NULL, nu = 0, tiesPred = FALSE, tiesObs = FALSE, verb = FALSE){
  if(is.null(w)) w <- rep(1,length(obs))
  if(tiesObs) tiesPred <- T
  
  if(!((length(w) == length(obs)) || (length(w) == length(pred)))) stop('obs, pred and w (if defined) should have the same length.')
  if(any(w < 0)) stop('all weights in w should be positive.')
  
  a <- Sys.time()
  order.o  <- order(pred)
  countInv <- countInversions(obs[order.o], w = w[order.o], nu = nu, verb = verb)
  nb.disc <- countInv[1]
  nb.conc  <- countInv[2]
  if(tiesObs){
    nb.ties.obs <- countInv[3]
  }else{
    nb.ties.obs <- 0
  }
  no.nb.conc <- 0; no.nb.disc <- 0; nb.ties.pred <- 0; nb.ties.pred <- 0; nb.ties.pred.obs <- 0;
  group <- vec_group_loc(pred)
  numbUniqTies <- which(lengths(group$loc)>1) #We know already the number of ties
  if(length(numbUniqTies)>0){
    for(i in 1:length(numbUniqTies)){
      # determine number of concordant pairs and number of discordant pairs that are not valid
      indTie <- group$loc[[numbUniqTies[i]]] #already ordered
      obsTie <- obs[indTie]
      wTie <- w[indTie]
      countInvTie <- countInversions(obsTie, w = wTie, nu = nu, verb = verb)
      no.nb.disc <- no.nb.disc + countInvTie[1]
      no.nb.conc  <- no.nb.conc + countInvTie[2]
      if(tiesObs) nb.ties.pred.obs <- nb.ties.pred.obs + countInvTie[3]
    }
    nb.conc <- nb.conc - no.nb.conc
    nb.disc <- nb.disc - no.nb.disc
    if(tiesObs){
      nb.ties.obs <- nb.ties.obs - nb.ties.pred.obs
    }
  }
  
  if(tiesPred){
    nb.ties.pred <- no.nb.disc + no.nb.conc
  }
  b <- Sys.time()
  # Using the formula of the C-index
  return(list(concProb = (2 * nb.conc + nb.ties.pred + nb.ties.obs + nb.ties.pred.obs) / (2 * (nb.conc + nb.disc + nb.ties.obs + nb.ties.pred + nb.ties.pred.obs)),
              time = b - a, nb.conc = nb.conc, nb.disc = nb.disc, nb.ties.pred = nb.ties.pred, nb.ties.obs = nb.ties.obs, nb.ties.pred.obs = nb.ties.pred.obs))
  #return(c(nb.conc, nb.disc, nb.ties.pred, nb.ties.obs, nb.ties.pred.obs))
}

concProb_cont_bruteForce <- function(obs, pred, w = NULL, nu = 0, tiesPred = FALSE, tiesObs = FALSE){
  if(is.null(w)) w <- rep(1,length(obs))
  if(tiesObs) tiesPred <- T
  
  if(!((length(w) == length(obs)) || (length(w) == length(pred)))) stop('obs, pred and w (if defined) should have the same length.')
  if(any(w < 0)) stop('all weights in w should be positive.')
  
  a <- Sys.time()
  nb.conc <- nb.disc <- nb.ties.pred <- nb.ties.obs <- nb.ties.pred.obs <- 0
  for (i in 1:length(obs)) {
    goodInd <- which(obs > obs[i] + nu)
    goodPred <- pred[goodInd]
    goodW <- w[goodInd]
    nb.conc <- nb.conc + sum(goodW[goodPred > pred[i]]) * w[i]
    nb.disc <- nb.disc + sum(goodW[goodPred < pred[i]]) * w[i]
    if(tiesPred){
      tiesInd <- which(goodPred == pred[i])
      nb.ties.pred <- nb.ties.pred + sum(goodW[tiesInd]) * w[i]
    } 
    if(tiesObs){
      tiesInd <- which(obs == obs[i] + nu)
      if(nu == 0) tiesInd <- tiesInd[tiesInd > i]
      if(length(tiesInd) > 0 ){
        tiesPredSel <- pred[tiesInd]
        tiesW <- w[tiesInd]
        nb.ties.pred.obsTemp <- 0
        predObsInd <- which(tiesPredSel == pred[i])
        if(length(predObsInd) > 0) nb.ties.pred.obsTemp <- sum(tiesW[predObsInd]) * w[i] 
        nb.ties.pred.obs <- nb.ties.pred.obs + nb.ties.pred.obsTemp
        nb.ties.obs <- nb.ties.obs + sum(tiesW) * w[i] - nb.ties.pred.obsTemp
      }
    }
  }
  b <- Sys.time()
  return(list(concProb = (2 * nb.conc + nb.ties.pred + nb.ties.obs + nb.ties.pred.obs) / (2 * (nb.conc + nb.disc + nb.ties.obs + nb.ties.pred + nb.ties.pred.obs)),
              time = b - a, nb.conc = nb.conc, nb.disc = nb.disc, nb.ties.pred = nb.ties.pred, nb.ties.obs = nb.ties.obs, nb.ties.pred.obs = nb.ties.pred.obs))
}
