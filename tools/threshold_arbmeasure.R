threshold_consistency <- function(Ws, prop, ...) {
  # %   Inputs: Ws,     N-by-N-by-M group of M weighted connectivity matrices
  # %           p,      proportion of weights to preserve
  # %                       range:  p=1 (all weights preserved) to
  # %                               p=0 (no weights removed)
  
  # Adapted from https://github.com/breakspear/threshold-consist and GB code at: 
  # https://github.com/PennBBL/baumNetCoupling/blob/master/GLB_threshold_consistency.m
  
  Wmean=apply(Ws, c(1,2), mean) # group mean connectivity matrix
  Wcv=apply(Ws,c(1,2), sd)/Wmean # coefficient of variation across the group
  W_thr=threshold_arbmeasure(Wmean,-Wcv,prop); # threshold by low CV = high consistency
  return(W_thr)
}

threshold_arbmeasure <- function(W,
                                 M,
                                 p, 
                                 copy=F){
  #adapted from threshold_proportional here: 
  #https://github.com/ncullen93/bctR/blob/b3c1c1bdb76b1a6d94020497b478a292241264a0/R/Utility.R
  
  #' This function "thresholds" the connectivity matrix by preserving a
  #' proportion p (0<p<1) of the strongest weights based on another matrix M. All other weights, and
  #' all weights on the main diagonal (self-self connections) are set to 0.
  #' 
  #' Note: For data w/ negative numbers, we consider Absolute Value of weights.
  #' 
  #' @param W : an R Matrix - Weighted Connectivity Matrix
  #' @param M : an R Matrix - Matrix with values to threshold by
  #' @param p : a float - Proportional Weight Threshold
  #' @param copy : a boolean - Whether to modify in place or not
  #' 
  #' @return W : Threshold Connectivity Matrix
  
  stopifnot(p > 0, p <= 1)
  #if (copy) W <- W
  n <- nrow(W) # number of nodes
  diag(W) <- 0 # clear the diagonal, modifies in place
  
  # if symmetric matrix, ensure symmetry is preserved
  s.flag <- F
  ud <- 1
  if (isSymmetric(W)) {
    s.flag <- T
    ud <- 2
    W[lower.tri(W)] <- 0
  }
  
  #I <- order(abs(W),decreasing=T) # sort indices by value magnitude
  I <- order(M, decreasing = T) #sort by M
  en <- as.integer(round((n*n-n) * p/ud))+1 # num. links to preserve
  W[I[en:length(I)]] <- 0 # remove the smallest links
  
  if (s.flag) W <- W + t(W) # add back the lower triangle
  
  return(W)
}
