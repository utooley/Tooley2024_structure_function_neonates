sc_fc_vect <- function(fc_matrix_or_vector, sc_matrix_or_vector, subj_list_fc_sc, subj_order_fc, subj_order_sc,nregions){
  # fc_matrix_or_vector: an n x n x m 3D array with n regions and m subjects, or a m x (nxn) vector of subjects by edges
  # sc_matrix_or_vector: an n x n x m 3D array with n regions and m subjects, or a m x (nxn) vector of subjects by edges
  # subj_order_fc: order of subjects in FC matrix dimension m
  # subj_order_sc: order of subjects in SC matrix dimension m
  # subj_list_fc_sc: a list of subjects to run sc-fc coupling for
  sc.fc.coupling <- matrix(nrow=length(subj_list_fc_sc), ncol=1)
  if(missing(nregions)) {
    nregions = 333
  } else {
    nregions = nregions
  }
  for (i in 1:length(subj_list_fc_sc)){
    subject <- subj_list_fc_sc[i]
    fc_index <- which(subj_order_fc==subject)
    sc_index <- which(subj_order_sc==subject)
    
    func_mat <- fc_matrix_or_vector[,,fc_index] 
    struc_mat <- sc_matrix_or_vector[,,sc_index]
    #make sure of how data is input if not a 3D array
    
    #make sure identity is NA in both
    for (j in 1:nregions){
      struc_mat[struc_mat ==0] <- NA
      #func_mat[struc_mat ==0] <- NA
      struc_mat[j,j] <- NA
      func_mat[j,j] <- NA #make identity NA
    } 
    sc.fc.coupling[i,1] <- cor(c(func_mat), c(struc_mat),use = "na.or.complete", method="spearman")
  }
  scfc <- list("mean.subject.sc.fc"= sc.fc.coupling)
  return(scfc)
}

sc_fc_coupl_func <- function(fc_matrix_or_vector, sc_matrix_or_vector, subj_list_fc_sc, subj_order_fc, subj_order_sc, nregions){
  # fc_matrix_or_vector: an n x n x m 3D array with n regions and m subjects, or a m x (nxn) vector of subjects by edges
  # sc_matrix_or_vector: an n x n x m 3D array with n regions and m subjects, or a m x (nxn) vector of subjects by edges
  # subj_order_fc: order of subjects in FC matrix dimension m
  # subj_order_sc: order of subjects in SC matrix dimension m
  # subj_list_fc_sc: a list of subjects to run sc-fc coupling for
  # nregions: number of regions n
  if(missing(nregions)) {
    nregions = 333
  } else {
    nregions = nregions
  }
  sc.fc.coupling <- matrix(nrow = length(subj_list_fc_sc), ncol = nregions)
  sc.non.zero.edges <-matrix(nrow = length(subj_list_fc_sc), ncol = nregions)
  
  for (i in 1:length(subj_list_fc_sc)){
    subject <- subj_list_fc_sc[i]
    fc_index <- which(subj_order_fc==subject)
    sc_index <- which(subj_order_sc==subject)
    
    func_mat <- fc_matrix_or_vector[,,fc_index] 
    struc_mat <- sc_matrix_or_vector[,,sc_index]
    #make sure of how data is input if not a 3D array
    
    #make sure identity is NA in both
    for (j in 1:nregions){
      struc_mat[struc_mat ==0] <- NA
      #func_mat[struc_mat ==0] <- NA
      struc_mat[j,j] <- NA
      func_mat[j,j] <- NA #make identity NA
    } 
    for (l in 1:nregions){
      sc.fc.coupling[i,l] <- cor(func_mat[,l], struc_mat[,l],use = "na.or.complete", method="spearman") #communicability
      sc.non.zero.edges[i,l] <- sum(!is.na(struc_mat[,l]))
    }
  }
  mean.regional.sc.fc <- colMeans(sc.fc.coupling, na.rm = T)
  mean.subject.sc.fc <- rowMeans(sc.fc.coupling, na.rm = T)
  scfc <- list("mean.regional.sc.fc" = mean.regional.sc.fc, "mean.subject.sc.fc"= mean.subject.sc.fc, " sc.fc.coupling" = sc.fc.coupling,
               "sc.non.zero.edges" = sc.non.zero.edges)
  return(scfc)
}

sc_fc_coupl_func_hemi <- function(fc_matrix_or_vector, sc_matrix_or_vector, subj_list_fc_sc, subj_order_fc, subj_order_sc, nregions, hemi_index){
  # fc_matrix_or_vector: an n x n x m 3D array with n regions and m subjects, or a m x (nxn) vector of subjects by edges
  # sc_matrix_or_vector: an n x n x m 3D array with n regions and m subjects, or a m x (nxn) vector of subjects by edges
  # subj_order_fc: order of subjects in FC matrix dimension m
  # subj_order_sc: order of subjects in SC matrix dimension m
  # subj_list_fc_sc: a list of subjects to run sc-fc coupling for
  # nregions: number of regions n
  # hemi_index: a vector of 1s and 2s indicating the two hemispheres
  if(missing(nregions)) {
    nregions = 333
  } else {
    nregions = nregions
  }
  sc.fc.coupling <- matrix(nrow = length(subj_list_fc_sc), ncol = nregions)
  sc.non.zero.edges <-matrix(nrow = length(subj_list_fc_sc), ncol = nregions)
  hemi.sc.fc.coupling<-matrix(nrow = length(subj_list_fc_sc), ncol = nregions)
  hemi.alt.sc.fc.coupling<-matrix(nrow = length(subj_list_fc_sc), ncol = nregions)
  
  for (i in 1:length(subj_list_fc_sc)){
    subject <- subj_list_fc_sc[i]
    fc_index <- which(subj_order_fc==subject)
    sc_index <- which(subj_order_sc==subject)
    
    func_mat <- fc_matrix_or_vector[,,fc_index] 
    struc_mat <- sc_matrix_or_vector[,,sc_index]
    #make sure of how data is input if not a 3D array
    
    #make sure identity is NA in both
    for (j in 1:nregions){
      struc_mat[struc_mat ==0] <- NA
      #func_mat[struc_mat ==0] <- NA
      struc_mat[j,j] <- NA
      func_mat[j,j] <- NA #make identity NA
    } 
    for (l in 1:nregions){
      sc.fc.coupling[i,l] <- cor(func_mat[,l], struc_mat[,l],use = "na.or.complete", method="spearman") #communicability
      sc.non.zero.edges[i,l] <- sum(!is.na(struc_mat[,l]))
      #within-hemisphere only structure-function coupling
      regions_system <- which(hemi_index == hemi_index[l])
      hemi.sc.fc.coupling[i,l] <- cor(func_mat[regions_system,l], struc_mat[regions_system ,l],use = "na.or.complete", method="spearman")
      #between-hemisphere
      regions_not_system <- which(hemi_index != hemi_index[l])
      hemi.alt.sc.fc.coupling[i,l] <- cor(func_mat[regions_not_system,l], struc_mat[regions_not_system ,l],use = "na.or.complete", method="spearman") 
    }
  }
  mean.regional.sc.fc <- colMeans(sc.fc.coupling, na.rm = T)
  mean.subject.sc.fc <- rowMeans(sc.fc.coupling, na.rm = T)
  mean.regional.hemi.sc.fc <- colMeans(hemi.sc.fc.coupling, na.rm = T)
  mean.subject.hemi.sc.fc <- rowMeans(hemi.sc.fc.coupling, na.rm = T)
  mean.regional.hemi.alt.sc.fc <- colMeans(hemi.alt.sc.fc.coupling, na.rm = T)
  mean.subject.hemi.alt.sc.fc <- rowMeans(hemi.alt.sc.fc.coupling, na.rm = T)
  scfc <- list("mean.regional.sc.fc" = mean.regional.sc.fc, "mean.subject.sc.fc"= mean.subject.sc.fc, " sc.fc.coupling" = sc.fc.coupling,
               "sc.non.zero.edges" = sc.non.zero.edges,
               "mean.regional.hemi.sc.fc" = mean.regional.hemi.sc.fc, "mean.subject.hemi.sc.fc"= mean.subject.hemi.sc.fc, " sc.fc.coupling.hemi" = hemi.sc.fc.coupling,
               "mean.regional.other.hemi.sc.fc" =  mean.regional.hemi.alt.sc.fc, "mean.subject.other.hemi.sc.fc"= mean.subject.hemi.alt.sc.fc, " sc.fc.coupling.other.hemi" = hemi.alt.sc.fc.coupling)
  return(scfc)
}
#testing
# fc_matrix_or_vector <- big.FC.mat
# sc_matrix_or_vector <- big.SC.prob.symm.norm.roi.thresh
# subj_list_fc_sc <- subj_list_fc_sc
# subj_order_fc <- fc_data$modid
# subj_order_sc <- sc_prob_subj_list_nonas
# nregions <- 324
# system_assignments <- gordon_networks_324

sc_fc_coupl_func_sys <- function(fc_matrix_or_vector, sc_matrix_or_vector, subj_list_fc_sc, subj_order_fc, subj_order_sc, nregions, system_assignments) {
  # fc_matrix_or_vector: an n x n x m 3D array with n regions and m subjects, or a m x (nxn) vector of subjects by edges
  # sc_matrix_or_vector: an n x n x m 3D array with n regions and m subjects, or a m x (nxn) vector of subjects by edges
  # subj_order_fc: order of subjects in FC matrix dimension m
  # subj_order_sc: order of subjects in SC matrix dimension m
  # subj_list_fc_sc: a list of subjects to run sc-fc coupling for
  # nregions: number of regions n
  # system_assignments: a vector of system assignments for each region
  
  if(missing(nregions)) {
    nregions = 333
  } else {
    nregions = nregions
  }
  sc.fc.coupling <- matrix(nrow = length(subj_list_fc_sc), ncol = nregions)
  sc.non.zero.edges <-matrix(nrow = length(subj_list_fc_sc), ncol = nregions)
  system.sc.fc.coupling <-matrix(nrow = length(subj_list_fc_sc), ncol = nregions)
  system.sc.fc.coupling.between <-matrix(nrow = length(subj_list_fc_sc), ncol = nregions)
  
  for (i in 1:length(subj_list_fc_sc)){
    subject <- subj_list_fc_sc[i]
    fc_index <- which(subj_order_fc==subject)
    sc_index <- which(subj_order_sc==subject)
    
    func_mat <- fc_matrix_or_vector[,,fc_index] 
    struc_mat <- sc_matrix_or_vector[,,sc_index]
    #make sure of how data is input if not a 3D array
    
    #make sure identity is NA in both
    for (j in 1:nregions){
      struc_mat[struc_mat ==0] <- NA
      #func_mat[struc_mat ==0] <- NA
      struc_mat[j,j] <- NA
      func_mat[j,j] <- NA #make identity NA
    } 
    for (l in 1:nregions){
      sc.fc.coupling[i,l] <- cor(func_mat[,l], struc_mat[,l],use = "na.or.complete", method="spearman") #communicability
      sc.non.zero.edges[i,l] <- sum(!is.na(struc_mat[,l]))
      #within-system only sc-fc coupling
      regions_system <- which(system_assignments == system_assignments[l])
      regions_not_system <- which(system_assignments != system_assignments[l])
      system.sc.fc.coupling[i,l] <- cor(func_mat[regions_system,l], struc_mat[regions_system ,l],use = "na.or.complete", method="spearman") #within-network only
      system.sc.fc.coupling.between[i,l] <- cor(func_mat[regions_not_system,l], struc_mat[regions_not_system,l],use = "na.or.complete", method="spearman") #between other networks only
    }
  }
  mean.regional.sc.fc <- colMeans(sc.fc.coupling, na.rm = T)
  mean.subject.sc.fc <- rowMeans(sc.fc.coupling, na.rm = T)
  mean.regional.within.system.sc.fc <- colMeans(system.sc.fc.coupling, na.rm = T)
  mean.subject.within.system.sc.fc <- rowMeans(system.sc.fc.coupling, na.rm = T)
  mean.regional.between.system.sc.fc <- colMeans(system.sc.fc.coupling.between, na.rm = T)
  mean.subject.between.system.sc.fc <- rowMeans(system.sc.fc.coupling.between, na.rm = T)
  scfc <- list("mean.regional.sc.fc" = mean.regional.sc.fc, "mean.subject.sc.fc"= mean.subject.sc.fc, " sc.fc.coupling" = sc.fc.coupling,
               "sc.non.zero.edges" = sc.non.zero.edges, 
               "mean.regional.within.sc.fc"=mean.regional.within.system.sc.fc, "mean.subject.within.system.sc.fc"=mean.subject.within.system.sc.fc, "sc.fc.within.system"=system.sc.fc.coupling,
               "mean.regional.between.sc.fc"=mean.regional.between.system.sc.fc, "mean.subject.between.system.sc.fc"=mean.subject.between.system.sc.fc,"sc.fc.between.system"=system.sc.fc.coupling.between)
  return(scfc)
}
