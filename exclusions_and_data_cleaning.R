# Setup -------------------------------------------------------------------
library(tidyverse)
library(stats)
library(parallel)
library(lm.beta)
library(summarytools)
library(psych)
library(GGally)
library(cifti)
library(ciftiTools)
library(ggseg)
library(ggsegGordon)
library(visreg)
library(mgcv)
library(NetworkToolbox)
library(qgraph)
library(reshape2)
source("~/Box/tools/threshold_arbmeasure.R")

# Load demo data ----------------------------------------------------------
exclusions <- read.csv("~/Box/Tooley 01_18_2023 eLABE Requested Data/subjects_all_inclusion_exclusion_01_25_2022.csv")
demographics <- read.csv("~/Box/Tooley 01_18_2023 eLABE Requested Data/Tooley_RC_20230118_updated08JUN23.csv")

demo_data_all <- left_join(demographics, exclusions, by="modid");dim(demo_data_all)

#exclude pre-term babies, MRI_injury exclusion, NICU exclusions, and birthweight exclusion and IRB exclusion
demo_data <- demo_data_all %>% filter(.,combined_exclusion_initial_analyses ==0 & irb_exclusion == 0);dim(demo_data)
#remame PMA
demo_data$PMA_scan <- demo_data$mri_test_pma_scan_dob
#make sex a factor
demo_data$child_sex <- factor(demo_data$child_sex, labels=c("Male","Female"))

#neighborhood decile is a character after birth
demo_data$neighborhood_statedecile_y2 <- as.numeric(demo_data$neighborhood_statedecile_y2)
demo_data$neighborhood_natlcentile_y2 <- as.numeric(demo_data$neighborhood_natlcentile_y2)
demo_data$neighborhood_statedecile_y3 <- as.numeric(demo_data$neighborhood_statedecile_y3)
demo_data$neighborhood_natlcentile_y3 <- as.numeric(demo_data$neighborhood_natlcentile_y3)

# Paths setup -------------------------------------------------------------
ciftiTools.setOption('wb_path', '/Applications/workbench/')
library(rgl) #to use ciftiTools graphics

gordon_parcel_matching <- read.csv("~/Box/tools/parcellations/Gordon_fs_LR/Gordon333Atlas.Parcels.LabelKey.csv")
fc_matrix_dir="~/data/smyser/smyser1/wunder/eLABe/gordon_pconns_plus_atlas_subcortical/full_mats/"
sc_prob_matrix_dir="~/CHPC/mnt/beegfs/scratch/tooley/connectivity/birth/"
sc_prob_matrix_dir="~/10.20.145.4/SMYSER03/smyser3/neonatal/ursula/FSL_probabilistic_connectivity/Y0/" #make it the one on the NIL instead

gordon_networks <- str_split_i(gordon_parcel_matching$label, "_", 3);gordon_networks <- as.factor(gordon_networks)
#medialparietal is cinguloparietal, parieto-occipital is parietooccipital, and #none includes retrosplenial temporal
scale_fill_gordon <- function(...){
  ggplot2:::manual_scale(
    'fill', 
    values = setNames(c(c("#F434F4", "#962996", "#FF2929","#28FE29", "#F6F629", "#A967FD","gray","#FEFED5", 
                                   "#282828", "#29D6D6","#F79029", "#299B9B", "#2929D0")), levels(gordon_networks)), 
                                   ...
  )
}
scale_color_gordon <- function(...){
  ggplot2:::manual_scale(
    'color', 
    values = setNames(c(c("#F434F4", "#962996", "#FF2929","#28FE29", "#F6F629", "#A967FD","gray","#FEFED5", 
                                   "#282828", "#29D6D6","#F79029", "#299B9B", "#2929D0")), levels(gordon_networks)), 
                                   ...
  )
}

#Load symmetric probabilistic SC matrices --------------------------------------------------------
#this sample includes the participants who had V1_b instead of V1_a
sc_subjects <- data.frame(list.files(sc_prob_matrix_dir)); colnames(sc_subjects) <- "modid"
sc_subjects <- unique(sc_subjects);dim(sc_subjects)
filenames <- paste0("~/CHPC/mnt/beegfs/scratch/tooley/connectivity/birth/", sc_subjects$modid, "/", sc_subjects$modid, "_connectome_xdist_20_DISCRETE_NONOVER/fdt_network_matrix")
filenames <- data.frame(filenames)

waytotal_filenames <- paste0("~/CHPC/mnt/beegfs/scratch/tooley/connectivity/birth/", sc_subjects$modid, "/", sc_subjects$modid, "_connectome_xdist_20_DISCRETE_NONOVER/waytotal")

#loop to load SC matrices
remove(big.SC.determ.mat,big.SC.determ.commun.mat)
sc_prob_subj_list <- vector()
missing_columns <-list()
subject_waytotals<- array(NA,c(length(sc_subjects$modid),333))
big.SC.prob.symm.waytotal.norm <- array(NA,c(333,333,length(sc_subjects$modid)))
big.SC.prob.symm.regional.waytotal.norm<- array(NA,c(333,333,length(sc_subjects$modid)))
big.SC.prob.symm.norm <- array(NA,c(333,333,length(sc_subjects$modid)))
big.SC.prob.symm <- array(NA,c(333,333,length(sc_subjects$modid)))
a=0
for (i in 1:dim(demo_data)[1]){ 
  if (file.exists(filenames$filenames[i]))
  {
    print(filenames$filenames[i])
    a=a+1 #if wanted to index without NAs, don't need otherwise
    sc_prob_subj_list[i] <- str_extract(filenames$filenames[i],"MOD....")
    SC.sub.probabilistic <- read.table(filenames$filenames[i], col.names = gordon_parcel_matching$label)
    
    missing_columns[[i]] <- which(colMeans(SC.sub.probabilistic)==0) #which columns are missing from this participant
    
    SC.sub.prob.symm <- data.frame(nrow = 333, ncol=333)
    #make a symmetric matrix, and make identity NA
    for (a in 1:333){
      for (b in 1:333){
        mean <- (SC.sub.probabilistic[a,b] + SC.sub.probabilistic[b,a])/2
        SC.sub.prob.symm[a,b] <-mean
        SC.sub.prob.symm[b,a] <- mean
        SC.sub.prob.symm[a,a] <- NA
      }
    }
    colnames(SC.sub.prob.symm) <- gordon_parcel_matching$label
    
    #normalize to total edge count
    total_edges <- sum(SC.sub.prob.symm, na.rm = T) 
    SC.sub.prob.symm.norm <- SC.sub.prob.symm/total_edges
    
    #normalize by waytotals
    subject_waytotals[i,] <- t(read.table(waytotal_filenames[i]))
    total_waytotal <- sum(subject_waytotals[i,])
    SC.sub.prob.symm.waytotal.norm <- SC.sub.prob.symm/total_waytotal
    
    #normalize each column by waytotals
    SC.sub.prob.symm.regional.waytotal.norm <- sweep(SC.sub.prob.symm, 2, subject_waytotals[i,], FUN = '/')

    big.SC.prob.symm.regional.waytotal.norm[,,i] <- as.matrix(SC.sub.prob.symm.regional.waytotal.norm)
    big.SC.prob.symm.waytotal.norm[,,i] <- as.matrix(SC.sub.prob.symm.waytotal.norm)
    big.SC.prob.symm.norm[,,i] <- as.matrix(SC.sub.prob.symm.norm)
    big.SC.prob.symm[,,i] <- as.matrix(SC.sub.prob.symm)
  }
}
#make a subject list without NAs
sc_prob_subj_list
length(sc_prob_subj_list)
sc_prob_subj_list_nonas <- as.character(na.exclude(sc_prob_subj_list))
length(sc_prob_subj_list_nonas)
#make data without NAs for missing participants
subject_waytotals <- as.data.frame(subject_waytotals)
subject_waytotals$total_waytotal <- rowSums(subject_waytotals)
subject_waytotals$modid <- sc_prob_subj_list
subject_waytotals<- subject_waytotals[!is.na(subject_waytotals$modid),]
big.SC.prob.symm.norm <- big.SC.prob.symm.norm[,,!is.na(sc_prob_subj_list)]
big.SC.prob.symm <- big.SC.prob.symm[,,!is.na(sc_prob_subj_list)]
big.SC.prob.symm.waytotal.norm <- big.SC.prob.symm.waytotal.norm[,,!is.na(sc_prob_subj_list)]
big.SC.prob.symm.regional.waytotal.norm <- big.SC.prob.symm.regional.waytotal.norm[,,!is.na(sc_prob_subj_list)]
dim(big.SC.prob.symm)
save(big.SC.prob.symm.norm,big.SC.prob.symm, big.SC.prob.symm.waytotal.norm,big.SC.prob.symm.regional.waytotal.norm, subject_waytotals, sc_prob_subj_list, sc_prob_subj_list_nonas, file="~/Box/projects/in_progress/struct_funct_neonates/data/clean_data_july_2023/n267_sc_prob_symm_norm_waytotals_matrix_array.RData")
load("~/Box/projects/in_progress/struct_funct_neonates/data/clean_data_july_2023/n267_sc_prob_symm_norm_waytotals_matrix_array.RData")

# Normalize matrix by ROI size (between parcels) -----------------------------------
sc_prob_matrix_dir="~/10.20.145.4/SMYSER03/smyser3/neonatal/ursula/FSL_probabilistic_connectivity/Y0/" #make it the one on the NIL instead
output_dir="~/Box/projects/in_progress/struct_funct_neonates/data/clean_data_july_2023/sc_matrices/birth/probabilistic/"

#loop for SC matrices
sc_prob_subj_list <- vector()
big.SC.prob.symm.norm.roi <- array(NA,c(333,333,length(sc_prob_subj_list_nonas)))
big.SC.prob.symm.norm.roi.norm.edges <- array(NA,c(333,333,length(sc_prob_subj_list_nonas)))
total_edges <- vector()

for (i in 1:length(sc_prob_subj_list_nonas)){
  subject <-  sc_prob_subj_list_nonas[i]
  print(subject)
  if(!file.exists(paste0(output_dir,"/SymmNormRoi/",subject,"_birth_gordon333_fsl_prob_pd_symm_norm_roi.tsv"))){
  long_subject <- sc_subjects$modid[startsWith(sc_subjects$modid, subject)]
  parcel_file <- read_table(paste0(sc_prob_matrix_dir, long_subject, "/",long_subject, "_Gordon_parcels_DISCRETE_NONOVER.txt"), col_names = F)
  SC.sub.prob.symm <- big.SC.prob.symm[,,i]
  SC.sub.prob.symm.norm.roi <- array(NA,c(333,333))
  parc_size <- vector()
  #get parcel volume for a participant
  for (a in 1:333){
    #print(a)
    parcel <- parcel_file$X1[a]
    #parcel_filename <-  gsub("/scratch/tooley/connectivity/",sc_prob_matrix_dir, parcel) #sub out the cluster path for local mounted path
    parcel_filename <-  gsub("/scratch/tooley/connectivity/birth/",sc_prob_matrix_dir, parcel) #for everyone run through after MOD2250_V1_a this works, for earlier subjects, use the line above.
    #cmd = sprintf('~/fsl/share/fsl/bin/fslstats %s -V | cut -d " " -f1', parcel_filename);
    cmd <- capture.output(cat('~/fsl/share/fsl/bin/fslstats', parcel_filename,' -V | cut -d \" \" -f1'))
    sizey <- system(cmd, intern=T);
    parc_size[a] = as.numeric(sizey);
  }
  #normalize
  for (a in 1:333){
    for (b in 1:333){
      if (a == b){
        SC.sub.prob.symm.norm.roi[a,a] <- NA
      } else {
        total_voxels <- parc_size[a] + parc_size[b]
        SC.sub.prob.symm.norm.roi[a,b] <- SC.sub.prob.symm[a,b]/total_voxels
        SC.sub.prob.symm.norm.roi[b,a] <- SC.sub.prob.symm[b,a]/total_voxels
      }
    }
  }
  big.SC.prob.symm.norm.roi[,,i] <- as.matrix(SC.sub.prob.symm.norm.roi)
  write.table(SC.sub.prob.symm.norm.roi, file=paste0(output_dir,"/SymmNormRoi/",subject,"_birth_gordon333_fsl_prob_pd_symm_norm_roi.tsv"), col.names = F, row.names=F)
  
  #normalize to total edge count, since we used the non-edge-normalized values here
  total_edges[i] <- sum(SC.sub.prob.symm, na.rm=T)
  SC.sub.prob.symm.norm.roi.norm.edges <- SC.sub.prob.symm.norm.roi/total_edges[i]
  write.table(SC.sub.prob.symm.norm.roi.norm.edges, file=paste0(output_dir,"/SymmNormRoiNormEdges/",subject,"_birth_gordon333_fsl_prob_pd_symm_norm_roi_norm_edges.tsv"), col.names = F, row.names=F)
  
  big.SC.prob.symm.norm.roi.norm.edges[,,i] <- as.matrix(SC.sub.prob.symm.norm.roi.norm.edges)
  } 
  SC.sub.prob.symm <- big.SC.prob.symm[,,i]
  total_edges[i] <- sum(SC.sub.prob.symm, na.rm=T)
  #write them out
  big.SC.prob.symm.norm.roi[,,i] <- as.matrix(read.table(file=paste0(output_dir,"/SymmNormRoi/",subject,"_birth_gordon333_fsl_prob_pd_symm_norm_roi.tsv")))
  big.SC.prob.symm.norm.roi.norm.edges[,,i] <- as.matrix(read.table(file=paste0(output_dir,"/SymmNormRoiNormEdges/",subject,"_birth_gordon333_fsl_prob_pd_symm_norm_roi_norm_edges.tsv")))
}

dim(big.SC.prob.symm.norm.roi)
dim(big.SC.prob.symm.norm.roi.norm.edges)
total_edges <- data.frame(cbind(sc_prob_subj_list_nonas,total_edges))
total_edges$total_edges <- as.numeric(total_edges$total_edges)
save(big.SC.prob.symm.norm.roi, big.SC.prob.symm.norm.roi.norm.edges,total_edges,sc_prob_subj_list_nonas,big.SC.prob.symm,big.SC.prob.symm.norm, file="~/Box/projects/in_progress/struct_funct_neonates/data/clean_data_july_2023/sc_matrices/birth/probabilistic/n267_SC_array_matrix_normed_ROI_norm_edges.RData" )
load("~/Box/projects/in_progress/struct_funct_neonates/data/clean_data_july_2023/sc_matrices/birth/probabilistic/n267_SC_array_matrix_normed_ROI_norm_edges.RData")

# Consistency-based thresholding of SC matrix ------------------------------------------
load("~/Box/projects/in_progress/struct_funct_neonates/data/clean_data_july_2023/sc_matrices/birth/probabilistic/n267_SC_array_matrix_normed_ROI_norm_edges.RData")
output_dir="~/Box/projects/in_progress/struct_funct_neonates/data/clean_data_july_2023/sc_matrices/birth/probabilistic/"
group.SC.prob.matrix.cv.thresh <- threshold_consistency(big.SC.prob.symm.norm.roi.norm.edges, 0.75)
group.SC.prob.matrix.cv.thresh[1:10,1:10] #look at it makes sure it works

#make all subject matrices CV-thresholded
myfun <- function(mat){
  mat[which(group.SC.prob.matrix.cv.thresh==0)] <- 0
  return(mat)
}
big.SC.prob.symm.thresh.vec <- apply(big.SC.prob.symm.norm.roi.norm.edges, MARGIN=3, FUN=myfun)
big.SC.prob.symm.norm.roi.thresh<- array(big.SC.prob.symm.thresh.vec, dim=dim(big.SC.prob.symm.norm))
dim(big.SC.prob.symm.norm.roi.thresh)
save(big.SC.prob.symm.norm.roi,big.SC.prob.symm.norm.roi.thresh, big.SC.prob.symm.norm.roi.norm.edges,total_edges,sc_prob_subj_list_nonas,big.SC.prob.symm,big.SC.prob.symm.norm, file="~/Box/projects/in_progress/struct_funct_neonates/data/clean_data_july_2023/sc_matrices/birth/probabilistic/n267_SC_array_normed_to_ROI_size_edges_norm_CV75.RData" )

#write out consistency-based thresholding matrices 
for (i in 1:dim(big.SC.prob.symm.norm.roi.thresh)[3]){
  subject=sc_prob_subj_list_nonas[i]
  print(subject)
  mat=big.SC.prob.symm.norm.roi.thresh[,,i]
  write.table(mat, file=paste0(output_dir,"/SymmNormRoiNormEdgesCV75/",subject,"_birth_gordon333_fsl_prob_pd_symm_norm_roi_cv75.csv"), sep = ",", col.names = F, row.names = F)
}

# Figure out exclusions for diffusion data ------------------------------
load("~/Box/projects/in_progress/struct_funct_neonates/data/clean_data_july_2023/sc_matrices/birth/probabilistic/n267_SC_array_normed_to_ROI_size_edges_norm_CV75.RData")
sc_prob_subj_list_nonas

#let's check how many of the diffusion people Jeanette has been excluding are still in my sample
jeanette_exclusion_reasons <- readxl::read_xlsx("~/Box/projects/in_progress/struct_funct_neonates/data/clean_data_july_2023/subj_lists/Ursala_data_check_usability.xlsx")
colnames(jeanette_exclusion_reasons)

intersect(sc_prob_subj_list_nonas, jeanette_exclusion_reasons$`Unique Exclusion IDs`) #I included 12 people that jeanette excluded
intersect(sc_prob_subj_list_nonas, jeanette_exclusion_reasons$`Removed SBREF`) #0 of the SBRef artifacts are in here
intersect(sc_prob_subj_list_nonas, jeanette_exclusion_reasons$`Removed IRB`) #I forgot to exclude the IRB person, but they don't have FC data so should still exclude them
intersect(sc_prob_subj_list_nonas, jeanette_exclusion_reasons$`Removed preterm (<37)`) # some preemies in my list, exclude them
# no NICU, injury, low birthweight, bad T2, BBR fail, or partial sequence in my list
intersect(sc_prob_subj_list_nonas, jeanette_exclusion_reasons$`Removed for Injury`) #No injury
intersect(sc_prob_subj_list_nonas, jeanette_exclusion_reasons$`Removed for NICU stay`) #No NICU
intersect(sc_prob_subj_list_nonas, jeanette_excluson_reasons$`Removed for Low Birthweight`) #No low birthweight
intersect(sc_prob_subj_list_nonas, jeanette_exclusion_reasons$`Removed for bad or no T2`) #No bad or no T2
intersect(sc_prob_subj_list_nonas, jeanette_exclusion_reasons$`Removed for partial Sequence`) #No partial sequence
intersect(sc_prob_subj_list_nonas, jeanette_exclusion_reasons$`Removed for SD = 5/6`) #2 participants who were removed for SD = 5/6 in my list

#exclude the IRB person, along with preemies and SD 5/6 from the sample, how many people with diffusion data do I have?
length(sc_prob_subj_list_nonas)
sc_prob_subj_list_nonas_filt <- sc_prob_subj_list_nonas[!sc_prob_subj_list_nonas %in% jeanette_exclusion_reasons$`Removed for SD = 5/6`];length(sc_prob_subj_list_nonas_filt)
sc_prob_subj_list_nonas_filt <- sc_prob_subj_list_nonas_filt[!sc_prob_subj_list_nonas_filt %in% jeanette_exclusion_reasons$`Removed preterm (<37)`];length(sc_prob_subj_list_nonas_filt)
sc_prob_subj_list_nonas_filt <- sc_prob_subj_list_nonas_filt[!sc_prob_subj_list_nonas_filt %in% jeanette_exclusion_reasons$`Removed IRB`];length(sc_prob_subj_list_nonas_filt)

#n=255 with usable SC data contained in the list 
sc_prob_subj_list_nonas_filt #this is what to use, though the SC arrays are in the order of sc_prob_subj_list_nonas
save(sc_prob_subj_list_nonas_filt,big.SC.prob.symm.norm.roi,big.SC.prob.symm.norm.roi.thresh, big.SC.prob.symm.norm.roi.norm.edges,total_edges,sc_prob_subj_list_nonas,big.SC.prob.symm,big.SC.prob.symm.norm, file="~/Box/projects/in_progress/struct_funct_neonates/data/clean_data_july_2023/sc_matrices/birth/probabilistic/n267_SC_array_normed_to_ROI_size_edges_norm_CV75.RData" )
load("~/Box/projects/in_progress/struct_funct_neonates/data/clean_data_july_2023/sc_matrices/birth/probabilistic/n267_SC_array_normed_to_ROI_size_edges_norm_CV75.RData")
write.csv(sc_prob_subj_list_nonas_filt, "~/Box/projects/in_progress/within_between_network_longitudinal/data/subjLists/n255_sc_birth_healthy_FT_good_data.csv")

# Read in FC data from pconns ---------------------------------------------
#read in FC matrices
files1 <- data.frame(list.files(fc_matrix_dir));files1 <-  filter(files1, grepl("_gordon_parcel_plus_term_N50_eLABe_atlas_subcort.pconn.nii",list.files.fc_matrix_dir.)); 
fc_subj_list <- str_remove(files1$list.files.fc_matrix_dir., "_gordon_parcel_plus_term_N50_eLABe_atlas_subcort.pconn.nii") %>% 
  str_remove(.,"~/data/smyser/smyser1/wunder/eLABe/gordon_pconns_plus_atlas_subcortical/full_mats/") %>% data.frame();colnames(fc_subj_list) <- c("long_modid")
dim(fc_subj_list)
fc_subj_list$modid <- gsub(pattern = "_V1_.", "", fc_subj_list$long_modid)
fc_data <- demo_data %>% filter(modid %in% fc_subj_list$modid) #this should be only the healthy full-term babies without IRB exclusion, n=261
dim(fc_data)

#read only the healthy FT pconns (n=261)
big.FC.mat <- array(NA,c(333,333,length(fc_data$modid)))
for (i in 1:dim(fc_data)[1]){
  subject <-fc_data$modid[i]
  fc_data$long_modid <- fc_subj_list$long_modid[startsWith(fc_subj_list$long_modid, subject)]
  FC.pconn <- cifti::read_cifti(paste0(fc_matrix_dir,fc_data$long_modid[i],"_gordon_parcel_plus_term_N50_eLABe_atlas_subcort.pconn.nii"))
  FC.pconn$data[FC.pconn$data>7] <- NA #eliminate the 7.25s on the diagonal
  big.FC.mat[,,i] <- as.matrix(FC.pconn$data[1:333,1:333])
  fc_data$avg_weight[i] <- mean(FC.pconn$data, na.rm=T)
}
dim(big.FC.mat)
save(big.FC.mat, fc_data, file="~/Box/projects/in_progress/struct_funct_neonates/data/clean_data_july_2023/fc_matrices/n261_pconns_matrix_array.RData")
load("~/Box/projects/in_progress/struct_funct_neonates/data/clean_data_july_2023/fc_matrices/n261_pconns_matrix_array.RData")

#read instead all the pconns that can possible read in (n=319), and save them to filter out afterwards
big.FC.mat <- array(NA,c(333,333,length(fc_subj_list$long_modid)))
for (i in 1:dim(fc_subj_list)[1]){
  FC.pconn <- cifti::read_cifti(paste0(fc_matrix_dir,fc_subj_list$long_modid[i],"_gordon_parcel_plus_term_N50_eLABe_atlas_subcort.pconn.nii"))
  FC.pconn$data[FC.pconn$data>7] <- NA #eliminate the 7.25s on the diagonal
  big.FC.mat[,,i] <- as.matrix(FC.pconn$data[1:333,1:333])
  fc_subj_list$avg_weight[i] <- mean(FC.pconn$data, na.rm=T)
}
dim(big.FC.mat)
big.FC.mat <- big.FC.mat[,,1:319]
save(big.FC.mat, fc_subj_list, file="~/Box/projects/in_progress/struct_funct_neonates/data/clean_data_july_2023/fc_matrices/n319_pconns_matrix_array.RData")

# Visualize avg FC matrix -------------------------------------------------
library(mdpeer)
library(viridis)
#average FC matrix for visualization
avg_fc_mat <- apply(big.FC.mat, c(1,2), FUN=mean, na.rm=T); avg_fc_mat[is.na(avg_fc_mat)] <- 0
#reorder by labels?
rownames(avg_fc_mat) <- gordon_networks; colnames(avg_fc_mat) <- gordon_networks
m <- avg_fc_mat
m <- m[, order(as.factor(colnames(m)))]
m <- m[order(as.factor(rownames(m))),]
avgmatrix <- as.matrix(m);rownames(avgmatrix) <- NULL;colnames(avgmatrix) <- NULL
heatmap(avgmatrix, Rowv = NA, Colv = NA, col = palf(20), labRow =  as.factor(rownames(avg_fc_mat))[order(as.factor(rownames(avg_fc_mat)))], revC = T)
palf <- colorRampPalette(turbo(10), bias = 1)
plot(vizu.mat(avgmatrix, title="Average FC Matrix", legend.title = "z", geom_tile.colour = NA,fill.scale.limits = c(min(avg_fc_mat), max(avg_fc_mat)),colors.palette =palf(20))) +
  theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_colorbar(ticks.colour = NA))

# Visualize avg SC matrix -------------------------------------------------
big.SC.prob.symm.norm.roi.norm.edges[big.SC.prob.symm.norm.roi.norm.edges==0] <- NA
avg_sc_mat <- apply(big.SC.prob.symm.norm.roi.norm.edges, c(1,2), FUN=mean, na.rm=T);#avg_sc_mat[is.na(avg_sc_mat)] <- 0
#reorder by labels?
rownames(avg_sc_mat) <- gordon_networks; colnames(avg_sc_mat) <- gordon_networks
m <- avg_sc_mat
m <- m[, order(as.factor(colnames(m)))]
m <- m[order(as.factor(rownames(m))),]
avgmatrix <- as.matrix(m);rownames(avgmatrix) <- NULL;colnames(avgmatrix) <- NULL
palf <- colorRampPalette(turbo(10), bias = 1.7)
avgmatrix[is.na(avgmatrix)] <- 0
plot(vizu.mat(avgmatrix, title="Average SC Matrix", legend.title = "z",colors.palette =palf(10), geom_tile.colour = NA)) +
  theme(plot.title = element_text(hjust = 0.5)) + guides(fill=guide_colorbar(ticks.colour = NA))

