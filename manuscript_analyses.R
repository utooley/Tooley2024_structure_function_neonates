# Setup ------------------------------------------------------------------
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
source("~/Box/tools/scfc_coupling_fxn.R")
source("~/Box/tools/raincloud.R") #from https://github.com/RainCloudPlots/RainCloudPlots/blob/master/tutorial_R/R_rainclouds.R

# Paths setup -------------------------------------------------------------
ciftiTools.setOption('wb_path', '/Applications/workbench/')
library(rgl) #to use ciftiTools graphics

gordon_parcel_matching <- read.csv("~/Box/tools/parcellations/Gordon_fs_LR/Gordon333Atlas.Parcels.LabelKey.csv")

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

# Load RData files with cleaned data --------------------------------------
load("~/Box/projects/in_progress/struct_funct_neonates/data/clean_data_july_2023/sc_matrices/birth/probabilistic/n267_SC_array_normed_to_ROI_size_edges_norm_CV75.RData")
load("~/Box/projects/in_progress/struct_funct_neonates/data/clean_data_july_2023/fc_matrices/n261_pconns_matrix_array.RData")

# Regional SC-FC controlling for ROI size ---------------------------------

subj_list_fc_sc <- intersect(fc_data$modid,sc_prob_subj_list_nonas_filt);length(subj_list_fc_sc)

#Controlling for ROI size, normalized to edge weight, CV-based thresholding at 75%
scfc.norm.roi.thresh75 <- sc_fc_coupl_func(big.FC.mat,big.SC.prob.symm.norm.roi.thresh, subj_list_fc_sc, fc_data$modid, sc_prob_subj_list_nonas)
data <- data.frame(scfc.norm.roi.thresh75$mean.regional.sc.fc); colnames(data) <- "mean.regional.sc.fc"

#plot it
data$label <- gordon_parcel_matching$label
ggplot(data=data,mapping = aes(fill=mean.regional.sc.fc)) +
  geom_brain(atlas = gordon, position = position_brain(hemi ~ side)) +
  theme_void() + 
  theme(legend.position = "bottom") + labs(fill="Spearman's rho")+
  scale_fill_gradient2(low= "lightblue", mid = "white", high = "orange",midpoint = mean(data$mean.regional.sc.fc), na.value = "white",guide = "colourbar", aesthetics = "fill")

# Take out small parcels where > 5 people had no parcel in the cortical ribbon --------------------
# As per Jeanette
# take out NA parcels from SC and FC matrices and redo
columns_with_na <- colSums(is.na(scfc.norm.roi.thresh75$` sc.fc.coupling`))
indices_to_remove <- which(columns_with_na>5)
indices_remaining <- setdiff(c(1:333), indices_to_remove)
big.SC.prob.symm.norm.roi.thresh.324 <- big.SC.prob.symm.norm.roi.thresh[-indices_to_remove,-indices_to_remove,]
gordon_parcel_matching_324 <- gordon_parcel_matching[-indices_to_remove,]
big.FC.mat.324 <- big.FC.mat[-indices_to_remove,-indices_to_remove,]
gordon_networks_324 <- gordon_networks[-indices_to_remove]
dim(big.FC.mat.324)
save(indices_to_remove, indices_remaining, gordon_networks_324, big.SC.prob.symm.norm.roi.thresh.324,gordon_parcel_matching_324, big.FC.mat.324,subj_list_fc_sc,sc_prob_subj_list_nonas,sc_prob_subj_list_nonas_filt, fc_data, scfc.norm.roi.thresh75.324, file = "~/Box/projects/in_progress/struct_funct_neonates/data/clean_data_july_2023/sc_fc_gordon_324/n239_data_for_324_parcel_scfc.RData")
load("~/Box/projects/in_progress/struct_funct_neonates/data/clean_data_july_2023/sc_fc_gordon_324/n239_data_for_324_parcel_scfc.RData")

#structure-function coupling without that have NAs in < 5 subjects
scfc.norm.roi.thresh75.324 <- sc_fc_coupl_func(big.FC.mat.324,big.SC.prob.symm.norm.roi.thresh.324, subj_list_fc_sc, fc_data$modid, sc_prob_subj_list_nonas, 324)
data <- data.frame(scfc.norm.roi.thresh75.324$mean.regional.sc.fc); colnames(data) <- "mean.regional.sc.fc"
#plot it
data$label <- gordon_parcel_matching_324$label
ggplot(data=data,mapping = aes(fill=mean.regional.sc.fc)) +
  geom_brain(atlas = gordon, position = position_brain(hemi ~ side)) +
  theme_void() + 
  theme(legend.position = "bottom") + labs(fill="Spearman's rho")+
  scale_fill_gradient2(low= "lightblue", mid = "white", high = "orange",midpoint = mean(data$mean.regional.sc.fc), na.value = "white",guide = "colourbar", aesthetics = "fill")

save(scfc.norm.roi.thresh75.324, scfc.norm.roi.thresh75, file = "~/Box/projects/in_progress/struct_funct_neonates/data/clean_data_july_2023/sc_fc_gordon_324/n239_scfc_objects_normalized_to_ROI_size.RData")
load("~/Box/projects/in_progress/struct_funct_neonates/data/clean_data_july_2023/sc_fc_gordon_324/n239_scfc_objects_normalized_to_ROI_size.RData")
load("~/Box/projects/in_progress/struct_funct_neonates/data/clean_data_july_2023/sc_fc_gordon_324/n239_data_for_324_parcel_scfc.RData")

# Does SC-FC variation correspond to the S-A axis? ----------------------------------
#Spin tests from https://github.com/frantisekvasa/rotate_parcellation
source("~/Box/tools/rotate_parcellation/R/perm.sphere.p.R") 
load("~/Box/projects/in_progress/struct_funct_neonates/data/rotations_gordon324_volume_MNI_10000x.Rdata")
sa_axis_gordon <- read.csv("~/Box/tools/parcellations/Gordon_fs_LR/SensorimotorAssociation_Axis_Gordon333.csv")
cbind(gordon_parcel_matching$label, sa_axis_gordon$label) #is it the same? Yay it is!
sa_axis_gordon_324 <- sa_axis_gordon[-indices_to_remove,]

#correlation with age effects
cor.test(data$mean.regional.sc.fc,sa_axis_gordon_324$SA.Axis.ranks, method = "spearman")
cor.test(data$mean.regional.sc.fc,sa_axis_gordon_324$SA.Axis.zscores,method = "spearman")

print(perm.sphere.p(data$mean.regional.sc.fc,sa_axis_gordon_324$SA.Axis.zscores, rotations.324, "spearman"))
print(perm.sphere.p(data$mean.regional.sc.fc,sa_axis_gordon_324$SA.Axis.ranks, rotations.324, "spearman"))

# Write out scfc coupling to workbench --------
values <- ifelse(gordon_networks=="None",0,scfc.norm.roi.thresh75$mean.regional.sc.fc) #don't write out to the None network
#dim is 347, because of subcortical, so pad it with 0s
new.347.values <- rep(0, 347)
#make 324 vector into 333 to write out
new.347.values[indices_remaining] <- values
output_file=paste0("~/Box/projects/in_progress/struct_funct_neonates/data/clean_data_july_2023/for_viz/scfc_for_workbench/scfc.norm.roi.thresh75_mean_regional_value.txt")
write.table(new.347.values, output_file, col.names = F, row.names = F)
output_pconn=paste0("~/Box/projects/in_progress/struct_funct_neonates/data/clean_data_july_2023/for_viz/scfc_for_workbench/scfc.norm.roi.thresh75_mean_regional_value.ptseries.nii")
command = sprintf("-cifti-convert -from-text % s  ~/Box/projects/in_progress/struct_funct_neonates/data/for_viz/MOD1022_V1_a_ALFF_std_parcellated.ptseries.nii % s", output_file, output_pconn) 
ciftiTools::run_wb_cmd(command, intern = FALSE, ignore.stdout = NULL, ignore.stderr = NULL)

# Age effects on structure-function ---------------------------------------
load("~/Box/projects/in_progress/struct_funct_neonates/data/clean_data_july_2023/sc_fc_gordon_324/n239_scfc_objects_normalized_to_ROI_size.RData")
load("~/Box/projects/in_progress/struct_funct_neonates/data/clean_data_july_2023/fc_matrices/n261_pconns_matrix_array.RData")

sc_fc_data <- filter(fc_data, modid %in% subj_list_fc_sc)
sc_fc_data <- data.frame(cbind(sc_fc_data,scfc.norm.roi.thresh75.324$mean.subject.sc.fc))
FD <- read.csv("~/Box/Tooley 01_18_2023 eLABE Requested Data/N319_avg_FD_retained_02_28_2023.csv") %>% rename(modid=MODID);
FD$modid <- str_remove(FD$modid,"_V1_a")
sc_fc_data <- left_join(sc_fc_data, FD, by="modid")
icv <- read.csv("~/Box/Tooley 01_18_2023 eLABE Requested Data/tooley_icv_tbv_elabe_neonate.csv") %>% rename(modid=MODID)
sc_fc_data <- left_join(sc_fc_data, icv, by="modid")
#read in diffusion motion data
eddy_motion <- read.csv("~/Box/projects/in_progress/struct_funct_neonates/data/elabe_FSL_604_eddy_quad_results_2023_01_10.csv");
sc_fc_data <- left_join(sc_fc_data, eddy_motion, by="modid")

#average SC-FC across cortex, controlling for covariates
age_scfc_model <- lm(scfc.norm.roi.thresh75.324.mean.subject.sc.fc~PMA_scan+child_sex+avg_FD_of_retained_frames+retained_frames+Relative_Motion, data=sc_fc_data)
summary(age_scfc_model);lm.beta(age_scfc_model)
visreg(age_scfc_model, "PMA_scan", xlab="Age (months)", ylab= "Structure-function coupling")
visreg((lm(scfc.norm.roi.thresh75.324.mean.subject.sc.fc~PMA_scan+child_sex+avg_FD_of_retained_frames+retained_frames+Relative_Motion, data=sc_fc_data)))

##plot age effect at the parcel level
subject_parcel_scfc_values <- scfc.norm.roi.thresh75.324$` sc.fc.coupling`; subject_parcel_scfc_values <- cbind(sc_fc_data, subject_parcel_scfc_values)
parcel_age_pvals<- lapply(colnames(select(subject_parcel_scfc_values,matches("^[[:digit:]]"))), function(x) { summary(lm(as.formula(paste0("`",x,"` ~ PMA_scan+child_sex+avg_FD_of_retained_frames+retained_frames+Relative_Motion")), data=subject_parcel_scfc_values))$coef[2,4]})
parcel_age_pvals <- unlist(parcel_age_pvals)
parcel_age_pvals_fdr <- cbind(parcel_age_pvals,p.adjust(parcel_age_pvals,method = "fdr"))
#get age betas
parcel_age_betas<- lapply(colnames(select(subject_parcel_scfc_values,matches("^[[:digit:]]"))), function(x) { lm.beta(lm(as.formula(paste0("`",x,"` ~ PMA_scan+child_sex+avg_FD_of_retained_frames+retained_frames+Relative_Motion")), data=subject_parcel_scfc_values))$standardized.coefficients[[2]]})
parcel_age_betas <- unlist(parcel_age_betas)
data_scfc <- data.frame(cbind(parcel_age_pvals_fdr,parcel_age_betas))
colnames(data_scfc) <- c("parcel_age_pvals","parcel_age_pvals_fdr","parcel_age_betas")

data_scfc$values <- ifelse(data_scfc$parcel_age_pvals_fdr < 0.05, data_scfc$parcel_age_betas, 0)
#data_scfc$values <- data_scfc$parcel_age_betas
data_scfc$label <- gordon_parcel_matching_324$label
ggplot(data=data_scfc,mapping = aes(fill=values)) +
  geom_brain(atlas = gordon, position = position_brain(hemi ~ side)) +
  theme_void() + 
  theme(legend.position = "bottom") + labs(fill=paste("Significant \u03B2s"))+
  scale_fill_gradient2(low= "blue", mid = "white", high = "red", guide = "colourbar", aesthetics = "fill")

# Write out age effects to workbench --------------------------------------
#remove 'None' network values
values <- ifelse(gordon_networks_324=="None",0,data_scfc$parcel_age_pvals_fdr)
#dim is 347, because of subcortical, so pad it with 0s
new.347.values <- rep(0, 347)
#make 324 vector into 333
new.347.values[indices_remaining] <- values
output_file=paste0("~/Box/projects/in_progress/struct_funct_neonates/data/clean_data_july_2023/for_viz/scfc_for_workbench/n239_scfc_age_effect_lin_pvalues_fdr.txt")
write.table(new.347.values, output_file, col.names = F, row.names = F)
output_pconn=paste0("~/Box/projects/in_progress/struct_funct_neonates/data/clean_data_july_2023/for_viz/scfc_for_workbench/n239_scfc_age_effect_lin_pvalues_fdr.ptseries.nii")
command = sprintf("-cifti-convert -from-text % s  ~/Box/projects/in_progress/struct_funct_neonates/data/for_viz/MOD1022_V1_a_ALFF_std_parcellated.ptseries.nii % s", output_file, output_pconn) 
ciftiTools::run_wb_cmd(command, intern = FALSE, ignore.stdout = NULL, ignore.stderr = NULL)

# Age effects in raincloud ------------------------------------------------
data_scfc$network <- gordon_networks_324
mean_effects <- data_scfc  %>% filter(., network !="None") %>% group_by(network) %>% 
  summarise(mean_sig_only=mean(values), mean_all_values=mean(parcel_age_betas), med_sig_only=median(values),
            med_all_values=median(parcel_age_betas)) %>% arrange(med_all_values)
mean_effects

#plot boxplot
g <- ggplot(data = filter(data_scfc, network !="None"), aes(y = reorder(network, -parcel_age_betas, median), x = parcel_age_betas, fill=network)) +
  geom_boxplot(width = .5, alpha = 0.8) +
  geom_point(aes(y = reorder(network, -parcel_age_betas, median), x = parcel_age_betas, color = network), alpha=0.5) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_fill_gordon() +
  scale_color_gordon() +
  theme_bw() +
  xlab("Age effects on structure-function coupling")+
  ylab("")+
  raincloud_theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) #fix the offset
g

#permutation test following Vàzquez-Rodriguez et al. (2019)
ci <- gordon_networks_324
# number of communities / classes
nci <- max(as.numeric(ci))
dumdum <- model.matrix(~factor(ci)-1)
# get mean R-square for each community
ci_mu <- t(data_scfc$parcel_age_betas) %*% dumdum / colSums(dumdum)

# label-permuting null model, set number of permutations
nperm <- 10000
# initialize community-wise mean R-square for each permutation
ci_perm <- matrix(0, nrow = nci, ncol = nperm)

for (ii in 1:nperm) {
  # permute community assignments (without replacement)
  p <- sample(length(ci))
  
  # dummy-code permuted community assignments
  dumdum_perm <- model.matrix(~factor(ci[p])-1)
  
  # get mean R-square for permuted community assignments
  ci_perm[, ii] <- t(data_scfc$parcel_age_betas) %*% dumdum_perm / colSums(dumdum_perm)
  
  cat("permutation", ii, "out of", nperm, "done\n")
}

ci_z <- numeric(nci)
for (ii in 1:nci) {
  # mean R-square for real community ii
  x <- ci_mu[ii]
  
  # mean R-square for permuted community ii
  mu <- mean(ci_perm[ii, ])
  
  # standard deviation of R-square for permuted community ii
  sigma <- sd(ci_perm[ii, ])
  
  # z-score
  ci_z[ii] <- (x - mu) / sigma
}

#ci_z <- ci_z[ci_z>-3]
i <- order(ci_z)
gordon_networks_324_no_none <- gordon_networks_324[gordon_networks_324!="None"];gordon_networks_324_no_none <- droplevels(gordon_networks_324_no_none)

df <- data.frame(Network = levels(gordon_networks_324)[i],
                 Zscore = ci_z[i])
df$p_values <- 2 * (1 - pnorm(abs(df$Zscore)))

# create the plot using ggplot
ggplot(df, aes(x = reorder(Network, Zscore, mean), y = Zscore, fill=ifelse(p_values<0.05,1,0))) +
  geom_bar(stat = "identity") +
  xlab("Network") +
  ylab("z-score") +
  scale_x_discrete(labels = levels(gordon_networks_324)[i]) +
  theme_classic()+
  theme(legend.position = "none")

# Do SC-FC age effects correspond to the S-A axis? ----------------------------------
source("~/Box/tools/rotate_parcellation/R/perm.sphere.p.R")
load("~/Box/projects/in_progress/struct_funct_neonates/data/rotations_gordon324_volume_MNI_10000x.Rdata")
sa_axis_gordon <- read.csv("~/Box/tools/parcellations/Gordon_fs_LR/SensorimotorAssociation_Axis_Gordon333.csv")
cbind(gordon_parcel_matching$label, sa_axis_gordon$label) #is it the same? Yay it is!
sa_axis_gordon_324 <- sa_axis_gordon[-indices_to_remove,]

#correlation with age effects
cor.test(data_scfc$parcel_age_betas,sa_axis_gordon_324$SA.Axis.zscores,method = "spearman")
cor.test(data_scfc$parcel_age_betas,sa_axis_gordon_324$SA.Axis.ranks,method = "spearman")

print(perm.sphere.p(data_scfc$parcel_age_betas,sa_axis_gordon_324$SA.Axis.zscores, rotations.324, "spearman"))
print(perm.sphere.p(data_scfc$parcel_age_betas,sa_axis_gordon_324$SA.Axis.ranks, rotations.324, "spearman"))
data_scfc$SA.Axis.ranks <- sa_axis_gordon_324$SA.Axis.ranks
data_scfc$SA.Axis.zscores <- sa_axis_gordon_324$SA.Axis.zscores

ggplot(data_scfc, aes(x = SA.Axis.ranks, y = parcel_age_betas)) + 
  geom_point(aes(color = parcel_age_betas)) +
  scale_color_gradientn(colours = c("#64CDF6","#597FC0","#4863AE","#ED1C24","#F47920"),values = scales::rescale(c(-0.37, -0.2,-0.1, 0.1, 0.2), to =c(0,1)), ) +
  geom_smooth(method = "lm")+
  theme_classic(base_size = 16)+
  theme(legend.position = "none")

# Load network metrics on SC matrices --------------------------------
#just use those calculated in MATLAB
sc_network_metrics_matlab <- read.csv("~/Box/projects/in_progress/struct_funct_neonates/data/clean_data_july_2023/sc_matrices/birth/probabilistic/SymmNormRoiNormEdgesCV75/n255_nreg324_birth_sc_probabilistic_norm_roi_norm_edge_CV_thresh75_within_between_gordon_withmodulpartcoef.csv") %>% rename(modid=subjs)

#load demographic data
exclusions <- read.csv("~/Box/Tooley 01_18_2023 eLABE Requested Data/subjects_all_inclusion_exclusion_01_25_2022.csv")
demographics <- read.csv("~/Box/Tooley 01_18_2023 eLABE Requested Data/Tooley_RC_20230118_updated08JUN23.csv")

demo_data_all <- left_join(demographics, exclusions, by="modid");dim(demo_data_all)

#exclude pre-term babies, MRI_injury exclusion, NICU exclusions, and birthweight exclusion and IRB exclusion
demo_data <- demo_data_all %>% filter(.,combined_exclusion_initial_analyses ==0 & irb_exclusion == 0);dim(demo_data)
#remame PMA
demo_data$PMA_scan <- demo_data$mri_test_pma_scan_dob
#make sex a factor
demo_data$child_sex <- factor(demo_data$child_sex, labels=c("Male","Female"))

# Examine age and SC network metrics ---------------------------
sc_network_metrics_matlab <- left_join(sc_network_metrics_matlab, demo_data, by="modid")
icv <- read.csv("~/Box/Tooley 01_18_2023 eLABE Requested Data/tooley_icv_tbv_elabe_neonate.csv") %>% rename(modid=MODID)
sc_network_metrics_matlab<- left_join(sc_network_metrics_matlab,icv,by="modid")
eddy_motion <- read.csv("~/Box/projects/in_progress/struct_funct_neonates/data/elabe_FSL_604_eddy_quad_results_2023_01_10.csv");
sc_network_metrics_matlab<- left_join(sc_network_metrics_matlab,eddy_motion,by="modid")
#make it only subjects who are in the SC-FC sample
sc_network_metrics_matlab <- filter(sc_network_metrics_matlab, modid %in% subj_list_fc_sc)

#average network metrics, controlling for covariates, total_edges changes it considerably
#participation coefficient
sc_part_coef_model <- lm(scale(part_coef_avg)~scale(PMA_scan)+child_sex+scale(Relative_Motion)+scale(intracranial_volume), data=sc_network_metrics_matlab) #total_edges, intracranial volume
sc_part_coef_model <- lm(part_coef_avg~PMA_scan+child_sex+Relative_Motion+intracranial_volume, data=sc_network_metrics_matlab) #total_edges, intracranial volume

summary(sc_part_coef_model);lm.beta(sc_part_coef_model);visreg(sc_part_coef_model, "PMA_scan", xlab="Age (months)", ylab= "Participation coefficient")
performance::check_model(sc_part_coef_model);visreg(sc_part_coef_model)

#modularity
sc_modul_model <- lm(modul~PMA_scan+child_sex+Relative_Motion+child_sex+intracranial_volume, data=sc_network_metrics_matlab) #total_edges, intracranial volume
summary(sc_modul_model);lm.beta(sc_modul_model);visreg(sc_modul_model, "PMA_scan", xlab="Age (months)", ylab= "Modularity")
performance::check_model(sc_modul_model)

# Load network metrics on FC matrices --------------------------------
#just use those calculated in matlab for the age x SES paper! Participation coefficient in R has issues, is reversed from in MATLAB for unsigned networks.
y0_network_metrics_matlab <- read.csv("~/Box/projects/in_progress/within_between_network_longitudinal/data/birth/fc/n319_nreg324_birth_within_between_gordon_withmodulpartcoef.csv") %>% rename(modid=Var1)
y0_network_metrics_matlab$modid <- str_remove(y0_network_metrics_matlab$modid,"_V1_a");y0_network_metrics_matlab$modid <- str_remove(y0_network_metrics_matlab$modid,"_V1_b")
FD <- read.csv("~/Box/Tooley 01_18_2023 eLABE Requested Data/N319_avg_FD_retained_02_28_2023.csv") %>% rename(modid=MODID);
FD$modid <- str_remove(FD$modid,"_V1_a");FD$modid <- str_remove(FD$modid,"_V1_b")

y0_network_metrics_matlab <- left_join(y0_network_metrics_matlab, FD, by="modid")
y0_network_metrics_matlab <- left_join(y0_network_metrics_matlab, demo_data, by="modid")

# Age and FC metrics------------------------------------------
y0_network_metrics_matlab$avg_part_coef <- (y0_network_metrics_matlab$part_coef_neg+y0_network_metrics_matlab$part_coef_pos)/2
#only the participants in the SC-FC sample
y0_network_metrics_matlab <- filter(y0_network_metrics_matlab, modid %in% subj_list_fc_sc)

#participation coefficient
fc_part_coef_model <- lm(avg_part_coef~PMA_scan+child_sex+avg_FD_of_retained_frames+retained_frames, data=y0_network_metrics_matlab) #total_edges, intracranial volume
summary(fc_part_coef_model);lm.beta(fc_part_coef_model);visreg(fc_part_coef_model, "PMA_scan", xlab="Age (months)", ylab= "Participation coefficient")
performance::check_model(fc_part_coef_model)
#gams don't change anything, not needed

#with avg weight
fc_part_coef_model <- lm(avg_part_coef~PMA_scan+child_sex+avg_FD_of_retained_frames+retained_frames+avgweight, data=y0_network_metrics_matlab) #total_edges, intracranial volume
summary(fc_part_coef_model);lm.beta(fc_part_coef_model);visreg(fc_part_coef_model, "PMA_scan", xlab="Age (months)", ylab= "Participation coefficient")
performance::check_model(fc_part_coef_model)

#modularity
fc_modul_model <- lm(modul~PMA_scan+child_sex+avg_FD_of_retained_frames+retained_frames, data=y0_network_metrics_matlab) #total_edges, intracranial volume
summary(fc_modul_model);lm.beta(fc_modul_model);visreg(fc_modul_model, "PMA_scan", xlab="Age (months)", ylab= "Modularity")
performance::check_model(fc_modul_model)

#with avg weight
fc_modul_model <- lm(modul~PMA_scan+child_sex+avg_FD_of_retained_frames+retained_frames+avgweight, data=y0_network_metrics_matlab) #total_edges, intracranial volume
summary(fc_modul_model);lm.beta(fc_modul_model);visreg(fc_modul_model, "PMA_scan", xlab="Age (months)", ylab= "Modularity")
performance::check_model(fc_modul_model)

# Average structure-function coupling within network per participant --------
#Load SC-FC data again
load("~/Box/projects/in_progress/struct_funct_neonates/data/clean_data_july_2023/sc_fc_gordon_324/n239_data_for_324_parcel_scfc.RData")
sc_fc_data <- filter(fc_data, modid %in% subj_list_fc_sc)
sc_fc_data <- data.frame(cbind(sc_fc_data,scfc.norm.roi.thresh75.324$mean.subject.sc.fc))
FD <- read.csv("~/Box/Tooley 01_18_2023 eLABE Requested Data/N319_avg_FD_retained_02_28_2023.csv") %>% rename(modid=MODID);
FD$modid <- str_remove(FD$modid,"_V1_a");FD$modid <- str_remove(FD$modid,"_V1_b")
sc_fc_data <- left_join(sc_fc_data, FD, by="modid")
icv <- read.csv("~/Box/Tooley 01_18_2023 eLABE Requested Data/tooley_icv_tbv_elabe_neonate.csv") %>% rename(modid=MODID)
sc_fc_data <- left_join(sc_fc_data, icv, by="modid")
#read in diffusion motion data
eddy_motion <- read.csv("~/Box/projects/in_progress/struct_funct_neonates/data/elabe_FSL_604_eddy_quad_results_2023_01_10.csv");
sc_fc_data <- left_join(sc_fc_data, eddy_motion, by="modid")

subject_parcel_scfc_values <- scfc.norm.roi.thresh75.324$` sc.fc.coupling`; subject_parcel_scfc_values <- cbind(sc_fc_data, subject_parcel_scfc_values)

subject_parcel_scfc_values_row <- subject_parcel_scfc_values %>% rowwise(modid)
oldnames <- colnames(subject_parcel_scfc_values_row[,108:431])
newnames <- paste0(as.character(gordon_networks_324),"_V",1:324)
#rename rows to be the name of the network
for(i in 1:324) names(subject_parcel_scfc_values_row)[names(subject_parcel_scfc_values_row) == oldnames[i]] = newnames[i]

subject_parcel_scfc_values_row <- subject_parcel_scfc_values_row %>% mutate(avg_default_scfc=mean(c_across(contains("Default")), na.rm = T), avg_fpn_scfc= mean(c_across(contains("FrontoParietal")),na.rm = T),
                                                                            avg_vis_scfc=mean(c_across(contains("Visual")),na.rm = T), avg_aud_scfc=mean(c_across(contains("Auditory")),na.rm = T),
                                                                            avg_medparietal_scfc=mean(c_across(contains("MedialParietal")),na.rm = T), avg_co_scfc=mean(c_across(contains("CinguloOperc")),na.rm = T),
                                                                            avg_van_scfc=mean(c_across(contains("VentralAttn")),na.rm = T),avg_dan_scfc=mean(c_across(contains("DorsalAttn")),na.rm = T),
                                                                            avg_smh_scfc=mean(c_across(contains("SMhand")),na.rm = T),avg_smm_scfc=mean(c_across(contains("SMmouth")),na.rm = T),
                                                                            avg_sal_scfc=mean(c_across(contains("Salience")),na.rm = T),avg_paroccip_scfc=mean(c_across(contains("ParietoOccip")),na.rm = T)) %>% 
  ungroup()

length(subject_parcel_scfc_values_row$avg_smm_scfc)
save(subject_parcel_scfc_values_row, file= "~/Box/projects/in_progress/struct_funct_neonates/data/clean_data_july_2023/sc_fc_gordon_324/n239_rowwise_scfc_average.RData")
load("~/Box/projects/in_progress/struct_funct_neonates/data/clean_data_july_2023/sc_fc_gordon_324/n239_rowwise_scfc_average.RData")

# Is within-network SC-FC related to eyetracking variables? -------------------------
#load instead the data from Aidan, which has a few more subjects
eyetracking_y1 <- read.csv("~/Box/Tooley 01_18_2023 eLABE Requested Data/tooley_y1_eyetracking_cb1_data.csv")
eyetracking_y1_2nd <- read.csv("~/Box/Tooley 01_18_2023 eLABE Requested Data/tooley_y1_eyetracking_cb2_data.csv")

#merge in calibration accuracy data
library(readxl)
calibration_accuracy <- read_excel("~/Box/ELABE_eyetrack_trajectories/DATA/Y1_Eyetracking_percent_valid_gaze_samples.xlsx")
colnames(calibration_accuracy);calibration_accuracy$modid <- calibration_accuracy$MODID

eyetracking_y1_data <- rbind(eyetracking_y1_2nd,eyetracking_y1)
eyetracking_y1_data$modid <- eyetracking_y1_data$MODID
#only have Avg.TTFF, average across center and peripheral
eyetracking_y1_data$Avg.Center.Peripheral.TTFF <- rowMeans(dplyr::select(eyetracking_y1_data, c("Center.Avg.TTFF.Face", "Peripheral.Avg.TTFF.Face")), na.rm = T)

sc_fc_eyetrack_y1 <- left_join(subject_parcel_scfc_values_row, eyetracking_y1_data, by= "modid")
sc_fc_eyetrack_y1 <- left_join(sc_fc_eyetrack_y1, calibration_accuracy, by= "modid")
describe(sc_fc_eyetrack_y1$Avg.Center.Peripheral.TTFF)#102 kids with eye-tracking at Y1 who are healthy full-term and have birth fMRI data
sc_fc_eyetrack_y1 %>% filter(!is.na(Avg.TTFF) & !is.na(avg_fpn_scfc)) %>% dplyr::select(modid) %>% dim()
#102 with FC and dMRI data at birth

#look at mean TTTF associations with calibration accuracy, precision, and % valid trials
cor.test(sc_fc_eyetrack_y1$Avg.Center.Peripheral.TTFF, sc_fc_eyetrack_y1$calibration_precision_deg)
cor.test(sc_fc_eyetrack_y1$Avg.Center.Peripheral.TTFF, sc_fc_eyetrack_y1$percent_valid)
cor.test(sc_fc_eyetrack_y1$Avg.Center.Peripheral.TTFF, sc_fc_eyetrack_y1$calibration_accuracy_deg)

#just look at average sc-fc overall predicing eye-tracking
TTFF.avg.model <- lm(scale(Avg.Center.Peripheral.TTFF)~ scale(child_age_y1_assessment)+scale(percent_valid)+scale(scfc.norm.roi.thresh75.324.mean.subject.sc.fc), data=sc_fc_eyetrack_y1)
TTFF.avg.model <- lm(Avg.Center.Peripheral.TTFF~ child_age_y1_assessment+percent_valid+scfc.norm.roi.thresh75.324.mean.subject.sc.fc, data=sc_fc_eyetrack_y1)
summary(TTFF.avg.model);lm.beta(TTFF.avg.model)
visreg(TTFF.avg.model, "scfc.norm.roi.thresh75.324.mean.subject.sc.fc", xlab="Average structure-function coupling", ylab="Reaction time (ms)")

visreg(lm(Avg.Center.Peripheral.TTFF~ child_age_y1_assessment+child_sex+percent_valid+scfc.norm.roi.thresh75.324.mean.subject.sc.fc, data=sc_fc_eyetrack_y1),
       "scfc.norm.roi.thresh75.324.mean.subject.sc.fc", xlab="Average SC-FC", main="Average Time to First Fixation (TTFF)")

#Motor and visual and attentional networks, a priori
par(mfrow=c(1,3))
networks <- c("avg_smh_scfc","avg_smm_scfc","avg_vis_scfc", "avg_dan_scfc","avg_co_scfc","avg_van_scfc" )
pvals <- list()
par(mfrow=c(2,3))
for (network in networks){
  l <- lm(as.formula(paste0("Avg.Center.Peripheral.TTFF~ child_age_y1_assessment+percent_valid+", network)), data=sc_fc_eyetrack_y1)
  print(summary(l));pvals[network] <- summary(l)$coefficients[4,4]
  visreg(l, network ,ylab= "Average TTFF", xlab = network, main = paste(network, "p =", round(summary(l)$coefficients[4,4],2)))
}
unlist(pvals)
pvals <- data.frame(cbind(networks,unlist(pvals),p.adjust(unlist(pvals),method = "fdr")))
pvals

#look at all networks
networks=sc_fc_eyetrack_y1 %>% select(ends_with("_scfc")) %>% colnames()
networks
pvals=list()
par(mfrow=c(3,4))
for (network in networks){
  l <- lm(as.formula(paste0("Avg.Center.Peripheral.TTFF ~ child_age_y1_assessment+percent_valid+", network)), data=sc_fc_eyetrack_y1)
  print(summary(l));pvals[network] <- summary(l)$coefficients[4,4]
  visreg(l, network ,ylab= "Average TTFF", xlab = network, main = paste(network, "p =", round(summary(l)$coefficients[4,4],2)))
}
unlist(pvals)
pvals <- data.frame(cbind(networks,unlist(pvals),p.adjust(unlist(pvals),method = "fdr")))
pvals

#Plot only networks that are sig after FDR
#DAN
summary(lm(Avg.Center.Peripheral.TTFF~ child_age_y1_assessment+percent_valid+avg_dan_scfc, data=sc_fc_eyetrack_y1))
lm.beta(lm(Avg.Center.Peripheral.TTFF~ child_age_y1_assessment+percent_valid+avg_dan_scfc, data=sc_fc_eyetrack_y1))
visreg(lm(Avg.Center.Peripheral.TTFF~ child_age_y1_assessment+percent_valid+avg_dan_scfc, data=sc_fc_eyetrack_y1),
       "avg_dan_scfc",  xlab="Dorsal attention SC-FC", main="", ylim=c(0,1000))

summary(lm(Center.Avg.TTFF.Face~ child_age_y1_assessment+percent_valid+avg_dan_scfc, data=sc_fc_eyetrack_y1))
summary(lm(Peripheral.Avg.TTFF.Face~ child_age_y1_assessment+percent_valid+avg_dan_scfc, data=sc_fc_eyetrack_y1))

#CO 
summary(lm(Avg.Center.Peripheral.TTFF~ child_age_y1_assessment+percent_valid+avg_co_scfc, data=sc_fc_eyetrack_y1))
lm.beta(lm(Avg.Center.Peripheral.TTFF~ child_age_y1_assessment+percent_valid+avg_co_scfc, data=sc_fc_eyetrack_y1))
visreg(lm(Avg.Center.Peripheral.TTFF~ child_age_y1_assessment+percent_valid+avg_co_scfc, data=sc_fc_eyetrack_y1),
       "avg_co_scfc",  xlab="Cingulo-opercular SC-FC", main="", ylim=c(0,1000))
summary(lm(Center.Avg.TTFF.Face~ child_age_y1_assessment+percent_valid+avg_co_scfc, data=sc_fc_eyetrack_y1))
summary(lm(Peripheral.Avg.TTFF.Face~ child_age_y1_assessment+percent_valid+avg_co_scfc, data=sc_fc_eyetrack_y1))

#VIS
summary(lm(Avg.Center.Peripheral.TTFF~ child_age_y1_assessment+percent_valid+avg_vis_scfc, data=sc_fc_eyetrack_y1))
lm.beta(lm(Avg.Center.Peripheral.TTFF~ child_age_y1_assessment+percent_valid+avg_vis_scfc, data=sc_fc_eyetrack_y1))
visreg(lm(Avg.Center.Peripheral.TTFF~ child_age_y1_assessment+percent_valid+avg_vis_scfc, data=sc_fc_eyetrack_y1),
       "avg_vis_scfc",  xlab="Cingulo-opercular SC-FC", main="", ylim=c(0,1000))
summary(lm(Center.Avg.TTFF.Face~ child_age_y1_assessment+percent_valid+avg_vis_scfc, data=sc_fc_eyetrack_y1))
summary(lm(Peripheral.Avg.TTFF.Face~ child_age_y1_assessment+percent_valid+avg_vis_scfc, data=sc_fc_eyetrack_y1))
