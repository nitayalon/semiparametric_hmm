### Source file for HiddenMArkovModel class
library(parallel)
library(magrittr)
library(Jmisc)
hmm_class_dir <- "~/MastersThesis/Main/Code/Mock_data_classes/hidden_markov_model_class"

fileName <- function(hmm_class_dir, file_name)
{
  paste0(hmm_class_dir,"/",file_name)
}

source(fileName(hmm_class_dir,"valid_hmm_class.R"))
source(fileName(hmm_class_dir,"hmm_model_class.R"))
source(fileName(hmm_class_dir,"generate_observation_from_hmm.R"))
source(fileName(hmm_class_dir,"create_transition_matrix.R"))
source(fileName(hmm_class_dir,"create_grid_parameters_for_hmm.R"))
source(fileName(hmm_class_dir,"compute_theta_method.R"))
source(fileName(hmm_class_dir,"compute_b.R"))
source(fileName(hmm_class_dir,"compute_starting_vector.R"))
source(fileName(hmm_class_dir,"compute_forward_vector.R"))
source(fileName(hmm_class_dir,"compute_backward_vector.R"))
source(fileName(hmm_class_dir,"compute_density_matrix.R"))
source(fileName(hmm_class_dir,"plot_hmm_observations.R"))
source(fileName(hmm_class_dir,"estimate_transition_matrix.R"))
source(fileName(hmm_class_dir,"create_grid_parameters_for_hmm.R"))
source(fileName(hmm_class_dir,"estimate_transition_matrix.R"))
sourceAll("~/MastersThesis/Main/Code/Mock_data_classes/Method_helper_for_HMM_class/")
sourceAll("~/MastersThesis/Main/Code/geneate_data/Double_exponential/")
sourceAll("~/MastersThesis/Main/Code/geneate_data/Normal/")
