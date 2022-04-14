# This file takes a matrix of data and generates 
# the design matrices necessary to fit the stan model.  
# 

# The data needs to look like this: 
# data (counts, rows = observations, columns = species)
# a column for "site", (rows)
# a column for "tech_rep", (technical replicates)
# a set of columns for 

# If you have mock community data it needs to look like this: 
# data (counts, rows = observations, columns = species).  Same column order as the sample data.
# a column for "site", (rows)
# a column for "tech_rep", (technical replicates)
# a set of columns for the true specified community for each species:
  # same species in the same order as the the sampling data.
  # The proportions should be transformed into "alr" space and 
  # include the letters "alr" in the title of each column. Note that all values need to be non-infitite.
  # 

library(tidyverse)
library(MCMCpack)
library(compositions)
library(rstan)
library(dplyr)

########################################################################
#### Create data frames that can be read into Stan model
########################################################################

NOM <- as.name(colnames(p_mock_all)[1])
   formula_a <- eval(NOM) ~ N_pcr_mock -1
   model_frame <- model.frame(formula_a, p_mock_all)
   model_vector_a_mock <- model.matrix(formula_a, model_frame) %>% as.numeric()
   N_pcr_mock_small <- cbind(N_pcr_mock, p_mock_all) %>%  filter(tech_rep == 1) %>% pull(N_pcr_mock)
   formula_b <- eval(NOM) ~ N_pcr_mock_small -1
   model_frame <- model.frame(formula_b, p_mock_all%>% filter(tech_rep==1))
   model_vector_a_mock_small <- model.matrix(formula_b, model_frame) %>% as.numeric()
   
  N_obs_mock       <- nrow(p_mock_all)
  
# unknown communities second
  # species compositions (betas)

NOM <- as.name(colnames(p_samp_all)[1])    
  
  p_samp_all$site <- as.factor(p_samp_all$site)
  N_site = length(unique(p_samp_all$site))
  p_samp_all$tech_rep <- as.factor(p_samp_all$tech_rep)
  if(N_site == 1){
    formula_b <- eval(NOM) ~ 1  
  } else {
    formula_b <- eval(NOM) ~ site
  }
  
  model_frame <- model.frame(formula_b, p_samp_all)
  model_matrix_b_samp <- model.matrix(formula_b, model_frame)

  # choose a single representative for each site to make predictions to
  model_frame <- model.frame(formula_b, p_samp_all %>% filter(tech_rep==1))
  model_matrix_b_samp_small <- model.matrix(formula_b, model_frame)
    
  # efficiencies (alpha)
  formula_a <- eval(NOM) ~ N_pcr_samp -1
  model_frame <- model.frame(formula_a, p_samp_all)
  model_vector_a_samp <- model.matrix(formula_a, model_frame) %>% as.numeric()
  N_pcr_samp_small <- cbind(N_pcr_samp, p_samp_all) %>% filter(tech_rep == 1) %>% pull(N_pcr_samp)
  formula_b <- eval(NOM) ~ N_pcr_samp_small -1
  
  model_frame <- model.frame(formula_b, p_samp_all %>% filter(tech_rep==1))
  model_vector_a_samp_small <- model.matrix(formula_b, model_frame) %>% as.numeric()
  
  #counters 
  N_obs_samp_small <- nrow(model_matrix_b_samp_small)
  N_obs_samp <- nrow(p_samp_all)
  N_b_samp_col <- ncol(model_matrix_b_samp)
  
  
#### Make Stan objects
  
  stan_data <- list(
    N_species = ncol(p_samp_all)-2,   # Number of species in data
    N_obs_samp = nrow(p_samp_all), # Number of observed community samples and tech replicates ; this will be Ncreek * Nt * Nbiol * Ntech * 2 [for upstream/downstream observations]
    N_obs_mock = nrow(p_mock_all), # Number of observed mock samples, including tech replicates
    N_obs_samp_small = nrow(p_samp_all %>% filter(tech_rep == 1)), # Number of unique observed community samples ; this will be Ncreek * Nt * Nbiol * 2 [for upstream/downstream observations]

    # Observed data of community matrices
    sample_data = p_samp_all %>% dplyr::select(contains("sp")),
    mock_data   = mock %>% dplyr::select(contains("sp")),
    
    # True proportions for mock community
    #mock_true_prop = p_mock_all %>% dplyr::select(contains("sp")),
    alr_mock_true_prop = p_mock_all %>% dplyr::select(contains("alr")),

    # vectors of PCR numbers
    N_pcr_samp = N_pcr_samp,
    N_pcr_mock = N_pcr_mock,
    
    # Design matrices: field samples
    N_b_samp_col = N_b_samp_col,
    model_matrix_b_samp = model_matrix_b_samp,
    model_matrix_b_samp_small = as.array(model_matrix_b_samp_small),
    model_vector_a_samp = model_vector_a_samp,
    model_vector_a_samp_small = as.array(model_vector_a_samp_small),
      
    # Design matrices: mock community samples
    model_vector_a_mock = as.array(model_vector_a_mock),

    # totalDNA = rep(10, nrow(p_samp_all %>% filter(tech_rep == 1))), #vector of total DNA concentrations for each community observed; this is arbitrary until we have field measurements
    # 
    # up_idx =  as.array(seq(2, nrow(p_samp_all %>% filter(tech_rep == 1)), 2)), #vector of length (N_obs_samp_small/2), giving index of upstream sites
    # down_idx =  as.array(seq(1, nrow(p_samp_all %>% filter(tech_rep == 1)), 2)), #vector of length (N_obs_samp_small/2), giving index of downstream sites
    
    # Priors
    alpha_prior = c(0,0.1),  # normal prior
    beta_prior = c(0,5),    # normal prior
    tau_prior = c(1,2)   # gamma prior
  )
  