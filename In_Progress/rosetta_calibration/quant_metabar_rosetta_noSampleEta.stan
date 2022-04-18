data { // 
  
  int N_species; // Number of species in data
  int N_obs_samp;  // Number of observed samples 
  int N_obs_samp_small;  // Number of observed samples for individual sites.
  int N_obs_mock; // Number of observed mock samples

  // Observed data of community matrices
  int sample_data[N_obs_samp,N_species];
  // Observed data of mock community matrices
  int mock_data[N_obs_mock,N_species];  
  // True proportions for mock community
  matrix[N_obs_mock,N_species] alr_mock_true_prop ;
 // matrix[N_obs_mock_small,N_species] alr_mock_true_prop_small ;
    
  // Design matrices: field samples
  int N_b_samp_col;
  matrix[N_obs_samp,N_b_samp_col]  model_matrix_b_samp;
  matrix[N_obs_samp_small,N_b_samp_col]  model_matrix_b_samp_small;

  vector[N_obs_samp]  model_vector_a_samp;
  vector[N_obs_samp_small]  model_vector_a_samp_small;

  // Design matrices: mock community samples
  //int N_b_mock_col;
  //matrix[N_obs_mock,N_b_mock_col]  model_matrix_b_mock;
  vector[N_obs_mock]  model_vector_a_mock;
 //vector[N_obs_mock_small]  model_vector_a_mock_small;

  // vector[N_obs_samp_small] totalDNA; //total DNA concentrations for each field sample
  // int up_idx[N_obs_samp_small/2]; //index of upstream sites
  // int down_idx[N_obs_samp_small/2]; //index of downstream sites

  // Priors
  real alpha_prior[2]; // Parameters of normal distribution for prior on alphas
  real beta_prior[2]; // Parameters of normal distribution for prior on betas
  real tau_prior[2]; // Parameters of gamma distribution for prior on tau (observation precision)
}

transformed data {
  matrix[N_obs_samp,N_b_samp_col+1] model_matrix_samp;
  
  model_matrix_samp = append_col(model_matrix_b_samp, model_vector_a_samp);
}

parameters {
  real<lower=0> tau; // single overdispersion sd for multinomial.
  // vector<lower=0>[N_species-1] tau; // single overdispersion sd for multinomial.
  vector[N_species-1] alpha_raw;
  vector[N_b_samp_col] beta_raw[N_species-1];
  // vector[N_obs_samp] eta_samp_raw[N_species-1];
  vector[N_obs_mock] eta_mock_raw[N_species-1];
}

transformed parameters {
  vector[N_species] alpha; // coefficients
  vector[N_b_samp_col] beta[N_species]; // coefficients
  // vector[N_obs_samp] eta_samp[N_species]; // overdispersion coefficients
  vector[N_obs_mock] eta_mock[N_species]; // overdispersion coefficients
  matrix[N_obs_samp,N_species] mu_samp; // estimates, in log space
  matrix[N_obs_mock,N_species] mu_mock; // estimates, in log space
  
  { // local variables delcaration
  matrix[N_obs_samp,N_species] logit_val_samp;
  matrix[N_obs_mock,N_species] logit_val_mock;
  matrix[N_species,N_obs_samp] prob_samp_t;
  matrix[N_species,N_obs_mock] prob_mock_t;
    
    
    // Fixed effects components.
    alpha[1:(N_species-1)] = alpha_prior[1] + alpha_raw * alpha_prior[2]; 
          // non-centered param beta ~ normal(alpha_prior[1], alpha_prior[2])
    alpha[N_species] = 0; // final species is zero (reference species)
  
    beta[N_species] = rep_vector(0.0,N_b_samp_col);// final species is zero (reference species)
    for (l in 1:(N_species-1)){
        beta[l] = beta_prior[1] + beta_raw[l] * beta_prior[2]; 
        // non-centered param beta ~ normal(beta_prior[1], beta_prior[2])
    }

    // add over dispersion by making the itercept term (beta) have a normal distribution.
    //tau = rep_vector(tau_base,N_species-1);
    eta_mock[N_species] = rep_vector(0.0,N_obs_mock); // final species is zero (reference species)
    // eta_samp[N_species] = rep_vector(0.0,N_obs_samp); // final species is zero (reference species)
    // // random effects vector of vectors
    for (l in 1:(N_species-1)) {
      eta_mock[l] = eta_mock_raw[l] * tau ; // non-centered param eta_mock ~ normal(0,tau)
      // eta_samp[l] = eta_samp_raw[l] * tau[l] ; // non-centered param eta_samp ~ normal(0,tau)
    }

  // from betas, alphas, and etas we can calculate sample-specific mu
    for (n in 1:N_species) {
      logit_val_samp[,n] = model_matrix_samp * append_row(beta[n],alpha[n]); 
                              // eta_samp[n];
      logit_val_mock[,n] = alr_mock_true_prop[,n] + 
                                model_vector_a_mock * alpha[n]
                                + eta_mock[n];
    }
    for(m in 1:N_obs_samp){
      prob_samp_t[,m] = softmax(transpose(logit_val_samp[m,]));
    }
    for(m in 1:N_obs_mock){
      prob_mock_t[,m] = softmax(transpose(logit_val_mock[m,]));
    }
    
    for (n in 1:N_species) {
      mu_samp[,n] = transpose(prob_samp_t)[,n] ;
      mu_mock[,n] = transpose(prob_mock_t)[,n] ;
    // if(n==1){print("log_prob 1 ",log_prob);}
  }

    // print("MU_SAMP_row",mu_samp[1,]);
    // print("SUM_MU_SAMP",sum(mu_samp[1,]));
    // 
    // print("MU_SAMP_col",mu_samp[,1]);
    // print("SUM_MU_SAMP_col",sum(mu_samp[,1]));

  } // end local variables declaration
  // print("mu 1 ",mu[1]);
  // print(N_pcr_samp*log(1+alpha[1]));
}

model {
  for(i in 1:N_obs_samp){ 
    sample_data[i,] ~  multinomial(transpose(mu_samp[i,]));  
  }
  for(i in 1:N_obs_mock){
    mock_data[i,]   ~  multinomial(transpose(mu_mock[i,]));  
  }
  // Priors
  for(i in 1:(N_species-1)){
    beta_raw[i] ~ std_normal();      //prior of normal(beta_prior[1],beta_prior[2]);
    // eta_samp_raw[i] ~ normal(0, tau[i]); // N(0,tau)
    eta_mock_raw[i] ~ normal(0, tau); // N(0,tau)
  }
  alpha_raw ~ std_normal(); // prior of normal(alpha_prior[1],alpha_prior[2]);
  tau ~ gamma(tau_prior[1],tau_prior[2]); //
}

generated quantities {
  matrix[N_obs_samp_small,N_species] int_samp_small;
  // matrix[N_obs_samp_small,N_species] abs_abundance;
  // matrix[N_obs_samp_small/2,N_species] culvert_diff;  //difference in abundance, upstream-downstream
  // matrix[N_obs_samp_small/2,N_species] est_culvert_effect; //estimated effect size as function of baseline downstream abundance
      // estimates, in proportion space one for each field sample.
  
    { // local variables delcaration
    matrix[N_obs_samp_small,N_species] logit_val_samp_ppd;
    matrix[N_species,N_obs_samp_small] prob_samp_t_ppd;
  
      for(n in 1:N_species) {
        logit_val_samp_ppd[,n] = model_matrix_b_samp_small * beta[n] ; //+
      }
      for(m in 1:N_obs_samp_small){
        prob_samp_t_ppd[,m] = softmax(transpose(logit_val_samp_ppd[m,]));
      }
      for(n in 1:N_species) {
        int_samp_small[,n] = transpose(prob_samp_t_ppd)[,n] ;
      }
    } // End local variable declaration.

    // for(m in 1:N_obs_samp_small){
    //   abs_abundance[m,] = int_samp_small[m,] * totalDNA[m];  
    // }
    // 
    // culvert_diff = (abs_abundance[up_idx,] - abs_abundance[down_idx,]);
    // 
    // for (j in 1:(N_obs_samp_small/2)){
    //   est_culvert_effect[j,] = culvert_diff[j,] ./ abs_abundance[down_idx[j],];
    // }
    
}
