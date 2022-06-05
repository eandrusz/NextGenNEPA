data {
  // qPCR calibration info
  int N_std_curve; //number of observations, across all dilutions, in standard curve. For example if 5 dilutions and 3 replicates of each, N_std_curve == 15.
  int N_std_samples; //number of unique biological entities being tested (flattening technical replicates)
  int std_sample_idx[N_std_curve]; 
  vector[N_std_samples] known_concentration; //known concentration of each of these observations; natural scale
  vector[N_std_curve] y_ct_std_curve; //observed Ct values for standard curve
  
  // qPCR info from environmental samples
  int N_qPCR_envir; //number of individual observations, across all envir samples
  int N_biol_estimates; //number of concentration estimates for individual biological samples (water bottles) * unique species  [e.g., 5 bottles with 1 species and 2 bottles with a second species would be N = 7]
  // int<lower = 0> sp_idx_qPCR[N_biol_estimates]; //species index linking these qPCR observations to the amplicon/metabarcoding observations in the larger dataset
  int<lower = 0> obs_samp_idx[N_qPCR_envir];
  // int<lower = 0> obs_samp_small_idx[N_biol_estimates]; //site index linking qPCR observations to the unique bottles of water in the larger database
  real y_ct_envir[N_qPCR_envir]; //observed Ct values from qPCR of environmental samples, for a selected set of species that will tie estimates of posterior proportions to absolute abundances

}

transformed data {
  vector[N_std_samples] conc;
  
  conc = log10(known_concentration);
}

parameters {
  real<lower=0> beta_std_curve_0;
  real<upper=0> beta_std_curve_1;
  real<lower=0> sigma_std_curve;
  
  //qPCR environmental samples
  vector[N_biol_estimates] envir_concentration; //log10 scale
  real<lower=0> sigma_envir_qPCR;
  // vector[N_biol_estimates] mu_qPCR_envir;
  
}

transformed parameters{
  vector[N_std_samples] mu_std_curve;
  vector[N_biol_estimates] mu_qPCR_envir;
  
  //calibration
  for(i in 1:N_std_samples){
    mu_std_curve[i] = beta_std_curve_0 + beta_std_curve_1 * conc[i];
  }
  //qPCR environmental
  for(i in 1:N_biol_estimates){
    mu_qPCR_envir[i] = beta_std_curve_0 + beta_std_curve_1 * envir_concentration[i];
  }

  
}

model {
   for(i in 1:N_std_curve){
      y_ct_std_curve[i] ~ normal(mu_std_curve[std_sample_idx[i]], sigma_std_curve);
    }
  
   // qPCR enviro calculation
   for(i in 1:N_qPCR_envir){
      y_ct_envir[i] ~ normal(mu_qPCR_envir[obs_samp_idx[i]], sqrt(sigma_std_curve^2 + sigma_envir_qPCR^2));
   }


  sigma_std_curve ~ gamma(1,10);
  sigma_envir_qPCR ~ gamma(1,10);
  beta_std_curve_0 ~ normal(40, 1);
  beta_std_curve_1 ~ normal(-3, .5);
  envir_concentration ~ normal(0, 5); //log10 scale
   
}

