####Functions to simulate communities, and the observation of those communities with eDNA


#### PROCESS MODEL
SimCommunityProcess <- function(Nspp = 4,  #number of species
                                Nt = 10,    #number of timepoints
                                Ncreek = 5, #number of independent creeks surveyed 
                                Nbiol = 3, #number of biological replicates
                                Nindiv_total = 1e5, #total population size at starting time
                                Seed = runif(1, 0,1e5)
){
  library(MCMCpack)
  library(tidyverse)
  seed <- Seed
  set.seed(Seed)
  
  
  x0_mean = rmultinom(1,  #mean starting number of individuals for each of the Nspp species, across all creeks
                  Nindiv_total, 
                  rdirichlet(Ncreek, rep(5, times = Nspp))) %>% 
    matrix(nrow = Nspp, byrow = T) #proportion of indiv for each species; creeks in columns, species in rows

  culvert_effect = rnorm(n = Nspp, sd = .3) #effect of culvert 

  trend = sample(c("random", "increasing", "decreasing", "sinusoidal"), Nspp, replace = T, prob = c(.35,.2,.35,.1))
  slope <- rep(NA, Nspp)  #slope parameter over series of timepoints
  slope[trend == "random"] <- 0
  slope[trend == "increasing"] <- runif(sum(trend == "increasing"), 0, 20)
  slope[trend == "decreasing"] <- runif(sum(trend == "decreasing"), -20, 0)
  #do sinusoidal in dataframe itself, for ease of coding
  
  #create dataframe
  Process <- expand.grid(
    i = 1:Nspp,  #species
    t = 1:Nt, #time
    j = 1:Ncreek, #creek
    k = 1:Nbiol, #replicates
    x_0 = NA #initial observations
  )  
  for (zz in 1:Nspp){
    for (ww in 1:Ncreek){
      Process$x_0[Process$i == zz & Process$j == ww] <- rep(rpois(Nbiol, x0_mean[zz, ww]), each = Nt)  #starting populations for each species in each creek, each species grouped around a common mean  
    }
  }
    #add in info on slope/trend  
    Process <- Process %>%
      left_join(data.frame(i = 1:Nspp, slope = slope, trend = trend))

  #for each sample, the expected number of individuals at time t is:
  Process$lambda <- (Process$x_0 +   
                       Process$slope*Process$t +
                       Process$x_0 * rnorm(nrow(Process), 0, .1)) #some process variability 
  
  Process$lambda[Process$trend == "sinusoidal"] <- (
    Process$x_0[Process$trend == "sinusoidal"] *
      (sin(Process$t[Process$trend == "sinusoidal"])) +
      rlnorm(sum(Process$trend == "sinusoidal"), 1, 1))
  
  Process$lambda[Process$lambda < 0] <- 0 #enforce lower limit of zero (consequence of using a linear model for temporal trend)                       
  
  #add in upstream expected values as deterministic function of downstream expected values; culvert effect is here a property of species
  Process <- Process %>% 
    left_join(data.frame(culvert_effect, "i" = 1:Nspp)) %>% 
    mutate(upstream = lambda + (lambda * culvert_effect))

  return(Process)
}

#Process <- SimCommunityProcess(Nt = 1)
#visualize
# Process %>%
#   filter(i == 3, j ==1) %>%
#   ggplot(aes(x= t, y = lambda, color = j)) +
#     geom_point() +
#     geom_point(aes(x = t, y = upstream), color = "grey")

#### OBSERVATION MODEL
SimCommunityObservation <- function(Npcr = 35, #PCR cycles
                                    processOutput = Process, #output from SimCommunityProcess.R
                                    NtechReps = 3, #number of biological replicates per sample
                                    #NlocationsPerCreek = 2, #e.g., upstream and downstream of a culvert
                                    #a = runif(Nspp, .5, .8),#rbeta(Nspp, 20, 5), #distribution of species-level amplification efficiencies 
                                    ReadsPerSample = 2e05,
                                    AmpliconDistrib = "NegativeBinomial", # "Poisson" or "NegativeBinomial"
                                    Seed = runif(1, 0, 10000),
                                    MOCK = TRUE
){
  library(tidyverse)
  set.seed(Seed)
  
  Nspp = length(unique(processOutput$i)) #this is dumb; assumes column named `i` for species
  a = rbeta(Nspp, 40, 10) #runif(Nspp, .6, .9)
  
  Observation <- expand_grid(Process, 
              m = 1:NtechReps) %>% 
    dplyr::select(i, t, j, k, m, lambda, upstream) %>% 
    rename(species = i, time = t, creek = j, biol = k, tech = m, downstream = lambda) %>% 
    pivot_longer(c(downstream, upstream), names_to = "site", values_to = "count") %>% 
    group_by(time,creek,biol,tech,site) %>%  #for each timepoint and creek and site and bottle and tech rep,  
    mutate(b_proportion = count/sum(count), #the proportion of DNA expected for each species i
           a = a) %>%  #amplification efficiency for each species
    ungroup() %>% 
    mutate(Namplicons = if(AmpliconDistrib == "Poisson"){
      rpois(nrow(.), lambda = b_proportion*(a + 1)^Npcr)} else 
        {if (AmpliconDistrib == "NegativeBinomial") rnbinom(nrow(.), mu = b_proportion*(a + 1)^Npcr, size = 50)} #here, asserting a precision term; higher values mean more precision. Real PCR is probably below 1, frustratingly.
      ) %>%   
    group_by(time,creek,biol,tech,site) %>%  ##for each timepoint and creek and site and bottle and tech rep,  
    mutate(proportion_amplicons = Namplicons/sum(Namplicons),
           Nreads = rmultinom(1, ReadsPerSample, proportion_amplicons), #observed reads, for fixed seq depth
           proportion_reads = Nreads/sum(Nreads))
  
  Observation$creek_time_site_idx <- match(paste(Observation$time, Observation$creek, Observation$site), 
                                           unique(paste(Observation$time, Observation$creek, Observation$site)))
  Observation$BiolSample_idx <- match(paste(Observation$time, Observation$creek, Observation$site, Observation$biol), 
                                           unique(paste(Observation$time, Observation$creek, Observation$site, Observation$biol)))
  Observation$PCRsample_idx <- match(paste(Observation$time, Observation$creek, Observation$site, Observation$biol, Observation$tech), 
                                      unique(paste(Observation$time, Observation$creek, Observation$site, Observation$biol, Observation$tech)))
  Observation$species_time_site_idx <- match(paste(Observation$time, Observation$species, Observation$site), 
                                             unique(paste(Observation$time, Observation$species, Observation$site)))  
  
  if (MOCK == TRUE){
    Mock <- expand_grid(species = 1:Nspp, 
                        tech = 1:NtechReps,
                        site = 1:3) %>% 
      group_by(tech, site) %>% 
      mutate(b_proportion = 1/Nspp,
             a = a) %>% 
      ungroup() %>% 
      mutate(Namplicons = if(AmpliconDistrib == "Poisson"){
        rpois(nrow(.), lambda = b_proportion*(a + 1)^Npcr)} else 
        {if (AmpliconDistrib == "NegativeBinomial") rnbinom(nrow(.), mu = b_proportion*(a + 1)^Npcr, size = 50)} #here, asserting a precision term; higher values mean more precision. Real PCR is probably lowish, frustratingly.
      ) %>%   
      group_by(tech) %>%  
      mutate(proportion_amplicons = Namplicons/sum(Namplicons),
             Nreads = rmultinom(1, ReadsPerSample, proportion_amplicons), #observed reads, for fixed seq depth
             proportion_reads = Nreads/sum(Nreads))
  }
  
  return(
    list(
      Observation = Observation %>% dplyr::select(-Namplicons, -proportion_amplicons),
      Mock = Mock %>% dplyr::select(-Namplicons, -proportion_amplicons),
      N_pcr_mock = Npcr)
    )
}

#obs <- SimCommunityObservation()
#visualize
# observed %>%
#   ungroup() %>%
#   filter(species == 3, creek ==1) %>%
#   dplyr::select(time, count, proportion_reads, site) %>%
#   pivot_longer(c(count, proportion_reads)) %>%
#   ggplot(aes(x= time, y = value, color = site)) +
#     geom_point() +
#     facet_wrap(~name, scales = "free_y")


makeDesign <- function(obs){
  require(tidyverse)
  
  alrTransform <- function(MOCK){
    require(tidyverse)
    require(compositions)
    
    p_mock <- MOCK %>% 
      dplyr::select(species, site, tech, b_proportion) %>% 
      rename(tech_rep = tech) %>% 
      pivot_wider(names_from = species, values_from = b_proportion, values_fill = 1e-9) %>% 
      ungroup() 
    colnames(p_mock)[3:(length(unique(MOCK$species))+2)] <- paste0("alr_", 1:length(unique(MOCK$species)))
    
    p_mock <- alr(p_mock[,3:ncol(p_mock)]) %>% as.matrix() %>% as.data.frame()
      p_mock[,length(unique(MOCK$species))] <- 0  #add reference zero expressly
      names(p_mock)[length(unique(MOCK$species))] <- paste0("alr_", length(unique(MOCK$species)))

    p_mock <-  cbind(MOCK %>% dplyr::select(tech, site) %>% distinct(),
            p_mock) %>% ungroup()
    names(p_mock)[1] <- "tech_rep"
   
            
    return(p_mock)
  }
  
  mock <- obs$Mock
  observed <- obs$Observation
  
  p_mock_all <- alrTransform(mock)
  
  mock <- mock %>% 
    dplyr::select(species, site, tech, Nreads) %>% 
    ungroup() %>% 
    #mutate(site = 1) %>% 
    mutate(species = paste0("sp_", species)) %>% 
    pivot_wider(names_from = species, values_from = Nreads, values_fill = 0)
  N_pcr_mock <- rep(obs$N_pcr_mock, nrow(p_mock_all)) #assumes all have the same Npcr
  
  
  p_samp_all <- observed %>% 
    ungroup() %>% 
    unite(time, creek, site, biol, col = "site") %>% 
    dplyr::select(site, tech, species, Nreads) %>% 
    rename(tech_rep = tech) %>% 
    mutate(species = paste0("sp_", species)) %>% 
    arrange(species) %>% 
    pivot_wider(names_from = species, values_from = Nreads, values_fill = 0)
  N_pcr_samp <- rep(43, nrow(p_samp_all))
  
  source(here("In_Progress","rosetta_calibration","Metabarcoding_make_design_files.R"), local = TRUE)
  
  return(stan_data)
}

#stan_data <- makeDesign(obs)




# makeInitBetaRaw <- function(obs){
#   require(tidyverse)
# 
#   alrTransform <- function(samp){
#     require(tidyverse)
#     require(compositions)
# 
#     p_samp <- obs %>%
#       ungroup() %>%
#       dplyr::select(species, time, creek, biol, tech, site, proportion_reads) %>%
#       filter(tech == 1) %>%
#       pivot_wider(names_from = species, values_from = proportion_reads)
#     colnames(p_samp)[6:(length(unique(obs$species))+5)] <- paste0("alr_", 1:length(unique(obs$species)))
#     p_samp <- p_samp %>% dplyr::select(starts_with("alr"))
#     p_samp <- alr(p_samp)
#     return(p_samp)
#   }
#   samp <- obs$Observation
#   p_samp <- alrTransform(samp)
#   return(p_samp)
# }






#initialization function for stan                  
stan_init_f3 <- function(N_species){
  A <- list(
    tau = runif(N_species-1,0.1,0.5),
    alpha_raw = runif(N_species-1,-0.1,0.1)
  )
  return(A)
}
stan_init_f4 <- function(N_species, Observation){
  
  A <- list(
    tau = runif(N_species-1,0.1,0.5),
    alpha_raw = runif(N_species-1,-0.1,0.1),
    beta_raw = t(as.matrix(makeInitBetaRaw(Observation$Observation)))
  )
  return(A)
}
