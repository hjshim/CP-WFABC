# Example: CP_WFABC_haploid_modelchoice - detects changing selection trajectories using model choice and estimates parameters of interest (s1, s2, CP) for haploid populations -> Influenza as an example with h fixed
source("CP_WFABC_haploid_modelchoice.R")

main <-function(){
  N_allele = read.table("example_haploid.txt", sep = ",", check.names=FALSE)
  rownames(N_allele) <- c("HA_48","HA_1395","NA_582","NA_823","M_147","NS_820","NP_159","PB1_33")
  N=176
  t=105
  sample_times=c(1, 14, 27, 40, 53, 66, 79, 92, 105)
  N_sample=rep(1000,9)
  min_freq=0.02 # ascertainment
  max_sims=1
  no_sim=100000
  best_sim=100
  post_2D_M1=FALSE # require MASS package
  post_graph=TRUE

  for (n_s in 1:length(N_allele[,1])){
    nonzero <- which(!N_allele[n_s,] == 0)
    t0=sample_times[nonzero[1]] # where allele appears with the same initial frequency as the data 
    s_start=sample_times[nonzero[1]] # where selection starts
    CP_WFABC_haploid_modelchoice(N,t,t0,s_start,sample_times,N_sample,N_allele[n_s,],min_freq,max_sims,no_sim,best_sim,set_seed=TRUE,post_graph,post_2D_M1)
  }
}

main()
