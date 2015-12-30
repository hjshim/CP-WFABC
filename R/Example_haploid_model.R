# Example: CP_WFABC_haploid_modelchoice - detects changing selection trajectories using model choice and estimates parameters of interest (s1, s2, CP) for haploid populations -> Influenza as an example with h fixed
source("CP_WFABC_haploid_modelchoice.R")

main <-function(){
  N_allele = read.table("example_haploid2.txt", sep = ",", check.names=FALSE)
  rownames(N_allele) <- c("PB1_1119","HA_1395","NP_1104","NP_1396")
  N=226
  sample_times=c(1, 14, 27, 40, 53, 66, 79, 92, 105, 118, 131, 144, 157)
  N_sample=rep(1000,13)
  min_freq=0.02 # ascertainment
  max_sims=1
  no_sim=10000000
  best_sim=1000
  post_2D_M1=FALSE # require MASS package
  post_graph=TRUE
  set_seed=TRUE

  for (n_s in 1:length(N_allele[,1])){
    CP_WFABC_haploid_modelchoice(N,sample_times,N_sample,N_allele[n_s,],min_freq,max_sims,no_sim,best_sim,set_seed,post_graph,post_2D_M1)
  }
}

main()
