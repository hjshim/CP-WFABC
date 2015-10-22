# Example: CP_WFABC_haploid_modelchoice - detects changing selection trajectories using model choice and estimates parameters of interest (s1, s2, CP) for haploid populations -> Influenza as an example with h fixed
source("CP_WFABC_haploid_modelchoice.R")

main <-function(){
  N=176
  t=105
  t0=1
  s_start=1
  sample_times=c(1, 14, 27, 40, 53, 66, 79, 92, 105)
  N_sample=rep(1000,9)
  N_allele = read.table("example_haploid.txt", sep = ",", check.names=FALSE)
  rownames(N_allele) <- c("HA_48","HA_1395","NA_582","NA_823","M_147","NS_820","NP_159","PB1_33")
  min_freq=0.02
  max_sims=1
  no_sim=1000000
  best_sim=1000
  post_2D_M1=FALSE # require MASS package
  post_graph=TRUE
  
  CP_WFABC_haploid_modelchoice(N,t,t0,s_start,sample_times,N_sample,N_allele,min_freq,max_sims,no_sim,best_sim,set_seed=TRUE,post_graph,post_2D_M1)
}

main()
