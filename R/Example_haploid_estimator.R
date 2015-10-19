# Example: CP_WFABC_haploid_estimator - estimates parameters of interest (s1, s2, CP) for haploid populations -> Influenza as an example
source("CP_WFABC_haploid_estimator.R")
N=176
t=105
t0=1
s_start=1
sample_times=c(1, 14, 27, 40, 53, 66, 79, 92, 105)
N_sample=rep(1000,9)
N_a1=c(0, 0, 1, 0, 1, 3, 96, 362, 772, 974, 990, 998, 995)
N_allele=rbind.data.frame(N_a1)
min_freq=0.02
max_sims=1
no_sim=100000
best_sim=1000
post_2D_M1=FALSE # require MASS package
post_graph=TRUE
CP_WFABC_haploid_estimator(N,t,t0,s_start,sample_times,N_sample,N_allele,min_freq,max_sims,no_sim,best_sim,set_seed=TRUE,post_graph,post_2D_M1)