# Example: CP_WFABC_diploid_estimator - estimates parameters of interest (s1, s2, CP and h) for diploid populations -> Panaxia dominula as an example with h estimated
source("CP_WFABC_diploid_estimator.R")
N=500
t=60
t0=1
h_fixed=FALSE
h_given=1
s_start=1
sample_times=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,49,50,51,52,53,54,55,56,57,58,59,60)
N_sample=c(234,922,410,538,992,744,1972,2682,1932,1016,2388,1162,3042,2178,2328,630,2616,3224,2766,960,378,344,46,118,62,162,74,100,262,1092,950,1292,680,462,1694,114,340,24,80,1244,1634,1050,4132,4580,966,452,284,664,910,828,838)
N_a1=c(26,63,22,30,45,48,84,100,69,29,88,29,108,56,67,7,78,148,102,21,7,7,1,1,0,2,0,0,3,38,31,9,5,1,11,2,3,0,1,6,19,7,34,45,23,8,8,10,17,5,11)
N_allele=rbind.data.frame(N_a1)
min_freq=0.02
max_sims=1
no_sim=100000
best_sim=1000
post_2D_M1=FALSE # require MASS package
post_graph=TRUE
CP_WFABC_diploid_estimator(N,t,t0,h_fixed,h_given,s_start,sample_times,N_sample,N_allele,min_freq,max_sims,no_sim,best_sim,set_seed=TRUE,post_graph,post_2D_M1)