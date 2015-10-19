# Example: CP_WFABC_haploid_modelchoice - detects changing selection trajectories using model choice and estimates parameters of interest (s1, s2, CP) for haploid populations -> Influenza as an example with h fixed
source("CP_WFABC_haploid_modelchoice.R")
N=176
t=105
t0=1
s_start=1
sample_times=c(1, 14, 27, 40, 53, 66, 79, 92, 105)
N_sample=rep(1000,9)
a1=c(1, 1, 5, 72, 380, 938, 947, 943, 923)
a2=c(0, 2, 16, 128, 520, 988, 999, 999, 999)
a3=c(0, 0, 0, 20, 199, 971, 988, 990, 983)
a4=c(1, 3, 96, 362, 772, 974, 990, 998, 995)
a5=c(0, 0, 46, 140, 136, 421, 771, 809, 842)
a6=c(0, 0, 0, 0, 0, 1, 134, 215, 517)
N_allele <- rbind.data.frame(a1,a2,a3,a4,a5,a6)
rownames(N_allele) <- c("HA_48","HA_1395","NA_582","NA_823","M_147","NS_820")
colnames(N_allele) <- c("1", "14", "27", "40", "53", "66", "79", "92", "105")
N_allele <- N_allele[order(N_allele[,1]),]
min_freq=0.02
max_sims=1
no_sim=10000
best_sim=100
post_2D_M1=FALSE # require MASS package
post_graph=TRUE
#de novo
CP_WFABC_haploid_modelchoice(N,t,t0,s_start,sample_times,N_sample,N_allele[which(N_allele[,1]<(min_freq*N_sample[1])),],min_freq,max_sims,no_sim,best_sim,set_seed=TRUE,post_graph,post_2D_M1)
#standing
if(max(N_allele[,1]) > (min_freq*N_sample[1])){
  for (n_s in (min_freq*N_sample[1]):max(N_allele[,1])){
    if(length(N_allele[which(N_allele[,1]==n_s),1]) != 0){
      CP_WFABC_haploid_modelchoice(N,t,t0,s_start,sample_times,N_sample,N_allele[which(N_allele[,1]==n_s),],min_freq,max_sims,no_sim,best_sim,set_seed=TRUE,post_graph,post_2D_M1)
    }
  }
} 

