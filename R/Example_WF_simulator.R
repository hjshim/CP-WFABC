# Example: WF_2s_simulator - simulates WF trajectories with a change in selection intensity from s1 to s2 at fluc_t
source("WF_2s_simulator.R")
N=1000
t=100
fluc_t=50
j=1
t0=1
s1=0.3
s2=-0.3
h=0.5
s_start=1
ploidy=1
N_sample=rep(100,10)
sample_times=c(1,12,23,34,45,56,67,78,89,100)
max_sims=1
WF_trajectory(N,t,fluc_t,j,t0,s1,s2,h,s_start,ploidy,N_sample,sample_times,max_sims) 