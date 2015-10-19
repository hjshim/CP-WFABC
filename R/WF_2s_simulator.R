# simulates the trajectory of a Wright-Fisher model with 2 selection coefficients

WF_trajectory<-function(N,t,fluc_t,j,t0,s1,s2,h=0.5,s_start,ploidy,N_sample=c(),sample_times=c(),max_sims) 
{ 
  # number of sampling time points
  nb_times = length(sample_times)
  N = round(N)
  
  # time-sampling check
  if(s1==s2){
    if(nb_times < 2) stop("Please provide length(sample_times) > 1 for time-sampling")
    if(t < 2) stop("Please provide t > 1 for time-sampling")
  } else {
    if(nb_times < 4) stop("Please provide length(sample_times) > 3 for change point detection")
    if(t < 4) stop("Please provide t > 3 for change point detection")
  }
  
  nsims=0
  repeat{
    nsims=nsims+1
    # initialize the trajectory
    x=numeric(t) 
    x[t0]=j
    # set 1st selection coefficient
    if(is.function(s1)) cur_s=s1() else cur_s=s1    
    # set Ne
    if(is.function(N)) cur_N=N() else cur_N=N    
    # loop over generations until fluc_t generations with s1
    for (i in (t0+1):fluc_t)
    {
      # set s1
      my_s=cur_s
      if(i<=s_start) my_s=0
      # calculate fitness
      wAA=1+my_s;wAa=1+my_s*h;waa=1
      wA=1+my_s;wa=1
      p=x[i-1]/(cur_N)
      # calculate sampling probabilities
      if (ploidy==2) {  
        prob=(wAA*p^2+wAa*p*(1-p))/(wAA*p^2+wAa*2*p*(1-p)+waa*(1-p)^2)
      }
      else {
        prob=wA*p/(wA*p+wa*(1-p))      
      }    
      # binomial sampling to the next generation
      x[i]=rbinom(1,cur_N,prob)
      # stop if allele is lost
      if (x[i]==0) break
    }
    # set 2nd selection coefficient
    if(is.function(s2)) cur_s=s2() else cur_s=s2    
    # loop over generations after fluc_t generations with s2
    for (m in (fluc_t+1):t)
    {
      # set s2
      my_s=cur_s
      # calculate fitness
      wAA=1+my_s;wAa=1+my_s*h;waa=1
      wA=1+my_s;wa=1
      p=x[m-1]/(cur_N)
      # calculate sampling probbilities
      if (ploidy==2) {  
        prob=(wAA*p^2+wAa*p*(1-p))/(wAA*p^2+wAa*2*p*(1-p)+waa*(1-p)^2)
      }
      else {
        prob=wA*p/(wA*p+wa*(1-p))      
      }    
      # binomial sampling to the next generation
      x[m]=rbinom(1,cur_N,prob)
      # stop if allele is lost
      if (x[m]==0) break
    }
    
    # do sampling
    sample_freq=rbinom(nb_times,N_sample,x[sample_times]/(cur_N))

    if (nsims>=max_sims) 
    {
      break
    }
  }
  
  # create result list
  res=list("s1"=s1,"s2"=s2,"fluc_t"=fluc_t,"N"=cur_N,"h"=h)
  res[["N_x"]]=x 
  res[["N_A2"]]=sample_freq
  res[["N_sample"]]=N_sample
  res[["times"]]=sample_times    
  res[["nsims"]]=nsims
  res
}

#WF_trajectory(N=100,t=100,fluc_t=50,j=1,t0=1,s1=0.5,s2=-0.5,,s_start=1,ploidy=1,nb_times=10,N_sample=100,sample_times=c(),max_sims=10)
#print(res)
