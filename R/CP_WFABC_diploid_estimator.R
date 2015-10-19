source("WF_2s_simulator.R")
source("mode.R")

CP_WFABC_diploid_estimator<-function(N=1000,t=100,t0=1,h_fixed=TRUE,h_given=0.5,s_start=1,sample_times=c(),N_sample=c(),N_allele=data.frame(),min_freq=0.02,max_sims=1,no_sim=1000000,best_sim=1000,set_seed=TRUE,post_graph=FALSE,post_2D_M1=FALSE)
{
  # ploidy
  ploidy <- 2
  N <- round(N)
  
	# reproducible random numbers
	if(set_seed==TRUE) { set.seed(123) }
  nb_times <- length(sample_times)

  # Initial allele frequency: de novo mutation or standing variation
  if(N_allele[1,1]/N_sample[1] < min_freq){
    j=1
  } else {
    j=round(N_allele[1,1])
  }
  
  # Vectors&Matrices for randomly generated inputs for simulated trajectories
  all_s1 <- numeric(no_sim)
  all_s2 <- numeric(no_sim)
  all_CP <- numeric(no_sim)
	if(h_fixed==FALSE){
    all_h <- numeric(no_sim)
	} else {
	  h <- h_given
	}
  
  # Vectors&Matrices for summary statistics from simulated trajectories
  all_max_CUSUM <- numeric(no_sim)
  all_Fs <- matrix(nrow=(nb_times-1), ncol=no_sim)
  all_Fs_sign <- matrix(nrow=(nb_times-1), ncol=no_sim)

  
  # Simulated trajectories with uniform prior for s1, s2, CP (and h for diploid case)
  for(m in 1:no_sim){
    repeat{
    s_1 <- runif(1, -1, 1)
    s_2 <- runif(1, -1, 1)
    CP <- floor(runif(1, 2, t-1))
    if(h_fixed==FALSE){
      h <- sample(c(0,0.5,1),1,replace=TRUE)
    } 
    res=WF_trajectory(N,t,CP,j,t0,s_1,s_2,h,s_start,ploidy,N_sample,sample_times,max_sims) 
    if((max(res$N_A2/N_sample)>=min_freq) & (res$N_x[CP]!=N) & (res$N_x[CP]!=0)) break
    }
    
    all_s1[m]=res$s1
    all_s2[m]=res$s2
    all_CP[m]=res$fluc_t
    if(h_fixed==FALSE){ 
      all_h[m]=res$h
    }
    
    # calculate simulated Fs stats
    Allele_f <- res$N_A2 / res$N_sample
    Fs <- numeric(nb_times-1)
    Fs_sign <- numeric(nb_times-1)
    for (i in 1:(nb_times-1))
    {
      z <- (Allele_f[i+1]+Allele_f[i])/2
      stat <- (Allele_f[i+1]-Allele_f[i])^2/(z*(1-z))
      #Fs'
      H_mean <- c(res$N_sample[i+1],res$N_sample[i])
      Har_mean <- 1/mean(1/H_mean)
      stat <- abs((stat*(1-1/(2*Har_mean))-2/Har_mean)/((res$times[i+1]-res$times[i])*(1+stat/4)*(1-1/res$N_sample[i+1])))
      if(stat=="NaN")
        Fs[i]=0
      else
        Fs[i]=stat
      
      if(Allele_f[i+1]>=Allele_f[i])
        Fs_sign[i]=1
      else
        Fs_sign[i]=-1
    }
  
    # compute simulated CUSUM
    Fs_mean <- mean(Fs*Fs_sign)
    sum <- 0
    sum_i <- c(0)
    for (i in 1:(nb_times-1))
    {
      sum <- (sum + Fs[i]*Fs_sign[i]-Fs_mean) 
      sum_i <- append(sum_i,sum)
    }    
    abs_sum_i <- abs(sum_i)
    max_CUSUM <- nnet::which.is.max(abs_sum_i)
    
    all_Fs[,m] <- Fs
    all_Fs_sign[,m] <- Fs_sign
    all_max_CUSUM[m] <- max_CUSUM
  }
  
  # Simulated: output files
  file_all_s1=paste("M1_sim_trajectory_all_s1",N,t,t0,j,s_start,ploidy,nb_times,min_freq,max_sims,no_sim,".txt",sep="_")
  cat(all_s1,sep=" ",file=file_all_s1)
  file_all_s2=paste("M1_sim_trajectory_all_s2",N,t,t0,j,s_start,ploidy,nb_times,min_freq,max_sims,no_sim,".txt",sep="_")
  cat(all_s2,sep=" ",file=file_all_s2)
  file_all_CP=paste("M1_sim_trajectory_all_CP",N,t,t0,j,s_start,ploidy,nb_times,min_freq,max_sims,no_sim,".txt",sep="_")
  cat(all_CP,sep=" ",file=file_all_CP)
  if(h_fixed==FALSE){
    file_all_h=paste("M1_sim_trajectory_all_h",N,t,t0,j,s_start,ploidy,nb_times,min_freq,max_sims,no_sim,".txt",sep="_")
    cat(all_h,sep=" ",file=file_all_h)    
  }
  file_all_max_CUSUM=paste("M1_sim_trajectory_max_CUSUM",N,t,t0,j,s_start,ploidy,nb_times,min_freq,max_sims,no_sim,".txt",sep="_")
  cat(all_max_CUSUM,sep=" ",file=file_all_max_CUSUM)
  file_all_Fs=paste("M1_sim_trajectory_all_Fs",N,t,t0,j,s_start,ploidy,nb_times,min_freq,max_sims,no_sim,".txt",sep="_")
  cat(all_Fs,sep=" ",file=file_all_Fs)
  file_all_Fs_sign=paste("M1_sim_trajectory_all_Fs_sign",N,t,t0,j,s_start,ploidy,nb_times,min_freq,max_sims,no_sim,".txt",sep="_")
  cat(all_Fs_sign,sep=" ",file=file_all_Fs_sign)
  
  pdf(paste("M1_prior_s1",N,t,t0,j,s_start,ploidy,nb_times,min_freq,max_sims,no_sim,".pdf",sep="_"),8,8)
  plot(density(all_s1),lwd=2,main=paste("Ne=",N,"; s1 prior"),xlab="s1")
  dev.off()
  pdf(paste("M1_prior_s2",N,t,t0,j,s_start,ploidy,nb_times,min_freq,max_sims,no_sim,".pdf",sep="_"),8,8)
  plot(density(all_s2),lwd=2,main=paste("Ne=",N,"; s2 prior"),xlab="s2")
  dev.off()
  pdf(paste("M1_prior_CP",N,t,t0,j,s_start,ploidy,nb_times,min_freq,max_sims,no_sim,".pdf",sep="_"),8,8)
  plot(density(all_CP),lwd=2,main=paste("Ne=",N,"; CP prior"),xlab="CP")
  dev.off()
  if(h_fixed==FALSE){
    pdf(paste("M1_prior_h",N,t,t0,j,s_start,ploidy,nb_times,min_freq,max_sims,no_sim,".pdf",sep="_"),8,8)
    barplot(c(length(which(all_h==0)),length(which(all_h==0.5)),length(which(all_h==1))),names.arg=c("0","0.5","1"),xlab="h",main=paste("Ne=",N,"; h prior"))
    dev.off()
  }
  
  
  # Summary of model choice & parameter inference result
  file_result=paste("Summary",N,t,t0,j,s_start,ploidy,nb_times,min_freq,max_sims,no_sim,".txt",sep="_")
  
  # Oberved values
  for (N_a in 1:nrow(N_allele)) {    
  	# observed
    Allele_f <- as.numeric(N_allele[N_a,] / N_sample)
  	Fs <- numeric(nb_times-1)
  	Fs_sign <- numeric(nb_times-1)
  	
  	# calculate observed Fs stats
  	for (i in 1:(nb_times-1))
  	{
  	  z <- (Allele_f[i+1]+Allele_f[i])/2
  	  stat <- (Allele_f[i+1]-Allele_f[i])^2/(z*(1-z))
  	  #Fs'
  	  H_mean <- c(N_sample[i+1],N_sample[i])
  	  Har_mean <- 1/mean(1/H_mean)
  	  stat <- abs((stat*(1-1/(2*Har_mean))-2/Har_mean)/((sample_times[i+1]-sample_times[i])*(1+stat/4)*(1-1/N_sample[i+1])))
  	  if(stat=="NaN")
  	    Fs[i]=0
  	  else
  	    Fs[i]=stat
  	  
  	  if(Allele_f[i+1]>=Allele_f[i])
  	    Fs_sign[i]=1
  	  else
  	    Fs_sign[i]=-1
  	}
  	
  	# compute observed CUSUM
  	Fs_mean <- mean(Fs*Fs_sign)
  	sum <- 0
  	sum_i <- c(0)
  	for (i in 1:(nb_times-1))
  	{
  	  sum <- (sum + Fs[i]*Fs_sign[i]-Fs_mean) 
  	  sum_i <- append(sum_i,sum)
  	}    
  	abs_sum_i <- abs(sum_i)
  	max_CUSUM <- nnet::which.is.max(abs_sum_i)
    
    # Best estimates
  	all_d_ss <- numeric(no_sim)
  
  	# Calculate d(D,Dobs)
  	for(m in 1:no_sim){
  	  d_ss=0
  	  if(all_max_CUSUM[m]==max_CUSUM){
  	    for (i in 1:(nb_times-1))
  	    {
  	      d_ss <- d_ss + (all_Fs[i,m]*all_Fs_sign[i,m]-Fs[i]*Fs_sign[i])^2
  	    }
  	  } else 
  	  {
  	    d_ss <- Inf
  	  }
  	  all_d_ss[m]=sqrt(d_ss)
  	}
    
  	# Parameter estimation
  	best_index <- head(order(abs(all_d_ss)),round(best_sim))
  	
  	best_s1=all_s1[best_index] 
  	best_s1_mode <- mode(best_s1)
  	best_s2=all_s2[best_index]
  	best_s2_mode <- mode(best_s2)  
  	best_CP=all_CP[best_index]
  	best_CP_mode <- mode(best_CP)
    
    cat(paste(round(best_s1_mode,3),"\t"),file=file_result,append=TRUE)
    cat(paste(round(best_s2_mode,3),"\t"),file=file_result,append=TRUE)
    cat(paste(round(best_CP_mode,3),"\t"),file=file_result,append=TRUE)
    
  	# Posterior graph results
    if(post_graph==TRUE){
    	pdf(paste("M1_posterior_s1",row.names(N_allele[N_a,]),N,t,t0,j,s_start,ploidy,nb_times,min_freq,max_sims,no_sim,".pdf",sep="_"),8,8)
    	plot(density(t(best_s1)),lwd=2,main=paste("Allele=",row.names(N_allele[N_a,]),"; s1 posterior mode=",signif(best_s1_mode,3)),xlab="s1")
    	dev.off()
    	pdf(paste("M1_posterior_s2",row.names(N_allele[N_a,]),N,t,t0,j,s_start,ploidy,nb_times,min_freq,max_sims,no_sim,".pdf",sep="_"),8,8)
    	plot(density(t(best_s2)),lwd=2,main=paste("Allele=",row.names(N_allele[N_a,]),"; s2 posterior mode=",signif(best_s2_mode,3)),xlab="s2")
    	dev.off()
    	pdf(paste("M1_posterior_CP",row.names(N_allele[N_a,]),N,t,t0,j,s_start,ploidy,nb_times,min_freq,max_sims,no_sim,".pdf",sep="_"),8,8)
    	plot(density(t(best_CP)),lwd=2,main=paste("Allele=",row.names(N_allele[N_a,]),"; CP posterior mode=",signif(best_CP_mode,3)),xlab="CP")
    	dev.off()
    }
    
  	if(h_fixed==FALSE){
  	  best_h=all_h[best_index]
  	  best_h_mode <- mode(best_h)
  	  cat(paste(round(best_h_mode,3),"\t"),file=file_result,append=TRUE)
  	  # Posterior graphs
  	  if(post_graph==TRUE){
    	  pdf(paste("M1_posterior_h",N,t,t0,j,s_start,ploidy,nb_times,min_freq,max_sims,no_sim,row.names(N_allele[N_a,]),".pdf",sep="_"),8,8)
    	  barplot(c(length(which(best_h==0)),length(which(best_h==0.5)),length(which(best_h==1))),names.arg=c("0","0.5","1"),xlab="h",main=paste("Allele=",row.names(N_allele[N_a,]),"; h posterior mode=",round(best_h_mode,1)))
    	  dev.off()
  	  }
  	}
    
  	if(post_2D_M1==TRUE){
  	  library(MASS)
  	  pdf(paste("M1_2D_posterior_s1_CP",row.names(N_allele[N_a,]),N,t,t0,j,s_start,ploidy,nb_times,min_freq,max_sims,no_sim,".pdf",sep="_"),8,8)
  	  p <- kde2d(best_s1,best_CP,n=500)
  	  image(p,xlab="s1",ylab="CP",main=paste("Allele=",row.names(N_allele[N_a,]),"; M1: 2D posterior for s1 and CP"))
  	  dev.off()
  	  pdf(paste("M1_2D_posterior_s2_CP",row.names(N_allele[N_a,]),N,t,t0,j,s_start,ploidy,nb_times,min_freq,max_sims,no_sim,".pdf",sep="_"),8,8)
  	  p <- kde2d(best_s2,best_CP,n=500)
  	  image(p,xlab="s2",ylab="CP",main=paste("Allele=",row.names(N_allele[N_a,]),"; M1: 2D posterior for s2 and CP"))
  	  dev.off()
  	  pdf(paste("M1_2D_posterior_s1_s2",row.names(N_allele[N_a,]),N,t,t0,j,s_start,ploidy,nb_times,min_freq,max_sims,no_sim,".pdf",sep="_"),8,8)
  	  z <- kde2d(best_s1,best_s2,n=500)
  	  image(z,xlab="s1",ylab="s2",main=paste("Allele=",row.names(N_allele[N_a,]),"; M1: 2D posterior for s1 and s2"))
  	  dev.off()
  	}
  }
  
  cat("\n",file=file_result,append=TRUE)
  
}


