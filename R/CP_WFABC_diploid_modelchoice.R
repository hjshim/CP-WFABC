source("WF_2s_simulator.R")
source("mode.R")

CP_WFABC_diploid_modelchoice<-function(N=1000,h_fixed=TRUE,h_given=0.5,sample_times=c(),N_sample=c(),N_allele=data.frame(),min_freq=0.02,max_sims=1,no_sim=1000000,best_sim=1000,set_seed=TRUE,post_graph=FALSE,post_2D_M1=FALSE)
{
  original_parameters <- list(max_sims=max_sims,no_sim=no_sim,set_seed=set_seed,h_fixed=h_fixed,h_given=h_given,min_freq=min_freq)
  data_parameters <- list(best_sim=best_sim,post_graph=post_graph,post_2D_M1=post_2D_M1)
  initialize_parameters <- initialize_N(2,N,sample_times,N_sample,N_allele)
  result_M0 <- simulate_M0(original_parameters, initialize_parameters)
  result_M1 <- simulate_M1(original_parameters, initialize_parameters)
  Data_result <- calculate_Data(original_parameters,initialize_parameters,result_M0,result_M1,data_parameters)  
}

initialize_N<-function(ploidy,N,sample_times,N_sample,N_allele)
{
  # initialize population
  ploidy <- ploidy
  N <- round(N)
  nb_times <- length(sample_times)
  t <- sample_times[nb_times]

  # Check when sample allele is non-zero
  nonzero <- which(!N_allele == 0)
  t0 <- sample_times[nonzero[1]] # where allele appears with the same initial frequency as the data 
  s_start <- sample_times[nonzero[1]] # where selection starts
  
  # Initial allele frequency: de novo mutation or standing variation
  j=round((N_allele[,nonzero[1]]/N_sample[nonzero[1]])*N)
  if(j < 1) { j=1 }
  
  list(ploidy=ploidy,N=N,sample_times=sample_times,N_sample=N_sample,N_allele=N_allele,nb_times=nb_times,j=j,t=t,t0=t0,s_start=s_start)
}

simulate_M0<-function(original_parameters, initialize_parameters)
{
  # name parameters for initialize_N
  ploidy=initialize_parameters$ploidy
  N=initialize_parameters$N
  nb_times=initialize_parameters$nb_times
  j=initialize_parameters$j
  sample_times=initialize_parameters$sample_times
  N_sample=initialize_parameters$N_sample
  N_allele=initialize_parameters$N_allele
  t=initialize_parameters$t
  t0=initialize_parameters$t0
  s_start=initialize_parameters$s_start
  # name original parameters
  max_sims=original_parameters$max_sims
  no_sim=original_parameters$no_sim
  set_seed=original_parameters$set_seed
  h_fixed=original_parameters$h_fixed
  h_given=original_parameters$h_given
  min_freq=original_parameters$min_freq
  
  # reproducible random numbers
  if(set_seed==TRUE) { set.seed(123) }
  
  # Vectors&Matrices for randomly generated inputs for M0 simulated trajectories
  M0_all_s1 <- numeric(no_sim)
  if(h_fixed==FALSE){
    M0_all_h <- numeric(no_sim)
  } else {
    h <- h_given
  }
  
  # Vectors&Matrices for summary statistics from M0 simulated trajectories
  M0_all_max_CUSUM <- numeric(no_sim)
  M0_all_Fs <- matrix(nrow=(nb_times-1), ncol=no_sim)
  M0_all_Fs_sign <- matrix(nrow=(nb_times-1), ncol=no_sim)
#  file_M0_traj=paste("M0_sim_trajectory",N,t,t0,j,s_start,ploidy,nb_times,min_freq,max_sims,no_sim,".txt",sep="_")
  
  # M0 Simulated trajectories with uniform prior for s1, s2, CP (and h for diploid case)
  for(m in 1:no_sim){
    repeat{
      s_1 <- runif(1, -1, 1)
      CP <- floor(runif(1, t0+1, t-1))
      if(h_fixed==FALSE){
        h <- sample(c(0,0.5,1),1,replace=TRUE)
      } 
      res=WF_trajectory(N,t,CP,j,t0,s_1,s_1,h,s_start,ploidy,N_sample,sample_times,max_sims) 
      if((max(res$N_A2/N_sample)>=min_freq)) break
    }
#    cat(res$N_A2,"\n",sep=" ",file=file_M0_traj, append=TRUE)
    
    M0_all_s1[m]=res$s1
    if(h_fixed==FALSE){ 
      M0_all_h[m]=res$h
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
    
    M0_all_Fs[,m] <- Fs
    M0_all_Fs_sign[,m] <- Fs_sign
    M0_all_max_CUSUM[m] <- max_CUSUM
  }
  
  # Simulated: output files
#   file_M0_all_s1=paste("M0_sim_trajectory_all_s1",N,t,t0,j,s_start,ploidy,nb_times,min_freq,max_sims,no_sim,".txt",sep="_")
#   cat(M0_all_s1,sep=" ",file=file_M0_all_s1)
#   if(h_fixed==FALSE){
#     file_M0_all_h=paste("M0_sim_trajectory_all_h",N,t,t0,j,s_start,ploidy,nb_times,min_freq,max_sims,no_sim,".txt",sep="_")
#     cat(M0_all_h,sep=" ",file=file_M0_all_h)    
#   }
#   file_M0_all_max_CUSUM=paste("M0_sim_trajectory_all_max_CUSUM",N,t,t0,j,s_start,ploidy,nb_times,min_freq,max_sims,no_sim,".txt",sep="_")
#   cat(M0_all_max_CUSUM,sep=" ",file=file_M0_all_max_CUSUM)
#   file_M0_all_Fs=paste("M0_sim_trajectory_all_Fs",N,t,t0,j,s_start,ploidy,nb_times,min_freq,max_sims,no_sim,".txt",sep="_")
#   cat(M0_all_Fs,sep=" ",file=file_M0_all_Fs)
#   file_M0_all_Fs_sign=paste("M0_sim_trajectory_all_Fs_sign",N,t,t0,j,s_start,ploidy,nb_times,min_freq,max_sims,no_sim,".txt",sep="_")
#   cat(M0_all_Fs_sign,sep=" ",file=file_M0_all_Fs_sign)  

  # M0 Simulated: prior graphs
  pdf(paste("M0_prior_s1",rownames(N_allele),".pdf",sep="_"),8,8)
  plot(density(M0_all_s1),lwd=2,main=paste("M0: Ne=",N,"; s1 prior"),xlab="s1")
  dev.off()
   if(h_fixed==FALSE){
    pdf(paste("M0_prior_h",rownames(N_allele),".pdf",sep="_"),8,8)
    barplot(c(length(which(M0_all_h==0)),length(which(M0_all_h==0.5)),length(which(M0_all_h==1))),names.arg=c("0","0.5","1"),xlab="h",main=paste("M0: Ne=",N,"; h prior"))
    dev.off()
  }

  if(h_fixed==FALSE){
    list(M0_all_Fs=M0_all_Fs,M0_all_Fs_sign=M0_all_Fs_sign,M0_all_max_CUSUM=M0_all_max_CUSUM,M0_all_s1=M0_all_s1,M0_all_h=M0_all_h)
  } else {
    list(M0_all_Fs=M0_all_Fs,M0_all_Fs_sign=M0_all_Fs_sign,M0_all_max_CUSUM=M0_all_max_CUSUM,M0_all_s1=M0_all_s1)
  }
}
  
simulate_M1<-function(original_parameters, initialize_parameters)
{
  # name parameters for initialize_N
  ploidy=initialize_parameters$ploidy
  N=initialize_parameters$N
  nb_times=initialize_parameters$nb_times
  j=initialize_parameters$j
  sample_times=initialize_parameters$sample_times
  N_sample=initialize_parameters$N_sample
  N_allele=initialize_parameters$N_allele
  t=initialize_parameters$t
  t0=initialize_parameters$t0
  s_start=initialize_parameters$s_start
  # name original parameters
  max_sims=original_parameters$max_sims
  no_sim=original_parameters$no_sim
  set_seed=original_parameters$set_seed
  h_fixed=original_parameters$h_fixed
  h_given=original_parameters$h_given
  min_freq=original_parameters$min_freq
  
  # reproducible random numbers
  if(set_seed==TRUE) { set.seed(123) }
  
  # Vectors&Matrices for randomly generated inputs for M1 simulated trajectories
  M1_all_s1 <- numeric(no_sim)
  M1_all_s2 <- numeric(no_sim)
  M1_all_CP <- numeric(no_sim)
	if(h_fixed==FALSE){
    M1_all_h <- numeric(no_sim)
	} else {
	  h <- h_given
	}
#  file_M1_traj=paste("M1_sim_trajectory",N,t,t0,j,s_start,ploidy,nb_times,min_freq,max_sims,no_sim,".txt",sep="_")
  
  # Vectors&Matrices for summary statistics from M1 simulated trajectories
  M1_all_max_CUSUM <- numeric(no_sim)
  M1_all_Fs <- matrix(nrow=(nb_times-1), ncol=no_sim)
  M1_all_Fs_sign <- matrix(nrow=(nb_times-1), ncol=no_sim)

  # M1 Simulated trajectories with uniform prior for s1, s2, CP (and h for diploid case)
  for(m in 1:no_sim){
    repeat{
    s_1 <- runif(1, -1, 1)
    s_2 <- runif(1, -1, 1)
    CP <- floor(runif(1, t0+1, t-1))
    if(h_fixed==FALSE){
      h <- sample(c(0,0.5,1),1,replace=TRUE)
    } 
    res=WF_trajectory(N,t,CP,j,t0,s_1,s_2,h,s_start,ploidy,N_sample,sample_times,max_sims) 
    if((max(res$N_A2/N_sample)>=min_freq) & (res$N_x[CP]!=N) & (res$N_x[CP]!=0)) break
    }
 #   cat(res$N_A2,"\n",sep=" ",file=file_M1_traj, append=TRUE)
    
    M1_all_s1[m]=res$s1
    M1_all_s2[m]=res$s2
    M1_all_CP[m]=res$fluc_t
    if(h_fixed==FALSE){ 
      M1_all_h[m]=res$h
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
    
    M1_all_Fs[,m] <- Fs
    M1_all_Fs_sign[,m] <- Fs_sign
    M1_all_max_CUSUM[m] <- max_CUSUM
  }
  
  # M1 Simulated: output files
#   file_M1_all_s1=paste("M1_sim_trajectory_all_s1",N,t,t0,j,s_start,ploidy,nb_times,min_freq,max_sims,no_sim,".txt",sep="_")
#   cat(M1_all_s1,sep=" ",file=file_M1_all_s1)
#   file_M1_all_s2=paste("M1_sim_trajectory_all_s2",N,t,t0,j,s_start,ploidy,nb_times,min_freq,max_sims,no_sim,".txt",sep="_")
#   cat(M1_all_s2,sep=" ",file=file_M1_all_s2)
#   file_M1_all_CP=paste("M1_sim_trajectory_all_CP",N,t,t0,j,s_start,ploidy,nb_times,min_freq,max_sims,no_sim,".txt",sep="_")
#   cat(M1_all_CP,sep=" ",file=file_M1_all_CP)
#   if(h_fixed==FALSE){
#     file_M1_all_h=paste("M1_sim_trajectory_all_h",N,t,t0,j,s_start,ploidy,nb_times,min_freq,max_sims,no_sim,".txt",sep="_")
#     cat(M1_all_h,sep=" ",file=file_M1_all_h)    
#   }
#   file_M1_all_max_CUSUM=paste("M1_sim_trajectory_all_max_CUSUM",N,t,t0,j,s_start,ploidy,nb_times,min_freq,max_sims,no_sim,".txt",sep="_")
#   cat(M1_all_max_CUSUM,sep=" ",file=file_M1_all_max_CUSUM)
#   file_M1_all_Fs=paste("M1_sim_trajectory_all_Fs",N,t,t0,j,s_start,ploidy,nb_times,min_freq,max_sims,no_sim,".txt",sep="_")
#   cat(M1_all_Fs,sep=" ",file=file_M1_all_Fs)
#   file_M1_all_Fs_sign=paste("M1_sim_trajectory_all_Fs_sign",N,t,t0,j,s_start,ploidy,nb_times,min_freq,max_sims,no_sim,".txt",sep="_")
#   cat(M1_all_Fs_sign,sep=" ",file=file_M1_all_Fs_sign)
  
  # M1 Simulated: prior graphs
  pdf(paste("M1_prior_s1",rownames(N_allele),".pdf",sep="_"),8,8)
  plot(density(M1_all_s1),lwd=2,main=paste("M1: Ne=",N,"; s1 prior"),xlab="s1")
  dev.off()
  pdf(paste("M1_prior_s2",rownames(N_allele),".pdf",sep="_"),8,8)
  plot(density(M1_all_s2),lwd=2,main=paste("M1: Ne=",N,"; s2 prior"),xlab="s2")
  dev.off()
  pdf(paste("M1_prior_CP",rownames(N_allele),".pdf",sep="_"),8,8)
  plot(density(M1_all_CP),lwd=2,main=paste("M1: Ne=",N,"; CP prior"),xlab="CP")
  dev.off()
  if(h_fixed==FALSE){
    pdf(paste("M1_prior_h",rownames(N_allele),".pdf",sep="_"),8,8)
    barplot(c(length(which(M1_all_h==0)),length(which(M1_all_h==0.5)),length(which(M1_all_h==1))),names.arg=c("0","0.5","1"),xlab="h",main=paste("M1: Ne=",N,"; h prior"))
    dev.off()
  }
  
  if(h_fixed==FALSE){
    list(M1_all_Fs=M1_all_Fs,M1_all_Fs_sign=M1_all_Fs_sign,M1_all_max_CUSUM=M1_all_max_CUSUM,M1_all_s1=M1_all_s1,M1_all_s2=M1_all_s2,M1_all_CP=M1_all_CP,M1_all_h=M1_all_h)
  } else {
    list(M1_all_Fs=M1_all_Fs,M1_all_Fs_sign=M1_all_Fs_sign,M1_all_max_CUSUM=M1_all_max_CUSUM,M1_all_s1=M1_all_s1,M1_all_s2=M1_all_s2,M1_all_CP=M1_all_CP)
  }
}
  
calculate_Data<-function(original_parameters,initialize_parameters,result_M0,result_M1,data_parameters)
{
  # name parameters for initialize_N
  ploidy=initialize_parameters$ploidy
  N=initialize_parameters$N
  nb_times=initialize_parameters$nb_times
  j=initialize_parameters$j
  sample_times=initialize_parameters$sample_times
  N_sample=initialize_parameters$N_sample
  N_allele=initialize_parameters$N_allele
  t=initialize_parameters$t
  t0=initialize_parameters$t0
  s_start=initialize_parameters$s_start
  # name original parameters
  max_sims=original_parameters$max_sims
  no_sim=original_parameters$no_sim
  set_seed=original_parameters$set_seed
  h_fixed=original_parameters$h_fixed
  h_given=original_parameters$h_given
  min_freq=original_parameters$min_freq
  # name Data parameters
  best_sim=data_parameters$best_sim
  post_graph=data_parameters$post_graph
  post_2D_M1=data_parameters$post_2D_M1
  
  # name results M0
  M0_all_Fs=result_M0$M0_all_Fs
  M0_all_Fs_sign=result_M0$M0_all_Fs_sign
  M0_all_max_CUSUM=result_M0$M0_all_max_CUSUM
  M0_all_s1=result_M0$M0_all_s1
  if(h_fixed==FALSE){
    M0_all_h=result_M0$M0_all_h
  }
  # name results M1
  M1_all_Fs=result_M1$M1_all_Fs
  M1_all_Fs_sign=result_M1$M1_all_Fs_sign
  M1_all_max_CUSUM=result_M1$M1_all_max_CUSUM
  M1_all_s1=result_M1$M1_all_s1
  M1_all_s2=result_M1$M1_all_s2
  M1_all_CP=result_M1$M1_all_CP
  if(h_fixed==FALSE){
    M1_all_h=result_M1$M1_all_h
  }
  
  # Summary of model choice & parameter inference result
  file_result=paste("Summary",N,ploidy,".txt",sep="_")
  
  # Oberved values
    Allele_f <- as.numeric(N_allele / N_sample)
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
    
    # ABC model choice
  	M0_all_d_ss <- numeric(no_sim)
  	M1_all_d_ss <- numeric(no_sim)
  	all_d_ss <- numeric(no_sim+no_sim)
  	all_model <- c(rep(0,no_sim),rep(1,no_sim))
  	best_model <- numeric(best_sim)
    
  	# Calculate d(D,Dobs)
  	for(m in 1:no_sim){
  	  d_ss=0
  	  if(M0_all_max_CUSUM[m]==max_CUSUM){
  	    for (i in 1:(nb_times-1))
  	    {
  	      d_ss <- d_ss + (M0_all_Fs[i,m]*M0_all_Fs_sign[i,m]-Fs[i]*Fs_sign[i])^2
  	    }
  	  } else 
  	  {
  	    d_ss <- Inf
  	  }
  	  M0_all_d_ss[m]=sqrt(d_ss)
  	}
  	
  	for(m in 1:no_sim){
  	  d_ss=0
  	  if(M1_all_max_CUSUM[m]==max_CUSUM){
  	    for (i in 1:(nb_times-1))
  	    {
  	      d_ss <- d_ss + (M1_all_Fs[i,m]*M1_all_Fs_sign[i,m]-Fs[i]*Fs_sign[i])^2
  	    }
  	  } else 
  	  {
  	    d_ss <- Inf
  	  }
  	  M1_all_d_ss[m]=sqrt(d_ss)
  	}
    
    # ABC model choice: Calculate Bayes factor M1/M0 = changing selection/constant selection
  	all_d_ss=c(M0_all_d_ss,M1_all_d_ss)
  	best_index <- head(order(abs(all_d_ss)),round(best_sim))
  	M0_best_index = best_index[best_index[]<=no_sim]
  	M1_best_index = (best_index[best_index[]>no_sim])-no_sim
  	best_model=all_model[best_index]
  	
  	posterior_model1=sum(best_model==1)/best_sim
  	BF=posterior_model1/(1-posterior_model1)
  	cat(paste(row.names(N_allele),"\t"),file=file_result,append=TRUE)
  	cat(paste(round(posterior_model1,3),"\t"),file=file_result,append=TRUE)
  	cat(paste(round(BF,3),"\t"),file=file_result,append=TRUE)
    
    # Parameter estimation
  	M0_best_s1=M0_all_s1[M0_best_index] 	
  	M1_best_s1=M1_all_s1[M1_best_index] 
  	M1_best_s2=M1_all_s2[M1_best_index]
  	M1_best_CP=M1_all_CP[M1_best_index]

  	if(length(M0_best_s1)==0 | length(M0_best_s1)==1) {
  	  M0_best_s1_mode <- mean(M0_best_s1)
  	  cat(paste(round(M0_best_s1_mode,3),"\t"),file=file_result,append=TRUE)
  	  if(h_fixed==FALSE){
  	    M0_best_h=M0_all_h[M0_best_index]
  	    M0_best_h_mode <- mean(M0_best_h)
  	    cat(paste(round(M0_best_h_mode,3),"\t"),file=file_result,append=TRUE)
  	  }
  	}	else {
  	  M0_best_s1_mode <- mode(M0_best_s1)
  	  cat(paste(round(M0_best_s1_mode,3),"\t"),file=file_result,append=TRUE)
  	  # Posterior graphs
  	  if(post_graph==TRUE){
    	  pdf(paste("M0_posterior_s1",rownames(N_allele),".pdf",sep="_"),8,8)
    	  plot(density(t(M0_best_s1)),lwd=2,main=paste("Allele=",rownames(N_allele),"; BF=",signif(BF,3),"; M0: s1 posterior mode=",signif(M0_best_s1_mode,3)),xlab="s1")
    	  dev.off()
  	  }
  	  if(h_fixed==FALSE){
  	    M0_best_h=M0_all_h[M0_best_index]
  	    M0_best_h_mode <- mode(M0_best_h)
  	    cat(paste(round(M0_best_h_mode,3),"\t"),file=file_result,append=TRUE)
        # Posterior graphs
  	    if(post_graph==TRUE){
    	    pdf(paste("M0_posterior_h",rownames(N_allele),".pdf",sep="_"),8,8)
    	    barplot(c(length(which(M0_best_h==0)),length(which(M0_best_h==0.5)),length(which(M0_best_h==1))),names.arg=c("0","0.5","1"),xlab="h",main=paste("Allele=",rownames(N_allele),"; BF=",signif(BF,3),"; M0: h posterior mode=",round(M0_best_h_mode,1)))
    	    dev.off()
  	    }
  	  }
  	}
    
  	if(length(M1_best_s1)==0 | length(M1_best_s1)==1) {
  	  M1_best_s1_mode <- mean(M1_best_s1)
  	  M1_best_s2_mode <- mean(M1_best_s2)
  	  M1_best_CP_mode <- mean(M1_best_CP)
  	  cat(paste(round(M1_best_s1_mode,3),"\t"),file=file_result,append=TRUE)
  	  cat(paste(round(M1_best_s2_mode,3),"\t"),file=file_result,append=TRUE)
  	  cat(paste(round(M1_best_CP_mode,3),"\t"),file=file_result,append=TRUE)
  	  if(h_fixed==FALSE){
  	    M1_best_h=M1_all_h[M1_best_index]
  	    M1_best_h_mode <- mean(M1_best_h)
  	    cat(paste(round(M1_best_h_mode,3),"\t"),file=file_result,append=TRUE)
  	  }
  	} else {
  	  M1_best_s1_mode <- mode(M1_best_s1)
  	  M1_best_s2_mode <- mode(M1_best_s2)
  	  M1_best_CP_mode <- mode(M1_best_CP)
  	  cat(paste(round(M1_best_s1_mode,3),"\t"),file=file_result,append=TRUE)
  	  cat(paste(round(M1_best_s2_mode,3),"\t"),file=file_result,append=TRUE)
  	  cat(paste(round(M1_best_CP_mode,3),"\t"),file=file_result,append=TRUE)
      # Posterior graphs
  	  if(post_graph==TRUE){
    	  pdf(paste("M1_posterior_s1",rownames(N_allele),".pdf",sep="_"),8,8)
    	  plot(density(t(M1_best_s1)),lwd=2,main=paste("Allele=",rownames(N_allele),"; BF=",signif(BF,3),"; M1: s1 posterior mode=",signif(M1_best_s1_mode,3)),xlab="s1")
    	  dev.off()
    	  pdf(paste("M1_posterior_s2",rownames(N_allele),".pdf",sep="_"),8,8)
    	  plot(density(t(M1_best_s2)),lwd=2,main=paste("Allele=",rownames(N_allele),"; BF=",signif(BF,3),"; M1: s2 posterior mode=",signif(M1_best_s2_mode,3)),xlab="s2")
    	  dev.off()
    	  pdf(paste("M1_posterior_CP",rownames(N_allele),".pdf",sep="_"),8,8)
    	  plot(density(t(M1_best_CP)),lwd=2,main=paste("Allele=",rownames(N_allele),"; BF=",signif(BF,3),"; M1: CP posterior mode=",signif(M1_best_CP_mode,3)),xlab="CP")
    	  dev.off()
  	  }
  	  if(h_fixed==FALSE){
  	    M1_best_h=M1_all_h[M1_best_index]
  	    M1_best_h_mode <- mode(M1_best_h)
  	    cat(paste(round(M1_best_h_mode,3),"\t"),file=file_result,append=TRUE)
  	    # Posterior graphs
  	    if(post_graph==TRUE){
    	    pdf(paste("M1_posterior_h",rownames(N_allele),".pdf",sep="_"),8,8)
    	    barplot(c(length(which(M1_best_h==0)),length(which(M1_best_h==0.5)),length(which(M1_best_h==1))),names.arg=c("0","0.5","1"),xlab="h",main=paste("Allele=",rownames(N_allele),"; BF=",signif(BF,3),"; M1: h posterior mode=",round(M1_best_h_mode,1)))
    	    dev.off()
  	    }
  	  }
  	  if(post_2D_M1==TRUE){
  	    library(MASS)
  	    pdf(paste("M1_2D_posterior_s1_CP",rownames(N_allele),".pdf",sep="_"),8,8)
  	    p <- kde2d(M1_best_s1,M1_best_CP,n=500)
  	    image(p,xlab="s1",ylab="CP",main=paste("Allele=",row.names(N_allele),"; BF=",signif(BF,3),"; M1: 2D posterior for s1 and CP"))
  	    dev.off()
  	    pdf(paste("M1_2D_posterior_s2_CP",rownames(N_allele),".pdf",sep="_"),8,8)
  	    p <- kde2d(M1_best_s2,M1_best_CP,n=500)
  	    image(p,xlab="s2",ylab="CP",main=paste("Allele=",row.names(N_allele),"; BF=",signif(BF,3),"; M1: 2D posterior for s2 and CP"))
  	    dev.off()
  	    pdf(paste("M1_2D_posterior_s1_s2",rownames(N_allele),".pdf",sep="_"),8,8)
  	    z <- kde2d(M1_best_s1,M1_best_s2,n=500)
  	    image(z,xlab="s1",ylab="s2",main=paste("Allele=",row.names(N_allele),"; BF=",signif(BF,3),"; M1: 2D posterior for s1 and s2"))
  	    dev.off()
  	  }
  	}
    
  	cat("\n",file=file_result,append=TRUE)
}


