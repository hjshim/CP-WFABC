# README #

### What is this repository for? ###

CP-WFABC v1.0
	detects allele trajectories of changing selection from those of constant selection using ABC model choice, and jointly estimates the position of a change point as well as the strength of both corresponding selection coefficients (and dominance for diploid cases) using Wright-Fisher ABC methods


### How do I get set up? ###

# variable inputs
	N				population size (number of chromosomes: N 
					individuals for haploids and N/2 individuals for 
					diploids) (1000 by default)
	sample_times	vector of exact sampling times in generations 
	N_sample		vector of sample sizes in number of chromosomes
	min_freq		data ascertainment of a minimum frequency 
					at one of the sampling time points (0 for no 
					condition and 1 to condition on fixation)
	N_allele		data frame of observed SNPs in row (with rownames 
					=SNP names, column=sampled numbers))
	max_sims		maximum number of simulations to do before giving 
					up (1 by default)
	no_sim			number of simulated datasets to be created (1e6 
					by default)
	best_sim		number of best simulations to be used for 
					estimation and model choice (1e3 by default)
	set_seed		reproducible numbers (TRUE by default)
	post_graph		Posterior densities of M0 and M1 (FALSE by 
					default)
	post_2D_M1		2D posteriors of M1 estimates (s1&s2, s1&CP, 
					s2&CP) (FALSE by default)
	h_fixed			h to be fixed in diploid populations (TRUE by 
					default)
	h_given			h to be used if fixed (0.5 by default)
	
# fixed inputs
	ploidy			1 for haploids, 2 for diploids

# assumptions
	j				number of SNP appearing at time t0 in the 
					population (given as observed initial frequency 
					in N_allele * population size, except if j<1 
					given as j=1)
	t0				time where SNP appears in generations (same as 
					data given in N_allele)
	s_start			time in generation when selection starts (s=0 
					before s_start) (same as t0)

# output
	PDFs of prior graphs for simulated parameters
	Text files of summary of results (SNP_name 	M1_posterior_BF	M0_estimate_s1	(M0_estimate_h)	M1_estimate_s1	M1_estimate_ s2	M1_estimate_CP	(M1_estimate_h)) 
	PDFs of posterior graphs of parameters of interest (if TRUE)
	PDFs of 2D posterior graphs of parameters of interest (if TRUE)


# WF_2s_simulator{
	usage 
	{simulates a Wright-Fisher trajectory with changing selection or constant selection}
	arguments {N,t,fluc_t,j,t0,s1,s2,h,s_start,ploidy,N_sample,sample_times,max_sims}
		} 

# CP_WFABC_diploid_modelchoice{
	usuage 
	{detects allele trajectories of changing selection from those of constant selection using ABC model choice, and jointly estimates the position of a change point as well as the strength of both corresponding selection coefficients and dominance for a diploid population using Change-Point Wright-Fisher ABC methods}
	arguments {N,h_fixed,h_given,sample_times,N_sample,N_allele,min_freq,max_sims,no_sim,best_sim,set_seed,post_graph,post_2D_M1}
		} 

# CP_WFABC_haploid_modelchoice{
	usuage 
	{detects allele trajectories of changing selection from those of constant selection using ABC model choice, and jointly estimates the position of a change point as well as the strength of both corresponding selection coefficients for a haploid population using Change-Point Wright-Fisher ABC methods}
	arguments {N,sample_times,N_sample,N_allele,min_freq,max_sims,no_sim,best_sim,set_seed,post_graph,post_2D_M1}
		} 

### Contribution guidelines ###

# references
	Shim, H., Laurent, S., Matuszewski, S., Foll, M., Jeffrey D Jensen (in review) Detecting and quantifying changing selection intensities from time-sampled polymorphism data. 
	Foll, M.*, Shim, H.*, & Jeffrey D Jensen (2014b) WFABC: A Wright-Fisher ABC-Based Approach for Inferring Effective Population Sizes and Selection Coefficients from Time-Sampled Data. Molecular Ecology Resources, 1, 87-98.
	Foll, M., Poh, Y.-P., Renzette, N., Ferrer-Admetlla, A., Bank, C., Shim, H., Malaspinas, A.S., Ewing, G., Liu, P., Wegmann, D., Caffrey, D.R., Zeldovich, K.B., Bolon, D.N., Wang, J.P., Kowalik, T.F., Schiffer, C.A., Finberg, R.W. & Jensen, J.D. (2014a) Influenza Virus Drug Resistance: A Time-Sampled Population Genetics Perspective. PLoS Genetics, 10, e1004185.


# keyword
	population genetics; fluctuating selection; change point analysis; time-sampled data; approximate Bayesian computation; Wright-Fisher model; experimental design 


# examples
	see Example_diploid_model.R
	see Example_haploid_model.R

### Who do I talk to? ###

# packageMaintainer
	Hyunjin Shim [jinenstar@gmail.com]