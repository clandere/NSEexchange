####### Calculate mean of eta distribution                               ########
####### WARNING: this method takes a very long time.                     ########
####### Runtime:  For a 14 codon sequence (runtime scales exponentially  ########
#######           with sequence length):                                 ########
#######           user   system elapsed                                  ########
#######           3.652   0.012   3.689                                  ########

#path.string <- strsplit(getwd(),'/')[[1]]

#if(path.string[length(path.string)-1]!='trunk'){
#  warning('Warning: you must load simulation_R-ext.so manually unless you are in a directory one below trunk.\nExample: trunk/R, trunk/preston, trunk/rundir\n')
#}else{ dyn.load('../C/eta_simulation/simulation_R-ext.so')}

#rm(path.string)

calc.mean.exact <- function(gene.list,codon.parms,A_1 = 4,A_2 = 4,B = 0.0025){
  # Purpose: Calculate the likelihood of a sequence gene.list given parameters codon.parms based on 
  #          the Boltzman-Potts Pr, which we calculate directly without estimation of denominator.
  # Inputs: gene.list - list containting $phi.value - expression rate
  #                                      $name - gene ID
  #                                      $gene.dat, which contains
  #                                             $codon - nucleotide char sequence
  #                                             $aa - amino acid char sequence
  #                                             $c_index - representation of codon sequence using integers
  #                                             $mut_rate - vector of mutation rates corresponding to c_index
  #                                             $elong_rat - vector of elong rates...
  #                                             $count - amino acid/codon count
  #         codon.parms - dataframe containing AA, c_index, codon key, mut rate, elong rate, and 
  #                    # of synonyms for every codon
  #         A1 - Cost of translation initiation (in ATP's)
  #               DEFAULT - 4
  #         A2 - Cost of peptide elongation (in ATP's)
  #               DEFAULT - 4
  #         B - background Nonsense Error Rate (in s^-1)
  #               DEFAULT - 0.0025
  # Outputs: runtime - summary of runtime information for .C call of class "proc_time"
  #          lik - likelihood of the sequence given parameters
  #          eta_mean - mean eta of all synonymous codon sequences
  #          eta_var - variance of eta for all synonymous codon sequences
  # Usage: Gene expression in units of "y"-
  #           likelihood_list <- calc.lik.exact(gene.list,codon.parms)
  #        Gene expression in units of "phi"
  #           likelihood_list <- calc.lik.exact(gene.list,codon.parms,Ne=1.36e7,q=4.19e-7)
  #        Change NSE rate and/or initiation cost
  #           likelihood_list <- calc.lik.exact(gene.list,codon.parms,A_1=2,B=0.00515)
  
  #########
  # Debug #
  #########
  
  #A_1 = 4;A_2 = 4;B = 0.0025
  #phi=1;gene.list=glf[[4]]
  
  ############################
  # SEQUENCE DATA/PARAMETERS #
  ############################
  
  CODON_INDEX = gene.list$gene.dat$c_index 
  CODON_INDEX <- CODON_INDEX[-length(CODON_INDEX)] #remove the last codon(stop codon) -- stop codon is not used in code
  AA_COUNT = length(CODON_INDEX)
  PHI = 0
  
  if(AA_COUNT > 50) {
    stop(paste('Sequence', gene.list$name, 'is too large to use exact likelihood calculation.'))
  }
  
  #####################
  # GENOME PARAMETERS #
  #####################
  
  ELONGATION_RATES = codon.parms$elong_rate[sort.list(codon.parms$c_index)]
  MUTATION_RATES = codon.parms$mut_rate[sort.list(codon.parms$c_index)]
  AT_BIAS=0.5 #at_bias parameter (currently unused)
  GAMMA = 1 #transition/transversion parameter (currently unused)
  IGNORE = 0 #Integer number to ignore last # amino acids
  AA.VEC = codon.parms$aa[sort.list(codon.parms$c_index)]
  CODON.VEC = codon.parms$codon[sort.list(codon.parms$c_index)]
 
  #####################################
  # Initialize output dummy variables #
  #####################################
  
  LIK = 0
  ETA_MEAN = 0
  ETA_VAR = 0
  ETA_OBS = 0
  ETA_MIN = 0
  ETA_MAX = 0
  NUM_BINS = 20
  BINS = rep(0,NUM_BINS)
  
  ###################
  # Call C function #
  ###################
  
  runtime <- system.time(OUT_DATA <- .C("calc_exact", as.integer(CODON_INDEX),       AA_COUNT = length(CODON_INDEX),   as.double(PHI),
                                                      as.double(ELONGATION_RATES),   as.double(MUTATION_RATES),        as.double(1),
                                                      as.double(A_1),                as.double(A_2),                   as.double(AT_BIAS),
                                                      as.double(B),                  as.integer(IGNORE),               as.double(GAMMA),
                                                      as.double(1),                  as.double(LIK),                   as.double(ETA_MEAN),
                                                      as.double(ETA_VAR),            as.double(ETA_OBS),               as.double(ETA_MIN),
                                                      as.double(ETA_MAX),            as.double(BINS),                  as.integer(NUM_BINS),
                                                      as.character(AA.VEC),          as.character(CODON.VEC),          as.integer(length(codon.parms$codon))))
  
  #llik <- log(OUT_DATA[[14]])
  eta.mean <- OUT_DATA[[15]]
  eta.var <- OUT_DATA[[16]]
  eta.obs <- OUT_DATA[[17]]
  eta.min <- OUT_DATA[[18]]
  eta.max <- OUT_DATA [[19]]
  #eta.bins <- OUT_DATA[[20]]
  
  #normalize eta bins  
  #eta.bins = eta.bins/sum(eta.bins)
  
  #bin.lims <- seq(eta.min,eta.max,length.out=NUM_BINS)
  
  c(eta.mean=eta.mean,eta.var=eta.var,eta.min=eta.min,eta.max=eta.max,eta.obs=eta.obs)#list(llik=c(llik=llik,eta.obs=eta.obs,eta.mean=eta.mean,eta.var=eta.var,eta.min=eta.min,eta.max=eta.max,runtime=runtime[[3]]),hist=list(bins=eta.bins,lims=bin.lims))
}
