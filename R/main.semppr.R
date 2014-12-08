main.semppr <- function(obs.genome,
                        codon.parms,
                        phi=NULL,
                        genome.parms=NULL,
                        aa.list = c('C'),
                        out.prefix=NULL,
                        num.MCMC.steps=3000,
                        BIS=4,GMT=0,MES=0,
                        phi.proposal.type = 'RN',
                        codon.proposal.type = 'LN',
                        aux.simulation.type = 'M',
                        parallel='lapply',
                        nseid='n-1',
                        n.cores = 2,
                        simulate.dataset=TRUE){
# Purpose: Fit an evolutionary model of nonsense errors described in Gilchrist et al. 2009 to genome data using an evolutionary
#          model described in Sella and Hirsh 2005 in order to estimate modl parameters. The enormous genotype space in this model 
#          renders the posterior distribution of model parameters 
#          doubly-intractible. Inference requires the MCMC Exchange method described in Murray et al. 2012 to circumvent calculation of parameter-dependent 
#          normalization constants in the parameter acceptance ratios. This method requires the simultaneous proposal of an auxiliary variable along with
#          model parameter proposal. The auxiliary variable appears in this model as a simulated gene or simulated genome, depending on whether the 
#          proposed parameter is gene-specific, such as expression rate, or genome-wide, such as codon mutation rate.
# Inputs: obs.genome - list of lists containing gene-specific information
#                                        gene.list - list containting 
#                                             $phi.value - expression rate
#                                             $name - gene ID
#                                             $gene.dat, which contains
#                                                 $codon - nucleotide char sequence
#                                                 $aa - amino acid char sequence
#                                                 $c_index - representation of codon sequence using integers
#                                                 $mut_rate - vector of mutation rates corresponding to c_index
#                                                 $elong_rat - vector of elong rates...
#                                                 $count - amino acid/codon count
#         genome.parms - list containing genome-wide parameters
#                                         $Q - Arbitrary scaling factor, default = 1
#                                         $Ne - Effective population size, default = 1
#                                         $A1 - Cost of translation initiation (in ATP's), default = 4
#                                         $A2 - Cost of peptide elongation (in ATP's), default = 4
#                                         $B - background Nonsense Error Rate (in s^-1), default = 0.0025 
#         codon.parms - dataframe with rows corresponding to codon data, including: 
#                                         [,1] - AA - Character corresponding to translated amino acid 
#                                         [,2] - codon - string containing NT sequence of codon 
#                                         [,3] - c_index - integer denoting codon index 
#         out.prefix - desired folder to write output files
#         min.num.steps - minimum number of MCMC steps for each gene
#         num.MCMC.steps - maximum number of MCMC steps for each gene
#         BIS - For Aux variable generation: Number of simulation steps per amino acid
#         GMT - For Aux variable generation: Amount of simulation time (only use with Evolution simulation method)
#         MES - For Aux variable generation: Number of simulation steps (you can use this or BIS or both)
#         proposal.type - Phi proposal type, options:
#                                      'RN' - reflecting normal proposal centered around current phi value. adapts to 
#                                             acceptance ratio >0.25 and <0.35
#                                      'LN' - Lognormal proposal centered around current phi value. adapts to accept
#                                             ratio >0.25 and <0.35
#         simulation.type - For Aux variable generation: simulation method, options:
#                                      'M' - MCMC
#                                      'E' - Evolution simulation
#         n.cores - number of cores used for mclapply
#         simulate.dataset - Flag to simulate datset 
#         phi.only - Flag to estimate phi only
# Outputs: Files located in out.prefix
# Usage: See 'R/test_script1.R' for a usage example

  #For Debugging 
  #out.prefix='output/';genome.parms=NULL;min.num.steps=500;num.MCMC.steps=3000;BIS=4;GMT=0;MES=0;phi.proposal.type = 'RN';aux.simulation.type = 'M';  n.cores = 2; simulate.dataset=TRUE;init.phi.dat=NULL;parallel='lapply'
  #setwd('~/WorkingCopies/ces3/branches/exchange_elpr/R/')
  #load('../data/subset_100genes.Rda')
  #load('../data/trna.CETS.flatmut.Rda')
  #obs.genome=list(glf.test.100[[1]],glf.test.100[[2]])
  #init.codon.parms=codon.parms
  #aa.list=c('C')
  
  #load required libraries and files
  #require(MASS)
  #source('../R/gen.seq.R')
  #source('../R/calc.eta.NSE.R')
  #source('../R/calcMu.R')
  #source('../R/gelman.conv.R')
  #source('../R/calcScuo.R')
  #source('../R/guess.phi.scuo.R')
  #source('../R/aa.accept.ratio.all.genes.R')
  #source('../R/aa.accept.ratio.one.gene.R')
  #for(file in c(dir(path='../R/',pattern='update*'),
  #              dir(path='../R/',pattern='initialize*'),
  #              dir(path='../R/',pattern='simulate*'),
  #              dir(path='../R/',pattern='construct*'),
  #              dir(path='../R/',pattern='initialize*'))) source(paste('../R/',file,sep=''))

  #dyn.load('../src/exchange_NSE.so')
    
  #begin timing 
  start.time = proc.time() 
  
  #### Initialization ####
  if(!is.null(aa.list)&&(aa.list=='All'||aa.list=='all')) {
    aa.list = c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 
                'K', 'L', 'N', 'P', 'Q', 'R', 'S', 'T', 
                'V', 'Y', 'Z')
  }else if(!is.null(aa.list)&&(aa.list=='None'||aa.list=='none')){
    aa.list=NULL
  }
  assign('nse.id', nseid, envir = ,.GlobalEnv)
  out.prefix = initialize.output(out.prefix)
  codon.parms = initialize(codon.parms, num.MCMC.steps)
  genome.parms = initialize.genome.parms(genome.parms)
  if(simulate.dataset){
    genome = simulate.data.all.genes(obs.genome,
                                     phi$obs,
                                     codon.parms$obs,
                                     codon.parms$init,
                                     genome.parms,
                                     MES,GMT,BIS,
                                     aux.simulation.type,
                                     n.cores,parallel)
  }
  phi <- initialize.phi(phi, genome, codon.parms, num.MCMC.steps, parallel)
  MCMC <- initialize.MCMC.info(phi$trace.dat, phi.proposal.type, aa.list, codon.parms$init, codon.proposal.type, BIS, MES, GMT, aux.simulation.type)
  genome <- initialize.genome(genome, phi$trace.dat, codon.parms$curr, genome.parms, MCMC, n.cores, parallel)
    
  #save initial codon.parms$curr and initialize phi vec
  
  step.phi = rep(0, length(genome))
  
  after.init.time = proc.time()
  
  #### End Initialization ####
  
  #### Begin MCMC ####
  i <- 0
  while(i < num.MCMC.steps){
    i <- i + 1
    #Update phi for each gene 
    ret <- update.phi.all.genes(genome, codon.parms$curr, genome.parms, MCMC, parallel=parallel, n.cores=n.cores)
    
    #Save phi values for this step **this could be costly
    for(j in 1:length(genome)){
      genome[[j]]$phi.curr = ret[[j]]$step.phi
      genome[[j]]$aux_c_index=ret[[j]]$aux_c_index
      step.phi[j] = ret[[j]]$step.phi
    }

    phi$trace.dat = rbind(phi$trace.dat,step.phi)
    #phi$trace.dat[i+1,] = step.phi
    
    #Update trna parms for each gene
    ret = update.codon.parms.all.aa(genome, phi$trace.dat[i+1,], codon.parms$curr, genome.parms, MCMC, aa.list, n.cores=n.cores, parallel=parallel)

    #Save codon parms and for this step **This could be costly
    for(j in 1:length(genome)) genome[[j]]$aux_c_index=ret$aux_c_index[[j]]
    
    codon.parms$trace[[i+1]] <- codon.parms$curr <- ret$curr.codon.parms
    
    #Update MCMC info
    MCMC = update.MCMC.info(i, MCMC, phi$trace.dat, codon.parms$trace, aa.list)
  }
  
  #### End MCMC ####
  
  #Stop clock time
  finish.time = proc.time()
  time.list = list(start.time=start.time,after.init.time=after.init.time,finish.time=finish.time)
  
  
  #Save files or return stuff
  Date = strsplit(as.character(Sys.time()),' ')[[1]][1]
  if(!is.null(out.prefix)){
    save(time.list,file=paste(out.prefix,Date,'_time.Rda',sep=''))
    save(codon.parms,file=paste(out.prefix,Date,'_codon.Rda',sep=''))
    save(phi,file=paste(out.prefix,Date,'_phi.Rda',sep=''))
    save(obs.genome,file=paste(out.prefix,Date,'_genome.Rda',sep=''))
  }else{
    return(list(obs.genome=obs.genome,
                phi=phi,
                codon.parms=codon.parms,
                time.list=time.list, mcmc=MCMC))
  }
}
