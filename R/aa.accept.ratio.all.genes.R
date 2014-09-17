aa.accept.ratio.all.genes <- function(genome, P.codon.parms, curr.codon.parms, genome.parms, MCMC, parallel='mclapply', n.cores=2){
  
  if(parallel=='mclapply'){
    if(length(genome)/n.cores < 5){
      mc.preschedule=FALSE
    }else{
      mc.preschedule=TRUE
    }
    ret = mclapply(X=genome,
                   FUN=aa.accept.ratio.one.gene,
                   P.codon.parms,
                   curr.codon.parms,
                   genome.parms,
                   MCMC,
                   mc.cores=n.cores,
                   mc.preschedule = mc.preschedule)
  }else if(parallel=='lapply'){
    ret = lapply(X=genome,
                 FUN=aa.accept.ratio.one.gene,
                 P.codon.parms,
                 curr.codon.parms,
                 genome.parms,
                 MCMC)
  }else{
    stop(paste('Function aa.accept.ratio.all.genes does not support parallelization of type ', '\"', parallel, '\"', sep=''))
  }

  aux_c_index <- lapply(ret, function(list)list$aux_c_index)
  
  #prior ratio
  prior.ratio <- 1 #change this if you want a prior on codon parameters
  #Metropolis-Hastings Ratio
  mh.ratio <- 1 #change this if using an asymmetric jumping distribution
  accept.ratio <- prod(unlist(lapply(ret, function(list)list$acceptance.ratio)))
  
  ret <- list(accept.ratio=prior.ratio*mh.ratio*accept.ratio, aux_c_index=aux_c_index)
}