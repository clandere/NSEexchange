initialize.genome <- function(genome,
                              phi.trace.dat,
                              codon.parms,
                              genome.parms,
                              MCMC,
                              n.cores,parallel='mclapply'){
  
  for(i in 1:length(genome)) genome[[i]]$phi.curr = phi.trace.dat[nrow(phi.trace.dat),i]
  #Burn in aux sequence to initial lambda and phi #
  if(parallel=='mclapply'){
    genome = mclapply(genome, 
                      initialize.gene,
                      codon.parms,
                      genome.parms,
                      MCMC,
                      mc.cores=n.cores,
                      mc.preschedule=TRUE)
    
  }else if(parallel=='lapply'){
    genome = lapply(genome, 
                    initialize.gene,
                    codon.parms,
                    genome.parms,
                    MCMC)
    
  }else{
    stop(paste('Function initialize.genome does not support parallelization of type ','\"',parallel,'\"',sep=''))
  }
  
  genome
}