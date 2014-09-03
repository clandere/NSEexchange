simulate.data.all.genes <- function(obs.genome,
                                    obs.phi,
                                    obs.codon.parms,
                                    init.codon.parms,
                                    obs.genome.parms,
                                    MES,GMT,BIS,
                                    aux.simulation.type,
                                    n.cores,parallel='mclapply'){
  
  #See if real phi value is included
  flag1 <- any(unlist(lapply(obs.genome,function(gene)(is.null(gene$phi.value) || is.na(gene$phi.value)))))
  flag2 <- is.null(obs.phi)
  if(flag1 && flag2){
    stop('Error: In initialize.obs.genome: Must pass observed phi values in argument obs.phi or
         contain observed phi values in obs.genome[[i]]$phi.value when simulating dataset. ')
  } 
  
  if(is.null(obs.codon.parms)){warning('obs.codon.parms for simulating dataset are not included. 
                                         Using init.codon.parms.');obs.codon.parms<-init.codon.parms}

  if(parallel=='mclapply'){
    obs.genome = mclapply(X=obs.genome, 
                      FUN=simulate.data.one.gene,
                      obs.codon.parms,
                      obs.genome.parms,
                      MES,GMT,BIS,
                      aux.simulation.type,
                      mc.cores=n.cores,
                      mc.preschedule=TRUE)
    
  }else if(parallel=='lapply'){
    obs.genome = lapply(X=obs.genome, 
                    FUN=simulate.data.one.gene,
                    codon.parms=obs.codon.parms,
                    genome.parms=obs.genome.parms,
                    MES,GMT,BIS,
                    aux.simulation.type)
    
  }else{
    stop(paste('Function initialize.obs.genome does not support parallelization of type ','\"',parallel,'\"',sep=''))
  }
  
}