update.phi.all.genes <- function(genome,codon.parms,pop.parms,MCMC,n.cores,parallel='mclapply'){
  
  if(parallel=='mclapply'){
    if(length(genome)/n.cores < 5){
      mc.preschedule=FALSE
    }else{
      mc.preschedule=TRUE
    }
    ret = mclapply(X=1:length(genome),
                   FUN = update.phi.one.gene,
                   genome,
                   codon.parms,
                   pop.parms,
                   MCMC,
                   mc.cores=n.cores,
                   mc.preschedule=mc.preschedule)
    #Note: mc.preschedule = TRUE may work better when the num.gene/n.cores ratio is large
    #                       
    #      mc.preschedule = FALSE may work better when there are only ~5 genes per core
  }else if(parallel=='lapply'){
    ret = lapply(X=1:length(genome),
                 FUN = update.phi.one.gene,
                 genome,
                 codon.parms,
                 pop.parms,
                 MCMC)
  }else{
    stop(paste('Function update.phi.all.genes does not support parallelization of type ','\"',parallel,'\"',sep=''))
  }
  
  ret
}
