initialize.phi <- function(phi, genome, codon.parms, mcmcsteps, parallel='mclapply')
{  
  if(is.null(phi))
  {
    phi <- list()
  }
  
  if(is.null(phi$init))
  {
    phi$init <- guess.phi.scuo(genome, codon.parms, parallel)  
    dim(phi$init) <- c(1,length(phi$init))
    colnames(phi$init) <- unlist(lapply(genome, function(list){list$name}))
  }
  
  if(is.null(phi$trace.dat))
  {
    phi$trace.dat <- phi$init
  }
    
  return(phi)
}