initialize.phi <- function(phi, genome, codon.parms, parallel='mclapply')
{  
  if(is.null(phi))
  {
    phi<-list()
  }
  
  if(is.null(phi$init))
  {
    phi$init <- guess.phi.scuo(genome, codon.parms, parallel)  
    dim(phi$init) <- c(1,length(phi$init))
  }
  
  if(is.null(phi$trace.dat))
  {
    phi$trace.dat = phi$init
  }
    
  return(phi)
}