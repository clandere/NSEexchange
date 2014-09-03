initialize.codon.parms <- function(codon.parms, tracelength)
{
  if(is.null(codon.parms$init))
  {
    if(!is.null(codon.parms$obs))
    {
      codon.parms$init <- codon.parms$obs
      warning('In initialize.codon.parms: codon.parms$init not included. Defaulting codon.parms$init = codon.parms$obs')
    }else{
      stop('In initialize.codon.parms: neither codon.parms$init or codon.parms$obs were found. Exiting...')
    }
  }
  if(is.null(codon.parms$curr))
  {
    codon.parms$curr <- codon.parms$init
  }
  if(is.null(codon.parms$trace))
  {
    ## initialize list with proper length and set first element of trace
    codon.parms$trace <- vector("list", tracelength)
    codon.parms$trace[[1]] <- codon.parms$curr
  }
  
  return(codon.parms)
}