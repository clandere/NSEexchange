update.codon.parms.one.aa <- function(genome, aa, curr.codon.parms, genome.parms, MCMC, parallel='mclapply', n.cores=2){
  #propose new parameters for this aa
  P.codon.parms <- propose.codon.parms.one.aa(curr.codon.parms, MCMC$codon[[aa]], aa)
  
#   cat(paste('mut', round(as.numeric(curr.codon.parms$mut_rate[6]), 4), 
#             'p.mut', round(as.numeric(P.codon.parms$mut_rate[6]), 4),
#             'elpr', round(as.numeric(curr.codon.parms$elong_pr[6]), 7),
#             'p.elpr',round(as.numeric(P.codon.parms$elong_pr[6]), 7), ' '))
  #calculate acceptance ratio for this set of parameters
  ret <- aa.accept.ratio.all.genes(genome, P.codon.parms, curr.codon.parms, genome.parms, MCMC, parallel, n.cores)
  accept.ratio <- ret$accept.ratio
  cat(paste('acceptance ratio:', round(accept.ratio,4), '\n'))
  
  #accept/reject
  if(!is.nan(accept.ratio) && (runif(1) < accept.ratio))
  {
    curr.codon.parms <- P.codon.parms
    cat("\t proposed parameters accepted\n")
  }
  
  ret=list(curr.codon.parms=curr.codon.parms, aux_c_index=ret$aux_c_index)
}

propose.codon.parms.one.aa <- function(curr.codon.parms, MCMC.aa, aa, mut.flag=TRUE){
  #Values 1:(n.codons-1) correspond to mutation rates
  #Values (n.codons):(2*n.codons-2) correspond to NSE pr
  P.values = exp(mvrnorm(1, MCMC.aa$mu, MCMC.aa$sigma))
  
  ind <- which(curr.codon.parms$aa == aa)
  n.codons <- length(ind)
  
  cat("not updating mutation rate! ")
  #Skip first codon mut rate for each aa
  #curr.codon.parms$mut_rate[ind[-1]] <- P.values[1:(n.codons-1)]
  
  #Don't skip first codon NSE pr
  if(nse.id=='n-1')
  {
    #curr.codon.parms$nse_pr[ind] <- P.values[n.codons:(2*n.codons-1)] #NSE id = n-1
    curr.codon.parms$nse_pr[ind[-1]] <- P.values[n.codons:(2*n.codons-2)] #NSE id = n-1
    curr.codon.parms$elong_pr <- 1-curr.codon.parms$nse_pr
  }else if(nse.id=='n')
  {
    curr.codon.parms$elong_pr[ind] <- 1-P.values[n.codons:(2*n.codons-1)] #NSE id = n
  }
  
  return(curr.codon.parms)
}
