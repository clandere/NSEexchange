initialize.covariance.mat <- function(codon.parm.subset){
  n.codons = length(codon.parm.subset$codon)
  #cat("initialize.covariance.mat -> n.codons");print(n.codons)
  if(n.codons==1){
    return(NULL)
  }else{
    if(nse.id=='n-1'){
      mu = c(as.numeric(codon.parm.subset$mut_rate[2:n.codons]),
             as.numeric(codon.parm.subset$nse_pr[2:n.codons])) #NSE id = n-1
      sigma = diag(2*n.codons-2)/100
    }else if(nse.id=='n'){
      mu = c(as.numeric(codon.parm.subset$mut_rate[2:n.codons]),
             as.numeric(codon.parm.subset$nse_pr[1:n.codons])) #NSE id = n
      sigma = diag(2*n.codons-1)/100
    }
    #cat("initialize.covariance.mat -> Mu");print(mu)
    #cat("initialize.covariance.mat -> Sigma");print(sigma)
    return(list(mu=mu,
                sigma=sigma))
  }
}
