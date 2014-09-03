initialize.covariance.mat <- function(codon.parm.subset){
  n.codons = length(codon.parm.subset$codon)
  if(n.codons==1){
    return(NULL)
  }else{
    if(nse.id=='n-1'){
      mu = c(log(as.numeric(codon.parm.subset$mut_rate[2:n.codons])),
             log(1-as.numeric(codon.parm.subset$elong_pr[2:n.codons])))#NSE id = n-1
      sigma = diag(2*n.codons-2)/100
    }else if(nse.id=='n'){
      mu = c(log(as.numeric(codon.parm.subset$mut_rate[2:n.codons])),
             log(1-as.numeric(codon.parm.subset$elong_pr[1:n.codons])))#NSE id = n
      sigma = diag(2*n.codons-1)/100
    }
    return(list(mu=mu,
                sigma=sigma))
  }
}
