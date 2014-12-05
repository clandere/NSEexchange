update.MCMC.info <- function(i, MCMC, phi.dat, codon.parms.trace, aa.list){
  
  curr.phi.dat = phi.dat[i+1,]
  prev.phi.dat = phi.dat[i,]
  curr.codon.parms = codon.parms.trace[[i+1]]
  prev.codon.parms = codon.parms.trace[[i]]
  
  #Update Phi info
  MCMC$phi$accept.history[i%%50,] = (curr.phi.dat!=prev.phi.dat) #current row is i+1 since we started with scuo guess
  MCMC$phi$curr.phi = curr.phi.dat
  
  #Adjust the phi proposal variance every 50 steps
  if(i%%50==0){
    accept.rate = sum(MCMC$phi$accept.history)/50
    
    #Gelman says optimal accept rate for one parameter
    #proposal is 0.44
    j=which(accept.rate > 0.49) #if accept rate is too big, make jumps bigger
    MCMC$phi$prop.var[j] = MCMC$phi$prop.var[j]*1.1
    
    j=which(accept.rate < 0.39) #if accept rate is too small, make jumps smaller
    MCMC$phi$prop.var[j] = MCMC$phi$prop.var[j]/1.1
  }
  
  #Update Codon info
  if(!is.null(aa.list)){ 
    for(aa in aa.list){
      ind = which(curr.codon.parms$aa==aa)
      if(nse.id=='n-1'){
        MCMC$codon[[aa]]$mu = 
          c(log(as.numeric(curr.codon.parms$mut_rate[ind[-1]])),
            log(1-as.numeric(curr.codon.parms$elong_pr[ind[-1]]))) #NSE id = n-1
      }else if(nse.id =='n'){
        MCMC$codon[[aa]]$mu = 
          c(log(as.numeric(curr.codon.parms$mut_rate[ind[-1]])),
            log(1-as.numeric(curr.codon.parms$elong_pr[ind]))) #NSE id = n
      }
      MCMC$codon[[aa]]$accept.history[i%%50] = 
        any(curr.codon.parms$elong_pr[ind]!=
               prev.codon.parms$elong_pr[ind])||
            any(curr.codon.parms$mut_rate[ind]!=
               prev.codon.parms$mut_rate[ind])
    }
    
    #Scale the covariance matrix
    #if(i%%50==0){
    #  for(aa in aa.list){
    #    accept.rate = sum(MCMC$codon[[aa]]$accept.history)/length(MCMC$codon[[aa]]$accept.history)
    #    cat(paste('AcceptRate:',accept.rate,MCMC$codon[[aa]]$sigma[1,1],'\n'))
    #    #Gelman says optimal acceptance rate for multi parameter
    #    #proposal is 0.23
    #    if(accept.rate > 0.28){
    #      MCMC$codon[[aa]]$sigma = MCMC$codon[[aa]]$sigma*1.2
    #    }else if(accept.rate < 0.18){
    #      MCMC$codon[[aa]]$sigma = MCMC$codon[[aa]]$sigma*0.8
    #    }
    #  }
    #}
    
    #Update covariance matrix
    if(i%%100==0){
      MCMC=update.covariance.matrix(codon.parms.trace,MCMC,aa.list,100,10)
    }
  }
  return(MCMC)
}
