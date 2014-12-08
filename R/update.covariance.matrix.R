update.covariance.matrix <- function(codon.parms.trace, MCMC, aa.list, window=100,factor=10){

  for(aa in aa.list){
    ind = which(codon.parms.trace[[1]]$aa==aa)
    
    #Make sure the sample has at least one acceptance. Otherwise
    #the covariance matrix will be 0 which is bad. If the sample
    #has no acceptances, the magnitude of the covariance matrix is
    #too large, so divide by factor
    accept.rate = sum(MCMC$codon[[aa]]$accept.history)/length(MCMC$codon[[aa]]$accept.history)
    cat(paste('\n update.covariance.matrix -> Accept Rate',accept.rate,'\n'))
    if(accept.rate > 0){
      #extract traces from trace list
      if(nse.id=='n-1'){
        outp=lapply(codon.parms.trace, function(cpt){ c(as.numeric(cpt$mut_rate[ind[-1]]) #log mut
                                                      ,as.numeric(cpt$nse_pr[ind[-1]])) } ) #NSE id = n-1        
#         outp=lapply(codon.parms.trace, function(cpt){ c(log(as.numeric(cpt$mut_rate[ind[-1]])) #log mut
#                                                       ,log(1-as.numeric(cpt$elong_pr[ind[-1]]))) } ) #NSE id = n-1
      }else if(nse.id=='n'){
        outp=lapply(codon.parms.trace, function(cpt){ c(as.numeric(cpt$mut_rate[ind[-1]]) #log mut
                                                        ,as.numeric(cpt$nse_pr[ind])) } ) #NSE id = n-1   
#         outp=lapply(codon.parms.trace, function(cpt){ c(log(as.numeric(cpt$mut_rate[ind[-1]])) #log mut
#                                                       ,log(1-as.numeric(cpt$elong_pr[ind]))) } ) #NSE id = n
      }

      #convert to matrix
      trace.mat <- do.call('rbind',outp)
      n.row <- nrow(trace.mat)
      
      #base covariance on most recent samples
      trace.mat <- trace.mat[(n.row-window):n.row,] 
      
      #output some shit
      cat(paste('Before covariance calculation:\n'))
      print(MCMC$codon[[aa]]$sigma)

      cat(paste('After covariance calculation:\n'))
      #calculate covariance
      covmat <- cov(trace.mat)
      if(nse.id == 'n') # setting max variance for n *** test!!!!
      {
        if(max(covmat) > 1){
          cat("setting max variance to 1+rand!")
          numel <- sum(covmat > 1)
          covmat[covmat > 1] <- 1 + rnorm(numel)   
        }
      }
      MCMC$codon[[aa]]$sigma <- covmat
      print(MCMC$codon[[aa]]$sigma)
    }else{
      MCMC$codon[[aa]]$sigma <- MCMC$codon[[aa]]$sigma/factor
    }
  }
  
  #return
  MCMC
}