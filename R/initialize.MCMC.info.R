initialize.MCMC.info <- function(phi.trace.dat,
                                 phi.proposal.type,
                                 aa.list,
                                 init.codon.parms,
                                 codon.proposal.type,
                                 BIS,MES,GMT,
                                 aux.simulation.type){
  #set up mcmc info
  MCMC = list(aux=list(),phi=list(),trna=list())
  
  #Auxiliary variable info
  MCMC$aux = list(BIS=BIS,MES=MES,GMT=GMT,simulation.type=aux.simulation.type)
  
  #Phi info
  MCMC$phi = list(curr.phi = phi.trace.dat[nrow(phi.trace.dat),],
                  prop.var = phi.trace.dat[nrow(phi.trace.dat),],
                  prop.type = phi.proposal.type,
                  accept.history = matrix(0, 50, ncol(phi.trace.dat)) #this records the last 50 accept/rejections\
  )
#   cat("initialize.MCMC.info -> phi cur: ");print(MCMC$phi$curr.phi)
#   cat("initialize.MCMC.info -> phi var: ");print(MCMC$phi$prop.var)
  
  #Codon parameter info  
  MCMC$codon = list()
  
  if(!is.null(aa.list)){
    levels.codon.parms = levels(as.factor(init.codon.parms$aa))
    exclude.levels=levels.codon.parms[!(levels.codon.parms %in% aa.list)]
    #MCMC$codon$mu = by(codon.parms,factor(codon.parms$aa,exclude=exclude.levels),)
    MCMC$codon = by(init.codon.parms,factor(init.codon.parms$aa,exclude=exclude.levels),initialize.covariance.mat)
    
    
    for(aa in aa.list){
      MCMC$codon[[aa]]$accept.history = rep(0,50)
    }
  }
  
  
  MCMC
}
