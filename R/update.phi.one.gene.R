update.phi.one.gene <- function(i,genome,codon.parms,pop.parms,MCMC){
    
  #1) Propose new phi 
  # Use random walk with lognormal distribution or reflecting normal. 
  # Note: reflecting normal seems to work better because when using a lognormal
  # distbn, the phi traces can get stuck in low-phi "wells" 
  if(MCMC$phi$prop.type == 'LN'){ #Lognormal Distribution has symmetric q'/q
    P.log.phi = rnorm(1,mean=log(MCMC$phi$curr.phi[i]),sd=MCMC$phi$prop.var[i])
    P.phi = exp(P.log.phi)
    mh.ratio = 1 #metropolis-hastings ratio -- symmetric jumping distribution
  }else if(MCMC$phi$prop.type == 'RN'){ #Reflecting Normal also has symmetric q'/q
    P.phi = rnorm(1,mean=MCMC$phi$curr.phi[i],sd=MCMC$phi$prop.var[i])
    if(P.phi<0) P.phi=-P.phi #reflecting normal distribution
    mh.ratio = 1 #metropolis-hastings ratio -- symmetric jumping distribution
  }else{
    error(paste('Error in SempprExchange.main.R: Invalid phi proposal type',proposal.type))
  }
  
  prior.ratio = 1 #assume flat prior - this will be replaced by s and k 
  #prior.ratio = phi.prior(P.phi,pop.parms$s,pop.parms$k)/phi.prior(MCMC$phi$curr.phi[i],pop.parms$s,pop.parms$k)
  
  #2) Generate auxiliary variable (simulate codon sequence)
  #      Russ thinks we need to start at a random sequence to make
  #      sure our samples are exact and independent like the 
  #      exchange algorithm requires, but it takes a very long time to burn in a random sequence
  #      for a high expression gene, so it isn't feasable to do so. Starting from the previous sequence
  #      assures that the aux variable is burned in and that we are pulling from the steady-state
  #      distribution. It also save a shitload of computational time. I compared the two approaches 
  #      and verified they produce the same results.
  
  #gene.list$aux_c_index <- rand.seq(gene.list,codon.parms) 
  #FIXME: I want to start from the previously accepted aux sequence, but I need 
  #       a way to keep track of it outside of this function. Maybe I need to 
  #       make a global structure and put it inside that.
  aux_c_index <- simulate.sequence(phi=P.phi,
                                   c_index=genome[[i]]$gene.dat$c_index,
                                   codon.parms=codon.parms,
                                   BIS=MCMC$aux$BIS,
                                   GMT=MCMC$aux$GMT,
                                   MES=MCMC$aux$MES,
                                   SIMULATION.METHOD=MCMC$aux$simulation.type)
  
  #3) Calculate acceptance rate
  #delta eta is defined as eta_obs - eta_aux
  eta_obs = calc.eta.NSE(codon_index=genome[[i]]$gene.dat$c_index,codon.parms,pop.parms)
  eta_aux = calc.eta.NSE(aux_c_index,codon.parms,pop.parms)
  delta.eta = eta_obs - eta_aux
  
  #NOTE: Mu cancels out of this acceptance ratio 
  accept.ratio = exp(pop.parms$Q*pop.parms$Ne*delta.eta*(MCMC$phi$curr.phi[i] - P.phi))*
    prior.ratio*
    mh.ratio
  
  #4) Draw random number and accept/reject
  #debug output
  #cat(paste('phi: ',gene.list$phi.value,' prop.phi: ',P.phi,' eta.obs: ',eta_obs,'eta.aux: ',eta_aux,' accept.rat: ',accept.ratio,'\n'))
  if(runif(1) < accept.ratio){
    step.phi = P.phi
  }else{
    step.phi = MCMC$phi$curr.phi[i]
  }

  ret=list(step.phi=step.phi,aux_c_index=aux_c_index)
}