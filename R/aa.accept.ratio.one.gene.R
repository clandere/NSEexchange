aa.accept.ratio.one.gene <- function(gene.list, P.codon.parms, curr.codon.parms, genome.parms, MCMC){
  
  #Generate aux_c_index
  aux_c_index = simulate.sequence(phi=gene.list$phi.curr,
                                  c_index=gene.list$aux_c_index,
                                  codon.parms=P.codon.parms,
                                  BIS=MCMC$aux$BIS,
                                  GMT=MCMC$aux$GMT,
                                  MES=MCMC$aux$MES,
                                  SIMULATION.METHOD=MCMC$aux$simulation.type)
  
  #cat("aa.accept.ratio.one.gene -> aux_c_index: ");print(aux_c_index)
#   cat("aa.accept.ratio.one.gene -> curr.codon.parms: ");print(curr.codon.parms)
#   cat("aa.accept.ratio.one.gene -> P.codon.parms: ");print(P.codon.parms)
  #Calculate acceptance ratio  
  eta_obs <- calc.eta.NSE(gene.list$gene.dat$c_index, curr.codon.parms, genome.parms)
  eta_obs_p <- calc.eta.NSE(gene.list$gene.dat$c_index, P.codon.parms, genome.parms)
  eta_sim <- calc.eta.NSE(aux_c_index, curr.codon.parms, genome.parms)
  eta_sim_p <- calc.eta.NSE(aux_c_index, P.codon.parms, genome.parms)
  
  mu_obs <- calcMu(gene.list$gene.dat$c_index, curr.codon.parms)
  mu_obs_p <- calcMu(gene.list$gene.dat$c_index, P.codon.parms)
  mu_sim <- calcMu(aux_c_index, curr.codon.parms)
  mu_sim_p <- calcMu(aux_c_index, P.codon.parms)
  
  acceptance.ratio <- ( (mu_obs_p * mu_sim)/(mu_obs * mu_sim_p) ) * exp((-genome.parms$Q * genome.parms$Ne * gene.list$phi.curr) * (eta_obs_p + eta_sim - eta_obs - eta_sim_p) )
#   cat("aa.accept.ratio.one.gene -> accept ratio: ");print(acceptance.ratio)
#   
#   cat("aa.accept.ratio.one.gene -> eta_obs: ");print(eta_obs)
#   cat("aa.accept.ratio.one.gene -> eta_obs_p: ");print(eta_obs_p)
#   cat("aa.accept.ratio.one.gene -> eta_sim: ");print(eta_sim)
#   cat("aa.accept.ratio.one.gene -> eta_sim_p: ");print(eta_sim_p)
# 
#   cat("aa.accept.ratio.one.gene -> mu_obs: ");print(mu_obs)
#   cat("aa.accept.ratio.one.gene -> mu_obs_p: ");print(mu_obs_p)
#   cat("aa.accept.ratio.one.gene -> mu_sim: ");print(mu_sim)
#   cat("aa.accept.ratio.one.gene -> mu_sim_p: ");print(mu_sim_p)
  
  
  ret <- list(aux_c_index=aux_c_index, acceptance.ratio=acceptance.ratio)
}