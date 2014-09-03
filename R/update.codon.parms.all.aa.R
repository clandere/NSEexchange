update.codon.parms.all.aa <- function(genome,phi.dat,curr.codon.parms,genome.parms,MCMC,aa.list,n.cores,parallel='mclapply'){
  
  if(!is.null(aa.list)){
    for(i in order(runif(length(aa.list)))){ #update aa in random order
      ret = update.codon.parms.one.aa(genome,aa.list[i],curr.codon.parms,genome.parms,MCMC,parallel,n.cores)
    }
    return(list(curr.codon.parms=ret$curr.codon.parms,aux_c_index=ret$aux_c_index))
  }else{
    return(list(curr.codon.parms=curr.codon.parms))
  }
  
  
}