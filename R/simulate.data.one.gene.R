simulate.data.one.gene <- function(gene.list,
                                   codon.parms,
                                   genome.parms,
                                   MES,GMT,BIS,
                                   simulation.type){
  
  cat(paste("\tSimulating", gene.list$name, ". . .\n"))
  
  gene.list$gene.dat$c_index[-length(gene.list$gene.dat$c_index)]=
    simulate.sequence(phi=gene.list$phi.value,
                      c_index=gene.list$gene.dat$c_index,
                      codon.parms=codon.parms,
                      BIS=5*BIS,
                      GMT=GMT,
                      MES=MES,
                      SIMULATION.METHOD=simulation.type,
                      A.1 = genome.parms$A.1,
                      A.2 = genome.parms$A.2,
                      Ne = genome.parms$Ne,
                      Q = genome.parms$Q)
  
  ## Update gene.dat$codon based on updated codon indices
  codon.vector <- character(length = length(gene.list$gene.dat$codon))
    #The stop codon has c_index 99, which doesn't correspond to codon.params$codon, so it is copied over directly.
    #This is bad, and I intend to change it
  codon.vector[length(codon.vector)] <- as.character(gene.list$gene.dat$codon[length(gene.list$gene.dat$codon)])
  
  for(ii in 1:(length(gene.list$gene.dat$codon)-1)){
    codon.vector[ii] <- codon.parms$codon[ which(codon.parms$c_index == gene.list$gene.dat$c_index[ii]) ]
  }
  gene.list$gene.dat$codon <- factor(codon.vector)    #$codon needs to be a factor, not a char vector
  gene.list

}
