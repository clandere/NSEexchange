simulate.data.one.gene <- function(gene.list,
                                   codon.parms,
                                   genome.parms,
                                   MES,GMT,BIS,
                                   simulation.type){
  
  gene.list$gene.dat$c_index[-length(gene.list$gene.dat$c_index)]=
    simulate.sequence(phi=gene.list$phi.value,
                      c_index=gene.list$gene.dat$c_index,
                      codon.parms=codon.parms,
                      BIS=5*BIS,
                      GMT=GMT,
                      MES=MES,
                      SIMULATION.METHOD=simulation.type,
                      A.1 = genome.parms$A1,
                      A.2 = genome.parms$A2,
                      Ne = genome.parms$Ne,
                      Q = genome.parms$Q)
  gene.list
}