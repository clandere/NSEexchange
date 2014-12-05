initialize.gene <- function(gene.list,codon.parms,genome.parms,MCMC){

  gene.list$aux_c_index = simulate.sequence(phi=gene.list$phi.curr,
                                            c_index=gene.list$gene.dat$c_index,
                                            codon.parms=codon.parms,
                                            BIS=5*MCMC$aux$BIS,
                                            GMT=MCMC$aux$GMT,
                                            MES=MCMC$aux$MES,
                                            SIMULATION.METHOD=MCMC$aux$simulation.type,
                                            A.1 = genome.parms$A1,
                                            A.2 = genome.parms$A2,
                                            Ne = genome.parms$Ne,
                                            Q = genome.parms$Q)

  gene.list$gene.dat$codon=NULL        # These next few lines reduce
  gene.list$gene.dat$aa=NULL           # the size of the genome object
  gene.list$gene.dat$elong_rate=NULL;  
  gene.list$gene.dat$mut_rate=NULL
  gene.list$gene.dat$count=NULL
  gene.list$gene.dat$c_index=as.integer(gene.list$gene.dat$c_index)
  gene.list$aux_c_index = as.integer(gene.list$aux_c_index)
  
  gene.list
}