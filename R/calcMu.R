calcMu <- function(c_index, codon.parms){
  n.codons=61
  #1) Calculate codon counts 
  codon.counts = unlist(lapply(1:n.codons,FUN=function(i) length(which(c_index==(i-1)))))
  #2) calc mu
  mut.rates = as.numeric(codon.parms$mut_rate[sort.list(codon.parms$c_index)])
  mut.rates <- exp(mut.rates)
  prod(mut.rates[1:n.codons]^codon.counts)
  #sum(mut.rates[1:n.codons]*codon.counts)
  
}
