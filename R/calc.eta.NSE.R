calc.eta.NSE <- function(codon_index,codon.parms,pop.parms = list(B = 0.0025,A1 = 4,A2 = 4)){
  #Purpose: Wrapper function for C code to calculate the eta value for a codon index based on 
  #         NSE model
  #Inputs: codon_index
  #        codon.parms - codon.parms structure. Can be created with SE.semppr.data.load
  #        pop.parms - parameters relevant to entire population of genes, list containing:
  #                      $B - background NSE rate
  #                      $A1 - Translation initiation
  #                      $A2 - Translation elongation
  #Outputs: eta - Eta value for codon index
  B = pop.parms$B
  A.1 = pop.parms$A1
  A.2 = pop.parms$A2

  if(codon_index[length(codon_index)]==99||codon_index[length(codon_index)]==61||codon_index[length(codon_index)]==62||codon_index[length(codon_index)]==63){
    codon_index <- codon_index[-length(codon_index)] #remove the last codon(stop codon) -- stop codon is not used in c
  }
  aa_count = length(codon_index)
  
  elong_pr = as.numeric(codon.parms$elong_pr[sort.list(codon.parms$c_index)])
  
  n_codons = length(elong_pr)
  
  eta = 0
  
  OUTPUT <- .C("calc_eta_NSE_R_wrapper",as.integer(codon_index),as.integer(aa_count),as.double(elong_pr),as.integer(n_codons),as.double(A.1),as.double(A.2),eta = as.double(eta))
  
  OUTPUT$eta
  
  
}


