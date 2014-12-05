opt.seq <- function(gene.list,codon.parms){
  aa.seq <- gene.list$gene.dat$aa
  
  if(is.na(aa.seq[length(aa.seq)])){
    aa.seq <- aa.seq[-length(aa.seq)]
  }
  
  #### This basically looks for the codon index with the best elongation prs 
  c.seq.opt <- unlist(lapply(aa.seq,FUN = function(aa) (codon.parms$c_index[which(codon.parms$aa==aa)])[which(codon.parms$elong_pr[which(codon.parms$aa==aa)] == max(codon.parms$elong_pr[which(codon.parms$aa==aa)]))[1]]))
  
  #return
  c.seq.opt
}

pess.seq <- function(gene.list,codon.parms){
  aa.seq <- gene.list$gene.dat$aa
  
  if(is.na(aa.seq[length(aa.seq)])){
    aa.seq <- aa.seq[-length(aa.seq)]
  }
  
  #### This basically looks for the codon index with the worst elongation prs 
  c.seq.pess <- unlist(lapply(aa.seq,FUN = function(aa) (codon.parms$c_index[which(codon.parms$aa==aa)])[which(codon.parms$elong_pr[which(codon.parms$aa==aa)] == min(codon.parms$elong_pr[which(codon.parms$aa==aa)]))[1]]))
  
  #return
  c.seq.pess
}

rand.seq <- function(gene.list,codon.parms){
  aa.seq <- gene.list$gene.dat$aa
  
  if(is.na(aa.seq[length(aa.seq)])){
    aa.seq <- aa.seq[-length(aa.seq)]
  }
  
  rand.seq <- unlist(lapply(aa.seq,FUN = function(aa){ candidates <- which(codon.parms$aa==aa)
                                                       (codon.parms$c_index[candidates])[ceiling(runif(n=1)*length(candidates))]
                                                       }))
  
  #return
  rand.seq
}