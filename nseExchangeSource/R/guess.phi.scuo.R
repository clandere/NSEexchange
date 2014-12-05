guess.phi.scuo <- function(genome,codon.parms,parallel){
  phi.dat = rlnorm(length(genome),0,1.5)
  
  #set max phi = 400. If phi gets too high, 
  #numerical precision becomes a problem because fitness=exp(-qNe*phi*eta)
  
  phi.dat[which(phi.dat>400)]=400
  
  scuo = calcScuo(genome,codon.parms,2,parallel)
  scuo.order = order(scuo,decreasing=TRUE)
  
  phi.dat=phi.dat[order(phi.dat,decreasing=TRUE)]
  phi.dat[scuo.order] = phi.dat
  
  phi.dat
}