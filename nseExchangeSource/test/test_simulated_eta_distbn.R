library(NSEexchange,quietly=TRUE)
data(subset_1gene)
data(S.cer.codon.parms)
load('data/benchmarked.CI.Rda')

test.gene=glf.test.1[[1]]
trna.dat=codon.parms$obs


num_steps=100
phi=1

#get eta distribution for evolution simulation
cat("\nGetting eta distribution for evolution simulation...    ")
t_c_index <- simulate.sequence(phi,test.gene$gene.dat$c_index,trna.dat,40,-10,SIMULATION.METHOD='E')
eta_evol = calc.eta.NSE(t_c_index,trna.dat)
for(i in 1:(num_steps-1)){
  t_c_index = simulate.sequence(phi,t_c_index,trna.dat,1,-1,SIMULATION.METHOD='E')
  eta_evol = c(eta_evol,calc.eta.NSE(t_c_index,trna.dat))
}
cat("Done\n")

#get eta distribution for mcmc simulation
cat("Getting eta distribution for mcmc simulation...         ")
t_c_index = simulate.sequence(phi,test.gene$gene.dat$c_index,trna.dat,40,0,SIMULATION.METHOD='M')
eta_mcmc = calc.eta.NSE(t_c_index,trna.dat)
for(i in 1:(num_steps-1)){
  t_c_index = simulate.sequence(phi,t_c_index,trna.dat,1,0,SIMULATION.METHOD='M')
  eta_mcmc = c(eta_mcmc,calc.eta.NSE(t_c_index,trna.dat))
}
cat("Done\n")

#test to see if confidence intervals overlap
cat("\nTesting to see if CI for sample means overlap with each other...        ")
ci_eta_evol = c(mean(eta_evol)-1.96*sqrt(var(eta_evol)/num_steps),mean(eta_evol)+1.96*var(eta_evol)/sqrt(num_steps))
ci_eta_mcmc = c(mean(eta_mcmc)-1.96*sqrt(var(eta_mcmc)/num_steps),mean(eta_mcmc)+1.96*var(eta_mcmc)/sqrt(num_steps))

if(ci_eta_evol[2]<ci_eta_mcmc[1]||ci_eta_evol[1]>ci_eta_mcmc[2]){cat("Error\n")}else{cat("Good\n")}

#test to see if distribution look like our benchmark distribution
cat("\nTesting to see if CI for sample means overlap with benchmarked CI...\n")
cat("Evol CI...     ")
if(ci_benchmark[2]<ci_eta_evol[1]||ci_benchmark[1]>ci_eta_evol[2]){cat("Error\n")}else{cat("Good\n")}
cat("MCMC CI...     ")
if(ci_benchmark[2]<ci_eta_mcmc[1]||ci_benchmark[1]>ci_eta_mcmc[2]){cat("Error\n")}else{cat("Good\n")}
