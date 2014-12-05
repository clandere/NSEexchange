library(NSEexchange,quietly=TRUE)
data(subset_1gene)
data(S.cer.codon.parms)
test.gene=glf.test.1[[1]]
codon.parms$init=codon.parms$obs

cat("\nObtaining samples from the posterior distribution of phi for one simulated gene...\n")
ret = main.semppr(glf.test.1,codon.parms,phi,aa.list='None',parallel='lapply',out.prefix=NULL)
phi.trace = ret$phi$trace.dat

cat("Testing if real phi lies within the 95% CR of estimated phi...      ")
CR95_phi = quantile(phi.trace[floor(length(phi.trace)/1.5):length(phi.trace)],probs=c(0.025,0.975))

if(test.gene$phi.value<CR95_phi[2]&&test.gene$phi.value>CR95_phi[1]){
  cat("Good\n")
}else{
  cat("\nReal phi did not lie within the 95% CR of estimated phi...\nTesting if real phi lies within the 99% CR of estimates phi...      ")
  CR99_phi = quantile(phi.trace[floor(length(phi.trace)/1.5):length(phi.trace)],probs=c(0.005,0.995))
  if(test.gene$phi.value<CR99_phi[2]&&test.gene$phi.value>CR99_phi[1]){
    cat("Good\n")
  }else{
    cat("Error\nReal phi did not lie within the 99% CR of estimated phi.\nSince this process is stochastic, this error does not necessarily\nmean that the model is not working properly. Run the test case again,\nand if this error repeatedly appears, the model isn't working properly.\n")
  }
}
