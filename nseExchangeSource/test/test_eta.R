library(NSEexchange,quietly=TRUE)
data(subset_1gene)
data(S.cer.codon.parms)
load('data/correct.moments.eta.Rda')

test.gene=glf.test.1[[1]]
trna.dat=codon.parms$obs

cat("Testing eta calculation, min/max eta, and eta mean and variance...\n")

outp<-as.numeric(calc.moments.eta(test.gene,trna.dat))

cat("Eta calculation...  ")
if(outp[5]==correct.moments.eta[5]){cat("Good\n")}else{"ERROR\n"}

cat("Eta min...          ")
if(outp[3]==correct.moments.eta[3]){cat("Good\n")}else{"ERROR\n"}

cat("Eta max...          ")
if(outp[4]==correct.moments.eta[4]){cat("Good\n")}else{"ERROR\n"}

cat("Eta mean...         ")
if(outp[1]==correct.moments.eta[1]){cat("Good\n")}else{"ERROR\n"}

cat("Eta var...          ")
if(outp[2]==correct.moments.eta[2]){cat("Good\n")}else{"ERROR\n"}


