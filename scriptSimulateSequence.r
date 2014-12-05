setwd("~/cubfits/preston")
rm(list=ls())
require(plyr, lib.loc = "~/cubfits/Dependencies/")
require(parallel)
library(NSEexchange, lib='~/cubfits/preston/')

data.folder <- "nseExchangeSource/data/data_files/"
dna.file <- paste(data.folder, "genome/S.cerevisiae.S288c.fasta", sep="");
mut.file <- paste(data.folder, "codon/Mut_rate/PNAS2011_S.cer.mut.tsv", sep="")
elong.rate.file <- paste(data.folder, "codon/Elong_rate/S.cerevisiae.2007.splitSer.tsv", sep="")
phi.file <- paste(data.folder, "phi/beyer.phi.tsv", sep="")

cat("Loading data . . .\n")
result <- NSEexchange:::semppr.data.load(dna.file = dna.file, mut.file = mut.file, elong.rate.file = elong.rate.file, phi.file = phi.file, elong.header=F, mut.header=F, phi.header=F)

obs.genome <- result$genome
codon.params <- result$codon.parms
phi <- result$phi

## clean up workspace for better overview during code development 
rm(list=c("dna.file", "mut.file", "elong.rate.file", "phi.file", "result", "data.folder"))
seq.ids <- mclapply(1:length(obs.genome), function(i) {paste(obs.genome[[i]]$name, collapse = "")}, mc.cores=4)

## simulate sequence
obs.genome.parms <- list(A.1 = 4, A.2 = 4, B = 0.0025, Ne = 1, Q = 1)
cat("Simulating Genome. . .\n")
sim.genome <- simulate.data.all.genes(obs.genome = obs.genome, obs.phi = phi, obs.codon.parms = codon.params, obs.genome.parms = obs.genome.parms, BIS=4, GMT=0, MES=0, n.cores=4, aux.simulation.type = 'M')
sim.seqs <- mclapply(1:length(sim.genome), function(i) {paste(sim.genome[[i]]$gene.dat$codon, collapse = "")}, mc.cores=8)
cat("Writing fasta. . .\n")
write.fasta(sequences = sim.seqs, names = seq.ids, file.out = "loganSimYeast.fasta")
