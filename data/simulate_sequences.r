#setwd("~/CodonUsageBias/R_nse")
rm(list=ls())
library(NSEexchange)

dna.file <- "data_files/genome/S.cerevisiae.S288c.fasta"
mut.file <- "data_files/codon/Mut_rate/PNAS2011_S.cer.mut.tsv"
elong.rate.file <- "data_files/codon/Elong_rate/S.cerevisiae.2007.splitSer.tsv"
phi.file <- "data_files/phi/beyer.phi.tsv"

result <- NSEexchange:::semppr.data.load(dna.file = dna.file, mut.file = mut.file, elong.rate.file = elong.rate.file, phi.file = phi.file, elong.header=F, mut.header=F, phi.header=F)

obs.genome <- result$genome
codon.params <- result$codon.parms
phi <- result$phi

## clean up workspace for better overview during code development 
rm(list=c("dna.file", "mut.file", "elong.rate.file", "phi.file", "result"))
seq.ids <- mclapply(1:length(obs.genome), function(i) {paste(obs.genome[[i]]$name, collapse = "")}, mc.cores=4)
## generate optimal an pesimal sequence
opt.seqs <- mclapply(1:length(obs.genome), function(i) {paste(codon.params$codon[opt.seq(obs.genome[[i]], codon.params)], collapse = "") }, mc.cores=4 )
write.fasta(sequences = opt.seqs, names = seq.ids, file.out = "yeast_opt_seqs.fasta")
pes.seqs <- mclapply(1:length(obs.genome), function(i) {paste(codon.params$codon[pess.seq(obs.genome[[i]], codon.params)], collapse = "") }, mc.cores=4 )
write.fasta(sequences = pes.seqs, names = seq.ids, file.out = "yeast_pes_seqs.fasta")

## simulate sequence
obs.genome.parms <- list(A.1 = 4, A.2 = 4, B = 0.0025, Ne = 1, Q = 1)
sim.genome <- simulate.data.all.genes(obs.genome = obs.genome, obs.phi = phi, obs.codon.parms = codon.params, obs.genome.parms = obs.genome.parms, BIS=4, GMT=0, MES=0, n.cores=4, aux.simulation.type = 'M')
sim.seqs <- mclapply(1:length(sim.genome), function(i) {paste(sim.genome[[i]]$gene.dat$codon, collapse = "")}, mc.cores=4)
write.fasta(sequences = sim.seqs, names = seq.ids, file.out = "yeast_sim_seqs.fasta")
