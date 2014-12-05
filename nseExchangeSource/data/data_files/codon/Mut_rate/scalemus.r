mu_unscaled <- read.table("PNAS2011_S.cer.mut.tsv")
colnames(mu_unscaled) <- c("AA", "Codon", "unscaled_M")

AAvec <- unique(mu_unscaled[[1]]);
mu_true <- mu_unscaled;
mu_vec <- numeric();
for(i in 1:length(AAvec)){
  tmpvec <- mu_unscaled[[length(mu_unscaled)]][grep(pattern = AAvec[i], x = mu_unscaled[[1]])];
  tmpvec <- tmpvec - tmpvec[length(tmpvec)]
  tmpvec <- tmpvec[-length(tmpvec)]
  mu_vec <- c(mu_vec, tmpvec)
  
  mu_true <- mu_true[-(length(mu_vec)+1),]
}

mu_true[[length(mu_true)+1]] <- mu_vec;
colnames(mu_true)[length(mu_true)] <- "scaled_M";
write.table(x=mu_true[,c(1,2,length(mu_true))], file="scaled_logMu.csv",quote=FALSE, sep="\t", row.names=FALSE);