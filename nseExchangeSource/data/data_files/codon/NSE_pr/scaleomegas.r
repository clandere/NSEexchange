omega_unscaled <- read.table("S.cerevisiae.2007.NSE.csv")
colnames(omega_unscaled) <- c("AA", "Codon", "Nse Prob")

omega_vec <- numeric(length(omega_unscaled[[length(omega_unscaled)]]));
for(i in 1:length(omega_unscaled[[length(omega_unscaled)]])){
  omega_vec[i] <- omega_unscaled[[length(omega_unscaled)]][i] / (1 - omega_unscaled[[length(omega_unscaled)]][i]);
}  

omega_unscaled[[length(omega_unscaled)+1]] <- omega_vec;
colnames(omega_unscaled)[length(omega_unscaled)] <- "unscaled omega"

AAvec <- unique(omega_unscaled[[1]]);
omega_true <- omega_unscaled;
omega_vec <- numeric();
for(i in 1:length(AAvec)){
  tmpvec <- omega_unscaled[[length(omega_unscaled)]][grep(pattern = AAvec[i], x = omega_unscaled[[1]])];
  tmpvec <- tmpvec - tmpvec[length(tmpvec)]
  tmpvec <- tmpvec[-length(tmpvec)]
  omega_vec <- c(omega_vec, tmpvec)
  
  omega_true <- omega_true[-(length(omega_vec)+1),]
}

omega_true[[length(omega_true)+1]] <- omega_vec;
colnames(omega_true)[length(omega_true)] <- "scaled omega";
write.table(x=omega_true[,c(1,2,length(omega_true))], file="scaled_omega.csv",quote=FALSE, sep="\t", row.names=FALSE);