plot.codon.parms <- function(codon.parms, aa.list, codon.list=NULL, n.col=4, burnin=0.5, ...){
  require(squash)
  
  if(!is.null(aa.list))
  {
    num.aa <- length(aa.list)
    ind.aa <- lapply(aa.list, function(aa) which(codon.parms$trace[[1]]$aa == aa))
    
    num.codon <- unlist(lapply(ind.aa, function(ind)length(ind)))
    mat <- NULL
    n.row <- rep(0,2)
    sum.n.codon <- 0
    for(i in 1:num.aa)
    {
      n.row[i] <- ceiling(num.codon[i]/n.col)
      if(num.codon[i] == n.col)
      {
        mat <- c(mat, 1:num.codon[i] + sum.n.codon)
      }else{
        fill <- n.col - (num.codon[i] %% n.col)
        if(fill == n.col) fill <- 0
        mat <- c(mat, 1:num.codon[i] + sum.n.codon, rep(0, fill))
      }
      sum.n.codon <- sum.n.codon + num.codon[i]
    }
    print(n.row)
    mat <- matrix(mat, sum(n.row), n.col, byrow=TRUE)
    
    par(oma=c(0, 3, 0, 0))
    layout(mat=mat)
    
    row.space <- 1/sum(n.row)
    
    for(i in 1:length(ind.aa))
    {
      for(j in 1:length(ind.aa[[i]]))
      {
        nse.trace <- unlist(lapply(codon.parms$trace,
                                  function(trace)log(1-as.numeric(trace$elong_pr[ind.aa[[i]][j]]))))
        mut.trace <- unlist(lapply(codon.parms$trace,
                                  function(trace) log(as.numeric(trace$mut_rate[ind.aa[[i]][j]]))))
        hist2(x=nse.trace[floor(burnin*length(nse.trace)+1):length(nse.trace)],
              y=mut.trace[floor(burnin*length(mut.trace)+1):length(mut.trace)],
              xlab='log(NSE Pr)\n',
              ylab='log(Mut Rate)',
              main=paste('Codon', codon.parms$trace[[1]]$codon[ind.aa[[i]][j]]),
              ...)
        if(!is.null(codon.parms$obs))
        {  
          mut.obs <- log(as.numeric(codon.parms$obs$mut_rate[ind.aa[[i]][j]]))
          nse.obs <- log(1-as.numeric(codon.parms$obs$elong_pr[ind.aa[[i]][j]]))
          
          abline(h=mut.obs,col='green')
          abline(v=nse.obs,col='green')
        }
        
      }
      #figure out where to put AA label
      if(i == 1)
      {
        space.above <- 0
        space.of <- sum(n.row[i]) * row.space
        pos.of <- 1-(space.above + space.of/2)
        print(pos.of)
        print(paste('Amino Acid', aa.list[1]))
        
        mtext(text = paste('Amino Acid', aa.list[i]), side = 2, outer=TRUE, at=c(pos.of), cex=1.3, line=0.5)
      }else{
        space.above <- sum(n.row[1:i-1])*row.space
        space.of <- sum(n.row[i])*row.space
        pos.of <- 1-(space.above + space.of/2)
        print(pos.of)
        mtext(text = paste('Amino Acid', aa.list[i]), side = 2, outer=TRUE, at=c(pos.of), cex=1.3, line=0.5)
      }
    }
  }
}
  