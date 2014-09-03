plot.phi <- function(phi, burnin = 0.5, main=NULL, lengths=NULL, ret=TRUE){
  require(ggplot2)
  if(is.null(phi$obs)){
    stop('In plot.phi: Must include observed phi values in phi$obs.')
  }
  real.phi = phi$obs
  
  n.samples = nrow(phi$trace.dat)
  num.genes = ncol(phi$trace.dat)
  
  est.phi = apply(phi$trace.dat,2,function(col)mean(col[ceiling(0.5*length(col))]))
  
  if(is.null(main)){
    main <- paste('Estimated vs Observed Log(phi) for',num.genes,'genes.')
  }
  
  x1 = median(log(real.phi))
  y1 = median(log(est.phi))
  plot.data = data.frame("y"=log(est.phi),"x"=log(real.phi))
  eqn.txt = lm_eqn(plot.data)
  if(!is.null(lengths)){
    plot.data = data.frame("Log_Estimated_Phi"=log(est.phi),"Log_Phi_obs"=log(real.phi),"length"=lengths)
    ggplot(plot.data,aes(x=Log_Phi_obs,y=Log_Estimated_Phi))+geom_point(aes(color=length))+
      geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)+
      geom_text(aes(x=-Inf,y=Inf),label=eqn.txt,parse=TRUE,hjust=-0.35,vjust=3.5)+
      ggtitle(label=main)+
      geom_abline(intercept=0,slope=1,color='green')
  }else{
    plot.data = data.frame("Log_Estimated_Phi"=log(est.phi),"Log_Phi_obs"=log(real.phi))
    ggplot(plot.data,aes(x=Log_Phi_obs,y=Log_Estimated_Phi))+geom_point()+
      geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)+
      geom_text(aes(x=-Inf,y=Inf),label=eqn.txt,parse=TRUE,hjust=-0.35,vjust=3.5)+
      ggtitle(label=main)+
      geom_abline(intercept=0,slope=1,color='green')
  }
  
  
}


lm_eqn = function(df){
  m = lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}