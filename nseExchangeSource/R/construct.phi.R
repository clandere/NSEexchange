construct.phi <- function(obs.genome=NULL,init.phi.file=NULL,init.phi.header=FALSE,obs.phi.file=NULL,obs.phi.header=FALSE,norm.expr=FALSE){
  #Check to see if we have the right information
  if(all(is.null(c(obs.genome,init.phi.file,obs.phi.file)))){
    stop('In construct.phi: Must include phi information.')
  }
  
  phi <- list()
  
  #$init
  if(!is.null(init.phi.file)){
    cat('Reading initial phi from file...\n')
    phi.dat = read.delim(init.phi.file,header=init.phi.header,stringsAsFactors=FALSE,'\t')
    phi$init = phi.dat[,2]
    if(norm.expr){
      phi.dat.norm = phi.dat[,2]/median(phi.dat[,2]) # transform data so that it has median = 1
      phi$init=phi.dat.norm
    }
    names(phi$init)<-phi.dat[,1]
  }else{
    cat('Initial phi values are not included...\n')
  }
  
  dim(phi$init) <- c(1,length(phi$init))
  
  #$obs
  if(!is.null(obs.phi.file)){
    cat('Reading observed phi from file...\n')
    phi.dat = read.delim(obs.phi.file,header=obs.phi.header,stringsAsFactors=FALSE,'\t')
    phi$obs = phi.dat[,2]
    if(norm.expr){
      phi.dat.norm = phi.dat[,2]/median(phi.dat[,2]) # transform data so that it has median = 1
      phi$obs=phi.dat.norm
    }
    names(phi$obs)<-phi.dat[,1]
  }else if(!is.null(obs.genome)){
    #check to see if obs phi is included in obs.genome
    obs.phi2 = unlist(lapply(obs.genome,function(list)list$phi.value))
    if(!is.null(obs.phi2)){
      cat('Setting observed phi values to phi values found in genome[[i]]$phi.value...\n')
      phi$obs <- obs.phi2
    }else{
      cat('Observed phi values are not included...\n')
      phi$obs <- NULL
    }
  }
  
  dim(phi$obs) <- c(1,length(phi$obs))
  
  class(phi) <- "phi"
  
  phi
}