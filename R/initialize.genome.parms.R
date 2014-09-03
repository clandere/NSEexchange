initialize.genome.parms<-function(init.genome.parms){
  if(!is.null(init.genome.parms)){
    #make sure genome parms are well-behaved
    if(any(!(names(init.genome.parms %in% c('Q','Ne','A1','A2'))))){
      warning('Warning: Input \'init.genome.parms\' is not well-behaved, reverting to default values')
      return(list(Q=1,Ne=1,A1=4,A2=4))
    }else{
      return(init.genome.parms)
    }
  }else{
    #Use default values for genome.parms
    return(list(Q=1,Ne=1,A1=4,A2=4))
  }
}
