
load.genome <- function(genome.fasta, codon.parms, phi=NULL, verbose=TRUE){  
  #recasting as a function so that we can compile
  
  sc.fasta = seqinr::read.fasta(file = genome.fasta, as.string = TRUE)
  
  num.genes=length(sc.fasta)
  gene.list=vector(mode="list",length=num.genes)
  
  if(is.null(phi))
  {
    phi <- rep(x = NA, num.genes)  
  }
  
  for(i in 1:num.genes){    
    #computing number of amino acids in a codon
    temp  = s2c(toupper(noquote(sc.fasta[[i]][1])))  #make base string uppercase and and then divide into separate words.
    ltemp = length(temp)
    numb.aa  = ltemp/3
    codon.id = as.data.frame(apply(matrix(temp,numb.aa,3,byrow=TRUE),1,codon.paste))
    names(codon.id) = "codon"
    ## TODO: try to replace by merge to not be dependend on plyr package!
    gene.dat = plyr::join(codon.id, codon.parms$obs, by ="codon", type="left", match='all')
    
    gene.dat$c_index[numb.aa] = 99  #last row is all missing values
    name = attr(sc.fasta[[i]],'name')
    gene.list[[i]] = list(name=name,phi.value=phi[i],gene.dat=gene.dat)
    if(i%%50 == 0) {if(verbose) print(i)}
  }
  return(gene.list)
}

codon.paste <- function(A){
  return( paste(A[1],A[2],A[3],sep="") )
}