semppr.data.load <- function(dna.file,
                             mut.file,
                             elong.rate.file=NULL,
                             elong.pr.file=NULL,
                             nse.pr.file=NULL,
                             phi.file=NULL,
                             norm.expr=FALSE,
                             Q=1,Ne=1,A1=4,A2=4,B=0.0025,
                             mut.header=TRUE,
                             elong.header=TRUE) {
  # Purpose: Easily load data in correct format for use with    
  #    calc.llik functions.                                     
  # Inputs: dna.file = name specifying dna file in fasta format. This file is required.
  #         nse.file = name specifying elongation rate file containing AA, codon nt's, 
  #           and elongation rates. NOTE: this is still elongation rate. 
  #         mut.file = name specifying mutation rate file containing AA, codon nt's, and
  #           elgonation rates. 
  #         norm.expr = flag to decide whether to scale         
  #           expression to have median 1.                        
  # Outputs: list(gene.list.full,codon.parms)                      
  # Usage: > Semppr_data_load(dna.name='/path/to/a.fasta',      
  #        + codon.parms='/path/to/b.tsv')      
  require(seqinr)
  require(stringr)
  require(plyr)
  
  #1) Make sure we have correct files 
  if(is.null(elong.rate.file)&&is.null(elong.pr.file)&&is.null(nse.pr.file)){
    stop("Error: In semppr.data.load - Include either a elong rate file, elong pr file, or NSE pr file.")
  }
  
  #1) Read in fasta files using mikes approach.
  sc.fasta = read.fasta(file = dna.file, as.string = TRUE)
  
  num.genes = length(sc.fasta)
  
  #2) Read in parameters and other shit.  Read tRNA Files.

    ##################################### 
    # Read from elongation and mutation #
    # files                             #
    #####################################
    
    if(!is.null(nse.pr.file)){
      elong.dat = read.delim(nse.pr.file,header=elong.header,sep = "\t",stringsAsFactors = FALSE)
      names(elong.dat) = c("AA","Codon","NsePr")
    }else if(!is.null(elong.pr.file)){
      elong.dat = read.delim(elong.pr.file,header=elong.header,sep = "\t",stringsAsFactors = FALSE)
      names(elong.dat) = c("AA","Codon","ElongPr")
    }else if(!is.null(elong.rate.file)){
      elong.dat = read.delim(elong.rate.file,header=elong.header,sep = "\t",stringsAsFactors = FALSE)
      names(elong.dat) = c("AA","Codon","ElongRate")
    }
    

    
    # Check for errors
    if(length(elong.dat[,1])!=64){
      if(length(elong.dat[,1])<61){
        stop("Error: Too few codons included in elongation/NSE file.\n")
      }else if(length(elong.dat[,1])==64&&any(elong.dat[62:64,1]!=c('X','X','X'))){
        stop("Error: Stop codons must occur at the end of the elongation file.\n")
      }else if(length(which(elong.dat[,1]=='X'))==0&&length(elong.dat[,1]==61)){
        elong.dat[62,] = c('X','TAA',0)
        elong.dat[63,] = c('X','TGA',0)
        elong.dat[64,] = c('X','TAG',0)
      }else{
        stop("Unknown Error in elongation file. Elongation file must contain 64 codons with stop codons or 61 codons without stop codons. Stop codons must appear at the end of the file.")
      }
    }
    
    mut.dat = read.delim(mut.file,header=mut.header,sep="\t",stringsAsFactors = FALSE)
    names(mut.dat) = c("AA","Codon","MutRate")
        
    # Check for errors
    if(length(mut.dat[,1]!=64)){
      if(length(mut.dat[,1])<61){
        stop("Error: Too few codons included in mutation file.\n")
      }else if(length(mut.dat[,1])==64&&any(mut.dat[62:64,1]!=c('X','X','X'))){
        stop("Error: Stop codons must occur at the end of the elongation file.\n")
      }else if(length(which(mut.dat[,1]=='X'))==0&&length(mut.dat[,1]==61)){
        mut.dat[62,] = c('X','TAA',0)
        mut.dat[63,] = c('X','TGA',0)
        mut.dat[64,] = c('X','TAG',0)
      }else{
        stop("Unknown Error in mutation file. Mutation file must contain 64 codons with stop codons or 61 codons without stop codons. Stop codons must appear at the end of the file.")
      }
    }
    
    #Assemble dataframe from the two dataframes
    codon.parms = elong.dat[,1:2]
    codon.parms[,3] = 0:63 #codon indices
    
    #Make sure codons are in the same order
    merged.dat = merge(elong.dat,mut.dat,by="Codon",sort=FALSE)

    if(!is.null(nse.pr.file)){
      codon.parms[,4] = 1-as.numeric(merged.dat$NsePr)
      codon.parms[,5] = as.numeric(merged.dat$NsePr)
      codon.parms[,6] = as.numeric(merged.dat$MutRate)
      names(codon.parms) = c("aa","codon","c_index","elong_pr","nse_pr","mut_rate")
    }else if(!is.null(elong.pr.file)){
      codon.parms[,4] = as.numeric(merged.dat$ElongPr)
      codon.parms[,5] = 1-as.numeric(merged.dat$ElongPr)
      codon.parms[,6] = as.numeric(merged.dat$MutRate)
      names(codon.parms) = c("aa","codon","c_index","elong_pr","nse_pr","mut_rate")
    }else if(!is.null(elong.rate.file)){
      codon.parms[,4] = as.numeric(merged.dat$ElongRate)/(as.numeric(merged.dat$ElongRate)+B) #elong_pr
      codon.parms[,5] = B/(as.numeric(merged.dat$ElongRate)+B) #nse_pr
      codon.parms[,6] = as.numeric(merged.dat$MutRate)
      codon.parms[,7] = as.numeric(merged.dat$ElongRate)
      names(codon.parms) = c("aa","codon","c_index","elong_pr","nse_pr","mut_rate","elong_rate")
    }
  
  #tabulate the number of codons associated with each amino acid.
  aa.count = as.data.frame(table(codon.parms$aa))
  names(aa.count)=c("aa","count")
  codon.parms = merge(codon.parms,aa.count,by.x="aa",by.y="aa")

  #Read in the z = phi data.
  
  if(!is.null(phi.file)){
    phi.dat = read.csv(phi.file,header=TRUE,stringsAsFactors=FALSE)
    phi = phi.dat[,2]
    if(norm.expr){
      phi.dat.norm = phi.dat[,2]/median(phi.dat[,2]) # transform data so that it has median = 1
      phi=phi.dat.norm
    }
  }else{
    phi.dat=NULL #initialize
    is.phi.included=TRUE
    for(i in 1:num.genes){
      name = attr(sc.fasta[[i]],'name')
      phi.dat[i] = as.numeric(str_split(name,'\t')[[1]][3])
      if(is.na(phi.dat[i])){
        is.phi.included=FALSE
      }
    }
    if(!is.phi.included){
      warning("Phi file not included and phi not contained or 
              formatted correctly in FASTA file. Setting all 
              phi values to NA.")
      for(i in 1:num.genes){
        phi.dat[i]=NA
      }
    }
    phi = phi.dat
    if(norm.expr&&is.phi.included){
      phi.dat.norm = phi.dat/median(phi.dat)
      phi=phi.dat.norm
    }
  }
  
  #3)Create a data frame to store codons in.  Within each gene it records each codon as a row 
  #with the coding amino acid, the codon, and the parameter values of the elongation and 
  #the translation rates.  At the optimization step will need to replace parameter values
  #with updated versions.
  
  
  #computing the size of the genome.
  lcodon = 0
  for(i in 1:num.genes){
    temp  = noquote(sc.fasta[[i]][1])
    lcodon = lcodon+nchar(temp)/3  
  }
  
  codon.paste = function(A){
    
    answer = paste(A[1],A[2],A[3],sep="")  
    
  }
  
  gene.matrix = function(sc.fasta,codon.parms,lcodon,phi){  
    #recasting as a function so that we can compile
    
    num.genes=length(sc.fasta)
    gene.list=vector(mode="list",length=num.genes)
    
    for(i in 1:num.genes){
      
      #computing number of amino acids in a codon
      temp  = s2c(toupper(noquote(sc.fasta[[i]][1])))  #make base string uppercase and and then divide into separate words.
      ltemp = length(temp)
      numb.aa  = ltemp/3
      codon.id= as.data.frame(apply(matrix(temp,numb.aa,3,byrow=TRUE),1,codon.paste))
      names(codon.id) = "codon"
      gene.dat = join(codon.id,codon.parms,by ="codon",type="left",match='all')
      
      gene.dat$c_index[numb.aa] = 99  #last row is all missing values
      name = attr(sc.fasta[[i]],'name')
      gene.list[[i]] = list(name=name,phi.value=phi[i],gene.dat=gene.dat)
      if(i%%50 ==0) print(i) 
      
    }
    answer = gene.list
  } 
  
  gene.list.full = gene.matrix(sc.fasta,codon.parms,lcodon,phi)
  
  phi.dat = unlist(lapply(gene.list.full,function(list) list$phi))
  
  
  answer = list(gene.list.full,codon.parms,phi.dat)
}