load.codon.parms <- function(init.nse.pr.file=NULL,init.elong.pr.file=NULL,
                                  init.elong.rate.file=NULL,init.mut.file=NULL,
                                  init.elong.header=FALSE,init.mut.header=FALSE,
                                  obs.nse.pr.file=NULL,obs.elong.pr.file=NULL,
                                  obs.elong.rate.file=NULL,obs.mut.file=NULL,
                                  obs.elong.header=FALSE,obs.mut.header=FALSE,
                                  B = NULL){
  #1) Make sure we have correct files 
  init.flag = any(!is.null(c(init.nse.pr.file,init.elong.pr.file,init.elong.rate.file)))&&!is.null(init.mut.file)
  obs.flag = any(!is.null(c(obs.nse.pr.file,obs.elong.pr.file,obs.elong.rate.file)))&&!is.null(obs.mut.file)
  
  if(!(init.flag||obs.flag)){
    stop("Error: In construct.codon.parms - Include a complete dataset for either initial or observed values.
              Complete datasets include: (Nonsense Error Pr file OR Elongation Pr file OR Elongation Rate file) AND
                                         Mutation Rate file.")
  }
  
  if(!is.null(init.elong.rate.file)||!is.null(obs.elong.rate.file)){
    if(is.null(B)){
      stop("Error: In construct.codon.parms - Must include background NSE rate when constructing codon parms
              using elongation rates. DEFAULT = 0.0025")
    }
  }
  
  codon.parms <- list()
  
  #### Initial ####
  
  if(init.flag){
  
    if(!is.null(init.nse.pr.file)){
      elong.dat = read.delim(init.nse.pr.file,header=init.elong.header,sep = "\t",stringsAsFactors = FALSE)
      names(elong.dat) = c("AA","Codon","NsePr")
    }else if(!is.null(init.elong.pr.file)){
      elong.dat = read.delim(init.elong.pr.file,header=init.elong.header,sep = "\t",stringsAsFactors = FALSE)
      names(elong.dat) = c("AA","Codon","ElongPr")
    }else if(!is.null(init.elong.rate.file)){
      elong.dat = read.delim(init.elong.rate.file,header=init.elong.header,sep = "\t",stringsAsFactors = FALSE)
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
    
    mut.dat = read.delim(init.mut.file,header=init.mut.header,sep="\t",stringsAsFactors = FALSE)
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
    init.codon.parms = elong.dat[,1:2]
    init.codon.parms[,3] = 0:63 #codon indices
    
    #Make sure codons are in the same order

    merged.dat = merge(elong.dat,mut.dat,by="Codon",sort=FALSE)
    
    
    if(!is.null(nse.pr.file)){
      #init.codon.parms[,4] = log(1-as.numeric(merged.dat$NsePr))
      init.codon.parms[,4] = log(as.numeric(merged.dat$NsePr))
      init.codon.parms[,5] = log(as.numeric(merged.dat$MutRate))
      #names(init.codon.parms) = c("aa","codon","c_index","elong_pr","nse_pr","mut_rate")
      names(init.codon.parms) = c("aa","codon","c_index","nse_pr","mut_rate")
    }else if(!is.null(elong.pr.file)){
      #init.codon.parms[,4] = as.numeric(merged.dat$ElongPr)
      init.codon.parms[,4] = log(1-as.numeric(merged.dat$ElongPr))
      init.codon.parms[,5] = log(as.numeric(merged.dat$MutRate))
      #names(init.codon.parms) = c("aa","codon","c_index","elong_pr","nse_pr","mut_rate")
      names(init.codon.parms) = c("aa","codon","c_index","nse_pr","mut_rate")
    }else if(!is.null(elong.rate.file)){
      #init.codon.parms[,4] = as.numeric(merged.dat$ElongRate)/(as.numeric(merged.dat$ElongRate)+B) #elong_pr
      init.codon.parms[,4] = log(B/(as.numeric(merged.dat$ElongRate)+B)) #nse_pr
      init.codon.parms[,5] = log(as.numeric(merged.dat$MutRate))
      init.codon.parms[,6] = log(as.numeric(merged.dat$ElongRate))
      #names(init.codon.parms) = c("aa","codon","c_index","elong_pr","nse_pr","mut_rate","elong_rate")
      names(init.codon.parms) = c("aa","codon","c_index","nse_pr","mut_rate","elong_rate")
    }    
    
    #tabulate the number of codons associated with each amino acid.
    aa.count = as.data.frame(table(init.codon.parms$aa))
    names(aa.count)=c("aa","count")
    init.codon.parms = merge(init.codon.parms,aa.count,by.x="aa",by.y="aa")
    
    codon.parms$init = init.codon.parms
  }
  
  if(obs.flag){
    ##### Observed ####
    
    if(!is.null(obs.nse.pr.file)){
      elong.dat = read.delim(obs.nse.pr.file,header=obs.elong.header,sep = "\t",stringsAsFactors = FALSE)
      names(elong.dat) = c("AA","Codon","NsePr")
    }else if(!is.null(obs.elong.pr.file)){
      elong.dat = read.delim(obs.elong.pr.file,header=obs.elong.header,sep = "\t",stringsAsFactors = FALSE)
      names(elong.dat) = c("AA","Codon","ElongPr")
    }else if(!is.null(obs.elong.rate.file)){
      elong.dat = read.delim(obs.elong.rate.file,header=obs.elong.header,sep = "\t",stringsAsFactors = FALSE)
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
    
    mut.dat = read.delim(obs.mut.file,header=obs.mut.header,sep="\t",stringsAsFactors = FALSE)
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
    obs.codon.parms = elong.dat[,1:2]
    obs.codon.parms[,3] = 0:63 #codon indices
    
    #Make sure codons are in the same order
    merged.dat = merge(elong.dat, mut.dat, by="Codon",sort=FALSE)
    
    if(!is.null(nse.pr.file)){
      #obs.codon.parms[,4] = log(1-as.numeric(merged.dat$NsePr))
      obs.codon.parms[,4] = log(as.numeric(merged.dat$NsePr))
      obs.codon.parms[,5] = log(as.numeric(merged.dat$MutRate))
      #names(obs.codon.parms) = c("aa","codon","c_index","elong_pr","nse_pr","mut_rate")
      names(obs.codon.parms) = c("aa","codon","c_index","nse_pr","mut_rate")
    }else if(!is.null(elong.pr.file)){
      #obs.codon.parms[,4] = as.numeric(merged.dat$ElongPr)
      obs.codon.parms[,4] = log(1-as.numeric(merged.dat$ElongPr))
      obs.codon.parms[,5] = log(as.numeric(merged.dat$MutRate))
      #names(obs.codon.parms) = c("aa","codon","c_index","elong_pr","nse_pr","mut_rate")
      names(obs.codon.parms) = c("aa","codon","c_index","nse_pr","mut_rate")
    }else if(!is.null(elong.rate.file)){
      #obs.codon.parms[,4] = as.numeric(merged.dat$ElongRate)/(as.numeric(merged.dat$ElongRate)+B) #elong_pr
      obs.codon.parms[,4] = log(B/(as.numeric(merged.dat$ElongRate)+B)) #nse_pr
      obs.codon.parms[,5] = log(as.numeric(merged.dat$MutRate))
      obs.codon.parms[,6] = log(as.numeric(merged.dat$ElongRate))
      #names(obs.codon.parms) = c("aa","codon","c_index","elong_pr","nse_pr","mut_rate","elong_rate")
      names(obs.codon.parms) = c("aa","codon","c_index","nse_pr","mut_rate","elong_rate")
    } 
    
    #tabulate the number of codons associated with each amino acid.
    aa.count = as.data.frame(table(obs.codon.parms$aa))
    names(aa.count)=c("aa","count")
    obs.codon.parms = merge(obs.codon.parms,aa.count,by.x="aa",by.y="aa")
  
    codon.parms$obs = obs.codon.parms
  }
  
  class(codon.parms) <-'codon.parms'
  
  codon.parms
}