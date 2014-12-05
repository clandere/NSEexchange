initialize.output <- function(out.prefix){

  if(!is.null(out.prefix)){
    #create output directory if it doesn't exist
    if(!file.exists(out.prefix))
    {
      dir.create(out.prefix)
    }
    #make sure out.prefix string ends in a '/'
    return(paste(paste(strsplit(out.prefix,'/')[[1]],collapse='/'),'/',sep=''))
  }else{
    return(out.prefix)
  }
}