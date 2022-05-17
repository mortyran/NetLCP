#' @title Data initialize
#' @description Initialize the data before using NetLCP.
#'
#'
#' @examples
#' dataInitialize()
#' @export

dataInitialize = function(){
  filePath = system.file("extdata", package = "NetLCP")
  # download
  if(file.exists(system.file("extdata", "NetLCPData.tar.gz", package = "NetLCP"))){
    print("NetLCP initializing, please be patient while we do something......")
    utils::untar(system.file("extdata", "NetLCPData.tar.gz", package = "NetLCP"), exdir = filePath)
    file.remove(paste0(filePath, "/NetLCPData.tar.gz"))
    if(file.exists(system.file("extdata", "extraction.rda", package = "NetLCP"))){
      return("NetLCP has been initialized!")
    }else{
      return("Some errors occur......")
    }
  }else if(file.exists(system.file("extdata", "extraction.rda", package = "NetLCP"))){
    return("Data initialization has been finished!")
  }else{
    print("NetLCP initializing, please be patient while we do something......")
    fileUrl = "http://hainmu-biobigdata.com/NetLCP/NetLCPData.tar.gz"
    utils::download.file(url = fileUrl, destfile = system.file("extdata", "NetLCPData.tar.gz", package = "NetLCP"))
    utils::untar(paste0(filePath, "/NetLCPData.tar.gz"), exdir = filePath)
    file.remove(paste0(filePath, "/NetLCPData.tar.gz"))
    if(file.exists(system.file("extdata", "extraction.rda", package = "NetLCP"))){
      return("NetLCP has been initialized!")
    }else{
      return("Some errors occur......")
    }
  }
  return("Some errors occur......")
}
