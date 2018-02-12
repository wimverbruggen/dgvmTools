#' Sensitivity analysis ED2 (one dimensional version)
#' Note that this will look for T-output files only! This assumes a certain directory structure of the ED2 output as well:
#'
#'                         /path/param.var_value/T-file-output.h5
#'
#' @param param.var The parameter that is being varied (prefix of the folders)
#' @param param.eval The parameter(s) of which we want to see the response (ED2 output names)
#' @param type Which kind of output files?
#' @param path Path to the sensitivity analysis
#'
#' @return A sensitivity analysis?
#' @export
sensitivity.ED2 <- function(param.var, param.eval, site, type, path, verbose=FALSE){

  # List all directories starting with param.var in "path" and get the parameter value
  sensdirs <- list.dirs(path = path, full.names = TRUE, recursive = FALSE)
  sensitivity <- list()

  # Go over all directories and get the ED2 output
  for(d in sensdirs){

    if(verbose) cat(" - Reading file",file.path(d),"\n")

    # Get the value of our parameter
    parval <- unlist(strsplit(d,split=paste0(param.var,"_")))[2] # watch out, this will only work if there's no other "param.val_" string in the path

    # Get the ED2 data
    edout <- readVar.ED2(path=file.path(d),sites=site,type=type,vars=param.eval)

    # Construct a sensitivity analysis object
    sensitivity[[paste(param.var,parval,sep="_")]] <- edout[[site]]
    sensitivity[[paste(param.var,parval,sep="_")]]$sens.param <- parval
  }

return(sensitivity)

}
