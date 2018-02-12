#' Read LPJ-GUESS output
#'
#' @param path.run Path to the LPJ-GUESS run, with all the output in a "output" subdirectory
#' @param sites Sites we are interested in
#' @param delete_nye [T/F] because there are some significant spikes on that day due to yearly allocation
#' @return Dataframe with all output
#'
#' @export
readOut.LPJG <- function(path.run,sites,delete_nye=FALSE){

  # TO DO
  # - Choose a good standard output format for ED! Preferably in CF format or so.

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  #                                                                                                 #
  #   P R E P A R E                                                                                 #
  #                                                                                                 #
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


  # Specify which datafiles we want to incorporate
  vars.daily  <- c("dmet","dflux","daily_aet")

  # The time and location related variables
  vars.postim <- c("Lon", "Lat", "Year", "Day")


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  #                                                                                                 #
  #   R E A D                                                                                       #
  #                                                                                                 #
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

  # Load packages
  suppressMessages(require(lubridate))

  # List all files we're interested in
  path.out <- paste0(path.run,"output")
  files  <- list.files(path = path.out, pattern = "\\.out$")

  # Read them all
  lpjg.output <- list()
  datalist = sub(".out","",files)
  for (i in 1:length(datalist)) lpjg.output[[datalist[i]]] <- read.table(file.path(path.out,files[i]),header=TRUE)

  # Read the gridlist as well
  gridfile  <- list.files(path = path.run, pattern = "gridlist.txt")
  griddata  <- read.table(file.path(path.run,gridfile[1]),header=FALSE,fill=TRUE,stringsAsFactors=FALSE)

  gridlist <- list()
  gridlist$lon  <- round(griddata$V1,digits=2)
  gridlist$lat  <- round(griddata$V2,digits=2)
  gridlist$site <- griddata$V3 # note that if the site label contains of more words, only the first one will be used (eg. "Wankama" instead of "Wankama Millet")

  # Rename sites to codes (by hand..) (TODO: change gridlist files so that they contain the site codes.. then delete below)
  gridlist$site[gridlist$site=="Dahra"]    <- "DAH"
  gridlist$site[gridlist$site=="Demokeya"] <- "DEM"
  gridlist$site[gridlist$site=="Agoufou"]  <- "AGG"
  gridlist$site[gridlist$site=="Kelma"]    <- "KEL"
  gridlist$site[gridlist$site=="Wankama"]  <- "WFM"

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  #                                                                                                 #
  #   C O N V E R T                                                                                 #
  #                                                                                                 #
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

  ##
  ## DAILY OUTPUT
  ##

  # Set up datalist
  lpjg.d <- list()
  postim <- lpjg.output[[vars.daily[[1]]]][vars.postim]
  for(i in 1:length(gridlist$site)){

    # Select datapoints corresponding to current site i
    sel <- which(postim$Lon == gridlist$lon[i] & postim$Lat == gridlist$lat[i])

    site <- gridlist$site[i]
    lpjg.d[[site]] <- list()

    # Set up POSIX
    lpjg.d[[site]]$posix <- as.POSIXct(paste(postim$Year[sel],postim$Day[sel],sep="-"),format="%Y-%j",tz=siteinfo[[site]]$tz)

    # Fill up with the data of interest
    for(vd in vars.daily){
      vars <- names(lpjg.output[[vd]]); vars <- vars[!(vars %in% vars.postim)]
      for(v in vars){
        lpjg.d[[site]][[v]] <- lpjg.output[[vd]][[v]][sel]
      }
    }
  }


  ##
  ## Rename variables
  ##

  var.rename <- list( Rainfall = "rain.rate",
                      Tair = "t.air",
                      SolIrrInc = "swinc",
                      mnee = "NEE",
                      mnpp = "NPP",
                      mgpp = "GPP",
                      mra  = "RA",
                      rh  = "RH"
  )

  ##
  ## Select sites (maybe we should do this earlier..)
  ##

  for(s in names(lpjg.d)){
    if(!s %in% sites) lpjg.d[[s]] <- NULL
  }


  # Actual renaming
  for(s in sites){
    for(v in names(lpjg.d[[s]])){
      if(v %in% names(var.rename)){
        lpjg.d[[s]][var.rename[[v]]] <- lpjg.d[[s]][v]
        lpjg.d[[s]][v] <- NULL
      }
    }
  }


  ## ----------------------------------------------------------------
  ## Change units and adopt to conventions
  ## ----------------------------------------------------------------

  sec_per_day <- 24*60*60


  for(s in sites){

    # Convert carbon fluxes from [kgC/m2/day] to [kgC/m2/sec]
    lpjg.d[[s]]$NEE <- -lpjg.d[[s]]$NEE/sec_per_day
    attr(lpjg.d[[s]]$NEE,"Units") <- "kgC m-2 s-1"

    lpjg.d[[s]]$GPP <- lpjg.d[[s]]$GPP/sec_per_day
    attr(lpjg.d[[s]]$GPP,"Units") <- "kgC m-2 s-1"

    lpjg.d[[s]]$NPP <- lpjg.d[[s]]$NPP/sec_per_day
    attr(lpjg.d[[s]]$NPP,"Units") <- "kgC m-2 s-1"

    lpjg.d[[s]]$RA  <- -lpjg.d[[s]]$RA/sec_per_day
    attr(lpjg.d[[s]]$RA,"Units") <- "kgC m-2 s-1"

    lpjg.d[[s]]$RH  <- -lpjg.d[[s]]$RH/sec_per_day
    attr(lpjg.d[[s]]$RH,"Units") <- "kgC m-2 s-1"

    # Convert rain.rate to mm/sec
    lpjg.d[[s]]$rain.rate  <- lpjg.d[[s]]$rain.rate/sec_per_day
    attr(lpjg.d[[s]]$rain.rate,"Units") <- "mm s-1"

    # Evapotranspiration
    lpjg.d[[s]]$ET <- lpjg.d[[s]]$ET/sec_per_day
    attr(lpjg.d[[s]]$ET,"Units") <- "mm s-1"
  }

  ## -------------------------------------------------------------------------------------------------------------
  ## Remove Dec 31st results, as there are usually some spikes resulting from yearly allocation.
  ## -------------------------------------------------------------------------------------------------------------

  if(delete_nye){
    for(s in sites){
      to_kill <- which(yday(lpjg.d[[s]]$posix)==1)-1 # 1st DOY minus one
      for(v in names(lpjg.d[[s]])){
        if(v != "posix"){
          lpjg.d[[s]][[v]][to_kill] <- NA
        }
      }
    }
  }


  return(lpjg.d)

}
