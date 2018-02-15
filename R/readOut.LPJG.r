#' Read LPJ-GUESS output
#'
#' @param path.run Path to the LPJ-GUESS run, with all the output in a "output" subdirectory
#' @param sites Sites we are interested in
#' @param freq Frequency of the output. Possible values: day, year
#' @param delete_nye [T/F] because there are some significant spikes on that day due to yearly allocation
#' @return Dataframe with all output
#'
#' @export
readOut.LPJG <- function(path.run,sites,freq="day",delete_nye=FALSE){

  # Input parameter control
  acceptable_freqs <- c("day","year")
  if(!freq %in% acceptable_freqs) stop("Frequency of output should be \"day\" or \"year\"")

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
  #                                                                                                 #
  #   P R E P A R E                                                                                 #
  #                                                                                                 #
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


  # Specify which datafiles we want to incorporate
  vars.daily  <- c("dmet","dflux","daily_aet")
  vars.yearly <- c("aaet","lai")


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
  ## Get output
  ##

  # Set up datalist
  lpjg.d <- list()
  lpjg.y <- list()

  postim.d <- lpjg.output[[vars.daily[[1]]]] [c("Lon","Lat","Year","Day")]
  postim.y <- lpjg.output[[vars.yearly[[1]]]][c("Lon","Lat","Year")]

  for(i in 1:length(gridlist$site)){

    if(gridlist$site[i] %in% sites){

      # Select datapoints corresponding to current site i
      sel.d <- which(postim.d$Lon == gridlist$lon[i] & postim.d$Lat == gridlist$lat[i])
      sel.y <- which(postim.y$Lon == gridlist$lon[i] & postim.y$Lat == gridlist$lat[i])

      site <- gridlist$site[i]
      lpjg.d[[site]] <- list()
      lpjg.y[[site]] <- list()

      # Set up POSIX
      lpjg.d[[site]]$posix <- as.POSIXct(paste(postim.d$Year[sel.d],postim.d$Day[sel.d],sep="-"),format="%Y-%j",tz=siteinfo[[site]]$tz)
      lpjg.y[[site]]$posix <- as.POSIXct(paste(postim.y$Year[sel.y],1,sep="-"),format="%Y-%j",tz=siteinfo[[site]]$tz)

      # Fill up with data: days
      for(vd in vars.daily){
        vars <- names(lpjg.output[[vd]]); vars <- vars[!(vars %in% vars.postim)]
        for(v in vars){
          lpjg.d[[site]][[v]] <- lpjg.output[[vd]][[v]][sel.d]
        }
      }

      # Fill up with data: years
      for(vy in vars.yearly){
        pfts <- names(lpjg.output[[vy]]); pfts <- pfts[!(pfts %in% vars.postim)]
        for(p in pfts){
          lpjg.y[[site]][[vy]][[p]] <- lpjg.output[[vy]][[p]][sel.y]
        }
        lpjg.y[[site]][[vy]] <- t(matrix(unlist(lpjg.y[[site]][[vy]]),ncol=13))
        rownames(lpjg.y[[site]][[vy]]) <- pfts
        colnames(lpjg.y[[site]][[vy]]) <- unique(year(lpjg.y[[site]]$posix))
      }

    }

  }


  ##
  ## Rename variables
  ##

  var.rename <- list( Rainfall = "rain",
                      Tair = "air.temp",
                      SolIrrInc = "rad.swinc",
                      mnee = "NEE",
                      mnpp = "NPP",
                      mgpp = "GPP",
                      mra  = "RA",
                      rh  = "RH",
                      aaet = "ET.pft",
                      lai = "LAI.pft"

  )

  # Actual renaming
  for(s in sites){
    for(v in names(lpjg.d[[s]])){
      if(v %in% names(var.rename)){
        lpjg.d[[s]][var.rename[[v]]] <- lpjg.d[[s]][v]
        lpjg.d[[s]][v] <- NULL
      }
    }
    for(v in names(lpjg.y[[s]])){
      if(v %in% names(var.rename)){
        lpjg.y[[s]][var.rename[[v]]] <- lpjg.y[[s]][v]
        lpjg.y[[s]][v] <- NULL
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

    lpjg.d[[s]]$RE  <- lpjg.d[[s]]$RA + lpjg.d[[s]]$RH
    attr(lpjg.d[[s]]$RE,"Units") <- "kgC m-2 s-1"

    # Convert rain to mm/sec
    lpjg.d[[s]]$rain  <- lpjg.d[[s]]$rain/sec_per_day
    attr(lpjg.d[[s]]$rain,"Units") <- "mm s-1"

    # Daily evapotranspiration
    lpjg.d[[s]]$ET <- lpjg.d[[s]]$ET/sec_per_day
    attr(lpjg.d[[s]]$ET,"Units") <- "mm s-1"

    # Annual output > units not changed
    attr(lpjg.y[[s]]$ET.pft, "Units") <- "mm year-1"
    attr(lpjg.y[[s]]$LAI.pft,"Units") <- "m2 m-2"

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


  if(freq=="day")  return(lpjg.d)
  if(freq=="year") return(lpjg.y)


}
