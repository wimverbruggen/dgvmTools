#' Make meteo driver .h5 files for the ED2 model. Based mainly on the code by Marcos Longo.
#' TODO: make this a bit more elegant looking.
#'
#' @param meteo list with variables following my naming conventions
#' @param height [m] height above ground level
#' @param name site name, will be used as prefix of all meteo files
#' @param path.header the header prefix, for the ED2 header file (eg. the directory on your HPC where these data will be located)
#' @param path.output output folder where all h5 files will be saved
#' @param verbose Give some process feedback to the user or not
#'
#' @author Wim Verbruggen (wim.verbruggen@ugent.be)
#' @author Marcos Longo
#'
#' @export
meteo.makeDriver.ED2 <- function(meteo,height,name,path.header,path.output,verbose=TRUE){

  # TODO: Check that all variables exist! For some values of height they may not exist yet!


  # Height as numeric and as string
  height.num <- height
  height.str <- paste0(".",height,"m")
  rm(height)

  # Convert to a format readable by the original "make_met_driver" script.
  datum <- list()
  datum$year      <- lubridate::year(  lubridate::with_tz(meteo$posix,tzone="UTC"))
  datum$month     <- lubridate::month( lubridate::with_tz(meteo$posix,tzone="UTC"))
  datum$day       <- lubridate::mday(  lubridate::with_tz(meteo$posix,tzone="UTC"))
  datum$hour      <- lubridate::hour(  lubridate::with_tz(meteo$posix,tzone="UTC"))
  datum$min       <- lubridate::minute(lubridate::with_tz(meteo$posix,tzone="UTC"))
  datum$sec       <- lubridate::second(lubridate::with_tz(meteo$posix,tzone="UTC"))

  datum$atm.prss  <- meteo[[paste0("air.press",height.str)]]            # [Pa]
  datum$atm.tmp   <- meteo[[paste0("air.temp", height.str)]]+273.15     # [K]
  datum$atm.shv   <- meteo[[paste0("air.q"   , height.str)]]            # [kg kg-1]
  datum$atm.vels  <- meteo[[paste0("wind.vel", height.str)]]            # [m s-1]
  datum$atm.vdir  <- datum$atm.vels*0
  datum$rshort.in <- meteo$rad.swinc                                    # [W m-2]
  datum$rlong.in  <- meteo$rad.lwinc                                    # [W m-2]
  datum$rain      <- meteo$rain                                         # [mm s-1]

  info = list( name     = name,
               longname = name,
               lon      = meteo$lon,
               lat      = meteo$lat,
               h5pref   = name,
               height   = height.num,
               dtdat    = unique(diff(lubridate::seconds(wfdei.3h$DAH$posix))) # eg. "1800" for 30 minutes
  )

  require(chron)
  require(h5)

  # Check/make output path
  if(!file.exists(path.output)) dir.create(path.output, showWarnings = TRUE, recursive = TRUE)


  if(verbose) cat(" + Processing data from ",info$longname,"...","\n",sep="")

  # Rename dates
  datum$when    = chron( paste(datum$month,datum$day,datum$year,sep="/")
                         , paste(datum$hour ,datum$min,datum$sec ,sep=":") )
  datum$today   = chron( paste(datum$month,datum$day,datum$year,sep="/") )
  datum$tomonth = chron( paste(datum$month,        1,datum$year,sep="/") )


  # Find the zenith angle, to split radiation into components.
  zen           = calc.zenith(lon=info$lon,lat=info$lat,when=datum$when,ed21=TRUE,zeronight=FALSE,meanval=TRUE,imetavg=1,nmean=60)
  datum$cosz    = zen$cosz

  # Split the incoming radiation into components, or estimate NIR from PAR if PAR measurements are available.                                                           #
  prss = mean(datum$atm.prss)
  rad  = rshort.bdown(rad.in=datum$rshort,atm.prss=datum$atm.prss,cosz=datum$cosz)
  datum$par.beam = rad$par.beam
  datum$par.diff = rad$par.diff
  datum$nir.beam = rad$nir.beam
  datum$nir.diff = rad$nir.diff

  # Decompose wind
  trigo          = (270. - datum$atm.vdir) * pi/180
  datum$atm.uspd = datum$atm.vels * cos(trigo)
  datum$atm.vspd = datum$atm.vels * sin(trigo)

  # Make ED2 output
  if(verbose) cat(" + Making ED-friendly output files... \n")

  # Make sure that the output directory exists, and if not, create it
  siteroot = file.path(path.output,info$h5pref)
  if(!file.exists(siteroot)) dir.create(siteroot)

  # List all possible unique month/year combinations
  unique.tomonth   = unique(datum$tomonth)
  n.unique.tomonth = length(unique.tomonth)

  for (um in sequence(n.unique.tomonth)){
      monyear.now = unique.tomonth[um]

      #----- Get current month and year, 3-letter month, and number of days in the month. -#
      month.now   = nummonths(monyear.now)
      year.now    = numyears (monyear.now)
      daymax.now  = daymax   (month=month.now,year=year.now)
      month.label = toupper(month.abb[month.now])
      year.label  = sprintf("%4.4i",year.now)
      #------------------------------------------------------------------------------------#


      #----- Print banner to entertain the user. ------------------------------------------#
      cat("   - Checking data from ",month.name[month.now]," ",year.now,"...","\n",sep="")
      #------------------------------------------------------------------------------------#


      #----- Find the indices of data that belong to this month and year. -----------------#
      sel    = datum$month == month.now & datum$year == year.now
      #------------------------------------------------------------------------------------#


      #----- Check that all data are there. -----------------------------------------------#
      nsel      = sum(sel)
      nexpected = daymax.now * (24*60*60) / info$dtdat
      #------------------------------------------------------------------------------------#

      #----- If the data are complete, make the output file. ------------------------------#
      if (nsel == nexpected){
        cat("     * Data series is complete, making the arrays...","\n")
        #---------------------------------------------------------------------------------#
        #      Create the matrices that will have the data.  These will use ED/NCEP name  #
        # convention, and will have a fake 2x2 matrix with the same data because the      #
        # tower may be needed in a site that is nearby but not with the same longitude    #
        # and latitude.                                                                   #
        #---------------------------------------------------------------------------------#
        dlwrf = aperm(a=array(datum$rlong.in[sel],dim=c(nsel,2,2)),perm=c(2,3,1))
        nbdsf = aperm(a=array(datum$nir.beam[sel],dim=c(nsel,2,2)),perm=c(2,3,1))
        nddsf = aperm(a=array(datum$nir.diff[sel],dim=c(nsel,2,2)),perm=c(2,3,1))
        prate = aperm(a=array(datum$rain    [sel],dim=c(nsel,2,2)),perm=c(2,3,1))
        pres  = aperm(a=array(datum$atm.prss[sel],dim=c(nsel,2,2)),perm=c(2,3,1))
        sh    = aperm(a=array(datum$atm.shv [sel],dim=c(nsel,2,2)),perm=c(2,3,1))
        tmp   = aperm(a=array(datum$atm.tmp [sel],dim=c(nsel,2,2)),perm=c(2,3,1))
        ugrd  = aperm(a=array(datum$atm.uspd[sel],dim=c(nsel,2,2)),perm=c(2,3,1))
        vbdsf = aperm(a=array(datum$par.beam[sel],dim=c(nsel,2,2)),perm=c(2,3,1))
        vddsf = aperm(a=array(datum$par.diff[sel],dim=c(nsel,2,2)),perm=c(2,3,1))
        vgrd  = aperm(a=array(datum$atm.vspd[sel],dim=c(nsel,2,2)),perm=c(2,3,1))
        #---------------------------------------------------------------------------------#

        #---------------------------------------------------------------------------------#
        #     Make sure that all data are there.  If not, crash!                          #
        #---------------------------------------------------------------------------------#
        if (any(is.na(dlwrf))) stop("dlwrf has missing values!")
        if (any(is.na(nbdsf))) stop("nbdsf has missing values!")
        if (any(is.na(nddsf))) stop("nddsf has missing values!")
        if (any(is.na(prate))) stop("prate has missing values!")
        if (any(is.na(pres ))) stop("pres  has missing values!")
        if (any(is.na(sh   ))) stop("sh    has missing values!")
        if (any(is.na(tmp  ))) stop("tmp   has missing values!")
        if (any(is.na(ugrd ))) stop("ugrd  has missing values!")
        if (any(is.na(vbdsf))) stop("vbdsf has missing values!")
        if (any(is.na(vddsf))) stop("vddsf has missing values!")
        if (any(is.na(vgrd ))) stop("vgrd  has missing values!")
        #---------------------------------------------------------------------------------#


        #---------------------------------------------------------------------------------#
        #     Make the dataset.                                                           #
        #---------------------------------------------------------------------------------#
        thisfile = file.path(siteroot
                             ,paste(info$h5pref,"_",year.label,month.label,".h5",sep=""))
        cat("     * Saving data to ",basename(thisfile),"...","\n",sep="")

        if(file.exists(thisfile)) file.remove(thisfile)
        fileh5 <- h5file(thisfile,"a")
        fileh5["dlwrf"] <- dlwrf
        fileh5["nbdsf"] <- nbdsf
        fileh5["nddsf"] <- nddsf
        fileh5["prate"] <- prate
        fileh5["pres"]  <- pres
        fileh5["sh"]    <- sh
        fileh5["tmp"]   <- tmp
        fileh5["ugrd"]  <- ugrd
        fileh5["vbdsf"] <- vbdsf
        fileh5["vddsf"] <- vddsf
        fileh5["vgrd"]  <- vgrd
        file.is.closed <- h5close(fileh5)
        if(!file.is.closed) stop("Error, file not closed!")
        rm(fileh5)
        memuse <- gc(); if((memuse[2,1]/memuse[2,5])>0.9) stop("C R I T I C A L : Process is taking up too much memory!")

        #---------------------------------------------------------------------------------#
        #     Make the dataset.                                                           #
        #---------------------------------------------------------------------------------#
        # thisfile = file.path(siteroot
        #                     ,paste(info$h5pref,"_",year.label,month.label,".h5",sep=""))
        # cat("     * Saving data to ",basename(thisfile),"...","\n",sep="")
        #
        #
        #
        # h5save(file=thisfile,dlwrf,nbdsf,nddsf,prate,pres,sh,tmp,ugrd,vbdsf,vddsf,vgrd)
        # H5close()
        #---------------------------------------------------------------------------------#
      }else{
        stop (" ---> Could not find all the data for this month!")
      }#end if
    }#end for
    #---------------------------------------------------------------------------------------#


    #---------------------------------------------------------------------------------------#
    #     Make the header for this experiment.                                              #
    #---------------------------------------------------------------------------------------#
    header  = file.path(siteroot,"ED_MET_DRIVER_HEADER")
    info    = c("# See README at the bottom of this file."
                ,"1"
                ,paste(path.header,info$h5pref,"_",sep="")
                ,paste(2,2,1.0,1.0,info$lon,info$lat)
                ,"12"
                ,paste("   'hgt'   'tmp'  'pres'    'sh'  'ugrd'  'vgrd' 'prate'"
                       ," 'dlwrf' 'nbdsf' 'nddsf' 'vbdsf' 'vddsf'",sep="")
                ,paste(c(info$height,rep(info$dtdat,times=11)),collapse=" ")
                ,paste(" 4 1 1 1 1 1 0 1 1 1 1 1",sep="")
                ,paste(" ")
                ,paste("!===========================================================!")
                ,paste("! README                                                    !")
                ,paste("!===========================================================!")
                ,paste("!     The header of the meteorological driver must contain  !")
                ,paste("! the following lines:                                      !")
                ,paste("!                                                           !")
                ,paste("! Line  1 : Banner, it will not be read;                    !")
                ,paste("! Line  2 : Number of file formats, hereafter N;            !")
                ,paste("! Lines 3+: For each of the N formats, add the following    !")
                ,paste("!           lines, going through a-f for the first format,  !")
                ,paste("!           then through a-f for the second format and so   !")
                ,paste("!            on:                                            !")
                ,paste("!    a. Prefixes of the file format;                        !")
                ,paste("!    b. nlon, nlat, deltalon, deltalat, lon0, lat0.  If     !")
                ,paste("!       lon and lat are also variables, only nlon and nlat  !")
                ,paste("!       will be used;                                       !")
                ,paste("!    c. Number of variables contained in this format;       !")
                ,paste("!    d. List of variables for each format (see Table 1);    !")
                ,paste("!    e. Frequency at which vares are updated, or the        !")
                ,paste("!       constant value if the variable type is 4;           !")
                ,paste("!    f. Variable type (see Table 2);                        !")
                ,paste("!                                                           !")
                ,paste("!===========================================================!")
                ,paste("! Table 1. Variable names recognized by ED.                 !")
                ,paste("!===========================================================!")
                ,paste("! -> lon    -  Longitude                        [    deg]   !")
                ,paste("! -> lat    -  Latitude                         [    deg]   !")
                ,paste("! -> hgt    -  Reference height                 [  m AGL]   !")
                ,paste("! -> tmp    -  Air temperature                  [      K]   !")
                ,paste("! -> pres   -  Pressure                         [     Pa]   !")
                ,paste("! -> sh     -  Specific humidity                [  kg/kg]   !")
                ,paste("! -> ugrd   -  Zonal wind                       [    m/s]   !")
                ,paste("! -> vgrd   -  Zonal wind                       [    m/s]   !")
                ,paste("! -> prate  -  Precipitation rate               [kg/m2/s]   !")
                ,paste("! -> dlwrf  -  Downward long wave radiation     [   W/m2]   !")
                ,paste("! -> nbdsf  -  Near-IR beam radiation           [   W/m2]   !")
                ,paste("! -> nddsf  -  Near-IR diffuse radiation        [   W/m2]   !")
                ,paste("! -> vbdsf  -  Visible beam radiation           [   W/m2]   !")
                ,paste("! -> vddsf  -  Visible beam radiation           [   W/m2]   !")
                ,paste("!===========================================================!")
                ,paste("!                                                           !")
                ,paste("!===========================================================!")
                ,paste("! Table 2. Variable types recognized by ED.                 !")
                ,paste("!===========================================================!")
                ,paste("!                                                           !")
                ,paste("! 0. Read gridded data - no time interpolation;             !")
                ,paste("! 1. Read gridded data - with time interpolatation;         !")
                ,paste("! 2. Read gridded data that is constant in time.            !")
                ,paste("!    If any of this is lon or lat, then deltalon, deltalat  !")
                ,paste("!    lon0, and lat0 will be ignored;                        !")
                ,paste("! 3. Read one value representing the whole grid, no time    !")
                ,paste("!   interpolation;                                          !")
                ,paste("! 4. Specify a constant for all polygons, constant in time. !")
                ,paste("!    In this case, give the constant value at line 'e'      !")
                ,paste("!    instead of the frequency.                              !")
                ,paste("!===========================================================!")
    )#end c
    #---------------------------------------------------------------------------------------#



    #---------------------------------------------------------------------------------------#
    #     Write the header.                                                                 #
    #---------------------------------------------------------------------------------------#
    write (x=info,file=header,ncolumns=1,append=FALSE,sep=" ")
    #---------------------------------------------------------------------------------------#
    #=======================================================================================#
    #=======================================================================================#

    #rm(datum  )
    #   H5close()

}
