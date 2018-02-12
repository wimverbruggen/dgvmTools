#==========================================================================================#
#==========================================================================================#
#      Function that finds the fraction of the day.                                        #
#------------------------------------------------------------------------------------------#
hms2frac = function(when){
  thishour  = chron::hours    (when)
  thismin   = chron::minutes  (when)
  thissec   = chron::seconds  (when)

  elapsed = thishour/24 + thismin/(24*60) + thissec/(24*60*60)
  return(elapsed)
}#end function
#==========================================================================================#
#==========================================================================================#

#==========================================================================================#
#==========================================================================================#
# Time conversion units                                                                    #
#------------------------------------------------------------------------------------------#
yr.day  <<- 365.2425         # # of days in a year                              [   day/yr]
day.sec <<- 86400.           # # of seconds in a day                            [    s/day]
day.min <<- 1440.            # # of minutes in a day                            [  min/day]
day.hr  <<- 24.              # # of hours in a day                              [   hr/day]
hr.sec  <<- 3600.            # # of seconds in an hour                          [     s/hr]
hr.min  <<- 60.              # # of minutes in an hour                          [   min/hr]
min.sec <<- 60.              # # of seconds in a minute                         [    s/min]
yr.sec  <<- yr.day * day.sec # # of seconds in a year                           [     s/yr]
#==========================================================================================#
#==========================================================================================#

#==========================================================================================#
#==========================================================================================#
#    This function finds the cosine of the zenith angle either for the right instant, or   #
# to the interval between two consecutive times.                                           #
#                                                                                          #
# Input variables:                                                                         #
#    - lon       - longitude of the point.  Mandatory, one value only.                     #
#    - lat       - latitude of the point.  Mandatory, one value only.                      #
#    - when      - time.  Mandatory, one point only or a vector.                           #
#    - ed21      - I shall use ED-2.1 method (TRUE/FALSE).  Default is TRUE                #
#    - zeronight - The cosine of zenith angle shall be set to zero at night.               #
#                  Default is FALSE                                                        #
#    - meanval   - I shall find the mean cosine of the integration time.  The beginning    #
#                  and the end are given by variable imetavg.  Default is FALSE.  In case  #
#                  it is TRUE but "when" has just one point, this flag will be ignored and #
#                  it will be solved as instantaneous.                                     #
#    - imetavg   - Which kind of time average was used?                                    #
#                  1 - averages ending at the reference time;                              #
#                  2 - averages beginning at the reference time;                           #
#                  3 - averages centred at the reference time.                             #
#    - nmean     - Number of intermediate points for the average                           #
#                                                                                          #
# The output is going to be a list with the following values:                              #
#    - cosz      - Cosine of zenith angle                                                  #
#    - zen       - The zenith angle in degrees                                             #
#    - height    - The sun height in degrees                                               #
#    - declin    - Declination in degrees                                                  #
#    - day       - Daytime (Sun above horizon)                                             #
#    - night     - Night time (Sun below 6 degrees below the horizon)                      #
#    (N.B. When both day and night are false, we consider it twilight.                     #
#------------------------------------------------------------------------------------------#
calc.zenith <- function (lon,lat,when,ed21=TRUE,zeronight=FALSE,meanval=FALSE,imetavg=1,nmean=120,...){

  #------------------------------------------------------------------------------------------#
  # General Radiation constants.                                                             #
  #------------------------------------------------------------------------------------------#
  pio180        <<- pi/ 180.
  capri         <<- -23.44 * pio180 # Tropic of Capricornium latitude             [      rad]
  shsummer      <<- -10             # Day of year of S.Hemisphere summer solstice [      day]
  solar         <<- 1.3533e3        # Solar constant                              [     W/m2]
  prefsea       <<- 101325.         # Reference sea level pressure                [       Pa]
  twothirds     <<- 2./3.          # 2/3                                          [      ---]
  cosz.min      <<- 0.03            # Minimum cosine of zenith angle
  cosz.highsun  <<- cos(84*pi/180)  # Zenith angle to not be called sunrise or sunset
  cosz.twilight <<- cos(96*pi/180)  # Cosine of the end of civil twilight
  fvis.beam.def <<- 0.43
  fnir.beam.def <<- 1.0 - fvis.beam.def
  fvis.diff.def <<- 0.52
  fnir.diff.def <<- 1.0 - fvis.diff.def
  #------------------------------------------------------------------------------------------#

  #------ Constants. ---------------------------------------------------------------------#
  dcoeff   = c( 0.006918, -0.399912,  0.070257, -0.006758,  0.000907, -0.002697,  0.001480)
  #---------------------------------------------------------------------------------------#


  #------ Find the number of elements. ---------------------------------------------------#
  ntimes  = length(when)
  if ((! meanval) | ntimes == 1) nmean = 1
  #---------------------------------------------------------------------------------------#


  #------ Make matrix of times to make the results time averages if needed be. -----------#
  if (nmean > 1){
    #------------------------------------------------------------------------------------#
    #     The minimum difference is safer than the mean in case the time series has      #
    # gaps.                                                                              #
    #------------------------------------------------------------------------------------#
    dwhen = diff(as.numeric(when))
    sel   = is.finite(dwhen)
    dwhen = dwhen[sel]
    dwhen = min(dwhen[dwhen > 0])
    #------------------------------------------------------------------------------------#


    #------------------------------------------------------------------------------------#
    #    Decide the beginning and ending times depending on imetavg.                     #
    #------------------------------------------------------------------------------------#
    if (imetavg == 1){
      #----- Averages ending at the reference time. ------------------------------------#
      na = 1 - nmean
      nz = 0
    }else if (imetavg == 2){
      #----- Averages starting at the reference time. ----------------------------------#
      na = 0
      nz = nmean - 1
    }else if (imetavg == 3){
      #---------------------------------------------------------------------------------#
      #     Averages centered at the reference time. The initial and ending times na    #
      # and nz will be slightly different depending on whether the number of mean       #
      # points is odd or even.                                                          #
      #---------------------------------------------------------------------------------#
      nz = floor(nmean/2) + 0.5 * ((nmean %% 2) - 1.0)
      na = - nz
    }else{
      print(paste(" ---> In function ed.zen: imetavg =",imetavg,".",sep=""))
      stop ("Invalid imetavg, it must be 1, 2, or 3!")
    }#end if
    #------------------------------------------------------------------------------------#



    #----- Averages ending at the reference time. ---------------------------------------#
    dtidx = seq(from=na,to=nz,by=1) / (nz - na + 1)
    WHEN  = chron( matrix(as.numeric(when),ncol=nmean,nrow=ntimes)
                   + matrix(dtidx,ncol=nmean,nrow=ntimes,byrow=TRUE) * dwhen)
    #------------------------------------------------------------------------------------#
  }else{
    #----- Single time, use only the instantaneous value. -------------------------------#
    WHEN  = matrix(as.numeric(when),ncol=nmean,nrow=ntimes)
  }#end if
  empty = as.numeric(WHEN) * NA
  #---------------------------------------------------------------------------------------#


  #------ Find the day of year, list of leap year times, and sun hour. -------------------#
  doy     = matrix(lubridate::yday(when)           ,ncol=nmean,nrow=ntimes)
  leap    = matrix(chron::leap.year(when)           ,ncol=nmean,nrow=ntimes)
  fracday = matrix(hms2frac (as.vector(WHEN)),ncol=nmean,nrow=ntimes)
  sunhr   = (fracday * 24 + lon / 15. + 24) %% 24
  #---------------------------------------------------------------------------------------#



  #------ Find the hour angle and its cosine. --------------------------------------------#
  hrangle = 15 * (sunhr - 12) * pi/180
  chra    = cos(hrangle)
  #---------------------------------------------------------------------------------------#

  #------ Find the declination
  if (ed21){
    doyfun = empty
    doyfun[!leap] = 2 * pi * (doy[!leap] - shsummer) / 365.
    doyfun[ leap] = 2 * pi * (doy[ leap] - shsummer) / 366.

    declin = capri * cos(doyfun)
  }else{
    doyfun = empty
    doyfun[!leap] = 2 * pi * (doy[!leap] - 1) / 365.
    doyfun[ leap] = 2 * pi * (doy[ leap] - 1) / 366.

    declin = ( dcoeff[1]
               + dcoeff[2] * cos(1.*doyfun) + dcoeff[3] * sin(1.*doyfun)
               + dcoeff[4] * cos(2.*doyfun) + dcoeff[5] * sin(2.*doyfun)
               + dcoeff[6] * cos(3.*doyfun) + dcoeff[7] * sin(3.*doyfun) )
  }#end if
  #---------------------------------------------------------------------------------------#

  #------ Find the cosine and sine of latitude and declination. --------------------------#
  clat = cos(pi/180*lat)
  slat = sin(pi/180*lat)
  cdec = matrix(cos(declin),ncol=nmean,nrow=ntimes)
  sdec = matrix(sin(declin),ncol=nmean,nrow=ntimes)

  #------ Find the cosine of the zenith angle, the zenith angle, and day/night flag. -----#
  cosz   = rowMeans(slat * sdec + clat * cdec * chra,...)
  zen    = acos(cosz) / pi/180
  hgt    = 90. - zen
  declin = declin / pi/180
  night  = cosz <  cosz.twilight
  day    = cosz >= cosz.min

  if (zeronight){
    cosz[night] =  0.
    hgt [night] =  0.
    zen [night] = 90.
  }#end if

  ans = list(cosz=cosz,zen=zen,hgt=hgt,declin=declin,day=day,night=night)
  return(ans)
}
#==========================================================================================#
#==========================================================================================#

#==========================================================================================#
#==========================================================================================#
#      This subroutine computes the split between direct and diffuse radiation, and        #
# between visible and near-infrared radiation using the method suggested by:               #
#                                                                                          #
# Weiss, A., J. M. Norman, 1985: Partitioning solar radiation into direct and diffuse,     #
#     visible and near-infrared components.  Agric. For. Meteorol., 34, 205-213. (WN85)    #
#                                                                                          #
# Input variables:                                                                         #
#                                                                                          #
#    * rad.in   - The incoming radiation at surface, in W/m2.  This can be either PAR,     #
#                 NIR, or the total shortwave radiation, but it must be in W/m2 in any of  #
#                 the cases.                                                               #
#    * atm.prss - The atmospheric pressure at the surface, in Pa.  An actual measurement   #
#                 is better, but if you don't have one just use some typical value given   #
#                 the place altitude (higher elevation sites get more radiation).          #
#    * cosz     - The cosine of zenith angle.  This can be estimated using function ed.zen #
#                 in file zen.r
#    * rad.type - The type of radiation provided in rad.in.  Default is total shortwave    #
#                 radiation, but the function also accepts PAR or NIR.  The value is case  #
#                 insensitive and only the first letter is checked.  "p" means PAR, "n"    #
#                 means NIR, and any other letter will be assumed shortwave.               #
#------------------------------------------------------------------------------------------#
rshort.bdown = function(rad.in,atm.prss,cosz,rad.type="rshort"){
  #---------------------------------------------------------------------------------------#
  #    Local constants.                                                                   #
  #---------------------------------------------------------------------------------------#
  #----- Extinction coefficient. (equations 1 and 4 of WN85) -----------------------------#
  par.beam.expext  = -0.185
  nir.beam.expext  = -0.060
  #----- This is the typical conversion of diffuse radiation in sunny days. --------------#
  par2diff.sun = 0.400
  nir2diff.sun = 0.600
  #----- Coefficients for various equations in WN85. -------------------------------------#
  wn85.06 = c( -1.1950, 0.4459, -0.0345 )
  wn85.11 = c(    0.90, 0.70  )
  wn85.12 = c(    0.88, 0.68  )
  #---------------------------------------------------------------------------------------#



  #---------------------------------------------------------------------------------------#
  #     Make rad type case insensitive, and retain only the first letter.                 #
  #---------------------------------------------------------------------------------------#
  rty = substring(tolower(rad.type),1,1)
  #---------------------------------------------------------------------------------------#


  #------ Initialise the radiation with NAs. ---------------------------------------------#
  par.beam    = NA * rad.in
  nir.beam    = NA * rad.in
  par.diff    = NA * rad.in
  nir.diff    = NA * rad.in
  par.full    = NA * rad.in
  nir.full    = NA * rad.in
  rshort.beam = NA * rad.in
  rshort.diff = NA * rad.in
  rshort.full = NA * rad.in
  par.max     = NA * rad.in
  nir.max     = NA * rad.in
  rshort.max  = NA * rad.in


  #------ Make day and night flags. ------------------------------------------------------#
  ntimes = length(cosz)
  night  = cosz <= cosz.min
  day    = ! night
  #---------------------------------------------------------------------------------------#



  #---------------------------------------------------------------------------------------#
  #     First thing to check is whether this is daytime or "night-time".  If the zenith   #
  # angle is too close to horizon, we assume it's dawn/dusk and all radiation goes to     #
  # diffuse.                                                                              #
  #---------------------------------------------------------------------------------------#
  par.beam    [night] = 0.0
  nir.beam    [night] = 0.0
  if (rty == "p"){
    par.diff    [night] = rad.in[night]
    nir.diff    [night] = fnir.diff.def * rad.in[night] / fvis.diff.def
  }else if(rty == "n"){
    par.diff    [night] = fvis.diff.def * rad.in[night] / fnir.diff.def
    nir.diff    [night] = rad.in[night]
  }else{
    par.diff    [night] = fvis.diff.def * rad.in[night]
    nir.diff    [night] = fnir.diff.def * rad.in[night]
  }#end if
  par.full    [night] = par.beam   [night] + par.diff   [night]
  nir.full    [night] = nir.beam   [night] + nir.diff   [night]
  rshort.beam [night] = par.beam   [night] + nir.beam   [night]
  rshort.diff [night] = par.diff   [night] + nir.diff   [night]
  rshort.full [night] = rshort.beam[night] + rshort.diff[night]
  par.max     [night] = 0.0
  nir.max     [night] = 0.0
  rshort.max  [night] = 0.0
  #---------------------------------------------------------------------------------------#



  #----- Save 1/cos(zen), which is the secant.  We will use this several times. ----------#
  secz      = 1. / cosz[day]
  log10secz = log10(secz)
  #---------------------------------------------------------------------------------------#


  #----- Total radiation at the top [  W/m2], using ED defaults. -------------------------#
  par.beam.top = fvis.beam.def * solar
  nir.beam.top = fnir.beam.def * solar
  #---------------------------------------------------------------------------------------#

  #---------------------------------------------------------------------------------------#
  #    Find the potential PAR components (beam, diffuse, total), using equations 1, 3,    #
  # and 9 of WN85.                                                                        #
  #---------------------------------------------------------------------------------------#
  par.beam.pot = ( par.beam.top
                   * exp ( par.beam.expext * (atm.prss[day] / prefsea) * secz) * cosz[day])
  par.diff.pot = par2diff.sun * (par.beam.top - par.beam.pot) * cosz[day]
  par.full.pot = par.beam.pot + par.diff.pot
  #------------------------------------------------------------------------------------#



  #---------------------------------------------------------------------------------------#
  #     Find the NIR absorption of 10 mm of precipitable water, using WN85 equation 6.    #
  #---------------------------------------------------------------------------------------#
  w10 = solar * 10 ** ((wn85.06[1]) + log10secz * (wn85.06[2] + wn85.06[3] * log10secz))
  #---------------------------------------------------------------------------------------#



  #---------------------------------------------------------------------------------------#
  #     Find the potential direct and diffuse near-infrared radiation, using equations    #
  # 4, 5, and 10 of WN85.                                                                 #
  #---------------------------------------------------------------------------------------#
  nir.beam.pot = ( ( nir.beam.top
                     * exp ( nir.beam.expext * (atm.prss[day] / prefsea) * secz) - w10 )
                   * cosz[day] )
  nir.diff.pot = nir2diff.sun * ( nir.beam.top - nir.beam.pot - w10 ) * cosz[day]
  nir.full.pot = nir.beam.pot + nir.diff.pot
  #---------------------------------------------------------------------------------------#



  #---------------------------------------------------------------------------------------#
  #     Total maximum radiation.                                                          #
  #---------------------------------------------------------------------------------------#
  par.max   [day] = par.full.pot
  nir.max   [day] = nir.full.pot
  rshort.max[day] = par.full.pot + nir.full.pot
  #---------------------------------------------------------------------------------------#



  #---------------------------------------------------------------------------------------#
  #     Find the actual total for PAR and NIR, using equations 7 and 8.                   #
  #---------------------------------------------------------------------------------------#
  if (rty == "p"){
    ratio      = rad.in[day] / par.full.pot
  }else if (rty == "n"){
    ratio      = rad.in[day] / nir.full.pot
  }else{
    ratio      = rad.in[day] / (par.full.pot + nir.full.pot)
  }#end if
  par.full[day] = ratio * par.full.pot
  nir.full[day] = ratio * nir.full.pot
  #---------------------------------------------------------------------------------------#



  #---------------------------------------------------------------------------------------#
  #     Find the fraction of PAR and NIR that stays as beam, using equations 11 and 12    #
  # of WN85.                                                                              #
  #---------------------------------------------------------------------------------------#
  #----- Make sure that the ratio is bounded. --------------------------------------------#
  aux.par  = pmin(wn85.11[1],pmax(0.,ratio))
  aux.nir  = pmin(wn85.12[1],pmax(0.,ratio))

  fvis.beam.act = ( par.beam.pot
                    * (1. - ((wn85.11[1] - aux.par)/wn85.11[2]) ^ twothirds)
                    / par.full.pot )
  fvis.beam.act = pmin(1.,pmax(0.,fvis.beam.act))

  fnir.beam.act = ( nir.beam.pot
                    * (1. - ((wn85.12[1] - aux.nir)/wn85.12[2]) ^ twothirds)
                    / nir.full.pot )
  fnir.beam.act = pmin(1.,pmax(0.,fvis.beam.act))

  fvis.diff.act = 1. - fvis.beam.act
  fnir.diff.act = 1. - fnir.beam.act
  #---------------------------------------------------------------------------------------#



  #---------------------------------------------------------------------------------------#
  #     Find the radiation components.                                                    #
  #---------------------------------------------------------------------------------------#
  par.beam    [day] = fvis.beam.act * par.full[day]
  par.diff    [day] = fvis.diff.act * par.full[day]
  nir.beam    [day] = fnir.beam.act * nir.full[day]
  nir.diff    [day] = fnir.diff.act * nir.full[day]
  rshort.beam [day] = par.beam   [day] + nir.beam   [day]
  rshort.diff [day] = par.diff   [day] + nir.diff   [day]
  rshort.full [day] = rshort.beam[day] + rshort.diff[day]
  #---------------------------------------------------------------------------------------#
  rshort.bdown = list( par.beam    = par.beam
                       , par.diff    = par.diff
                       , par.full    = par.full
                       , nir.beam    = nir.beam
                       , nir.diff    = nir.diff
                       , nir.full    = nir.full
                       , rshort.beam = rshort.beam
                       , rshort.diff = rshort.diff
                       , rshort.full = rshort.full
                       , par.max     = par.max
                       , nir.max     = nir.max
                       , rshort.max  = rshort.max
  )#end list
  return(rshort.bdown)
}#end function rshort.bdown
#==========================================================================================#
#==========================================================================================#









#==========================================================================================#
#==========================================================================================#
#      Function that determines whether the year is leap or not.                           #
#------------------------------------------------------------------------------------------#
is.leap = function(when){

  if (is.chron(when)){
    year = numyears(when)
  }else{
    year = when
  }#end if

  leaptf = year %% 400 == 0 | (year %% 4 == 0 & year %% 100 != 0)
  return(leaptf)
} #end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      Function that determines the number of days in a given month.                       #
#------------------------------------------------------------------------------------------#
daymax = function(month,year){
  mmm  = c(31,28,31,30,31,30,31,31,30,31,30,31)
  mday = mmm[month]

  addone       = month == 2 & is.leap(year)
  mday[addone] = mday[addone] + 1

  return(mday)
} #end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      Function that determines the number of the month given the character name.          #
#------------------------------------------------------------------------------------------#
mmm2mon = function(mmm,lang="English"){
  lang = tolower(lang)
  if (lang == "english"){
    m3l  = c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
  }else if(lang == "portuguese"){
    m3l  = c("jan","fev","mar","abr","mai","jun","jul","ago","set","out","nov","dez")
  }#end if

  mmmloc = tolower(substring(as.character(mmm),1,3))
  monout = match(mmmloc,m3l)
  return(monout)
} #end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      Function that determines 3-letter name of the month given their number.             #
#------------------------------------------------------------------------------------------#
mon2mmm = function(mon,lang="English"){
  lang = tolower(lang)
  if (lang == "english"){
    m3l  = c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
  }else if(lang == "portuguese"){
    m3l  = c("jan","fev","mar","abr","mai","jun","jul","ago","set","out","nov","dez")
  }#end if

  monout = m3l[mon]
  return(monout)
} #end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      Function that converts a chron object to numeric years.                             #
#------------------------------------------------------------------------------------------#
numyears = function(when){
  yrs    = years(when)
  lyrs   = levels(yrs)
  yrout  = as.numeric(lyrs[match(yrs,lyrs)])
  return(yrout)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      Function that converts a chron object to numeric months.                            #
#------------------------------------------------------------------------------------------#
nummonths = function(when){
  mos    = months(when)
  lmos   = levels(mos)
  moout  = match(mos,lmos)
  return(moout)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      Function that converts a chron object to numeric days.                              #
#------------------------------------------------------------------------------------------#
numdays = function(when){
  dys    = days(when)
  ldys   = levels(dys)
  dyout  = match(dys,ldys)
  return(dyout)
} #end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      Function that returns the dates as characters.                                      #
#------------------------------------------------------------------------------------------#
chardates = function(when){
  mymonth = substring(100   + nummonths(when),2,3)
  myday   = substring(100   + numdays  (when),2,3)
  myyear  = substring(10000 + numyears (when),2,5)
  mydate  = paste(mymonth,myday,myyear,sep="/")
  return(mydate)
} #end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      Function that returns the dates as characters.                                      #
#------------------------------------------------------------------------------------------#
label.dates = function(when,add.hours=TRUE){
  mymonth = substring(100   + nummonths(when),2,3)
  myday   = substring(100   + numdays  (when),2,3)
  myyear  = substring(10000 + numyears (when),2,5)
  mydate  = paste(myyear,mymonth,myday,sep="-")

  if (add.hours){
    mytime  = paste(substring(100 + hours  (when),2,3)
                    ,substring(100 + minutes(when),2,3)
                    ,substring(100 + seconds(when),2,3)
                    ,sep="")
    mylabel = paste(mydate,mytime,sep="-")
  }else{
    mylabel = mydate
  }#end if

  return(mylabel)
} #end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      Function that returns the times as characters.                                      #
#------------------------------------------------------------------------------------------#
chartimes = function(when){
  myhour = substring(100 + hours  (when),2,3)
  myminu = substring(100 + minutes(when),2,3)
  myseco = substring(100 + seconds(when),2,3)
  mytime = paste(myhour,myminu,myseco,sep=":")
  return(mytime)
} #end function
#==========================================================================================#
#==========================================================================================#



#==========================================================================================#
#==========================================================================================#
#      Function that finds the numeric version of the days.                                #
#------------------------------------------------------------------------------------------#
dayofyear = function(when){
  offdays   = c(0, 31,59,90,120,151,181,212,243,273,304,334,365)

  thisday   = numdays  (when)
  thismonth = nummonths(when)
  thisyear  = numyears (when)
  thisfrac  = hms2frac (when)

  addone    = as.integer(thismonth > 2 & is.leap(when))

  doy =  thisday + offdays[thismonth] + addone + thisfrac
  return(doy)
} #end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     This function appends several time-related variables for a given data frame.         #
#------------------------------------------------------------------------------------------#
alltimes = function(datin,lon,lat,ed21=TRUE,zeronight=FALSE,meanval=FALSE,imetavg=1
                    ,nmean=120,...){
  #------ Copy the input data frame, and call the other functions. -----------------------#
  datout = datin
  datout$year       = numyears (datout$when)
  datout$month      = nummonths(datout$when)
  datout$day        = numdays  (datout$when)
  datout$hour       = hours    (datout$when)
  datout$minu       = minutes  (datout$when)
  datout$today      = dates    (datout$when)
  datout$tomonth    = chron(paste(datout$month,1,datout$year,sep="/"))
  datout$doy        = dayofyear(datout$when)
  zenith            = ed.zen   (when=datout$when,lon=lon,lat=lat,ed21=ed21
                                ,zeronight=zeronight,meanval=meanval,imetavg=imetavg
                                ,nmean=nmean,...)
  datout$cosz       =   zenith$cosz
  datout$sunhgt     =   zenith$hgt
  datout$nighttime  =   zenith$night
  datout$daytime    =   zenith$day
  datout$twilight   = (! zenith$night) & (! zenith$day)
  datout$notdaytime = ! zenith$day

  return(datout)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      List of trimestral seasons.                                                         #
#------------------------------------------------------------------------------------------#
season <<- function(when,add.year=FALSE){


  #----- Get the year and month. ---------------------------------------------------------#
  year = numyears (when)
  mon  = nummonths(when)
  #---------------------------------------------------------------------------------------#



  #----- We don't give summer/winter, instead we make generic season names. --------------#
  sidx = c( 4, 4, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4)
  #---------------------------------------------------------------------------------------#



  #---------------------------------------------------------------------------------------#
  #     Assign the season names depending on the month and year.                          #
  #---------------------------------------------------------------------------------------#
  if (add.year){
    #----- Add year before the season. --------------------------------------------------#
    seasout      = paste(year       ,substring(100+sidx[mon     ],2,3),sep="")
    #------------------------------------------------------------------------------------#


    #----- December, January, and February have two years. ------------------------------#
    mm1          = mon %in%  c(1,2)
    seasout[mm1] = paste(year[mm1]-1,substring(100+sidx[mon[mm1]],2,3),sep="")
    #------------------------------------------------------------------------------------#
  }else{
    #----- No year to be tagged. --------------------------------------------------------#
    seasout = sidx[mon]
    #------------------------------------------------------------------------------------#
  }#end if
  #---------------------------------------------------------------------------------------#



  #----- Return variable. ----------------------------------------------------------------#
  return(seasout)
  #---------------------------------------------------------------------------------------#
}#end for
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#     This function creates a pretty time scale.  It is loosely based on pretty, but here  #
# we make extensive use of the chron functions, and define the suitable scales in a        #
# different way as time has a non-decimal scale.                                           #
#  The result is a list containing the levels, and nice labels for the plots.              #
#------------------------------------------------------------------------------------------#
pretty.time = function(when,n=10,...){

  #----- Find the 1st and last time. -----------------------------------------------------#
  whena = min(when,na.rm=TRUE)
  whenz = max(when,na.rm=TRUE)

  #----- Table for accepted bases for months, hours, minutes, and seconds. ---------------#
  base.months = c(1,2,3,4,6)
  base.days   = c(1,2,4,7,15)
  base.hours  = c(1,2,3,4,6,12)
  base.minsec = c(1,2,5,10,15,20,30)

  #----- Convert time to seconds. --------------------------------------------------------#
  when.sec = as.numeric(when) * day.sec

  #---------------------------------------------------------------------------------------#
  #    Find a first guess of the step size, so we decide whether to use years, months,    #
  # days, hours, minutes, or seconds.                                                     #
  #---------------------------------------------------------------------------------------#
  wstep.1st  = mean(diff(pretty(when.sec,n)))

  if (wstep.1st == 0){
    myunit=NA
    #---- Whatever, use just what comes out from the regular pretty. --------------------#
    vlevels = chron(pretty(when,n))
    vlabels = as.character(vlevels)
    padj    = rep(0,times=length(vlabels))
  }else if(wstep.1st / yr.sec > 0.8){
    myunit="years"
    #------------------------------------------------------------------------------------#
    #     Years are the best scale for the plots.                                        #
    #------------------------------------------------------------------------------------#
    yrrange = numyears(when)
    vlevels = pretty(yrrange,n)
    vlevels = dates(x=paste(1,1,vlevels,sep="/"))

    vlabels = paste(months(vlevels),years(vlevels),sep="-")
    padj    = rep(0,times=length(vlabels))
  }else if(wstep.1st / (30. * day.sec) > 0.8){
    myunit="months"
    #------------------------------------------------------------------------------------#
    #     Months are the best scale for the plots.                                       #
    #------------------------------------------------------------------------------------#
    #----- Find the time step that is the closest to the base. --------------------------#
    wstep     = wstep.1st / (30. * day.sec)
    whichbase = base.months[which.min(abs(wstep-base.months))]

    #----- Find the list of years to plot. ----------------------------------------------#
    allyears = numyears(when)
    yeara    = min(allyears,na.rm=TRUE)
    yearz    = max(allyears,na.rm=TRUE)+1
    vlevels  = seq.dates(from = paste(1,1,yeara,sep="/")
                         ,to   = paste(1,1,yearz,sep="/")
                         ,by   = "months")
    mon1st   = nummonths(vlevels)
    monlevs  = seq(from=1,to=12,by=whichbase)

    #----- Find the limits that will keep the labels not too far from the data. ---------#
    wlaba    = dates(paste(nummonths(whena),1,numyears(whena),sep="/"))
    monz     = nummonths(whenz) %% 12 + 1
    yearz    = numyears(whenz) + as.integer(monz == 1)
    wlabz    = dates(paste(monz,1,yearz,sep="/"))
    sel      = ( mon1st %in% monlevs
                 & vlevels >= min(wlaba,na.rm=TRUE)
                 & vlevels <= max(wlabz,na.rm=TRUE) )
    vlevels = dates(vlevels[sel],out.format="m/d/year")
    vlabels = paste(months(vlevels),years(vlevels),sep="-")
    padj    = rep(0,times=length(vlabels))

  }else if(wstep.1st / day.sec > 0.8){
    myunit="days"
    #------------------------------------------------------------------------------------#
    #     Days are the best scale for the plots, but we keep them tethered to months,    #
    # even if the grid becomes slightly irregular.                                       #
    #------------------------------------------------------------------------------------#
    #----- Find the time step that is the closest to the base. --------------------------#
    wstep     = wstep.1st / day.sec
    whichbase = base.days[which.min(abs(wstep-base.days))]
    #----- Find the list of years to plot. ----------------------------------------------#
    allyears = numyears(when)
    yeara    = min(allyears,na.rm=TRUE)
    yearz    = max(allyears,na.rm=TRUE)+1
    #------------------------------------------------------------------------------------#
    #     Impose the list of months to be from January to December, we will trim the     #
    # numbers later.                                                                     #
    #------------------------------------------------------------------------------------#
    montha   = 1
    monthz   = 12
    #------------------------------------------------------------------------------------#
    #     Impose the list of months to be from January to December, we will trim the     #
    # numbers later.                                                                     #
    #------------------------------------------------------------------------------------#
    daylevs=seq(from=1,to=31-whichbase+1,by=whichbase)
    #----- First guess for the levels. --------------------------------------------------#
    vlevels  = seq.dates(from = paste(1,1,yeara,sep="/")
                         ,to   = paste(1,1,yearz,sep="/")
                         ,by   = "days")
    day1st   = numdays(vlevels)

    #----- Find the limits that will keep the labels not too far from the data. ---------#
    wlaba    = dates(whena)
    dayz     = numdays(whenz) %% daymax(nummonths(whenz),numyears(whenz)) + 1
    monz     = 1 + (nummonths(whenz) - 1 + as.integer(dayz==1)) %% 12
    yearz    = numyears(whenz) + as.integer(monz == 1)
    wlabz    = dates(paste(monz,dayz,yearz,sep="/"))
    sel      = ( day1st %in% daylevs
                 & vlevels >= min(wlaba,na.rm=TRUE)
                 & vlevels <= max(wlabz,na.rm=TRUE) )
    vlevels = dates(vlevels[sel],out.format="m/d/y")
    vlabels    = paste(months(vlevels),days(vlevels),sep="/")

    padj    = rep(0,times=length(vlabels))

    sel        = vlevels == vlevels[1] | (numdays(vlevels) == 1 & nummonths(vlevels) == 1)
    vlabels[sel] = paste(months(vlevels[sel]),"/",days(vlevels[sel]),"\n"
                         ,years(vlevels[sel]),sep="")
    padj[sel]    = 0.5
  }else if(wstep.1st / hr.sec > 0.8){
    myunit="hours"
    #------------------------------------------------------------------------------------#
    #     Hours are the best scale for the plots.                                        #
    #------------------------------------------------------------------------------------#
    #----- Find the time step that is the closest to the base. --------------------------#
    wstep     = wstep.1st / hr.sec
    whichbase = base.hours[which.min(abs(wstep-base.hours))]
    #----- Find the list of days to plot. -----------------------------------------------#
    when1st  = dates(min(when  ,na.rm=TRUE))
    whenlast = dates(max(when+1,na.rm=TRUE))
    mydates  = seq.dates(from=when1st,to=whenlast,by="days")
    mytimes  = times(seq(from=0,to=day.sec-1,by=whichbase*hr.sec)) / day.sec
    ndates   = length(mydates)
    ntimes   = length(mytimes)
    #----- First guess for the levels. --------------------------------------------------#
    vlevels  = chron(dates=rep(x=mydates,each=ntimes),times=rep(x=mytimes,times=ndates))
    wlaba    = chron(dates=paste(nummonths(whena),numdays(whena),numyears(whena),sep="/"),
                     times=paste(hours(whena),0,0,sep=":"))
    hourz    = (hours(whenz) + 1) %% 24
    d2831    = daymax(nummonths(whenz),numyears(whenz))
    dayz     = (numdays(whenz) - 1 + as.integer(hourz == 0)) %% d2831 + 1
    monz     = (nummonths(whenz) - 1 + as.integer(dayz == 1)) %% 12 + 1
    yearz    = numyears(whenz) + as.integer(monz == 1)
    wlabz    = chron(dates=paste(monz,dayz,yearz,sep="/"),times=paste(hourz,0,0,sep=":"))
    sel      = ( vlevels >= min(wlaba,na.rm=TRUE)
                 & vlevels <= max(wlabz,na.rm=TRUE) )

    #------------------------------------------------------------------------------------#
    #     Make the labels, and put day and month information only on the first time of   #
    # the day, and all information in the second time, and the first time of the year.   #
    #------------------------------------------------------------------------------------#
    vlabels      = paste(substring(100+hours(vlevels),2,3)
                         ,substring(100+minutes(vlevels),2,3),sep=":")
    padj         = rep(0,times=length(vlabels))
    #----- First time of the day. -------------------------------------------------------#
    sel          = hours(vlevels) == 0
    vlabels[sel] = paste(substring(100+hours(vlevels[sel]),2,3),":"
                         ,substring(100+minutes(vlevels[sel]),2,3),"\n"
                         ,months(vlevels[sel]),"-",days(vlevels[sel]),sep="")
    padj[sel]    = 0.5
    #----- First time of the year. ------------------------------------------------------#
    sel          = ( vlevels == vlevels[1]
                     |  ( nummonths(vlevels) == 1 & numdays(vlevels) == 1
                          & hours(vlevels) == 0 ))
    vlabels[sel] = paste(substring(100+hours(vlevels[sel]),2,3),":"
                         ,substring(100+minutes(vlevels[sel]),2,3),"\n"
                         ,months(vlevels[sel]),"-",days(vlevels[sel]),"\n"
                         ,years(vlevels[sel]),sep="")
    padj[sel]    = 0.5
    #------------------------------------------------------------------------------------#


  }else if(wstep.1st / min.sec > 0.8){
    myunit="minutes"
    #------------------------------------------------------------------------------------#
    #     Minutes are the best scale for the plots.                                      #
    #------------------------------------------------------------------------------------#
    #----- Find the time step that is the closest to the base. --------------------------#
    wstep     = wstep.1st / min.sec
    whichbase = base.minsec[which.min(abs(wstep-base.minsec))]
    #----- Find the list of days to plot. -----------------------------------------------#
    when1st  = dates(min(when  ,na.rm=TRUE))
    whenlast = dates(max(when+1,na.rm=TRUE))
    mydates  = seq.dates(from=when1st,to=whenlast,by="days")
    mytimes  = times(seq(from=0,to=day.sec-1,by=whichbase*min.sec)) / day.sec
    ndates   = length(mydates)
    ntimes   = length(mytimes)
    #----- First guess for the levels. --------------------------------------------------#
    vlevels  = chron(dates=rep(x=mydates,each=ntimes),times=rep(x=mytimes,times=ndates))

    wlaba    = chron(dates=paste(nummonths(whena),numdays(whena),numyears(whena),sep="/"),
                     times=paste(hours(whena),minutes(whena),0,sep=":"))
    minz     = (minutes(whenz) + 1) %% 60
    hourz    = (hours(whenz) + as.integer(minz == 0)) %% 24
    d2831    = daymax(nummonths(whenz),numyears(whenz))
    dayz     = (numdays(whenz) - 1 + as.integer(hourz == 0)) %% d2831 + 1
    monz     = (nummonths(whenz) - 1 + as.integer(dayz == 1)) %% 12 + 1
    yearz    = numyears(whenz) + as.integer(monz == 1)
    wlabz    = chron(dates=paste(monz,dayz,yearz,sep="/")
                     ,times=paste(hourz,minz,0,sep=":"))
    sel      = ( vlevels >= min(wlaba,na.rm=TRUE)
                 & vlevels <= max(wlabz,na.rm=TRUE) )

    #------------------------------------------------------------------------------------#
    #     Make the labels, and put day and month information only on the first time of   #
    # the day, and all information in the second time, and the first time of the year.   #
    #------------------------------------------------------------------------------------#
    vlabels      = paste(substring(100+hours(vlevels),2,3)
                         ,substring(100+minutes(vlevels),2,3),sep=":")
    padj         = rep(0,times=length(vlabels))
    #----- First time of the day. -------------------------------------------------------#
    sel          = hours(vlevels) == 0 & minutes(vlevels) == 0
    vlabels[sel] = paste(substring(100+hours(vlevels[sel]),2,3),":"
                         ,substring(100+minutes(vlevels[sel]),2,3),"\n"
                         ,months(vlevels[sel]),"-",days(vlevels[sel]),sep="")
    padj[sel]    = 0.5
    #----- First time of the year. ------------------------------------------------------#
    sel          = ( vlevels == vlevels[1]
                     |  ( nummonths(vlevels) == 1 & numdays(vlevels) == 1
                          & hours(vlevels) == 0 & minutes(vlevels) == 0))
    vlabels[sel] = paste(substring(100+hours(vlevels[sel]),2,3),":"
                         ,substring(100+minutes(vlevels[sel]),2,3),"\n"
                         ,months(vlevels[sel]),"-",days(vlevels[sel]),"\n"
                         ,years(vlevels[sel]),sep="")
    padj[sel]    = 0.5
    #------------------------------------------------------------------------------------#



  }else{
    myunit="seconds"
    #------------------------------------------------------------------------------------#
    #     Minutes are the best scale for the plots.                                      #
    #------------------------------------------------------------------------------------#
    #----- Find the time step that is the closest to the base. --------------------------#
    wstep     = wstep.1st
    whichbase = base.minsec[which.min(abs(wstep-base.minsec))]
    #----- Find the list of days to plot. -----------------------------------------------#
    when1st  = dates(min(when  ,na.rm=TRUE))
    whenlast = dates(max(when+1,na.rm=TRUE))
    mydates  = seq.dates(from=when1st,to=whenlast,by="days")
    mytimes  = times(seq(from=0,to=day.sec-1,by=whichbase)) / day.sec
    ndates   = length(mydates)
    ntimes   = length(mytimes)
    #----- First guess for the levels. --------------------------------------------------#
    vlevels  = chron(dates=rep(x=mydates,each=ntimes),times=rep(x=mytimes,times=ndates))

    wlaba    = chron(dates=paste(nummonths(whena),numdays(whena),numyears(whena),sep="/"),
                     times=paste(hours(whena),minutes(whena),seconds(whena),sep=":"))
    secz     = (seconds(whenz) + 1) %% 60
    minz     = (minutes(whenz) + as.integer(secz == 0)) %% 60
    hourz    = (hours(whenz) + as.integer(minz == 0)) %% 24
    d2831    = daymax(nummonths(whenz),numyears(whenz))
    dayz     = (numdays(whenz) - 1 + as.integer(hourz == 0)) %% d2831 + 1
    monz     = (nummonths(whenz) - 1 + as.integer(dayz == 1)) %% 12 + 1
    yearz    = numyears(whenz) + as.integer(monz == 1)
    wlabz    = chron(dates=paste(monz,dayz,yearz,sep="/")
                     ,times=paste(hourz,minz,secz,sep=":"))
    sel      = ( vlevels >= min(wlaba,na.rm=TRUE)
                 & vlevels <= max(wlabz,na.rm=TRUE) )

    #------------------------------------------------------------------------------------#
    #     Make the labels, and put day and month information only on the first time of   #
    # the day, and all information in the second time, and the first time of the year.   #
    #------------------------------------------------------------------------------------#
    vlabels      = paste(substring(100+hours(vlevels),2,3)
                         ,substring(100+minutes(vlevels),2,3)
                         ,substring(100+seconds(vlevels),2,3),sep=":")
    padj         = rep(0,times=length(vlabels))
    #----- First time of the day. -------------------------------------------------------#
    sel          = hours(vlevels) == 0 & minutes(vlevels) == 0 & seconds(vlevels) == 0
    vlabels[sel] = paste(substring(100+hours(vlevels[sel]),2,3),":"
                         ,substring(100+minutes(vlevels[sel]),2,3),":"
                         ,substring(100+seconds(vlevels[sel]),2,3),"\n"
                         ,months(vlevels[sel]),"-",days(vlevels[sel]),sep="")
    padj[sel]    = 0.5
    #----- First time of the year. ------------------------------------------------------#
    sel          = ( vlevels == vlevels[1]
                     |  ( nummonths(vlevels) == 1 & numdays(vlevels) == 1
                          & hours(vlevels) == 0  & minutes(vlevels) == 0
                          & seconds(vlevels) == 0))
    vlabels[sel] = paste(substring(100+hours(vlevels[sel]),2,3),":"
                         ,substring(100+minutes(vlevels[sel]),2,3),":"
                         ,substring(100+seconds(vlevels[sel]),2,3),"\n"
                         ,months(vlevels[sel]),"-",days(vlevels[sel]),"\n"
                         ,years(vlevels[sel]),sep="")
    padj[sel]    = 0.5
    #------------------------------------------------------------------------------------#
  }#end if

  vresult=list(levels=vlevels,labels=vlabels,n=length(vlevels),scale=myunit,padj=padj)
  return(vresult)
}#end function
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#      List with the names of months and seasons.                                          #
#------------------------------------------------------------------------------------------#
season.list  <<- c("MAM","JJA","SON","DJF")
#==========================================================================================#
#==========================================================================================#
