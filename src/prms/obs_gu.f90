!***********************************************************************
! Reads and stores observed data from all specified measurement stations
!***********************************************************************
      MODULE PRMS_OBS
      USE PRMS_CONSTANTS, ONLY: MONTHS_PER_YEAR
      IMPLICIT NONE
!   Local Variables
      character(len=*), parameter :: MODDESC = 'Time Series Data'
      character(len=*), parameter :: MODNAME = 'obs'
      character(len=*), parameter :: Version_obs = '2021-08-13'
      INTEGER, SAVE :: Nlakeelev, Nwind, Nhumid, Rain_flag
!   Declared Variables
      INTEGER, SAVE :: Rain_day
      REAL, SAVE, ALLOCATABLE :: Pan_evap(:), Runoff(:), Precip(:)
      REAL, SAVE, ALLOCATABLE :: SubDailyPrecip(:, :), SubDailyObsTimes(:)                      ! PJT - 2018Jan03 - Sub-daily precip inputs
      INTEGER, SAVE :: NumPrecipObs                                                             ! PJT - 2018Jan03 - Sub-daily precip inputs
      REAL, SAVE, ALLOCATABLE :: Humidity(:), Wind_speed(:)
      REAL, SAVE, ALLOCATABLE :: Tmax(:), Tmin(:), Solrad(:), Snowdepth(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Streamflow_cfs(:), Streamflow_cms(:)
      ! Lake Module Variables
      REAL, SAVE, ALLOCATABLE :: Gate_ht(:), Lake_elev(:)
!   Declared Parameters
      INTEGER, SAVE :: Runoff_units, Rain_code(MONTHS_PER_YEAR)
      END MODULE PRMS_OBS

!***********************************************************************
!     main obs routine
!***********************************************************************
      INTEGER FUNCTION obs()
      USE PRMS_CONSTANTS, ONLY: RUN, SETDIMENS, DECL, INIT
      USE PRMS_MODULE, ONLY: Process_flag
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: obsdecl, obsinit, obsrun, obssetdims
!***********************************************************************
      obs = 0

      IF ( Process_flag==RUN ) THEN
        obs = obsrun()
      ELSEIF ( Process_flag==SETDIMENS ) THEN
        obs = obssetdims()
      ELSEIF ( Process_flag==DECL ) THEN
        obs = obsdecl()
      ELSEIF ( Process_flag==INIT ) THEN
        obs = obsinit()
      ENDIF

      END FUNCTION obs

!***********************************************************************
!     obssetdims - declares obs module specific dimensions
!***********************************************************************
      INTEGER FUNCTION obssetdims()
      USE PRMS_CONSTANTS, ONLY: MAXDIM
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: decldim
      EXTERNAL read_error
!***********************************************************************
      obssetdims = 0

      IF ( decldim('nlakeelev', 0, MAXDIM, &
     &     'Maximum number of lake elevations for any rating table data set')/=0 ) CALL read_error(7, 'nlakeelev')
      IF ( decldim('nwind', 0, MAXDIM, 'Number of wind-speed measurement stations')/=0 ) CALL read_error(7, 'nwind')
      IF ( decldim('nhumid', 0, MAXDIM, 'Number of relative humidity measurement stations')/=0 ) CALL read_error(7, 'nhumid')

      END FUNCTION obssetdims

!***********************************************************************
!     obsdecl - makes public variable declarations for the obs module
!   Declared Parameters
!     rain_code
!***********************************************************************
      INTEGER FUNCTION obsdecl()
      USE PRMS_CONSTANTS, ONLY: DOCUMENTATION, ACTIVE, OFF, xyz_dist_module
      USE PRMS_MODULE, ONLY: Model, Nratetbl, Ntemp, Nrain, Nsol, Nobs, Nevap, Nsnow, Precip_flag
      USE PRMS_MODULE, ONLY: NmaxPrecipObs                                                             ! PJT - 2018Jan03 - Sub-daily precip inputs
      USE PRMS_OBS
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: declvar, getdim, declparam
      EXTERNAL :: read_error, print_module
!***********************************************************************
      obsdecl = 0

      CALL print_module(MODDESC, MODNAME, Version_obs)

!   Declared Variables
      IF ( Nobs>0 ) THEN
        ALLOCATE ( Runoff(Nobs) )
        IF ( declvar(MODNAME, 'runoff', 'nobs', Nobs, 'real', &
     &       'Streamflow at each measurement station', &
     &       'runoff_units', Runoff)/=0 ) CALL read_error(8, 'runoff')
        ALLOCATE ( Streamflow_cfs(Nobs) )
        IF ( declvar(MODNAME, 'streamflow_cfs', 'nobs', Nobs, 'double', &
     &       'Streamflow at each measurement station', &
     &       'cfs', Streamflow_cfs)/=0 ) CALL read_error(8, 'streamflow_cfs')
        ALLOCATE ( Streamflow_cms(Nobs) )
        IF ( declvar(MODNAME, 'streamflow_cms', 'nobs', Nobs, 'double', &
     &       'Streamflow at each measurement station', &
     &       'cms', Streamflow_cms)/=0 ) CALL read_error(8, 'streamflow_cms')
        IF ( declparam(MODNAME, 'runoff_units', 'one', 'integer', &
     &       '0', '0', '1', &
     &       'Measured streamflow units', 'Measured streamflow units (0=cfs; 1=cms)', &
     &       'none')/=0 ) CALL read_error(1, 'runoff_units')
      ENDIF

      IF ( Nrain>0 ) THEN
        ALLOCATE ( Precip(Nrain) )
        IF ( declvar(MODNAME, 'precip', 'nrain', Nrain, 'real', &
     &       'Precipitation at each measurement station', &
     &       'precip_units', Precip)/=0 ) CALL read_error(8, 'precip')
      ENDIF

      !If required, create storage space for sub-daily precipitation values                             ! PJT - 2018Jan03 - begin Sub-daily precip inputs      
      IF ( NmaxPrecipObs>0 ) THEN                                                                       ! PJT - 2018Jan03    
        ALLOCATE ( SubDailyPrecip(Nrain, NmaxPrecipObs) )
        IF ( declvar(MODNAME, 'subdailyprecip', 'nrain,nmaxprecipobs', Nrain*NmaxPrecipObs, 'real', & 
     &       'Sub-daily precipitation values, read and populated for' //  &   
     &       ' each timestep sub-daily values are available, requires' // &  
     &       ' Green&Ampt module to be active.', &  
     &       'precip_units', SubDailyPrecip)/=0 ) CALL read_error(3, 'subdailyprecip')   

        ALLOCATE ( SubDailyObsTimes( NmaxPrecipObs) )
        IF ( declvar(MODNAME, 'subdailyobstimes', 'nmaxprecipobs', NmaxPrecipObs, 'real', &  
     &       'Timestep associated with each sub-daily precipitation value, ' //  &
     &       ' in decimal hours, requires Green&Ampt module to be active.', &    
     &       'hours', SubDailyObsTimes)/=0 ) CALL read_error(3, 'subdailyobstimes')

        IF ( declvar(MODNAME, 'numprecipobs', 'one', 1, 'integer', &
     &       'Number of sub-daily precip observations in this timestep', &
     &       'none', NumPrecipObs)/=0 ) CALL read_error(8, 'numprecipobs')                              ! PJT - 2018Jan03 
      ENDIF                                                                                             ! PJT - 2018Jan03 - end Sub-daily precip inputs 
  
      IF ( Ntemp>0 ) THEN
        ALLOCATE ( Tmin(Ntemp) )
        IF ( declvar(MODNAME, 'tmin', 'ntemp', Ntemp, 'real', &
     &       'Minimum air temperature at each measurement station', &
     &       'temp_units', Tmin)/=0 ) CALL read_error(8, 'tmin')
        ALLOCATE ( Tmax(Ntemp) )
        IF ( declvar(MODNAME, 'tmax', 'ntemp', Ntemp, 'real', &
     &       'Maximum air temperature at each measurement station', &
     &       'temp_units', Tmax)/=0 ) CALL read_error(8, 'tmax')
      ENDIF

      IF ( Nsol>0 ) THEN
        ALLOCATE ( Solrad(Nsol) )
        IF ( declvar(MODNAME, 'solrad', 'nsol', Nsol, 'real', &
     &       'Solar radiation at each measurement station', &
     &       'Langleys', Solrad)/=0 ) CALL read_error(8, 'solrad')
      ENDIF

      Nsnow = getdim('nsnow')
      IF ( Nsnow==-1 ) CALL read_error(6, 'nsnow')
      Nhumid = getdim('nhumid')
      IF ( Nhumid==-1 ) CALL read_error(6, 'nhumid')
      Nwind = getdim('nwind')
      IF ( Nwind==-1 ) CALL read_error(6, 'nwind')
      Nlakeelev = getdim('nlakeelev')
      IF ( Nlakeelev==-1 ) CALL read_error(6, 'nlakeelev')

      IF ( Model==DOCUMENTATION ) THEN
        IF ( Nsnow==0 ) Nsnow = 1
        IF ( Nhumid==0 ) Nhumid = 1
        IF ( Nwind==0 ) Nwind = 1
        IF ( Nlakeelev==0 ) Nlakeelev = 1
      ENDIF

      IF ( Nsnow>0 ) THEN
        ALLOCATE ( Snowdepth(Nsnow) )
        IF ( declvar(MODNAME, 'snowdepth', 'nsnow', Nsnow, 'real', &
     &       'Snow depth at each measurement station', &
     &       'inches', Snowdepth)/=0 ) CALL read_error(8, 'snowdepth')
      ENDIF

      IF ( Nevap>0 ) THEN
        ALLOCATE ( Pan_evap(Nevap) )
        IF ( declvar(MODNAME, 'pan_evap', 'nevap', Nevap, 'real', &
     &       'Pan evaporation at each measurement station', &
     &       'inches', Pan_evap)/=0 ) CALL read_error(8, 'pan_evap')
      ENDIF

      IF ( Nhumid>0 ) THEN
        ALLOCATE ( Humidity(Nhumid) )
        IF ( declvar(MODNAME, 'humidity', 'nhumid', Nhumid, 'real', &
     &       'Relative humidity at each measurement station', &
     &       'percentage', Humidity)/=0 ) CALL read_error(8, 'humidity')
      ENDIF

      IF ( Nwind>0 ) THEN
        ALLOCATE ( Wind_speed(Nwind) )
        IF ( declvar(MODNAME, 'wind_speed', 'nwind', Nwind, 'real', &
     &       'Wind speed at each measurement station', &
     &       'meters per second', Wind_speed)/=0 ) CALL read_error(8, 'wind_speed')
      ENDIF

!   Declared Parameters
      Rain_flag = OFF
      IF ( Precip_flag==xyz_dist_module ) Rain_flag = ACTIVE
      IF ( Rain_flag==ACTIVE .OR. Model==DOCUMENTATION ) THEN
        IF ( declvar(MODNAME, 'rain_day', 'one', 1, 'integer', &
     &       'Flag to set the form of any precipitation to rain (0=determine form; 1=rain)', &
     &       'none', Rain_day)/=0 ) CALL read_error(8, 'rain_day')
        IF ( declparam(MODNAME, 'rain_code', 'nmonths', 'integer', &
     &       '2', '1', '5', &
     &       'Flag indicating rule for precipitation station use', &
     &       'Monthly (January to December) flag indicating rule for'// &
     &       ' precipitation measurement station use (1=only'// &
     &       ' precipitation if the regression stations have'// &
     &       ' precipitation; 2=only precipitation'// &
     &       ' if any station in the basin has precipitation;'// &
     &       ' 3=precipitation if xyz says so; 4=only'// &
     &       ' precipitation if rain_day variable is set to 1; 5=only'// &
     &       ' precipitation if psta_freq_nuse stations have precipitation)', &
     &       'none')/=0 ) CALL read_error(1, 'rain_code')
      ENDIF

! Lake Variables
      IF ( Nratetbl>0 ) THEN
        ALLOCATE ( Gate_ht(Nratetbl) )
        IF ( declvar(MODNAME, 'gate_ht', 'nratetbl', Nratetbl, 'real', &
     &       'Height of the gate opening at each dam with a gate', &
     &       'inches', Gate_ht)/=0 ) CALL read_error(8, 'gate_ht')
      ENDIF

      IF ( Nlakeelev>0 ) THEN
        ALLOCATE ( Lake_elev(Nlakeelev) )
        IF ( declvar(MODNAME, 'lake_elev', 'nlakeelev', Nlakeelev, 'real', &
     &       'Elevation of each simulated lake surface', &
     &       'feet', Lake_elev)/=0 ) CALL read_error(8, 'lake_elev')
      ENDIF

      END FUNCTION obsdecl

!***********************************************************************
!     obsinit - initializes obs module
!***********************************************************************
      INTEGER FUNCTION obsinit()
      USE PRMS_CONSTANTS, ONLY: ACTIVE, OFF, MONTHS_PER_YEAR, CFS
      USE PRMS_MODULE, ONLY: Nratetbl, Ntemp, Nrain, Nsol, Nobs, Nevap, Nsnow
      USE PRMS_OBS
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: getparam
      EXTERNAL :: read_error
!***********************************************************************
      obsinit = 0

      Runoff_units = CFS
      IF ( Nobs>0 ) THEN
        IF ( getparam(MODNAME, 'runoff_units', 1, 'integer', Runoff_units)/=0 ) CALL read_error(2, 'runoff_units')
      ENDIF

      IF ( Rain_flag==ACTIVE ) THEN
        IF ( getparam(MODNAME, 'rain_code', MONTHS_PER_YEAR, 'integer', Rain_code)/=0 ) CALL read_error(2, 'rain_code')
      ENDIF

      IF ( Nobs>0 ) THEN
        Runoff = 0.0
        Streamflow_cfs = 0.0D0
        Streamflow_cms = 0.0D0
      ENDIF
      IF ( Nrain>0 ) Precip = 0.0
      Rain_day = OFF
      IF ( Ntemp>0 ) THEN
        Tmax = 0.0
        Tmin = 0.0
      ENDIF
      IF ( Nsol>0 ) Solrad = 0.0
      IF ( Nevap>0 ) Pan_evap = 0.0
      IF ( Nsnow>0 ) Snowdepth = 0.0
      IF ( Nlakeelev>0 ) Lake_elev = 0.0
      IF ( Nratetbl>0 ) Gate_ht = 0.0
      IF ( Nhumid>0 ) Humidity = 0.0
      IF ( Nwind>0 ) Wind_speed = 0.0

      END FUNCTION obsinit

! **********************************************************************
!     obsrun - runs obs module
! **********************************************************************
      INTEGER FUNCTION obsrun()
      USE PRMS_CONSTANTS, ONLY: CMS, CFS2CMS_CONV, ACTIVE
      USE PRMS_MODULE, ONLY: Nratetbl, Ntemp, Nrain, Nsol, Nobs, Nevap, Nsnow, &
                             Nowyear, Nowmonth, Nowday, forcing_check_flag
      USE PRMS_MODULE, ONLY: NmaxPrecipObs                                                      ! PJT - 2018Jan03 - Sub-daily precip inputs
      USE PRMS_OBS
      USE PRMS_CLIMATEVARS, ONLY: Ppt_zero_thresh
      IMPLICIT NONE
! Functions
      INTRINSIC DBLE
      INTEGER, EXTERNAL :: readvar
      EXTERNAL :: read_error
      DOUBLE PRECISION, EXTERNAL :: deltim, delnex                                              ! PJT - 2018Jan03 - Sub-daily precip inputs
      EXTERNAL :: dattim                                                                        ! PJT - 2018Jan03 - Sub-daily precip inputs
      INTEGER, EXTERNAL ::  read_line                                                           ! PJT - 2018Jan03 - Sub-daily precip inputs
! Local Variables
      INTEGER :: i
      INTEGER :: j, CheckTime(6), CallReadLine                                                  ! PJT - 2018Jan03 - Sub-daily precip inputs
      DOUBLE PRECISION :: dtn                                                                   ! PJT - 2018Jan03 - Sub-daily precip inputs
! **********************************************************************
      obsrun = 0

      IF ( Nobs>0 ) THEN
        IF ( readvar(MODNAME, 'runoff')/=0 ) CALL read_error(9, 'runoff')
        IF ( Runoff_units==CMS ) THEN
          DO i = 1, Nobs
            Streamflow_cms(i) = DBLE( Runoff(i) )
            Streamflow_cfs(i) = Streamflow_cms(i)/CFS2CMS_CONV
          ENDDO
        ELSE
          DO i = 1, Nobs
            Streamflow_cfs(i) = DBLE( Runoff(i) )
            Streamflow_cms(i) = Streamflow_cfs(i)*CFS2CMS_CONV
          ENDDO
        ENDIF
      ENDIF

      IF ( Nrain>0 ) THEN
        IF ( readvar(MODNAME, 'precip')/=0 ) CALL read_error(9, 'precip')
        IF ( Ppt_zero_thresh>0.0 ) THEN
          DO i = 1, Nrain
            IF ( Precip(i)<Ppt_zero_thresh ) Precip(i) = 0.0
          ENDDO
        ENDIF
      ENDIF

      IF ( Ntemp>0 ) THEN
        IF ( readvar(MODNAME, 'tmax')/=0 ) CALL read_error(9, 'tmax')
        IF ( readvar(MODNAME, 'tmin')/=0 ) CALL read_error(9, 'tmin')
        IF ( forcing_check_flag == ACTIVE ) THEN
          DO i = 1, Ntemp
            IF ( Tmin(i) < -98 .OR. Tmax(i) < -98 ) CYCLE
            IF ( Tmin(i) > Tmax(i) ) THEN
              PRINT *, 'WARNING, observed tmin > tmax: HRU, date, tmin, tmax:', &
                       i, Nowyear, Nowmonth, Nowday, Tmin(i), Tmax(i)
!              WRITE ( 862, * ) 'WARNING, observed tmin > tmax, swapped fort_866 and fort_867: HRU, date, tmin, tmax:', &
!                               i, Nowyear, Nowmonth, Nowday, Tmin(i), Tmax(i)
!              foo = Tmax(i)
!              Tmax(i) = Tmin(i)
!              Tmin(i) = foo
            ENDIF
          ENDDO
!          WRITE ( 866,  '(I4,2I3,3I2,22F8.2)' ) Nowyear, Nowmonth, Nowday, 0, 0, 0, (Tmax(i), i=1,Ntemp)
!          WRITE ( 867, '(I4,2I3,3I2,22F8.2)' ) Nowyear, Nowmonth, Nowday, 0, 0, 0, (Tmin(i), i=1,Ntemp)
        ENDIF
      ENDIF

      IF ( Nsol>0 ) THEN
        IF ( readvar(MODNAME, 'solrad')/=0 ) CALL read_error(9, 'solrad')
      ENDIF

      IF ( Nevap>0 ) THEN
        IF ( readvar(MODNAME, 'pan_evap')/=0 ) CALL read_error(9, 'pan_evap')
      ENDIF

      IF ( Nsnow>0 ) THEN
        IF ( readvar(MODNAME, 'snowdepth')/=0 ) CALL read_error(9, 'snowdepth')
      ENDIF

      IF ( Rain_flag==1 ) THEN
        IF ( Rain_code(Nowmonth)==4 ) THEN
          IF ( readvar(MODNAME, 'rain_day')/=0 ) CALL read_error(9, 'rain_day')
        ENDIF
      ENDIF

      IF ( Nlakeelev>0 ) THEN
        IF ( readvar(MODNAME, 'lake_elev')/=0 ) CALL read_error(9, 'lake_elev')
      ENDIF

      IF ( Nratetbl>0 ) THEN
        IF ( readvar(MODNAME, 'gate_ht')/=0 ) CALL read_error(9, 'gate_ht')
      ENDIF

      IF ( Nhumid>0 ) THEN
        IF ( readvar(MODNAME, 'humidity')/=0 ) CALL read_error(9, 'humidity')
      ENDIF

      IF ( Nwind>0 ) THEN
        IF ( readvar(MODNAME, 'wind_speed')/=0 ) CALL read_error(9, 'wind_speed')
      ENDIF

      !If specified, check for sub-daily precip values in the DATA file(s)                  ! PJT begin - 2018Jan03 - Sub-daily precip inputs  (Start)
      IF ( NmaxPrecipObs>0 ) THEN                                                           ! PJT - 2018Jan03
          
          !Clear the subdaily arrays
          SubDailyPrecip = 0.0          
          SubDailyObsTimes = 0.0
          !Zero the counter
          NumPrecipObs=0

          !Check if there are intensity values provided for this time step (returns the timespan between the current and next observation in hours)
          dtn = delnex()

          !Get the time stamp of the intensity observation
          CALL dattim('now', CheckTime)
          
          !Are there instantaneous records provided? 
          IF (dtn < 24)  THEN
                !A series of values have been provided, assume these are precip observations (must have activated the Green&Ampt package to get here anyways)
              
                !Is the first observation a daily value which we should ignore?
                IF (CheckTime(4) + CheckTime(5)+ CheckTime(6) == 0) THEN
                    !Call MMS to advance to the first sub-daily observation within the DATA file(s)
                    CallReadLine= read_line()
                ENDIF
                
            DO
                !MMS has been advanced to the first set of sub-daily observations, increment through all sub-daily values while storing 
                
                !Increment the counter
                NumPrecipObs=NumPrecipObs + 1
                
                !Get the timestamp of the intensity observation
                CALL dattim('now', CheckTime)
                
                !Store timestamp (decimal hours)
                SubDailyObsTimes(NumPrecipObs)=CheckTime(4)+CheckTime(5)/SNGL(60)+CheckTime(6)/SNGL(3600)
                !Sanity Check, all time values should be less than 24hr
                IF (SubDailyObsTimes(NumPrecipObs)>24) THEN
                    PRINT *, 'ERROR (obs.f90) - There are sub-daily timestamp values that exceed 24hrs.   Timestamp (hrs) =',SubDailyObsTimes(NumPrecipObs)
                ENDIF
                
                !Read the intesity observations for this timestamp
                IF ( readvar(MODNAME, 'precip')/=0 ) CALL read_error(9, 'precip')

                !Store the precip values
                DO i = 1, Nrain
                    SubDailyPrecip(i,NumPrecipObs)=Precip(i)
                ENDDO
               
                !Call for the duration of intensity interval
                dtn = delnex()

                !Does the next intensity value lie outside this computational day? (i.e., is the sum of the found intervals longer than 24hrs)
                IF (dtn + DBLE(CheckTime(4))+DBLE(CheckTime(5))/DBLE(60)+DBLE(CheckTime(6))/DBLE(3600)>24.0 - 1/DBLE(360000))  THEN  !1/DBLE(360000) is a quick and dirty check for round off error
                    !This value is from another computational day, EXIT the loop
                    EXIT
                ENDIF
                
                !Call MMS to advance to the next observation within the DATA file(s), will only stay in this loop if there are additional sub-daily values associated with this Julian day
                CallReadLine= read_line()

            END DO
           
            !Anytime sub-daily values are provided, the daily values are ignored completely (i.e., the total daily precip will need to be summed from the sub-daily values)
            Precip=0.0
            DO j = 1, NumPrecipObs
                DO i = 1, Nrain
                    Precip(i) =Precip(i)+ SubDailyPrecip(i,j)
                ENDDO
            ENDDO
                                                                                            ! PJT - 2018Jan03
          ENDIF                                                                             ! PJT end - 2018Jan03 - Sub-daily precip inputs  (End)

      ENDIF
 
      END FUNCTION obsrun
