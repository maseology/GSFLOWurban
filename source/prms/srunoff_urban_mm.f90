!***********************************************************************
! Computes surface runoff and infiltration for each HRU using a
! design specific to the analysis waterbalance computations for
! urban development.
! (Intended to be used in a (grid-based) distributed fashion.)
!
! Then main difference of this module is the incorperation of built-
! up structures (i.e., buildings) that is disconnected from overland
! flow processes.
! The module includes LID functionality, and storm drainage systems that 
! will interact with the moving groundater table.
!    
! mm, 12/08/2015 initial development     
!***********************************************************************
      MODULE PRMS_URBAN
        IMPLICIT NONE
        ! Local Variables
        CHARACTER(LEN=13), SAVE :: MODNAME
        ! Declared Variables
        DOUBLE PRECISION, SAVE :: Basin_dscn_stor, Basin_dscn_evap
        REAL, SAVE, ALLOCATABLE :: Dscn_stor(:), Dscn_stor_evap(:), Hru_dscnstorevap(:)
        REAL, SAVE, ALLOCATABLE :: Hru_frac_dscn(:), Hru_dscn(:), Dscn_frac_hru(:)
        REAL, SAVE, ALLOCATABLE :: Hru_infstor(:), Infstor_frac_hru(:)
        DOUBLE PRECISION, SAVE :: Basin_infstor, Basin_infstor_seep
        REAL, SAVE, ALLOCATABLE :: Inf_stor(:), Infstor_seep(:)
        DOUBLE PRECISION, SAVE :: Basin_stdrn_infil, Basin_stdrn_seep, Basin_urbanfarflow
        DOUBLE PRECISION, SAVE :: Basin_urban_to_ssr, Basin_urban_to_infstor
        REAL, SAVE, ALLOCATABLE :: Urban_to_ssr(:), Urban_to_infstor(:), Urban_to_stdrn(:), Stdrn_seep(:)
        REAL, SAVE, ALLOCATABLE :: Urban_to_farflow(:), Urban_to_strm_seg(:), Urban_to_lake(:)
        
        ! Declared Parameters
        INTEGER, SAVE, ALLOCATABLE :: Dscn_hru_id(:), Infstor_hru_id(:), Stdrn_hru_id(:)
        REAL, SAVE, ALLOCATABLE :: Imperv_frac_dscn(:), Hru_frac_infstor(:)
        REAL, SAVE, ALLOCATABLE :: Dscn_stor_max(:), Dscn_evap_frac(:)
        REAL, SAVE, ALLOCATABLE :: Dscn_to_stdrn(:), Dscn_to_infstor(:)
        REAL, SAVE, ALLOCATABLE :: Imperv_stor_seep(:), Imperv_stor_to_stdrn(:), Imperv_stor_to_infstor(:)
        REAL, SAVE, ALLOCATABLE :: Infstor_max(:), Infstor_invert(:)
        REAL, SAVE, ALLOCATABLE :: Infstor_seep_coef(:), Infstor_to_stdrn_coef(:)
        REAL, SAVE, ALLOCATABLE :: Stdrn_cond(:), Stdrn_invert(:)
        
      END MODULE PRMS_URBAN
    
!***********************************************************************
!     Main urban waterbalance routine
!***********************************************************************
      INTEGER FUNCTION srunoff_urban()
      USE PRMS_MODULE, ONLY: Process, Save_vars_to_file, Init_vars_from_file
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: urbandecl, urbaninit
      EXTERNAL :: urban_restart
!***********************************************************************
      srunoff_urban = 0

      IF ( Process(:4)=='decl' ) THEN
        srunoff_urban = urbandecl()
      ELSEIF ( Process(:4)=='init' ) THEN
        IF ( Init_vars_from_file==1 ) CALL urban_restart(1)
        srunoff_urban = urbaninit()
      ELSEIF ( Process(:5)=='clean' ) THEN
        IF ( Save_vars_to_file==1 ) CALL urban_restart(0)
      ENDIF

      END FUNCTION srunoff_urban    
    
!***********************************************************************
!     urbandecl - set up parameters for urban waterbalance computations
!   Declared Parameters
!     smidx_coef, smidx_exp, carea_max, imperv_stor_max, snowinfil_max
!     hru_area, soil_moist_max, soil_rechr_max, carea_min
!***********************************************************************
      INTEGER FUNCTION urbandecl()
      USE PRMS_URBAN
      USE PRMS_MODULE, ONLY: Nhru, Ndscn, Ninfstor, MAXDIM
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: decldim, declvar, declparam
      EXTERNAL read_error
! Local Variables
      CHARACTER(LEN=80), SAVE :: Version_srunoff_urban
!***********************************************************************
      urbandecl = 0
        
      IF ( decldim('ndscn', 0, MAXDIM, 'Number of disconnected reservoirs')/=0 ) CALL read_error(7, 'ndscn')
      IF ( decldim('ninfstor', 0, MAXDIM, 'Number of infiltration storage reservoirs')/=0 ) CALL read_error(7, 'ninfstor')
      
! Disconnected storage variables
      IF ( declvar(MODNAME, 'basin_dscn_evap', 'one', 1, 'double', &
     &     'Basin area-weighted average evaporation from disconnected storage', &
     &     'inches', Basin_dscn_evap)/=0 ) CALL read_error(3, 'basin_dscn_evap')

      IF ( declvar(MODNAME, 'basin_dscn_stor', 'one', 1, 'double', &
     &     'Basin area-weighted average storage in disconnected storage', &
     &     'inches', Basin_dscn_stor)/=0 ) CALL read_error(3, 'basin_dscn_stor')      

      IF ( declvar(MODNAME, 'basin_urbanfarflow', 'one', 1, 'double', &
     &     'Basin area-weighted average urban to farflow contirbution', &
     &     'inches', Basin_urbanfarflow)/=0 ) CALL read_error(3, 'basin_urbanfarflow')   
      
      IF ( declvar(MODNAME, 'basin_urban_to_ssr', 'one', 1, 'double', &
     &     'Basin area-weighted average urban runoff redirected to the subsurface', &
     &     'inches', Basin_urban_to_ssr)/=0 ) CALL read_error(3, 'basin_urban_to_ssr')   
      
      IF ( declvar(MODNAME, 'basin_urban_to_infstor', 'one', 1, 'double', &
     &     'Basin area-weighted average urban runoff redirected to infiltration reservoirs', &
     &     'inches', Basin_urban_to_infstor)/=0 ) CALL read_error(3, 'basin_urban_to_infstor')
      
      ALLOCATE ( Hru_dscn(Nhru) )
      IF ( declvar(MODNAME, 'hru_dscn', 'nhru', Nhru, 'real', &
     &     'Area of HRU that is disconnected', &
     &     'acres', Hru_dscn)/=0 ) CALL read_error(3, 'hru_dscn')
 
      ALLOCATE ( Dscn_frac_hru(Nhru) )
      IF ( declvar(MODNAME, 'dscn_frac_hru', 'nhru', Nhru, 'real', &
     &     'Fraction of disconnected reservoir contained within each HRU', &
     &     'decimal fraction', Dscn_frac_hru)/=0 ) CALL read_error(3, 'dscn_frac_hru')      
 
      ALLOCATE ( Hru_frac_dscn(Nhru) )
      IF ( declvar(MODNAME, 'hru_frac_dscn', 'nhru', Nhru, 'real', &
     &     'Fraction of HRU that is disconnected', &
     &     'decimal fraction', Hru_frac_dscn)/=0 ) CALL read_error(3, 'hru_frac_dscn')             
      
      ALLOCATE ( Dscn_stor(Ndscn) )
      IF ( declvar(MODNAME, 'dscn_stor', 'ndscn', Ndscn, 'real', &
     &     'Area-weighted average storage from disconnected reservoir', &
     &     'inches', Dscn_stor)/=0 ) CALL read_error(3, 'dscn_stor') 

      ALLOCATE ( Dscn_stor_evap(Ndscn) )
      IF ( declvar(MODNAME, 'dscn_stor_evap', 'ndscn', Ndscn, 'real', &
     &     'Area-weighted average evaporation from disconnected storage reservoir', &
     &     'inches', Dscn_stor_evap)/=0 ) CALL read_error(3, 'dscn_stor_evap')   
      
      ALLOCATE ( Hru_dscnstorevap(Nhru) )
      IF ( declvar(MODNAME, 'hru_dscnstorevap', 'nhru', Nhru, 'real', &
     &     'HRU area-weighted average evaporation from disconnected storage reservoir', &
     &     'inches', Hru_dscnstorevap)/=0 ) CALL read_error(3, 'hru_dscnstorevap') 
      
      ALLOCATE ( Urban_to_infstor(Ninfstor) )
      IF ( declvar(MODNAME, 'urban_to_infstor', 'ninfstor', Ninfstor, 'real', &
     &     'Quantity of water being directed to infiltration storage systems', &
     &     'inches', Urban_to_infstor)/=0 ) CALL read_error(3, 'urban_to_infstor') 
      
      ALLOCATE ( Urban_to_stdrn(Nhru) )
      IF ( declvar(MODNAME, 'urban_to_stdrn', 'nhru', Nhru, 'real', &
     &     'Quantity of water being directed from urban areas to storm drains', &
     &     'inches', Urban_to_stdrn)/=0 ) CALL read_error(3, 'urban_to_stdrn') 
      
      ALLOCATE ( Urban_to_farflow(Nhru) )
      IF ( declvar(MODNAME, 'urban_to_farflow', 'nhru', Nhru, 'real', &
     &     'Quantity of water being directed from urban areas to farfield', &
     &     'inches', Urban_to_farflow)/=0 ) CALL read_error(3, 'urban_to_farflow') 

      ALLOCATE ( Urban_to_strm_seg(Nhru) )
      IF ( declvar(MODNAME, 'urban_to_strm_seg', 'nhru', Nhru, 'real', &
     &     'Quantity of water being directed from urban areas to streams', &
     &     'inches', Urban_to_strm_seg)/=0 ) CALL read_error(3, 'urban_to_strm_seg') 
      
      ALLOCATE ( Urban_to_lake(Nhru) )
      IF ( declvar(MODNAME, 'urban_to_lake', 'nhru', Nhru, 'real', &
     &     'Quantity of water being directed from urban areas to lakes', &
     &     'inches', Urban_to_lake)/=0 ) CALL read_error(3, 'urban_to_lake') 
      
      ALLOCATE ( Urban_to_ssr(Nhru) )
      IF ( declvar(MODNAME, 'urban_to_ssr', 'nhru', Nhru, 'real', &
     &     'Quantity of water being directed from urban surfaces to ssr', &
     &     'inches', Urban_to_ssr)/=0 ) CALL read_error(3, 'urban_to_ssr') 
      
! Disconnected storage parameters
      ALLOCATE ( Imperv_frac_dscn(Nhru) )
      IF ( declparam(MODNAME, 'imperv_frac_dscn', 'nhru', 'real', &
     &     '0.0', '0.0', '1.0', &
     &     'Fraction of impervious area that is disconnected', &
     &     'Fraction of HRU impervious area that is disconnected', &
     &     'decimal fraction')/=0 ) CALL read_error(1, 'imperv_frac_dscn')    

      ALLOCATE ( Hru_frac_infstor(Nhru) )
      IF ( declparam(MODNAME, 'hru_frac_infstor', 'nhru', 'real', &
     &     '0.0', '0.0', '1.0', &
     &     'Fraction of HRU that overlays an infiltration storage reservoir', &
     &     'Fraction of HRU that overlays an infiltration storage reservoir', &
     &     'decimal fraction')/=0 ) CALL read_error(1, 'hru_frac_infstor')       
      
      ALLOCATE ( Dscn_stor_max(Ndscn) )
      IF ( declparam(MODNAME, 'dscn_stor_max', 'ndscn', 'real', &
     &     '0.05', '0.0', '40.0', &
     &     'Maximum disconnected retention storage', &
     &     'Maximum disconnected retention storage for each reservoir', &
     &     'inches')/=0 ) CALL read_error(1, 'dscn_stor_max')
     
      ALLOCATE ( Dscn_to_stdrn(Nhru) )
      IF ( declparam(MODNAME, 'dscn_to_stdrn', 'nhru', 'real', &
     &     '0.0', '0.0', '1.0', &
     &     'Fraction of disconnected storage that is directed to storm drainage', &
     &     'Fraction of water reaching the disconnected surface'// &
     &     ' that is directed to storm drainage systems', &
     &     'decimal fraction')/=0 ) CALL read_error(1, 'dscn_to_stdrn')           
      
      ALLOCATE ( Dscn_to_infstor(Nhru) )
      IF ( declparam(MODNAME, 'dscn_to_infstor', 'nhru', 'real', &
     &     '0.0', '0.0', '1.0', &
     &     'Fraction of disconnected storage that is directed into infiltration storage', &
     &     'Fraction of water reaching the disconnected surface'// &
     &     ' that is directed into on-site infiltration storage systems', &
     &     'decimal fraction')/=0 ) CALL read_error(1, 'dscn_to_infstor')  
      
      ALLOCATE ( Dscn_evap_frac(Nhru) )
      IF ( declparam(MODNAME, 'dscn_evap_frac', 'nhru', 'real', &
     &     '1.0', '0.0', '1.0', &
     &     'Fraction of unsatisfied potential evapotranspiration to apply to disconnected storage', &
     &     'Fraction of unsatisfied potential evapotranspiration to apply to disconnected storage', &
     &     'decimal fraction')/=0 ) CALL read_error(1, 'dscn_evap_frac')      
      
      ALLOCATE ( Dscn_hru_id(Nhru) )
      IF ( declparam(MODNAME, 'dscn_hru_id', 'nhru', 'integer', &
     &     '0', 'bounded', 'ndscn', &
     &     'Identification number of the disconnected reservoir associated with an HRU', &
     &     'Identification number of the disconnected reservoir associated with an HRU;'// &
     &     ' more than one HRU can be associated with each disconnected reservoir', &
     &     'none')/=0 ) CALL read_error(1, 'dscn_hru_id')  
      
! Storm drain variables
      IF ( declvar(MODNAME, 'basin_stdrn_seep', 'one', 1, 'double', &
     &     'Basin area-weighted average seepage from storm drainage network', &
     &     'inches', Basin_stdrn_seep)/=0 ) CALL read_error(3, 'basin_stdrn_seep')  
      
      ALLOCATE ( Stdrn_seep(Nhru) )
      IF ( declvar(MODNAME, 'stdrn_seep', 'nhru', Nhru, 'real', &
     &     'Seepage from HRU storm drains', &
     &     'inches', Stdrn_seep)/=0 ) CALL read_error(3, 'stdrn_seep') 
      
! Infiltration storage variables
      IF ( declvar(MODNAME, 'basin_infstor_seep', 'one', 1, 'double', &
     &     'Basin area-weighted average seepage from infiltration storage systems', &
     &     'inches', Basin_infstor_seep)/=0 ) CALL read_error(3, 'basin_infstor_seep')      
      
      IF ( declvar(MODNAME, 'basin_infstor', 'one', 1, 'double', &
     &     'Basin area-weighted average storage in infiltration storage systems', &
     &     'inches', Basin_infstor)/=0 ) CALL read_error(3, 'basin_infstor')        

      ALLOCATE ( Hru_infstor(Nhru) )
      IF ( declvar(MODNAME, 'hru_infstor', 'nhru', Nhru, 'real', &
     &     'Area of HRU that is underlain by an infiltration storage facility', &
     &     'acres', Hru_infstor)/=0 ) CALL read_error(3, 'hru_infstor')
 
      ALLOCATE ( Infstor_frac_hru(Nhru) )
      IF ( declvar(MODNAME, 'infstor_frac_hru', 'nhru', Nhru, 'real', &
     &     'Fraction of infiltration storage reservoir contained within each HRU', &
     &     'decimal fraction', Infstor_frac_hru)/=0 ) CALL read_error(3, 'infstor_frac_hru')
      
      ALLOCATE ( Inf_stor(Ninfstor) )
      IF ( declvar(MODNAME, 'inf_stor', 'ninfstor', Ninfstor, 'real', &
     &     'Area-weighted average storage from infiltration storage reservoir', &
     &     'inches', Inf_stor)/=0 ) CALL read_error(3, 'inf_stor') 

      ALLOCATE ( Infstor_seep(Ninfstor) )
      IF ( declvar(MODNAME, 'infstor_seep', 'ninfstor', Ninfstor, 'real', &
     &     'Seepage from infiltration storage reservoir', &
     &     'inches', Infstor_seep)/=0 ) CALL read_error(3, 'infstor_seep')        
    
! Infiltration storage parameters
      ALLOCATE ( Infstor_max(Ninfstor) )
      IF ( declparam(MODNAME, 'infstor_max', 'ninfstor', 'real', &
     &     '10.0', '0.0', '100.0', &
     &     'Maximum infiltration detention storage', &
     &     'Maximum infiltration detention storage for each reservoir', &
     &     'inches')/=0 ) CALL read_error(1, 'infstor_max')      
      
      ALLOCATE ( Infstor_invert(Ninfstor) )
      IF ( declparam(MODNAME, 'infstor_invert', 'ninfstor', 'real', &
     &     '0.0', '-1000.0', '30000.0', &
     &     'Invert elevation of infiltration detention storage reservoir', &
     &     'Invert elevation of infiltration detention storage reservoir', &
     &     'elev_units')/=0 ) CALL read_error(1, 'infstor_invert')  
      
      ALLOCATE ( Infstor_seep_coef(Nhru) )
      IF ( declparam(MODNAME, 'infstor_seep_coef', 'nhru', 'real', &
     &     '0.02', '0.0', '100.0', &
     &     'Basal conductance if infiltration storage reservoirs', &
     &     'Seepage conductance at the base of infiltration storage reservoirs for each HRU', &
     &     'inches/day')/=0 ) CALL read_error(1, 'infstor_seep_coef')     

      ALLOCATE ( Infstor_to_stdrn_coef(Nhru) )
      IF ( declparam(MODNAME, 'infstor_to_stdrn_coef', 'nhru', 'real', &
     &     '0.0', '0.0', '1.0', &
     &     'Coefficient used in linear drainage flow from infiltration storage to storm drainage', &
     &     'Coefficient used in linear drainage flow from infiltration'// &
     &     ' storage to storm drainage system', &
     &     'fraction/day')/=0 ) CALL read_error(1, 'infstor_to_stdrn_coef')      

      ALLOCATE ( Infstor_hru_id(Nhru) )
      IF ( declparam(MODNAME, 'infstor_hru_id', 'nhru', 'integer', &
     &       '0', 'bounded', 'ninfstor', &
     &       'Identification number of the infiltration reservoir associated with an HRU', &
     &       'Identification number of the infiltration reservoir associated with an HRU;'// &
     &       ' more than one HRU can be associated with each disconnected reservoir', &
     &       'none')/=0 ) CALL read_error(1, 'infstor_hru_id')     

! Storm drainge variables
      IF ( declvar(MODNAME, 'basin_stdrn_infil', 'one', 1, 'double', &
     &     'Basin area-weighted average infiltration loss into storm drainage systems', &
     &     'inches', Basin_stdrn_infil)/=0 ) CALL read_error(3, 'basin_stdrn_infil') 
      
! Storm drainge parameters
      ALLOCATE ( Stdrn_invert(Nhru) )
      IF ( declparam(MODNAME, 'stdrn_invert', 'nhru', 'real', &
     &     '0.0', '-1000.0', '30000.0', &
     &     'Invert elevation of HRU storm drain', &
     &     'Invert elevation of HRU storm drain', &
     &     'elev_units')/=0 ) CALL read_error(1, 'stdrn_invert')        
      
      ALLOCATE ( Stdrn_cond(Nhru) )
      IF ( declparam(MODNAME, 'stdrn_cond', 'nhru', 'real', &
     &     '35.0', '0.0', '100.0', &
     &     'Conductance of storm drain walls used to compute groundwater infiltration and leakage', &
     &     'Conductance of storm drain walls used to compute groundwater infiltration and leakage', &
     &     'inches/day')/=0 ) CALL read_error(1, 'stdrn_cond')       
      
      ALLOCATE ( Stdrn_hru_id(Nhru) )
      IF ( declparam(MODNAME, 'stdrn_hru_id', 'nhru', 'integer', &
     &       '0', '-1000', '1000', &
     &       'Identification number of the storm sewershed associated with an HRU', &
     &       'Identification number of the storm sewershed associated with an HRU;'// &
     &       ' more than one HRU can be associated with each storm sewershed;'// &
     &       ' <0 lake ID; >0 stream segment ID; =0 to farfield', &
     &       'none')/=0 ) CALL read_error(1, 'stdrn_hru_id')      
      
! Additional impervious area parameters
      ALLOCATE ( Imperv_stor_seep(Nhru) )
      IF ( declparam(MODNAME, 'imperv_stor_seep', 'nhru', 'real', &
     &     '0.0', '0.0', '1.0', &
     &     'Fraction of impervious storage that is allowed to infiltrate', &
     &     'Fraction of impervious storage that is allowed to infiltrate', &
     &     'fraction/day')/=0 ) CALL read_error(1, 'imperv_stor_seep')     

      ALLOCATE ( Imperv_stor_to_stdrn(Nhru) )
      IF ( declparam(MODNAME, 'imperv_stor_to_stdrn', 'nhru', 'real', &
     &     '1.0', '0.0', '1.0', &
     &     'Fraction of impervious storage that is directed to storm drainage system', &
     &     'Fraction of impervious storage that is directed to storm drainage system', &
     &     'fraction/day')/=0 ) CALL read_error(1, 'imperv_stor_to_stdrn')      

      ALLOCATE ( Imperv_stor_to_infstor(Nhru) )
      IF ( declparam(MODNAME, 'imperv_stor_to_infstor', 'nhru', 'real', &
     &     '0.0', '0.0', '1.0', &
     &     'Fraction of impervious storage that is directed to infiltration storage system', &
     &     'Fraction of impervious storage that is directed to infiltration storage system', &
     &     'fraction/day')/=0 ) CALL read_error(1, 'imperv_stor_to_infstor')
      
      END FUNCTION urbandecl     
    
!***********************************************************************
!     urbaninit - Initialize urban module - get parameter values
!***********************************************************************
      INTEGER FUNCTION urbaninit()
      USE PRMS_URBAN
      USE PRMS_MODULE, ONLY: Model, Nhru, Ndscn, Ninfstor, Print_debug, &
     &    Inputerror_flag, Parameter_check_flag, Init_vars_from_file
      USE PRMS_BASIN, ONLY: Hru_area, Hru_elev, Hru_imperv, Hru_percent_imperv, &
     &    NEARZERO
      IMPLICIT NONE
      INTEGER, EXTERNAL :: getdim, getparam
      EXTERNAL read_error, check_param_value, check_param_limits
! Local Variables
      INTEGER :: i, d, f, num_hrus
      REAL :: frac_sum
      INTEGER, ALLOCATABLE :: dscn_id_cnt(:), infstor_id_cnt(:)
      REAL, ALLOCATABLE :: dscn_area(:), infstor_area(:)
!***********************************************************************
      urbaninit = 0
        
      IF ( Init_vars_from_file==0 ) THEN
        Basin_dscn_stor = 0.0D0
        Basin_dscn_evap = 0.0D0
        Basin_infstor = 0.0D0
        Basin_infstor_seep = 0.0D0
        Basin_stdrn_seep = 0.0D0
        Basin_stdrn_infil = 0.0D0
        Basin_urbanfarflow = 0.0D0
        Basin_urban_to_ssr = 0.0D0
        Basin_urban_to_infstor = 0.0D0
        Dscn_stor = 0.0
        Dscn_stor_evap = 0.0
        Hru_dscnstorevap = 0.0
        Inf_stor = 0.0
        Infstor_seep = 0.0
        Stdrn_seep = 0.0
        Urban_to_infstor = 0.0
        Urban_to_ssr = 0.0
        Urban_to_stdrn = 0.0
        Urban_to_farflow = 0.0
        Urban_to_strm_seg = 0.0
        Urban_to_lake = 0.0
      ENDIF

      IF ( getparam(MODNAME, 'imperv_stor_seep', Nhru, 'real', Imperv_stor_seep)/=0 ) CALL read_error(2, 'imperv_stor_seep')
      IF ( getparam(MODNAME, 'stdrn_hru_id', Nhru, 'integer', Stdrn_hru_id)/=0 ) CALL read_error(2, 'stdrn_hru_id')
      IF ( getparam(MODNAME, 'imperv_stor_to_stdrn', Nhru, 'real', Imperv_stor_to_stdrn)/=0 ) CALL read_error(2, 'imperv_stor_to_stdrn')
      IF ( getparam(MODNAME, 'stdrn_invert', Nhru, 'real', Stdrn_invert)/=0 ) CALL read_error(2, 'stdrn_invert')
      IF ( getparam(MODNAME, 'stdrn_cond', Nhru, 'real', Stdrn_cond)/=0 ) CALL read_error(2, 'stdrn_cond')      
      
      IF ( Ndscn>0 ) THEN
        IF ( getparam(MODNAME, 'imperv_frac_dscn', Nhru, 'real', Imperv_frac_dscn)/=0 ) CALL read_error(2, 'imperv_frac_dscn')    
        IF ( getparam(MODNAME, 'dscn_hru_id', Nhru, 'integer', Dscn_hru_id)/=0 ) CALL read_error(2, 'dscn_hru_id')
        IF ( getparam(MODNAME, 'dscn_stor_max', Ndscn, 'real', Dscn_stor_max)/=0 ) CALL read_error(2, 'dscn_stor_max')
        IF ( getparam(MODNAME, 'dscn_evap_frac', Nhru, 'real', Dscn_evap_frac)/=0 ) CALL read_error(2, 'dscn_evap_frac')
        IF ( getparam(MODNAME, 'dscn_to_stdrn', Nhru, 'real', Dscn_to_stdrn)/=0 ) CALL read_error(2, 'dscn_to_stdrn')
        IF ( Ninfstor>0 ) THEN
          IF ( getparam(MODNAME, 'dscn_to_infstor', Nhru, 'real', Dscn_to_infstor)/=0 ) CALL read_error(2, 'dscn_to_infstor')          
        ENDIF
      ELSE
        Imperv_frac_dscn = 0.0
      ENDIF

      IF ( Ninfstor>0 ) THEN
        IF ( getparam(MODNAME, 'hru_frac_infstor', Nhru, 'real', Hru_frac_infstor)/=0 ) CALL read_error(2, 'hru_frac_infstor')      
        IF ( getparam(MODNAME, 'infstor_hru_id', Nhru, 'integer', Infstor_hru_id)/=0 ) CALL read_error(2, 'infstor_hru_id')
        IF ( getparam(MODNAME, 'infstor_max', Ninfstor, 'real', Infstor_max)/=0 ) CALL read_error(2, 'infstor_max')
        IF ( getparam(MODNAME, 'infstor_invert', Ninfstor, 'real', Infstor_invert)/=0 ) CALL read_error(2, 'infstor_invert')
        IF ( getparam(MODNAME, 'infstor_seep_coef', Nhru, 'real', Infstor_seep_coef)/=0 ) CALL read_error(2, 'infstor_seep_coef')
        IF ( getparam(MODNAME, 'infstor_to_stdrn_coef', Nhru, 'real', Infstor_to_stdrn_coef)/=0 ) CALL read_error(2, 'infstor_to_stdrn_coef')
        IF ( getparam(MODNAME, 'imperv_stor_to_infstor', Nhru, 'real', Imperv_stor_to_infstor)/=0 ) CALL read_error(2, 'imperv_stor_to_infstor')
      ELSE
        Hru_frac_infstor = 0.0
      ENDIF
      
      ! check values
      ALLOCATE ( dscn_id_cnt(Ndscn) )
      dscn_id_cnt = 0
      num_hrus = 0
      DO i = 1, Nhru
        d = Dscn_hru_id(i)
        IF ( d>0 ) THEN
          dscn_id_cnt(d) = dscn_id_cnt(d) + 1
          IF ( Imperv_frac_dscn(i)==0.0 ) THEN
            num_hrus = num_hrus + 1
            Imperv_frac_dscn(i)=0.001 
          ENDIF            
        ENDIF
      ENDDO
      IF ( num_hrus>0 .AND. Print_debug>-1 ) THEN
        WRITE (*, '(/,A,/,9X,A,/,9X,A,I7,/,9X,A,/)') &
     &         'WARNING: fraction of disconnected cover equals 0.0, but the', &
     &         'HRU has been assigned a disconnected storage reservoir ID.', &
     &         'number of HRUs for which this condition exists:', num_hrus, &
     &         'These disconnected cover fractions have been reset to 0.001.'
      ENDIF
      num_hrus = 0
      DO i = 1, Ndscn
        IF ( dscn_id_cnt(i)==0 ) num_hrus = num_hrus + 1
      ENDDO
      IF ( num_hrus>0 .AND. Print_debug>-1 ) THEN
        WRITE (*, '(/,A,I4,A,/,9X,A,/)') &
     &         'WARNING: ', num_hrus, ' disconnected reservoirs were not assigned ', &
     &         'to any HRU. These reservoirs will be ignored.'
      ENDIF
      !num_hrus = 0
      !DO i = 1, Nhru
      !  d = Dscn_hru_id(i)
      !  IF ( d==0 ) THEN
      !    IF ( Imperv_frac_dscn(i)>0.0 ) THEN
      !      Ndscn = Ndscn + 1
      !      Dscn_hru_id(i) = Ndscn
      !    ENDIF
      !  ENDIF
      !ENDDO
      
      ALLOCATE ( infstor_id_cnt(Ninfstor) )
      infstor_id_cnt = 0
      num_hrus = 0
      DO i = 1, Nhru
        f = Infstor_hru_id(i)
        IF ( f>0 ) THEN
          infstor_id_cnt(f) = infstor_id_cnt(f) + 1
          IF ( Hru_frac_infstor(i)==0.0 ) THEN
            num_hrus = num_hrus + 1
            Hru_frac_infstor(i) = 1.0 
          ENDIF            
        ENDIF
      ENDDO
      IF ( num_hrus>0 .AND. Print_debug>-1 ) THEN
        WRITE (*, '(/,A,/,9X,A,/,9X,A,I7,/,9X,A,/)') &
     &         'WARNING: fraction of HRU covering an infiltration reservoir', &
     &         'equals 0.0, but the HRU has been assigned a reservoir ID.', &
     &         'number of HRUs for which this condition exists:', num_hrus, &
     &         'These infiltration reservoir fractions have been reset to 1.0.'
      ENDIF
      num_hrus = 0
      DO i = 1, Ninfstor
        IF ( infstor_id_cnt(i)==0 ) num_hrus = num_hrus + 1
      ENDDO
      IF ( num_hrus>0 .AND. Print_debug>-1 ) THEN
        WRITE (*, '(/,A,I4,A,/,9X,A,/)') &
     &         'WARNING: ', num_hrus, ' infiltration reservoirs were not assigned ', &
     &         'to any HRU. These reservoirs will be ignored.'
      ENDIF
      !num_hrus = 0
      !DO i = 1, Nhru
      !  f = Infstor_hru_id(i)
      !  IF ( f==0 ) THEN
      !    IF ( Imperv_frac_dscn(i)>0.0 ) THEN
      !      Ninfstor = Ninfstor + 1
      !      Infstor_hru_id(i) = Ninfstor
      !    ENDIF
      !  ENDIF
      !ENDDO
      
      ! compute areas of disconnected storage
      ALLOCATE ( dscn_area(Ndscn), infstor_area(Ninfstor) )
      dscn_area = 0.0
      infstor_area = 0.0
      num_hrus = 0
      DO i = 1, Nhru
        IF ( Imperv_frac_dscn(i)>0.0 .AND. Hru_percent_imperv(i)==0.0 ) THEN
          num_hrus = num_hrus + 1
          Imperv_frac_dscn(i) = 0.0
        ENDIF
        IF ( Imperv_frac_dscn(i)>0.999 ) Imperv_frac_dscn(i) = 1.0
        IF ( Imperv_frac_dscn(i)<0.001 ) Imperv_frac_dscn(i) = 0.0
        Hru_frac_dscn(i) = Imperv_frac_dscn(i)*Hru_percent_imperv(i) 
        Hru_dscn(i) = Hru_area(i)*Hru_frac_dscn(i)
        Hru_infstor(i) = Hru_area(i)*Hru_frac_infstor(i)
        Hru_percent_imperv(i) = (1.0 - Imperv_frac_dscn(i))*Hru_percent_imperv(i) ! correct (connected) impervious areas
        Hru_imperv(i) = Hru_area(i)*Hru_percent_imperv(i)                         ! same
        d = Dscn_hru_id(i)
        f = Infstor_hru_id(i)
        IF ( d>0 ) dscn_area(d) = dscn_area(d) + Hru_dscn(i)
        IF ( f>0 ) infstor_area(f) = infstor_area(f) + Hru_infstor(i)
      ENDDO
      IF ( num_hrus>0 .AND. Print_debug>-1 ) THEN
        WRITE (*, '(/,A,/,9X,A,/,9X,A,I7,/,9X,A,/)') &
     &         'WARNING: fraction of impervious cover equals 0.0 and the', &
     &         'specified disconnected fraction is greater than 0.0', &
     &         'number of HRUs for which this condition exists:', num_hrus, &
     &         'These disconnected fractions will be reset to zero.'
      ENDIF  
      
      DO i = 1, Nhru
        d = Dscn_hru_id(i)
        IF ( d>0 ) THEN
          Dscn_frac_hru(i) = Hru_dscn(i)/dscn_area(d)
          IF ( Dscn_frac_hru(i)<NEARZERO ) THEN
            Dscn_frac_hru(i) = 0.0
            Dscn_hru_id(i) = 0
          ENDIF
          IF ( Dscn_frac_hru(i)>1.0 ) Dscn_frac_hru(i)=1.0
        ENDIF
        f = Infstor_hru_id(i)
        IF ( f>0 ) THEN
          Infstor_frac_hru(i) = Hru_infstor(i)/infstor_area(f)
          IF ( Infstor_frac_hru(i)<NEARZERO ) THEN
            Infstor_frac_hru(i) = 0.0
            Infstor_hru_id(i) = 0
          ENDIF
          IF ( Infstor_frac_hru(i)>1.0 ) Infstor_frac_hru(i)=1.0
        ELSE
          Imperv_stor_to_infstor(i) = 0.0
          Dscn_to_infstor(i) = 0.0
        ENDIF        
      ENDDO    
      
      IF ( Parameter_check_flag>0 ) THEN
        DO i = 1, Nhru        
          CALL check_param_value(i, 'imperv_stor_seep', Imperv_stor_seep(i), Inputerror_flag)
          CALL check_param_value(i, 'imperv_stor_to_stdrn', Imperv_stor_to_stdrn(i), Inputerror_flag)
          CALL check_param_limits(i, 'stdrn_cond', Stdrn_cond(i), 0.0, 1000.0, Inputerror_flag)
          IF ( Ndscn>0 ) THEN
            CALL check_param_value(i, 'dscn_evap_frac', Dscn_evap_frac(i), Inputerror_flag)
            CALL check_param_value(i, 'dscn_to_stdrn', Dscn_to_stdrn(i), Inputerror_flag)          
            IF ( Ninfstor>0 ) THEN
              CALL check_param_value(i, 'dscn_to_infstor', Dscn_to_infstor(i), Inputerror_flag)
            ENDIF
          ENDIF
          IF ( Ninfstor>0 ) THEN
            CALL check_param_value(i, 'imperv_stor_to_infstor', Imperv_stor_to_infstor(i), Inputerror_flag)
            CALL check_param_value(i, 'infstor_seep_coef', Infstor_seep_coef(i), Inputerror_flag)
            CALL check_param_value(i, 'infstor_to_stdrn_coef', Infstor_to_stdrn_coef(i), Inputerror_flag)              
          ENDIF
        ENDDO
        DO i = 1, Ndscn
          CALL check_param_limits(i, 'dscn_stor_max', Dscn_stor_max(i), 0.0, 40.0, Inputerror_flag)  
        ENDDO
        DO i = 1, Ninfstor
          CALL check_param_limits(i, 'infstor_max', Infstor_max(i), 0.0, 100.0, Inputerror_flag)  
        ENDDO        
      ENDIF
                  
      num_hrus = 0
      DO i = 1, Nhru
        frac_sum = Dscn_to_stdrn(i) + Dscn_to_infstor(i)
        IF ( frac_sum>1.0 ) num_hrus = num_hrus + 1
        IF ( Dscn_to_infstor(i)>=1.0 ) THEN
          Dscn_to_infstor(i) = 1.0
          Dscn_to_stdrn(i) = 0.0
        ELSEIF ( (Dscn_to_stdrn(i)+Dscn_to_infstor(i))>=1.0 ) THEN
          Dscn_to_stdrn(i) = 1.0-Dscn_to_infstor(i)
        ENDIF        
      ENDDO          
      IF ( num_hrus>0 .AND. Print_debug>-1 ) THEN
        WRITE (*, '(/,A,/,9X,A,/,9X,A,I7,/,9X,A,/,9X,A,/)') &
     &         'WARNING, fractions partitioning disconnected area inputs is over', &
     &         'specified: dscn_to_stdrn+dscn_to_infstor>1.0', &
     &         'number of HRUs for which this condition exists:', num_hrus, &
     &         'Specified fractions will be maintained in the priority:', &
     &         'dscn_to_infstor, dscn_to_stdrn'
      ENDIF      
      
      num_hrus = 0
      DO i = 1, Nhru
        frac_sum = Imperv_stor_to_stdrn(i) + Imperv_stor_to_infstor(i)
        IF ( frac_sum>1.0 ) THEN
          num_hrus = num_hrus + 1
          IF ( Imperv_stor_to_infstor(i)>=1.0 ) THEN
            Imperv_stor_to_infstor(i) = 1.0
            Imperv_stor_to_stdrn(i) = 0.0
          ELSE
            Imperv_stor_to_stdrn(i) = 1.0-Imperv_stor_to_infstor(i)
          ENDIF
        ENDIF       
      ENDDO          
      IF ( num_hrus>0 .AND. Print_debug>-1 ) THEN
        WRITE (*, '(/,A,/,9X,A,/,9X,A,I7,/,9X,A,/,9X,A,/)') &
     &         'WARNING, fractions partitioning storm drain inputs is over', &
     &         'specified: imperv_stor_to_stdrn+imperv_stor_to_infstor>1.0', &
     &         'number of HRUs for which this condition exists:', num_hrus, &
     &         'Specified fractions will be maintained in the priority:', &
     &         'imperv_stor_to_infstor, imperv_stor_to_stdrn'
      ENDIF 
      
      num_hrus = 0
      DO i = 1, Nhru
        IF ( Imperv_stor_seep(i)>1.0 ) THEN
          num_hrus = num_hrus + 1
          Imperv_stor_seep(i)=1.0
        ENDIF
      ENDDO
      IF ( num_hrus>0 .AND. Print_debug>-1 ) THEN
        WRITE (*, '(/,A,/,9X,A,I7,/,9X,A,/)') &
     &         'WARNING, impervious seepage coefficient is great than 1.0.', &
     &         'number of HRUs for which this condition exists:', num_hrus, &
     &         'These disconnected fractions will be set to 1.0.'
      ENDIF
      
      num_hrus = 0
      DO i = 1, Nhru
        IF ( Stdrn_invert(i)>Hru_elev(i) ) THEN
          num_hrus = num_hrus + 1
          Stdrn_invert(i) = Hru_elev(i)
        ENDIF
      ENDDO
      IF ( num_hrus>0 .AND. Print_debug>-1 ) THEN
        WRITE (*, '(/,A,/,9X,A,I7,/,9X,A,/)') &
     &         'WARNING: storm drainage invert greater than HRU elevation.', &
     &         'number of HRUs for which this condition exists:', num_hrus, &
     &         'These invert elevations will be set to their HRU elevations.'
      ENDIF
      
      END FUNCTION urbaninit
    
!***********************************************************************
!     urban_restart - write or read srunoff restart file
!***********************************************************************
      SUBROUTINE urban_restart(In_out)
      USE PRMS_MODULE, ONLY: Ndscn, Ninfstor, Restart_outunit, Restart_inunit
      USE PRMS_URBAN
      IMPLICIT NONE
      ! Argument
      INTEGER, INTENT(IN) :: In_out
      EXTERNAL check_restart, check_restart_dimen
      ! Local Variable
      CHARACTER(LEN=13) :: module_name
      INTEGER :: Ndscn_test, Ninfstor_test, ierr
!***********************************************************************
      IF ( In_out==0 ) THEN
        WRITE ( Restart_outunit ) MODNAME
        WRITE ( Restart_outunit ) Ndscn, Ninfstor, Basin_dscn_stor, Basin_dscn_evap, Basin_infstor, &
     &                            Basin_infstor_seep, Basin_stdrn_seep,  Basin_stdrn_infil, Basin_urbanfarflow, &
     &                            Basin_urban_to_ssr, Basin_urban_to_infstor
        WRITE ( Restart_outunit ) Dscn_stor
        WRITE ( Restart_outunit ) Dscn_stor_evap
        WRITE ( Restart_outunit ) Inf_stor
        WRITE ( Restart_outunit ) Infstor_seep
        WRITE ( Restart_outunit ) Stdrn_seep
        WRITE ( Restart_outunit ) Urban_to_infstor
        WRITE ( Restart_outunit ) Urban_to_ssr
        WRITE ( Restart_outunit ) Urban_to_stdrn
        WRITE ( Restart_outunit ) Urban_to_farflow
        WRITE ( Restart_outunit ) Urban_to_strm_seg
        WRITE ( Restart_outunit ) Urban_to_lake
      ELSE
        ierr = 0
        READ ( Restart_inunit ) module_name
        CALL check_restart(MODNAME, module_name)
        READ ( Restart_inunit ) Ndscn, Ninfstor, Basin_dscn_stor, Basin_dscn_evap, Basin_infstor, &
     &                          Basin_infstor_seep, Basin_stdrn_seep, Basin_stdrn_infil, Basin_urbanfarflow, &
     &                          Basin_urban_to_ssr, Basin_urban_to_infstor
        CALL check_restart_dimen('ndscn', Ndscn_test, Ndscn, ierr)
        CALL check_restart_dimen('ninfstor', Ninfstor_test, Ninfstor, ierr)
        IF ( ierr==1 ) STOP 
        READ ( Restart_inunit ) Dscn_stor
        READ ( Restart_inunit ) Dscn_stor_evap
        READ ( Restart_inunit ) Inf_stor
        READ ( Restart_inunit ) Infstor_seep
        READ ( Restart_inunit ) Stdrn_seep
        READ ( Restart_inunit ) Urban_to_infstor
        READ ( Restart_inunit ) Urban_to_ssr
        READ ( Restart_inunit ) Urban_to_stdrn
        READ ( Restart_inunit ) Urban_to_farflow
        READ ( Restart_inunit ) Urban_to_strm_seg
        READ ( Restart_inunit ) Urban_to_lake
      ENDIF
      END SUBROUTINE urban_restart