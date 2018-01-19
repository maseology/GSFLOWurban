!***********************************************************************
! Sums inflow to and outflow from PRMS ground-water reservoirs; outflow
! can be routed to downslope ground-water reservoirs and stream
! segments
!
! Can be used for depression storage
!
! mm modfications: Added an option to employ TOPMODEL. Intended to be
!                  used only in fully-distributed (ie, grid-based) mode
!***********************************************************************
! Modified 7/1997 J. Vaccaro to set a minimum value for groundwater flow
! by reading in a minimum ground-water storage value for each groundwater
! reservoir, if this value is set=0, then standard PRMS routine module.
! A minimum may represent an injection well, intrabasin transfer,
! contribution from larger regional gw system, or past residual storage
! modified 10/1/2008 rsregan to include Vaccaro code
!***********************************************************************
      MODULE PRMS_GWFLOW
      IMPLICIT NONE
!   Local Variables
      CHARACTER(LEN=11), SAVE :: MODNAME                                                                                    ! mm
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Gwstor_minarea(:), Gwin_dprst(:)
      DOUBLE PRECISION, SAVE :: Basin_gw_upslope, Basin_gwfarflow
      INTEGER, SAVE :: Gwminarea_flag, Ntop                                                                                 ! mm
      DOUBLE PRECISION, SAVE :: Basin_dnflow
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Zmean(:), Gwr_z(:), Topbasin_area(:)                                           ! mm
      REAL, SAVE, ALLOCATABLE :: g(:)                                                                                       ! mm
      INTEGER, SAVE, ALLOCATABLE :: Topbasin_gwrs(:)                                                                        ! mm
!   Declared Variables
      DOUBLE PRECISION, SAVE :: Basin_gwstor, Basin_gwflow, Basin_gwsink
      DOUBLE PRECISION, SAVE :: Basin_gwin
      DOUBLE PRECISION, SAVE :: Basin_gwstor_minarea_wb
      REAL, SAVE, ALLOCATABLE :: Gwres_flow(:), Gwres_sink(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Gw_upslope(:), Gwres_in(:)
      REAL, SAVE, ALLOCATABLE :: Hru_gw_cascadeflow(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Gw_in_soil(:), Gw_in_ssr(:), Hru_storage(:), Hru_lateral_flow(:)
      DOUBLE PRECISION, SAVE, ALLOCATABLE :: Gwstor_minarea_wb(:), Hru_streamflow_out(:)
      REAL, SAVE, ALLOCATABLE :: Topindex(:)                                                                                ! mm
!   Declared Parameters
      REAL, SAVE, ALLOCATABLE :: Gwflow_coef(:), Gwsink_coef(:)
      REAL, SAVE, ALLOCATABLE :: Gwstor_init(:), Gwstor_min(:)
      REAL, SAVE, ALLOCATABLE :: Lake_seep_elev(:), Elevlake_init(:), Gw_seep_coef(:)
      REAL, SAVE, ALLOCATABLE :: Gwtop_k(:), Gwtop_n(:)                                                                     ! mm
      INTEGER, SAVE, ALLOCATABLE :: Hru_topbasin(:)                                                                         ! mm
      END MODULE PRMS_GWFLOW

!***********************************************************************
!     Main gwflow routine
!***********************************************************************
      INTEGER FUNCTION gwflow()
      USE PRMS_MODULE, ONLY: Process, Save_vars_to_file, Init_vars_from_file
      IMPLICIT NONE
! Functions
      INTEGER, EXTERNAL :: gwflowdecl, gwflowinit, gwflowrun
      EXTERNAL gwflow_restart
!***********************************************************************
      gwflow = 0

      IF ( Process(:3)=='run' ) THEN
        gwflow = gwflowrun()
      ELSEIF ( Process(:4)=='decl' ) THEN
        gwflow = gwflowdecl()
      ELSEIF ( Process(:4)=='init' ) THEN
        IF ( Init_vars_from_file==1 ) CALL gwflow_restart(1)
        gwflow = gwflowinit()
      ELSEIF ( Process(:5)=='clean' ) THEN
        IF ( Save_vars_to_file==1 ) CALL gwflow_restart(0)
      ENDIF

      END FUNCTION gwflow

!***********************************************************************
!     gwflowdecl - set up parameters for groundwater computations
!   Declared Parameters
!     gwstor_init, gwflow_coef, gwsink_coef
!***********************************************************************
      INTEGER FUNCTION gwflowdecl()
      USE PRMS_GWFLOW
      USE PRMS_MODULE, ONLY: Nhru, Ngw, Model, Dprst_flag, Gw_flag, &                                                       ! mm
     &    Cascadegw_flag, Init_vars_from_file
      IMPLICIT NONE
! Functions
      INTRINSIC INDEX
      INTEGER, EXTERNAL :: declparam, declvar
      EXTERNAL read_error, print_module
! Local Variables
      CHARACTER(LEN=80), SAVE :: Version_gwflow
!***********************************************************************
      gwflowdecl = 0

      Version_gwflow = 'gwflow.f90 2016-06-01 10:10:00Z'
      IF ( Gw_flag==1 ) THEN                                                                                                ! mm begin
        MODNAME = 'gwflow'
      ELSEIF ( Gw_flag==2 ) THEN
        Version_gwflow = Version_gwflow(:11)//'_topmodel 2016-12-05Z mm $'
        MODNAME = 'gw_topmodel'
      ENDIF
      CALL print_module(Version_gwflow, 'Groundwater                 ', 90)                                                 ! mm end
  
! cascading variables and parameters
      IF ( Cascadegw_flag>0 .OR. Model==99 ) THEN
        ALLOCATE ( Gw_upslope(Ngw) )
        IF ( declvar(MODNAME, 'gw_upslope', 'ngw', Ngw, 'double', &
     &       'Groundwater flow received from upslope GWRs for each GWR', &
     &       'acre-inches', Gw_upslope)/=0 ) CALL read_error(3, 'gw_upslope')

        ALLOCATE ( Hru_gw_cascadeflow(Ngw) )
        IF ( declvar(MODNAME, 'hru_gw_cascadeflow', 'ngw', Ngw, 'real', &
     &       'Cascading groundwater flow from each GWR', &
     &       'inches', Hru_gw_cascadeflow)/=0 ) CALL read_error(3, 'hru_gw_cascadeflow')
      ENDIF

      ALLOCATE ( Gwres_flow(Ngw) )
      IF ( declvar(MODNAME, 'gwres_flow', 'ngw', Ngw, 'real', &
     &     'Groundwater discharge from each GWR to the stream network', &
     &     'inches', Gwres_flow)/=0 ) CALL read_error(3, 'gwres_flow')

      ALLOCATE ( Gwres_in(Ngw) )
      IF ( declvar(MODNAME, 'gwres_in', 'ngw', Ngw, 'double', &
     &     'Total inflow to each GWR from associated capillary and gravity reservoirs', &
     &     'acre-inches', Gwres_in)/=0 ) CALL read_error(3, 'gwres_in')

      ALLOCATE ( Gwres_sink(Ngw) )
      IF ( declvar(MODNAME, 'gwres_sink', 'ngw', Ngw, 'real', &
     &     'Outflow from GWRs to the groundwater sink; water is'// &
     &     ' considered underflow or flow to deep aquifers and does'// &
     &     ' not flow to the stream network', &
     &     'inches', Gwres_sink)/=0 ) CALL read_error(3, 'gwres_sink')

      ALLOCATE ( Gw_in_soil(Ngw) )
      IF ( declvar(MODNAME, 'gw_in_soil', 'ngw', Ngw, 'double', &
     &     'Drainage from capillary reservoir excess water for each GWR', &
     &     'acre-inches', Gw_in_soil)/=0 ) CALL read_error(3, 'gw_in_soil')

      ALLOCATE ( Gw_in_ssr(Ngw) )
      IF ( declvar(MODNAME, 'gw_in_ssr', 'ngw', Ngw, 'double', &
     &     'Drainage from gravity reservoir excess water for each GWR', &
     &     'acre-inches', Gw_in_ssr)/=0 ) CALL read_error(3, 'gw_in_ssr')

      IF ( declvar(MODNAME, 'basin_gwstor', 'one', 1, 'double', &
     &     'Basin area-weighted average of storage in GWRs', &
     &     'inches', Basin_gwstor)/=0 ) CALL read_error(3, 'basin_gwstor')

      IF ( declvar(MODNAME, 'basin_gwin', 'one', 1, 'double', &
     &     'Basin area-weighted average of inflow to GWRs', &
     &     'inches', Basin_gwin)/=0 ) CALL read_error(3, 'basin_gwin')

      IF ( declvar(MODNAME, 'basin_gwflow', 'one', 1, 'double', &
     &     'Basin area-weighted average of groundwater flow to the stream network', &
     &     'inches', Basin_gwflow)/=0 ) CALL read_error(3, 'basin_gwflow')

      IF ( declvar(MODNAME, 'basin_gwsink', 'one', 1, 'double', &
     &     'Basin area-weighted average of GWR outflow to the groundwater sink', &
     &     'inches', Basin_gwsink)/=0 ) CALL read_error(3, 'basin_gwsink')

      ALLOCATE ( Hru_streamflow_out(Nhru) )
      IF ( declvar(MODNAME, 'hru_streamflow_out', 'nhru', Nhru, 'double', &
     &     'Total flow to stream network from each HRU', &
     &     'cfs', Hru_streamflow_out)/=0 ) CALL read_error(3, 'Hru_streamflow_out')

      ALLOCATE ( Hru_lateral_flow(Nhru) )
      IF ( declvar(MODNAME, 'hru_lateral_flow', 'nhru', Nhru, 'double', &
     &     'Lateral flow to stream network from each HRU', &
     &     'inches', Hru_lateral_flow)/=0 ) CALL read_error(3, 'Hru_lateral_flow')

      ALLOCATE ( Hru_storage(Nhru) )
      IF ( declvar(MODNAME, 'hru_storage', 'nhru', Nhru, 'double', &
     &     'Storage for each HRU', &
     &     'inches', Hru_storage)/=0 ) CALL read_error(3, 'hru_storage')

      ALLOCATE ( Gwstor_minarea(Ngw) )
      IF ( Dprst_flag==1 ) ALLOCATE ( Gwin_dprst(Ngw) )

      IF ( Gw_flag==1 ) THEN                                                                                                ! mm begin
        IF ( Init_vars_from_file==0 ) THEN
          ALLOCATE ( Gwstor_init(Ngw) )
          IF ( declparam(MODNAME, 'gwstor_init', 'ngw', 'real', &
     &         '2.0', '0.0', '10.0', &
     &         'Initial storage in each GWR', &
     &         'Storage in each GWR at the beginning of a simulation', &
     &         'inches')/=0 ) CALL read_error(1, 'gwstor_init')
        ENDIF

        ALLOCATE ( Gwflow_coef(Ngw) )
        IF ( declparam(MODNAME, 'gwflow_coef', 'ngw', 'real', &
     &       '0.015', '0.0', '0.5', &
     &       'Groundwater routing coefficient', &
     &       'Linear coefficient in the equation to compute groundwater discharge for each GWR', &
     &       'fraction/day')/=0 ) CALL read_error(1, 'gwflow_coef')
      ELSEIF ( Gw_flag==2 ) THEN                 
        ALLOCATE ( Gwr_z(Ngw) )
        IF ( declvar(MODNAME, 'hru_z', 'ngw', Ngw, 'double', &
     &       'Depth to groundwater table at each GWR', &
     &       'inches', Gwr_z)/=0 ) CALL read_error(3, 'hru_z')

        ALLOCATE ( Topbasin_area(Ntop) )
        IF ( declvar(MODNAME, 'topbasin_area', 'ntop', Ngw, 'double', &
     &       'Area of TOPMODEL subbasins', &
     &       'acres', Topbasin_area)/=0 ) CALL read_error(3, 'topbasin_area')

        ALLOCATE ( Topindex(Ngw) )
        IF ( declvar(MODNAME, 'topindex', 'ngw', Ngw, 'real', &
     &       'Topographic index at each GWR', &
     &       'none', Topindex)/=0 ) CALL read_error(3, 'topindex')        
        
        IF ( Init_vars_from_file==0 ) THEN
          ALLOCATE ( Gwstor_init(Ngw) )
          IF ( declparam(MODNAME, 'hru_z_init', 'ngw', 'real', &
     &         '1.0', '0.0', '1000.0', &
     &         'Initial depth to watertable at each GWR', &
     &         'Watertable depth at each GWR at the beginning of a simulation', &
     &         'elev_units')/=0 ) CALL read_error(1, 'hru_z_init')
        ENDIF         
        
        ALLOCATE ( Gwflow_coef(Ntop) )
        IF ( declparam(MODNAME, 'topmodel_f', 'ntop', 'real', &
     &       '1.0', '0.0', '1.0', &
     &       'TOPMODEL parameter f', &
     &       'TOPMODEL parameter f', &
     &       'fraction')/=0 ) CALL read_error(1, 'topmodel_f')
        
        ALLOCATE ( Gwtop_k(Ngw) )
        IF ( declparam(MODNAME, 'topmodel_k', 'ngw', 'real', &
     &       '35.0', '0.0000001', '5000000.0', &
     &       'Groundwater conductance', &
     &       'Groundwater conductance for each GWR', &
     &       'inches/day')/=0 ) CALL read_error(1, 'topmodel_k')
        
        ALLOCATE ( Gwtop_n(Ngw) )
        IF ( declparam(MODNAME, 'topmodel_n', 'ngw', 'real', &
     &       '0.3', '0.01', '0.9', &
     &       'Near-surface porosity', &
     &       'Near-surface porosity for each GWR', &
     &       'none')/=0 ) CALL read_error(1, 'topmodel_n')       
        
        ALLOCATE ( Hru_topbasin(Ngw) )
        IF ( declparam(MODNAME, 'hru_topbasin', 'nhru', 'integer', &
     &       '1', 'bounded', 'ntop', &
     &       'Index of TOPMODEL basin assigned to each HRU', &
     &       'Index of TOPMODEL basin assigned to each HRU', &
     &       'none')/=0 ) CALL read_error(1, 'hru_topbasin')        
      ENDIF                                                                                                                 ! mm end

      ALLOCATE ( Gwsink_coef(Ngw) )
      IF ( declparam(MODNAME, 'gwsink_coef', 'ngw', 'real', &
     &     '0.0', '0.0', '0.05', &
     &     'Groundwater sink coefficient', &
     &     'Linear coefficient in the equation to compute outflow'// &
     &     ' to the groundwater sink for each GWR', &
     &     'fraction/day')/=0 ) CALL read_error(1, 'gwsink_coef')

      ALLOCATE ( Gwstor_min(Ngw) )
      IF ( declparam(MODNAME, 'gwstor_min', 'ngw', 'real', &
     &     '0.0', '0.0', '1.0', &
     &     'Minimum storage in each GWR', &
     &     'Minimum storage in each GWR to ensure storage is greater'// &
     &     ' than specified value to account for inflow from deep'// &
     &     ' aquifers or injection wells with the water source outside the basin', &
     &     'inches')/=0 ) CALL read_error(1, 'gwstor_min')

      ALLOCATE ( Gwstor_minarea_wb(Ngw) )
      IF ( declvar(MODNAME, 'gwstor_minarea_wb', 'ngw', Ngw, 'double', &
     &     'Storage added to each GWR when storage is less than gwstor_min', &
     &     'inches', Gwstor_minarea_wb)/=0 ) CALL read_error(3, 'gwstor_minarea_wb')

      IF ( declvar(MODNAME, 'basin_gwstor_minarea_wb', 'one', 1, 'double', &
     &     'Basin area-weighted average storage added to each GWR when storage is less than gwstor_min', &
     &     'inches', Basin_gwstor_minarea_wb)/=0 ) CALL read_error(3, 'basin_gwstor_minarea_wbr_minarea_wb')

      END FUNCTION gwflowdecl

!***********************************************************************
!     gwflowinit - Initialize gwflow module - get parameter values,
!                  compute initial values.
!***********************************************************************
      INTEGER FUNCTION gwflowinit()
      USE PRMS_GWFLOW
      USE PRMS_MODULE, ONLY: Ngw, Dprst_flag, Print_debug, Inputerror_flag, &
     &                       Cascadegw_flag, Init_vars_from_file, Gwr_swale_flag, Gw_flag                                   ! mm
      USE PRMS_BASIN, ONLY: Gwr_type, Hru_area, Basin_area_inv, Active_gwrs, Gwr_route_order, &
     &                      Elev_units, METERS2FEET                                                                         ! mm
      USE PRMS_FLOWVARS, ONLY: Gwres_stor, Pkwater_equiv
      USE PRMS_CASCADE, ONLY: Ncascade_gwr, Gwr_down, Gwr_down_frac                                                         ! mm
      USE PRMS_INTCP, ONLY: Hru_intcpstor
      USE PRMS_SRUNOFF, ONLY: Hru_impervstor, Dprst_stor_hru
      USE PRMS_SOILZONE, ONLY: Soil_moist_tot
      IMPLICIT NONE
      INTEGER, EXTERNAL :: getparam
      EXTERNAL read_error
      INTRINSIC DBLE, SNGL, IABS                                                                                            ! mm
! Local Variables
      REAL, ALLOCATABLE :: ca(:), slope(:)                                                                                  ! mm
      REAL :: cw2                                                                                                           ! mm
      INTEGER :: i, j, jjj, k                                                                                               ! mm
!***********************************************************************
      gwflowinit = 0

      IF ( getparam(MODNAME, 'gwsink_coef', Ngw, 'real', Gwsink_coef)/=0 ) CALL read_error(2, 'gwsink_coef')                ! mm begin
      IF ( Gw_flag==1 ) THEN
        IF ( getparam(MODNAME, 'gwflow_coef', Ngw, 'real', Gwflow_coef)/=0 ) CALL read_error(2, 'gwflow_coef')
        IF ( getparam(MODNAME, 'gwstor_min', Ngw, 'real', Gwstor_min)/=0 ) CALL read_error(2, 'gwstor_min')        
      ELSEIF ( Gw_flag==2 ) THEN
        ALLOCATE ( slope(Ngw), ca(Ngw) )
        IF ( getparam(MODNAME, 'topmodel_k', Ngw, 'real', Gwtop_k)/=0 ) CALL read_error(2, 'topmodel_k') 
        IF ( getparam(MODNAME, 'topmodel_n', Ngw, 'real', Gwtop_n)/=0 ) CALL read_error(2, 'topmodel_n') 
        IF ( getparam(MODNAME, 'hru_slope', Ngw, 'real', slope)/=0 ) CALL read_error(2, 'hru_slope')
        IF ( getparam(MODNAME, 'hru_topbasin', Ngw, 'integer', Hru_topbasin)/=0 ) CALL read_error(2, 'hru_topbasin')
        Ntop = 1
        DO j = 1, Active_gwrs
          IF ( Hru_topbasin(j)<1 ) THEN
            IF ( Print_debug>-1 ) PRINT *, 'WARNING, GWR: ', i, ' not set to a valid TOPMODEL basin. Now set to 1'
            Hru_topbasin(j) = 1              
          ELSEIF ( Hru_topbasin(j)>Ntop ) THEN
            Ntop = Hru_topbasin(j)
          ENDIF
        ENDDO
        ALLOCATE ( Zmean(Ntop), g(Ntop), Topbasin_gwrs(Ntop) )
        IF ( getparam(MODNAME, 'topmodel_f', Ntop, 'real', Gwflow_coef)/=0 ) CALL read_error(2, 'topmodel_f')
      ENDIF                                                                                                                 ! mm end

      Gwminarea_flag = 0
      Gwstor_minarea = 0.0D0
      Gwstor_minarea_wb = 0.0D0
      Basin_gwstor_minarea_wb = 0.0D0
      Hru_storage = 0.0D0                                                                                                   ! mm begin
      Basin_gwstor = 0.0D0
      DO i = 1, Active_gwrs
        ! GWR's cannot be swales unless gwr_swale_flag > 0
        IF ( Gwr_type(i)==3 ) THEN ! rsr, may need to add gwr_type and gwr_segment
          IF ( Gwr_swale_flag==0 ) THEN
            PRINT *, 'ERROR, GWRs cannot be swales when gwr_swale_flag = 0, GWR:', i
            Inputerror_flag = 1
          ELSEIF ( Gwr_swale_flag==1 ) THEN
            IF ( Print_debug>-1 ) PRINT *, 'WARNING, GWR:', i, ' is treated as a swale, flow sent to sink'
          ELSEIF ( Gwr_swale_flag==2 ) THEN
            IF ( Print_debug>-1 ) PRINT *, 'WARNING, GWR:', i, &
       &                                   ' is treated as a swale, flow sent to basin_cfs and hru_segment if > 0'
          ELSE
! maybe gwr_swale_flag = 3 abs(hru_segment) so hru_segment could be changed from 0 to allow HRU swales
            PRINT *, 'ERROR, invalid gwr_swale_flag value, specified as:', gwr_swale_flag
            Inputerror_flag = 1
          ENDIF
        ENDIF
      ENDDO
      IF ( Gw_flag==1 ) THEN
        IF ( Init_vars_from_file==0 ) THEN
          IF ( getparam(MODNAME, 'gwstor_init', Ngw, 'real', Gwstor_init)/=0 ) CALL read_error(2, 'gwstor_init')
          DO i = 1, Ngw
            Gwres_stor(i) = DBLE( Gwstor_init(i) )
          ENDDO
          DEALLOCATE ( Gwstor_init )
        ENDIF
        DO j = 1, Active_gwrs
          i = Gwr_route_order(j)
          Basin_gwstor = Basin_gwstor + Gwres_stor(i)*DBLE(Hru_area(i))
          IF ( Gwstor_min(i)>0.0 ) THEN
            Gwminarea_flag = 1
            Gwstor_minarea(i) = DBLE( Gwstor_min(i)*Hru_area(i) )
          ENDIF
          IF ( Gwflow_coef(i)>1.0 ) THEN
            IF ( Print_debug>-1 ) PRINT *, 'WARNING, gwflow_coef value > 1.0 for GWR:', i, Gwflow_coef(i)
          ENDIF

          Hru_storage(i) = DBLE( Soil_moist_tot(i) + Hru_intcpstor(i) + Hru_impervstor(i) ) + Gwres_stor(i) + Pkwater_equiv(i)
          IF ( Dprst_flag==1 ) Hru_storage(i) = Hru_storage(i) + Dprst_stor_hru(i)
        ENDDO
        IF ( Gwminarea_flag==0 ) DEALLOCATE ( Gwstor_minarea )
      ELSEIF ( Gw_flag==2 ) THEN
        IF ( Init_vars_from_file==0 ) THEN
          IF ( getparam(MODNAME, 'hru_z_init', Ngw, 'real', Gwstor_init)/=0 ) CALL read_error(2, 'hru_z_init')
          DO i = 1, Ngw
            cw2 = Gwstor_init(i)*12.0  
            IF ( Elev_units==1 ) cw2 = cw2*METERS2FEET
            Gwr_z(i) = DBLE( cw2 ) ! inches
          ENDDO
          DEALLOCATE ( Gwstor_init )
        ENDIF          
        DO j = 1, Active_gwrs
          ca(j) = Hru_area(j)
        ENDDO
        DO j = 1, Active_gwrs
          i = Gwr_route_order(j)
          DO k = 1, Ncascade_gwr(i)
            jjj = Gwr_down(k, i)
            IF ( jjj>0 ) THEN ! if gwr_down(k, i) > 0, cascade contributes to a downslope GWR
              ca(jjj) = ca(jjj) + Gwr_down_frac(k, i)*ca(i)
            ELSEIF ( jjj<0 ) THEN ! if gwr_down(k, i) < 0, cascade contributes to a stream
              jjj = IABS( jjj )
              ! in this current implimentation, TOPMODEL's groundwater reservoir doesn't contribute directly to streams, only to replenish soilzone reservoirs, which, in turn, will contribute to streams, lakes, and farfields
            ENDIF
          ENDDO 
          Hru_storage(i) = DBLE( Soil_moist_tot(i) + Hru_intcpstor(i) + Hru_impervstor(i) ) + Pkwater_equiv(i)
          IF ( Dprst_flag==1 ) Hru_storage(i) = Hru_storage(i) + Dprst_stor_hru(i)          
        ENDDO
        g = 0.0
        Zmean = 0.0D0
        Topbasin_area = 0.0D0
        Topbasin_gwrs = 0
        DO j = 1, Active_gwrs
          cw2 = Hru_area(j) ! cell width squared: assumes a grid-based HRU distribution (would still work otherwise, but code may require modification)
          Topindex(j) = ca(j) / Gwtop_k(j) / cw2 / slope(j)
          i = Hru_topbasin(j)
          g(i) = g(i) + Topindex(j)
          Topbasin_area(i) = Topbasin_area(i) + Hru_area(j)
          Topbasin_gwrs(i) = Topbasin_gwrs(i) + 1
        ENDDO
        DO j = 1, Ntop
          IF ( Topbasin_area(j)<=0.0 ) THEN
            PRINT *, 'ERROR: TOPMODEL basin indices must be sequential. No GWRs are given the index: ', j
            PRINT *, '       Highest index assigned: ', Ntop
            Inputerror_flag = 1
            CYCLE
          ENDIF
          g(j) = g(j) / Topbasin_area(j)         
        ENDDO
        DO j = 1, Active_gwrs
          i = Hru_topbasin(j)
          Zmean(i) = Zmean(i) + Gwr_z(j) - DBLE(g(i) - Topindex(j) / Gwflow_coef(i))
        ENDDO
        DO i = 1, Ntop
          Zmean(i) = Zmean(i) / Topbasin_gwrs(i)  
          Basin_gwstor = Basin_gwstor + Zmean(i)*Topbasin_area(i)       
        ENDDO
      ENDIF                                                                                                                 ! mm end
      Basin_gwstor = Basin_gwstor*Basin_area_inv

      IF ( Dprst_flag==1 ) Gwin_dprst = 0.0D0

      IF ( Init_vars_from_file==1 ) RETURN

! do only once, so restart uses saved values
      IF ( Cascadegw_flag>0 ) THEN
        Gw_upslope = 0.0D0
        Hru_gw_cascadeflow = 0.0
      ENDIF
      Gwres_flow = 0.0
      Gwres_in = 0.0
      Gwres_sink = 0.0
      Gw_in_ssr = 0.0D0
      Gw_in_soil = 0.0D0
      Basin_gwflow = 0.0D0
      Basin_gwsink = 0.0D0
      Basin_gwin = 0.0D0
      Basin_gwfarflow = 0.0D0
      Basin_gw_upslope = 0.0D0
      Basin_dnflow = 0.0D0
      Hru_streamflow_out = 0.0D0
      Hru_lateral_flow = 0.0D0

      END FUNCTION gwflowinit

!***********************************************************************
!     gwflowrun - Computes groundwater flow to streamflow and to
!                 groundwater sink
!***********************************************************************
      INTEGER FUNCTION gwflowrun()
      USE PRMS_GWFLOW
      USE PRMS_MODULE, ONLY: Dprst_flag, Print_debug, Cascadegw_flag, Gwr_swale_flag, Gw_flag                               ! mm
      USE PRMS_BASIN, ONLY: Active_gwrs, Gwr_route_order, &
     &    Basin_area_inv, Hru_area, Gwr_type, Hru_area_dble
      USE PRMS_FLOWVARS, ONLY: Soil_to_gw, Ssr_to_gw, Sroff, Ssres_flow, Gwres_stor, Pkwater_equiv
      USE PRMS_CASCADE, ONLY: Ncascade_gwr
      USE PRMS_SET_TIME, ONLY: Cfs_conv
      USE PRMS_SRUNOFF, ONLY: Dprst_seep_hru, Hru_impervstor, Dprst_stor_hru
      USE PRMS_INTCP, ONLY: Hru_intcpstor
      USE PRMS_SOILZONE, ONLY: Soil_moist_tot, Gwr_to_ssr                                                                   ! mm
      IMPLICIT NONE
! Functions
      EXTERNAL rungw_cascade, print_date
      INTRINSIC ABS, DBLE, DABS, SNGL
! Local Variables
      INTEGER :: i, j, jj                                                                                                   ! mm
      REAL :: dnflow
      DOUBLE PRECISION :: gwin, gwstor, gwsink, gwflow, gwstor_last, gwarea, far_gwflow, seepage                            ! mm
!***********************************************************************
      gwflowrun = 0

      IF ( Cascadegw_flag>0 ) THEN
        Gw_upslope = 0.0D0
        Basin_dnflow = 0.0D0
        Basin_gw_upslope = 0.0D0
      ENDIF

      Basin_gwstor_minarea_wb = 0.0D0
      Basin_gwflow = 0.0D0
      Basin_gwstor = 0.0D0
      Basin_gwsink = 0.0D0
      Basin_gwin = 0.0D0
      Basin_gwfarflow = 0.0D0
      IF ( Gw_flag==1 ) THEN                                                                                                ! mm begin
        DO j = 1, Active_gwrs
          i = Gwr_route_order(j)
          gwarea = Hru_area_dble(i)
          gwstor = Gwres_stor(i)*gwarea ! acre-inches
          ! soil_to_gw is for whole HRU, not just perv
          Gw_in_soil(i) = Soil_to_gw(i)*Hru_area(i)
          Gw_in_ssr(i) = Ssr_to_gw(i)*Hru_area(i)
          gwin = Gw_in_soil(i) + Gw_in_ssr(i)
          IF ( Cascadegw_flag>0 ) THEN
            gwin = gwin + Gw_upslope(i)
            Basin_gw_upslope = Basin_gw_upslope + Gw_upslope(i)
          ENDIF
          IF ( Dprst_flag==1 ) THEN
            !rsr, need basin variable for WB
            Gwin_dprst(i) = Dprst_seep_hru(i)*gwarea
            gwin = gwin + Gwin_dprst(i)
          ENDIF
          IF ( Gwminarea_flag==1 ) THEN
            ! check to be sure gwres_stor >= gwstor_minarea before computing outflows
            IF ( gwstor<Gwstor_minarea(i) ) THEN
              gwstor_last = gwstor
              gwstor = Gwstor_minarea(i)
              !rsr, keep track of change in storage for WB
              Gwstor_minarea_wb(i) = gwstor - gwstor_last
              Basin_gwstor_minarea_wb = Basin_gwstor_minarea_wb + Gwstor_minarea_wb(i)
              Gwstor_minarea_wb(i) = Gwstor_minarea_wb(i)/gwarea
              IF ( Print_debug>-1 ) PRINT *, 'Added to gwres_stor as storage < gwstor_min to GWR:', i, &
     &                                       ' amount:', Gwstor_minarea_wb(i)
            ELSE
              Gwstor_minarea_wb(i) = 0.0D0
            ENDIF
          ENDIF
          gwstor = gwstor + gwin
          Basin_gwin = Basin_gwin + gwin

! Compute groundwater discharge
          gwflow = gwstor*DBLE( Gwflow_coef(i) )

! Reduce storage by outflow
          gwstor = gwstor - gwflow

          gwsink = 0.0D0
          IF ( Gwsink_coef(i)>0.0 ) THEN
            gwsink = gwstor*DBLE( Gwsink_coef(i) )
            gwstor = gwstor - gwsink
          ENDIF
! if gwr_swale_flag = 1 swale GWR flow goes to sink, 2 included in stream network and cascades
! maybe gwr_swale_flag = 3 abs(hru_segment) so hru_segment could be changed from 0 to allow HRU swales
          IF ( Gwr_swale_flag==1 ) THEN
            IF ( Gwr_type(i)==3 ) THEN
              gwsink = gwsink + gwflow
              gwflow = 0.0D0
            ENDIF
          ENDIF
          Gwres_sink(i) = SNGL( gwsink/gwarea )
          Basin_gwsink = Basin_gwsink + gwsink
          Basin_gwstor = Basin_gwstor + gwstor

          dnflow = 0.0
          Gwres_flow(i) = SNGL( gwflow/gwarea )
          IF ( Cascadegw_flag>0 ) THEN
            IF ( Ncascade_gwr(i)>0 ) THEN
              CALL rungw_cascade(i, Ncascade_gwr(i), Gwres_flow(i), dnflow, far_gwflow)
              Hru_gw_cascadeflow(i) = dnflow + SNGL( far_gwflow )
              Basin_dnflow = Basin_dnflow + dnflow*gwarea
              Basin_gwfarflow = Basin_gwfarflow + far_gwflow*gwarea
            ENDIF
          ENDIF
          Basin_gwflow = Basin_gwflow + DBLE(Gwres_flow(i))*gwarea

          ! leave gwin in inch-acres
          Gwres_in(i) = gwin
          Gwres_stor(i) = gwstor/gwarea
          Hru_lateral_flow(i) = DBLE( Gwres_flow(i) + Sroff(i) + Ssres_flow(i) )
          ! Cfs_conv converts acre-inches per timestep to cfs
          Hru_streamflow_out(i) = gwarea*Cfs_conv*Hru_lateral_flow(i)
          Hru_storage(i) = DBLE( Soil_moist_tot(i) + Hru_intcpstor(i) + Hru_impervstor(i) ) + Gwres_stor(i) &
     &                     + Pkwater_equiv(i)
          IF ( Dprst_flag==1 ) Hru_storage(i) = Hru_storage(i) + Dprst_stor_hru(i)
        ENDDO
      ELSEIF ( Gw_flag==2 ) THEN
        DO j = 1, Active_gwrs
          i = Gwr_route_order(j)
          Gwr_to_ssr(i) = 0.0
          gwarea = Hru_area_dble(i)
          ! soil_to_gw is for whole HRU, not just perv
          Gw_in_soil(i) = Soil_to_gw(i)*gwarea
          Gw_in_ssr(i) = Ssr_to_gw(i)*gwarea
          gwin = Gw_in_soil(i) + Gw_in_ssr(i)
          IF ( Dprst_flag==1 ) THEN
            !rsr, need basin variable for WB
            Gwin_dprst(i) = Dprst_seep_hru(i)*gwarea
            gwin = gwin + Gwin_dprst(i)
          ENDIF
          Basin_gwin = Basin_gwin + gwin
          
! Compute groundwater discharge
          gwflow = 0.0D0
          gwsink = 0.0D0
          IF ( Gwsink_coef(i)>0.0 ) THEN
            gwsink = gwin*DBLE( Gwsink_coef(i) )
            gwin = gwin - gwsink
            IF ( gwin<0.0D0 ) THEN
              gwsink = gwsink + gwin
              gwin = 0.0D0
            ENDIF
            Gwres_sink(i) = SNGL( gwsink/gwarea )
            Basin_gwsink = Basin_gwsink + gwsink
          ENDIF          
          
          Gwres_flow(i) = SNGL( gwflow/gwarea )
          Basin_gwflow = Basin_gwflow + Gwres_flow(i)*gwarea
          
! Compute local watertable increases (depth to watertable decrease)
          jj = Hru_topbasin(i)
          Zmean(jj) = Zmean(jj) - gwin/Topbasin_area(jj)/DBLE(Gwtop_n(i))

          ! leave gwin in inch-acres
          Gwres_in(i) = gwin
          Hru_lateral_flow(i) = gwflow + DBLE( Sroff(i) + Ssres_flow(i) )
          ! Cfs_conv converts acre-inches per timestep to cfs
          Hru_streamflow_out(i) = gwarea*Cfs_conv*Hru_lateral_flow(i)
          Hru_storage(i) = DBLE( Soil_moist_tot(i) + Hru_intcpstor(i) + Hru_impervstor(i) ) + Pkwater_equiv(i)
          IF ( Dprst_flag==1 ) Hru_storage(i) = Hru_storage(i) + Dprst_stor_hru(i)
        ENDDO
        DO j = 1, Active_gwrs
          gwarea = Hru_area_dble(j)
          jj = Hru_topbasin(j)
          Gwr_z(j) = Zmean(jj) + DBLE(g(jj) - Topindex(j) / Gwflow_coef(jj))
          IF ( Gwr_z(j)<0.0 ) THEN ! watertable above surface, use excess to replenish soilzone
            seepage = -Gwr_z(j)*DBLE(Gwtop_n(j)) ! inches
            Gwr_to_ssr(j) = SNGL(seepage)
            Gwr_z(j) = 0.0D0
            Gwres_flow(i) = Gwres_flow(i) + SNGL( seepage/gwarea ) ! gwflow is being used to save groundwater seepage in TOPMODEL mode
            Basin_gwflow = Basin_gwflow + Gwres_flow(i)*gwarea
          ENDIF
        ENDDO
        ! reset average depth to watertables
        Zmean = 0.0D0
        DO j = 1, Active_gwrs
          i = Hru_topbasin(j)
          Zmean(i) = Zmean(i) + Gwr_z(j) - DBLE(g(i) - Topindex(j) / Gwflow_coef(i))
        ENDDO
        DO i = 1, Ntop
          Zmean(i) = Zmean(i) / Topbasin_gwrs(i)  
          Basin_gwstor = Basin_gwstor + Zmean(i)*Topbasin_area(i)       
        ENDDO        
      ENDIF                                                                                                                 ! mm end

      Basin_gwflow = Basin_gwflow*Basin_area_inv
      Basin_gwstor = Basin_gwstor*Basin_area_inv
      Basin_gwsink = Basin_gwsink*Basin_area_inv
      Basin_gwin = Basin_gwin*Basin_area_inv
      Basin_gwfarflow = Basin_gwfarflow*Basin_area_inv
      Basin_gw_upslope = Basin_gw_upslope*Basin_area_inv
      Basin_gwstor_minarea_wb = Basin_gwstor_minarea_wb*Basin_area_inv
      Basin_dnflow = Basin_dnflow*Basin_area_inv

      END FUNCTION gwflowrun

!***********************************************************************
!     Compute cascading GW flow
!***********************************************************************
      SUBROUTINE rungw_cascade(Igwr, Ncascade_gwr, Gwres_flow, Dnflow, Far_gwflow)
      USE PRMS_SRUNOFF, ONLY: Strm_seg_in, Strm_farfield
      USE PRMS_GWFLOW, ONLY: Gw_upslope
      USE PRMS_CASCADE, ONLY: Gwr_down, Gwr_down_frac, Cascade_gwr_area, Nsegmentp1
      ! Cfs_conv converts acre-inches per timestep to cfs
      USE PRMS_SET_TIME, ONLY: Cfs_conv
      IMPLICIT NONE
      INTRINSIC IABS, DBLE, SNGL
! Arguments
      INTEGER, INTENT(IN) :: Igwr, Ncascade_gwr
      REAL, INTENT(INOUT) :: Gwres_flow, Dnflow
      DOUBLE PRECISION, INTENT(OUT) :: Far_gwflow
! Local variables
      INTEGER :: j, k
!***********************************************************************
      Far_gwflow = 0.0D0
      DO k = 1, Ncascade_gwr
        j = Gwr_down(k, Igwr)
        ! Gwres_flow is in inches
! if gwr_down(k, Igwr) > 0, cascade contributes to a downslope GWR
        IF ( j>0 ) THEN
          Gw_upslope(j) = Gw_upslope(j) + Gwres_flow*Cascade_gwr_area(k, Igwr)
          Dnflow = Dnflow + Gwres_flow*Gwr_down_frac(k, Igwr)
! if gwr_down(k, Igwr) < 0, cascade contributes to a stream
        ELSEIF ( j<0 ) THEN
          j = IABS( j )
          IF ( j/=Nsegmentp1 ) THEN
            Strm_seg_in(j) = Strm_seg_in(j) + DBLE( Gwres_flow*Cascade_gwr_area(k, Igwr) )*Cfs_conv
          ELSE
            Strm_farfield = Strm_farfield + DBLE( Gwres_flow*Cascade_gwr_area(k, Igwr) )*Cfs_conv
            Far_gwflow = Far_gwflow + DBLE( Gwres_flow*Gwr_down_frac(k, Igwr) )
          ENDIF
        ENDIF
      ENDDO

      ! gwres_flow reduced by cascading flow to HRUs or farfield
      Gwres_flow = Gwres_flow - Dnflow - SNGL( Far_gwflow )
      IF ( Gwres_flow<0.0 ) Gwres_flow = 0.0

      END SUBROUTINE rungw_cascade

!***********************************************************************
!     gwflow_restart - write or read gwflow restart file
!***********************************************************************
      SUBROUTINE gwflow_restart(In_out)
      USE PRMS_MODULE, ONLY: Restart_outunit, Restart_inunit, Cascadegw_flag
      USE PRMS_GWFLOW
      ! Argument
      INTEGER, INTENT(IN) :: In_out
      EXTERNAL check_restart
      ! Local Variable
      CHARACTER(LEN=6) :: module_name
!***********************************************************************
      IF ( In_out==0 ) THEN
        WRITE ( Restart_outunit ) MODNAME
        WRITE ( Restart_outunit ) Basin_gwstor, Basin_gwflow, Basin_gwsink, Basin_gwin, Basin_gwstor_minarea_wb, &
     &          Gwminarea_flag, Basin_dnflow, Basin_gw_upslope, Basin_gwfarflow
        IF ( Gwminarea_flag==1 ) THEN
          WRITE ( Restart_outunit ) Gwstor_minarea_wb
          WRITE ( Restart_outunit ) Gwstor_minarea
        ENDIF
        WRITE ( Restart_outunit ) Hru_streamflow_out
        WRITE ( Restart_outunit ) Hru_lateral_flow
        WRITE ( Restart_outunit ) Hru_storage
        WRITE ( Restart_outunit ) Gwres_flow
        WRITE ( Restart_outunit ) Gwres_sink
        WRITE ( Restart_outunit ) Gwres_in
        WRITE ( Restart_outunit ) Gw_in_soil
        WRITE ( Restart_outunit ) Gw_in_ssr
        IF ( Cascadegw_flag>0 ) THEN
          WRITE ( Restart_outunit ) Gw_upslope
          WRITE ( Restart_outunit ) Hru_gw_cascadeflow
        ENDIF
        IF ( Gw_flag==2 ) THEN                                                                                              ! mm
          WRITE ( Restart_outunit ) Zmean
          WRITE ( Restart_outunit ) Gwr_z
        ENDIF                                                                                                               ! mm
      ELSE
        READ ( Restart_inunit ) module_name
        CALL check_restart(MODNAME, module_name)
        READ ( Restart_inunit ) Basin_gwstor, Basin_gwflow, Basin_gwsink, Basin_gwin, Basin_gwstor_minarea_wb, &
     &         Gwminarea_flag, Basin_dnflow, Basin_gw_upslope, Basin_gwfarflow
        IF ( Gwminarea_flag==1 ) THEN ! could be error if someone turns off gwstor_min for restart
          READ ( Restart_inunit ) Gwstor_minarea_wb
          READ ( Restart_inunit ) Gwstor_minarea
        ENDIF
        READ ( Restart_inunit ) Hru_streamflow_out
        READ ( Restart_inunit ) Hru_lateral_flow
        READ ( Restart_inunit ) Hru_storage
        READ ( Restart_inunit ) Gwres_flow
        READ ( Restart_inunit ) Gwres_sink
        READ ( Restart_inunit ) Gwres_in
        READ ( Restart_inunit ) Gw_in_soil
        READ ( Restart_inunit ) Gw_in_ssr
        IF ( Cascadegw_flag>0 ) THEN
          READ ( Restart_inunit ) Gw_upslope
          READ ( Restart_inunit ) Hru_gw_cascadeflow
        ENDIF
        IF ( Gw_flag==2 ) THEN                                                                                              ! mm
          WRITE ( Restart_inunit ) Zmean
          WRITE ( Restart_inunit ) Gwr_z
        ENDIF                                                                                                               ! mm
      ENDIF
      END SUBROUTINE gwflow_restart
