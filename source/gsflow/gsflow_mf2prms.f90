!***********************************************************************
! Transfer MODFLOW data to PRMS
!***********************************************************************
!     ******************************************************************
!     Mapping module to convert MODFLOW to PRMS states for use by GSFLOW
!   Declared Parameters
!     gvr_hru_id, gvr_cell_id
!     ******************************************************************
      INTEGER FUNCTION gsflow_mf2prms()
      USE GSFMODFLOW, ONLY: Mfq2inch_conv, Gwc_col, Gwc_row
      USE PRMS_SOILZONE, ONLY: Hrucheck, Gvr_hru_id, Gw2sm_grav, Gw2sm_grav_save
      USE GWFUZFMODULE, ONLY: SEEPOUT
      USE PRMS_MODULE, ONLY: Process, Nhrucell, Gvr_cell_id
      IMPLICIT NONE
! Functions
      EXTERNAL print_module
! Local Variables
      INTEGER :: i, irow, icol, icell
!      CHARACTER(LEN=14) :: MODNAME
! Save Variables
      CHARACTER(LEN=80), SAVE :: Version_gsflow_mf2prms
!***********************************************************************
      gsflow_mf2prms = 0

      IF ( Process(:3)=='run' ) THEN
        DO i = 1, Nhrucell
          IF ( Hrucheck(Gvr_hru_id(i))==1 ) THEN
            icell = Gvr_cell_id(i)
            irow = Gwc_row(icell)
            icol = Gwc_col(icell)
            Gw2sm_grav_save(i) = Gw2sm_grav(i)
            Gw2sm_grav(i) = SEEPOUT(icol, irow)*Mfq2inch_conv(i)
          ENDIF
        ENDDO

      ELSEIF ( Process(:4)=='decl' ) THEN
        Version_gsflow_mf2prms = 'gsflow_mf2prms.f90 2016-03-11 20:42:32Z'
        CALL print_module(Version_gsflow_mf2prms, 'GSFLOW MODFLOW to PRMS      ', 90)
!        MODNAME = 'gsflow_mf2prms'
      ENDIF

      END FUNCTION gsflow_mf2prms
