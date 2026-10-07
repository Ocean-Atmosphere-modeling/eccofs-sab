      SUBROUTINE checkdefs
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This subroutine checks activated C-preprocessing options for        !
!  consistency.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
      USE mod_strings
!
      USE strings_mod, ONLY : uppercase
!
      implicit none
!
!  Local variable declarations.
!
      integer :: iatms, ibbl, ibiology, idriver
      integer :: ivelHadv, ivelVadv, ivmix, nearshore
      integer :: is, lstr, ng
!
!-----------------------------------------------------------------------
!  Initialze counters
!-----------------------------------------------------------------------
!
      iatms = 0
      ibbl = 0
      ibiology = 0
      idriver = 0
      ivelHadv = 0
      ivelVadv = 0
      ivmix = 0
      nearshore = 0
!
!-----------------------------------------------------------------------
!  Report activated C-preprocessing options.
!-----------------------------------------------------------------------
!
      Coptions=' '
      IF (Master) WRITE (stdout,10)
  10  FORMAT (/,' Activated C-preprocessing Options:',/)
  20  FORMAT (1x,a,t36,a)
!
      IF (Master) THEN
        WRITE (stdout,20) TRIM(ADJUSTL(MyAppCPP)), TRIM(ADJUSTL(title))
      END IF
      is=LEN_TRIM(Coptions)+1
      lstr=LEN_TRIM(MyAppCPP)
      Coptions(is:is+lstr)=TRIM(ADJUSTL(MyAppCPP))
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is)=','
!
      IF (Master) WRITE (stdout,20) 'ANA_BSFLUX',                       &
     &   'Analytical kinematic bottom salinity flux'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_BSFLUX,'
!
      IF (Master) WRITE (stdout,20) 'ANA_BTFLUX',                       &
     &   'Analytical kinematic bottom temperature flux'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_BTFLUX,'
!
      IF (Master) WRITE (stdout,20) 'ASSUMED_SHAPE',                    &
     &   'Using assumed-shape arrays'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' ASSUMED_SHAPE,'
!
      IF (Master) WRITE (stdout,20) 'ATM_PRESS',                        &
     &   'Impose inverse barometer effect in pressure gradient term'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' ATM_PRESS,'
!
      IF (Master) WRITE (stdout,20) 'AVERAGES',                         &
     &   'Writing out time-averaged nonlinear model fields'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' AVERAGES,'
      IF (Master) WRITE (stdout,20) 'BOUNDARY_ALLREDUCE',               &
     &   'Using mpi_allreduce in mp_boundary routine'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+20)=' BOUNDARY_ALLREDUCE,'
!
      IF (Master) WRITE (stdout,20) 'BULK_FLUXES',                      &
     &   'Surface bulk fluxes parameterization'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' BULK_FLUXES,'
!
      IF (Master) WRITE (stdout,20) 'COLLECT_ALLREDUCE',                &
     &   'Using mpi_allreduce in mp_collect routine'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+19)=' COLLECT_ALLGATHER,'
!
      IF (Master) WRITE (stdout,20) 'CHARNOK',                          &
     &   'Charnock surface roughness from wind stress'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' CHARNOK,'
!
      IF (Master) WRITE (stdout,20) 'CRAIG_BANNER',                     &
     &   'Craig and Banner wave breaking surface flux'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' CRAIG_BANNER,'
!
      IF (Master) WRITE (stdout,20) 'COOL_SKIN',                        &
     &   'Surface cool skin correction'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' COOL_SKIN,'
!
      IF (Master) WRITE (stdout,20) 'CURVGRID',                         &
     &   'Orthogonal curvilinear grid'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' CURVGRID,'
!
      IF (Master) WRITE (stdout,20) 'DJ_GRADPS',                        &
     &   'Parabolic Splines density Jacobian (Shchepetkin, 2002)'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' DJ_GRADPS,'
!
      IF (Master) WRITE (stdout,20) 'EMINUSP',                          &
     &   'Compute Salt Flux using E-P'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' EMINUSP,'
      IF (Master) WRITE (stdout,20) '!GATHER_SENDRECV',                 &
     &   'Using mpi_gatherv in mp_gather2d/mp_gather3d routines'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' !GATHER_SENDRECV,'
!
      IF (Master) WRITE (stdout,20) 'GLS_MIXING',                       &
     &   'Generic Length-Scale turbulence closure'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' GLS_MIXING,'
      ivmix=ivmix+1
!
      IF (Master) WRITE (stdout,20) 'LIMIT_STFLX_COOLING',              &
     &   'Suppress further cooling if SST is at freezing point'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+21)=' LIMIT_STFLX_COOLING,'
!
      IF (Master) WRITE (stdout,20) 'KANTHA_CLAYSON',                   &
     &   'Kantha and Clayson stability function formulation'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' KANTHA_CLAYSON,'
!
      IF (Master) WRITE (stdout,20) 'LONGWAVE_OUT',                     &
     &   'Compute outgoing longwave radiation internally'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' LONGWAVE_OUT,'
!
      IF (Master) WRITE (stdout,20) 'MASKING',                          &
     &   'Land/Sea masking'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' MASKING,'
!
      IF (Master) WRITE (stdout,20) 'MIX_GEO_TS',                       &
     &   'Mixing of tracers along geopotential surfaces'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' MIX_GEO_TS,'
!
      IF (Master) WRITE (stdout,20) 'MIX_S_UV',                         &
     &   'Mixing of momentum along constant S-surfaces'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' MIX_S_UV,'
!
      IF (Master) WRITE (stdout,20) 'MPI',                              &
     &   'MPI distributed-memory configuration'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+4)=' MPI,'
!
      IF (Master) WRITE (stdout,20) '!MULTI_THREAD',                    &
     &   'Using mpi_init with single thread level processing'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' !MULTIPLE_THREAD,'
!
      IF (Master) WRITE (stdout,20) 'NONLINEAR',                        &
     &   'Nonlinear Model'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' NONLINEAR,'
!
      IF (Master) WRITE (stdout,20) 'NONLIN_EOS',                       &
     &   'Nonlinear Equation of State for seawater'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' NONLIN_EOS,'
!
      IF (Master) WRITE (stdout,20) 'NO_LBC_ATT',                       &
     &   'Not checking NetCDF global attribute NLM_LBC during restart'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' NO_LBC_ATT,'
!
      IF (Master) WRITE (stdout,20) 'N2S2_HORAVG',                      &
     &   'Horizontal smoothing of buoyancy and shear'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' N2S2_HORAVG,'
!
      IF (Master) WRITE (stdout,20) 'OMEGA_IMPLICIT',                   &
     &   'Adaptive, Courant-number based implicit vertical advection'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' OMEGA_IMPLICIT,'
!
      IF (Master) WRITE (stdout,20) 'PIO_LIB',                          &
     &   'Using Parallel-IO from the PIO library'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' PIO_LIB,'
!
      IF (Master) WRITE (stdout,20) 'PRESS_COMPENSATE',                 &
     &   'Compensate for boundary data without ATM pressure effect'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' PRESS_COMPENSATE,'
!
      IF (Master) WRITE (stdout,20) 'POWER_LAW',                        &
     &   'Power-law shape time-averaging barotropic filter'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' POWER_LAW,'
!
      IF (Master) WRITE (stdout,20) 'PROFILE',                          &
     &   'Time profiling activated'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' PROFILE,'
!
      IF (Master) WRITE (stdout,20) 'K_GSCHEME',                        &
     &   'Third-order upstream advection of TKE fields'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' K_GSCHEME,'
!
      IF (Master) WRITE (stdout,20) 'REDUCE_ALLREDUCE',                 &
     &   'Using mpi_allreduce in mp_reduce routine'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' REDUCE_ALLREDUCE,'
!
      IF (Master) WRITE (stdout,20) 'RST_SINGLE',                       &
     &   'Single precision fields in restart NetCDF file'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' RST_SINGLE,'
!
      IF (Master) WRITE (stdout,20) 'SALINITY',                         &
     &   'Using salinity'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' SALINITY,'
      IF (Master) WRITE (stdout,20) '!SCATTER_BCAST',                   &
     &   'Using mpi_scatterv in mp_scatter2d/mp_scatter3d routines'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' !SCATTER_BCAST,'
!
      IF (Master) WRITE (stdout,20) 'STEP2D_LF_AM3',                    &
     &   'Predictor/Corrector LF-AM3 stepping scheme'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' STEP2D_LF_AM3,'
!
      IF (Master) WRITE (stdout,20) 'SINGLE_PRECISION',                 &
     &   'Single precision arithmetic numerical kernel'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' SINGLE_PRECISION,'
!
      IF (Master) WRITE (stdout,20) 'SOLAR_SOURCE',                     &
     &   'Solar Radiation Source Term'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' SOLAR_SOURCE,'
!
      IF (Master) WRITE (stdout,20) 'SOLVE3D',                          &
     &   'Solving 3D Primitive Equations'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' SOLVE3D,'
!
      IF (Master) WRITE (stdout,20) 'TS_DIF2',                          &
     &   'Harmonic mixing of tracers'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' TS_DIF2,'
!
      IF (Master) WRITE (stdout,20) 'UV_ADV',                           &
     &   'Advection of momentum'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' UV_ADV,'
!
      IF (Master) WRITE (stdout,20) 'UV_COR',                           &
     &   'Coriolis term'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' UV_COR,'
!
      IF (Master) WRITE (stdout,20) 'UV_U3HADVECTION',                  &
     &   'Third-order upstream horizontal advection of 3D momentum'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' UV_U3HADVECTION,'
      ivelHadv=ivelHadv+1
!
      IF (Master) WRITE (stdout,20) 'UV_C4VADVECTION',                  &
     &   'Fourth-order centered vertical advection of momentum'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' UV_C4VADVECTION,'
      ivelVadv=ivelVadv+1
!
      IF (Master) WRITE (stdout,20) 'UV_QDRAG',                         &
     &   'Quadratic bottom stress'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' UV_QDRAG,'
      ibbl=ibbl+1
!
      IF (Master) WRITE (stdout,20) 'UV_VIS2',                          &
     &   'Harmonic mixing of momentum'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' UV_VIS2,'
!
      IF (Master) WRITE (stdout,20) 'VAR_RHO_2D',                       &
     &   'Variable density barotropic mode'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' VAR_RHO_2D,'
!
      IF (Master) WRITE (stdout,20) 'WIND_MINUS_CURRENT',               &
     &   'Compute effective wind by removing surface ocean current'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+20)=' WIND_MINUS_CURRENT,'
!
      IF (Master) WRITE (stdout,20) 'WTYPE_GRID',                       &
     &   'Spatially varying Jerlov water type index'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' WTYPE_GRID,'
!
!-----------------------------------------------------------------------
!  Stop if unsupported C-preprocessing options or report issues with
!  particular options.
!-----------------------------------------------------------------------
!
      CALL checkadj
!
!-----------------------------------------------------------------------
!  Check C-preprocessing options.
!-----------------------------------------------------------------------
!
!  Stop if more than one vertical closure scheme is selected.
!
      IF (Master.and.(ivmix.gt.1)) THEN
        WRITE (stdout,30)
  30    FORMAT (/,' CHECKDEFS - only one vertical closure scheme',      &
     &            ' is allowed.')
        exit_flag=5
      END IF
!
!  Stop if more that one bottom stress formulation is selected.
!
      IF (Master.and.(ibbl.gt.1)) THEN
        WRITE (stdout,40)
  40    FORMAT (/,' CHECKDEFS - only one bottom stress formulation is', &
     &            ' allowed.')
        exit_flag=5
      END IF
!
!  Stop if no bottom stress formulation is selected.
!
      IF (Master.and.(ibbl.eq.0)) THEN
        WRITE (stdout,50)
  50    FORMAT (/,' CHECKDEFS - no bottom stress formulation is',       &
     &            ' selected.')
        exit_flag=5
      END IF
!
!  Stop if more than one biological module is selected.
!
      IF (Master.and.(ibiology.gt.1)) THEN
        WRITE (stdout,60)
  60    FORMAT (/,' CHECKDEFS - only one biology MODULE is allowed.')
        exit_flag=5
      END IF
!
!  Stop if more that one model driver is selected.
!
      IF (Master.and.(idriver.gt.1)) THEN
        WRITE (stdout,70)
  70    FORMAT (/,' CHECKDEFS - only one model example is allowed.')
        exit_flag=5
      END IF
!
!  Stop if more than one advection scheme has been activated.
!
      IF (Master.and.(ivelHadv.gt.1)) THEN
        WRITE (stdout,140) 'horizontal','momentum','ivelHadv =',ivelHadv
        exit_flag=5
      END IF
      IF (Master.and.(ivelVadv.gt.1)) THEN
        WRITE (stdout,140) 'vertical','momentum','ivelVadv =',ivelVadv
        exit_flag=5
      END IF
 140  FORMAT (/,' CHECKDEFS - only one ',a,' advection scheme',         &
     &        /,13x,'is allowed for ',a,', ',a,1x,i1)
!
      RETURN
      END SUBROUTINE checkdefs
