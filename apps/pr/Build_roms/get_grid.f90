      MODULE get_grid_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module reads grid information from input file using either     !
!  the standard NetCDF library or the Parallel-IO (PIO) library.       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_netcdf
      USE mod_pio_netcdf
      USE mod_scalars
!
      USE exchange_2d_mod
      USE mp_exchange_mod, ONLY : mp_exchange2d
      USE nf_fread2d_mod,  ONLY : nf_fread2d
      USE strings_mod,     ONLY : FoundError, find_string
!
      implicit none
!
      PUBLIC  :: get_grid
      PRIVATE :: get_grid_nf90
      PRIVATE :: get_grid_pio
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE get_grid (ng, tile, model)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/get_grid.F"
!
!-----------------------------------------------------------------------
!  Read in GRID NetCDF file according to IO type.
!-----------------------------------------------------------------------
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
      SELECT CASE (GRD(ng)%IOtype)
        CASE (io_nf90)
          CALL get_grid_nf90 (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj)
        CASE (io_pio)
          CALL get_grid_pio (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj)
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) GRD(ng)%IOtype
          exit_flag=2
      END SELECT
      IF (FoundError(exit_flag, NoError, 91, MyFile)) RETURN
!
  10  FORMAT (' GET_GRID - Illegal input file type, io_type = ',i0,     &
     &        /,12x,'Check KeyWord ''INP_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE get_grid
!
!***********************************************************************
      SUBROUTINE get_grid_nf90 (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: cr, gtype, i, status, vindex
      integer :: Vsize(4)
!
      real(dp), parameter :: Fscl = 1.0_dp
      real(r8) :: Fmax, Fmin
      real(r8) :: wtype(LBi:UBi,LBj:UBj)
!
      character (len=256) :: ncname
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/get_grid.F"//", get_grid_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Inquire about the contents of grid NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, 153, MyFile)) RETURN
      ncname=GRD(ng)%name
!
!  Open grid NetCDF file for reading.
!
      IF (GRD(ng)%ncid.eq.-1) THEN
        CALL netcdf_open (ng, model, ncname, 0, GRD(ng)%ncid)
        IF (FoundError(exit_flag, NoError, 160, MyFile)) THEN
          WRITE (stdout,10) TRIM(ncname)
          RETURN
        END IF
      END IF
!
!  Check grid file dimensions for consitency.
!
      CALL netcdf_check_dim (ng, model, ncname,                         &
     &                       ncid = GRD(ng)%ncid)
      IF (FoundError(exit_flag, NoError, 170, MyFile)) RETURN
!
!  Inquire about the variables.
!
      CALL netcdf_inq_var (ng, model, ncname,                           &
     &                     ncid = GRD(ng)%ncid)
      IF (FoundError(exit_flag, NoError, 176, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Check if required variables are available.
!-----------------------------------------------------------------------
!
      IF (.not.find_string(var_name,n_var,'xl',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'xl', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'el',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'el', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'spherical',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'spherical', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'h',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'h', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'f',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'f', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'pm',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'pm', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'pn',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'pn', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'dndx',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'dndx', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'dmde',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'dmde', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'angle',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'angle', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'mask_rho',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'mask_rho', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'mask_u',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'mask_u', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'mask_v',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'mask_v', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'mask_psi',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'mask_psi', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'wtype_grid',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'wtype_grid', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (LuvSponge(ng)) THEN
        IF (.not.find_string(var_name,n_var,'visc_factor',vindex)) THEN
          IF (Master) WRITE (stdout,20) 'visc_factor', TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
      END IF
      IF (ANY(LtracerSponge(:,ng))) THEN
        IF (.not.find_string(var_name,n_var,'diff_factor',vindex)) THEN
          IF (Master) WRITE (stdout,20) 'diff_factor', TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
      END IF
!
!  Read in logical switch for spherical grid configuration.
!
      spherical=.FALSE.
      IF (find_string(var_name,n_var,'spherical',vindex)) THEN
        CALL netcdf_get_lvar (ng, model, ncname, 'spherical',           &
     &                        spherical,                                &
     &                        ncid = GRD(ng)%ncid)
        IF (FoundError(exit_flag, NoError, 342, MyFile)) RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Read in grid variables.
!-----------------------------------------------------------------------
!
!  Set Vsize to zero to deativate interpolation of input data to model
!  grid in "nf_fread2d".
!
      DO i=1,4
        Vsize(i)=0
      END DO
!
!  Scan the variable list and read in needed variables.
!
      IF (Master) WRITE (stdout,'(1x)')
!
      DO i=1,n_var
        SELECT CASE (TRIM(ADJUSTL(var_name(i))))
!
!  Read in basin X-length.
!
          CASE ('xl')
            CALL netcdf_get_fvar (ng, model, ncname,                    &
     &                            'xl', xl(ng),                         &
     &                            ncid = GRD(ng)%ncid)
            IF (FoundError(exit_flag, NoError, 370, MyFile)) EXIT
!
!  Read in basin Y-length.
!
          CASE ('el')
            CALL netcdf_get_fvar (ng, model, ncname,                    &
     &                            'el', el(ng),                         &
     &                            ncid = GRD(ng)%ncid)
            IF (FoundError(exit_flag, NoError, 378, MyFile)) EXIT
!
!  Read in bathymetry.
!
          CASE ('h')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % rmask,                         &
     &                        GRID(ng) % h)
            IF (FoundError(status, nf90_noerr, 398, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              Hmin(ng)=REAL(Fmin,dp)
              Hmax(ng)=REAL(Fmax,dp)
              IF (Master) THEN
                WRITE (stdout,30) 'bathymetry at RHO-points: h',        &
     &                            ng, TRIM(ncname), hmin(ng), hmax(ng)
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % h)
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % h)
!
!  Read in Land/Sea masking at RHO-points.
!
          CASE ('mask_rho')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % rmask,                         &
     &                        GRID(ng) % rmask)
            IF (FoundError(status, nf90_noerr, 457, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'mask on RHO-points: mask_rho',       &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % rmask)
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % rmask)
!
!  Read in Land/Sea masking at U-points.
!
          CASE ('mask_u')
            gtype=u2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % umask,                         &
     &                        GRID(ng) % umask)
            IF (FoundError(status, nf90_noerr, 510, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'mask on U-points: mask_u',           &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_u2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % umask)
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % umask)
!
!  Read in Land/Sea masking at V-points.
!
          CASE ('mask_v')
            gtype=v2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % vmask,                         &
     &                        GRID(ng) % vmask)
            IF (FoundError(status, nf90_noerr, 563, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'mask on V-points: mask_v',           &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_v2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % vmask)
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % vmask)
!
!  Read in Land/Sea masking at PSI-points.
!
          CASE ('mask_psi')
            gtype=p2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % pmask,                         &
     &                        GRID(ng) % pmask)
            IF (FoundError(status, nf90_noerr, 616, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'mask on PSI-points: mask_psi',       &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_p2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % pmask)
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % pmask)
!
!  Read in Jerlov water type.
!
          CASE ('wtype_grid')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % rmask,                         &
     &                        wtype)
            IF (FoundError(status, nf90_noerr, 710, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'Jerlov water type: wtype_grid',      &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                wtype)
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          wtype)
            MIXING(ng) % Jwtype = NINT(wtype)   ! nearest integer
!
!  Read in horizontal, spatially varying factor to increase/decrease
!  viscosity (nondimensional) in specific areas of the domain.
!
          CASE ('visc_factor')
            IF (LuvSponge(ng)) THEN
              gtype=r2dvar
              status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,        &
     &                          var_name(i), var_id(i),                 &
     &                          0, gtype, Vsize,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
     &                          GRID(ng) % rmask,                       &
     &                          MIXING(ng) % visc_factor)
              IF (FoundError(status, nf90_noerr, 759, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'horizontal viscosity sponge '//    &
     &                              'factor: visc_factor',              &
     &                              ng, TRIM(ncname), Fmin, Fmax
                END IF
              END IF
              IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
                CALL exchange_r2d_tile (ng, tile,                       &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  MIXING(ng) % visc_factor)
              END IF
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            MIXING(ng) % visc_factor)
            END IF
!
!  Read in horizontal, spatially varying factor to increase/decrease
!  diffusivity (nondimensional) in specific areas of the domain.
!
          CASE ('diff_factor')
            IF (ANY(LtracerSponge(:,ng))) THEN
              gtype=r2dvar
              status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,        &
     &                          var_name(i), var_id(i),                 &
     &                          0, gtype, Vsize,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
     &                          GRID(ng) % rmask,                       &
     &                          MIXING(ng) % diff_factor)
              IF (FoundError(status, nf90_noerr, 809, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'horizontal diffusivity sponge '//  &
     &                              'factor: diff_factor',              &
     &                              ng, TRIM(ncname), Fmin, Fmax
                END IF
              END IF
              IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
                CALL exchange_r2d_tile (ng, tile,                       &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  MIXING(ng) % diff_factor)
              END IF
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            MIXING(ng) % diff_factor)
            END IF
!
!  Read in Coriolis parameter.
!
          CASE ('f')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % rmask,                         &
     &                        GRID(ng) % f)
            IF (FoundError(status, nf90_noerr, 857, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'Coriolis parameter at RHO-points: f',&
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % f)
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % f)
!
!  Read in coordinate transfomation metrics (m) associated with the
!  differential distances in XI.
!
          CASE ('pm')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % rmask,                         &
     &                        GRID(ng) % pm)
            IF (FoundError(status, nf90_noerr, 911, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'reciprocal XI-grid spacing: pm',     &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % pm)
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % pm)
!
!  Read in coordinate transfomation metrics (n) associated with the
!  differential distances in ETA.
!
          CASE ('pn')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % rmask,                         &
     &                        GRID(ng) % pn)
            IF (FoundError(status, nf90_noerr, 965, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'reciprocal ETA-grid spacing: pn',    &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % pn)
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % pn)
!
!  Read in derivatives of inverse metrics factors: d(m)/d(eta).
!
          CASE ('dmde')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % rmask,                         &
     &                        GRID(ng) % dmde)
            IF (FoundError(status, nf90_noerr, 1019, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'ETA-derivative of inverse metric '// &
     &                            'factor pm: dmde',                    &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % dmde)
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % dmde)
!
!  Read in derivatives of inverse metrics factors: d(n)/d(xi).
!
          CASE ('dndx')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % rmask,                         &
     &                        GRID(ng) % dndx)
            IF (FoundError(status, nf90_noerr, 1073, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'XI-derivative of inverse metric '//  &
     &                            'factor pn: dndx',                    &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % dndx)
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % dndx)
!
!  Read in X-coordinates at PSI-points.
!
          CASE ('x_psi')
            gtype=p2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % pmask,                         &
     &                        GRID(ng) % xp)
            IF (FoundError(status, nf90_noerr, 1128, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'x-location of PSI-points: x_psi',    &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          .FALSE., .FALSE.,                       &
     &                          GRID(ng) % xp)
!
!  Read in Y-coordinates at PSI-points.
!
          CASE ('y_psi')
            gtype=p2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % pmask,                         &
     &                        GRID(ng) % yp)
            IF (FoundError(status, nf90_noerr, 1167, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'y-location of PSI-points: y-psi',    &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          .FALSE., .FALSE.,                       &
     &                          GRID(ng) % yp)
!
!  Read in X-coordinates at RHO-points.
!
          CASE ('x_rho')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % rmask,                         &
     &                        GRID(ng) % xr)
            IF (FoundError(status, nf90_noerr, 1206, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'x-location of RHO-points: x-rho',    &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            IF (.not.spherical) THEN
              LonMin(ng)=Fmin
              LonMax(ng)=Fmax
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          .FALSE., .FALSE.,                       &
     &                          GRID(ng) % xr)
!
!  Read in Y-coordinates at RHO-points.
!
          CASE ('y_rho')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % rmask,                         &
     &                        GRID(ng) % yr)
            IF (FoundError(status, nf90_noerr, 1261, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'y-location of RHO-points: y_rho',    &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            IF (.not.spherical) THEN
              LatMin(ng)=Fmin
              LatMax(ng)=Fmax
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          .FALSE., .FALSE.,                       &
     &                          GRID(ng) % yr)
!
!  Read in X-coordinates at U-points.
!
          CASE ('x_u')
            gtype=u2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % umask,                         &
     &                        GRID(ng) % xu)
            IF (FoundError(status, nf90_noerr, 1316, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'x-location of U-points: x_u',        &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          .FALSE., .FALSE.,                       &
     &                          GRID(ng) % xu)
!
!  Read in Y-coordinates at U-points.
!
          CASE ('y_u')
            gtype=u2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % umask,                         &
     &                        GRID(ng) % yu)
            IF (FoundError(status, nf90_noerr, 1367, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'y-location of U-points: y_u',        &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          .FALSE., .FALSE.,                       &
     &                          GRID(ng) % yu)
!
!  Read in X-coordinates at V-points.
!
          CASE ('x_v')
            gtype=v2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % vmask,                         &
     &                        GRID(ng) % xv)
            IF (FoundError(status, nf90_noerr, 1418, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'x-location of V-points: x_v',        &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          .FALSE., .FALSE.,                       &
     &                          GRID(ng) % xv)
!
!  Read in Y-coordinates at V-points.
!
          CASE ('y_v')
            gtype=v2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % vmask,                         &
     &                        GRID(ng) % yv)
            IF (FoundError(status, nf90_noerr, 1469, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'y-location of V-points: y_v',        &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          .FALSE., .FALSE.,                       &
     &                          GRID(ng) % yv)
!
!  Read in longitude at PSI-points.
!
          CASE ('lon_psi')
            IF (spherical) THEN
              gtype=p2dvar
              status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,        &
     &                          var_name(i), var_id(i),                 &
     &                          0, gtype, Vsize,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
     &                          GRID(ng) % pmask,                       &
     &                          GRID(ng) % lonp)
              IF (FoundError(status, nf90_noerr, 1521, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'longitude of PSI-points: lon_psi', &
     &                              ng, TRIM(ncname), Fmin, Fmax
                END IF
              END IF
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            .FALSE., .FALSE.,                     &
     &                            GRID(ng) % lonp)
            END IF
!
!  Read in latitude at PSI-points.
!
          CASE ('lat_psi')
            IF (spherical) THEN
              gtype=p2dvar
              status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,        &
     &                          var_name(i), var_id(i),                 &
     &                          0, gtype, Vsize,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
     &                          GRID(ng) % pmask,                       &
     &                          GRID(ng) % latp)
              IF (FoundError(status, nf90_noerr, 1562, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'latitude of PSI-points lat_psi',   &
     &                              ng, TRIM(ncname), Fmin, Fmax
                END IF
              END IF
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            .FALSE., .FALSE.,                     &
     &                            GRID(ng) % latp)
            END IF
!
!  Read in longitude at RHO-points.
!
          CASE ('lon_rho')
            IF (spherical) THEN
              gtype=r2dvar
              status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,        &
     &                          var_name(i), var_id(i),                 &
     &                          0, gtype, Vsize,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, LonMin(ng), LonMax(ng),           &
     &                          GRID(ng) % rmask,                       &
     &                          GRID(ng) % lonr)
              IF (FoundError(status, nf90_noerr, 1603, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'longitude of RHO-points: lon_rho', &
     &                              ng, TRIM(ncname),                   &
     &                              LonMin(ng), LonMax(ng)
                END IF
              END IF
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            .FALSE., .FALSE.,                     &
     &                            GRID(ng) % lonr)
            END IF
!
!  Read in latitude at RHO-points.
!
          CASE ('lat_rho')
            IF (spherical) THEN
              gtype=r2dvar
              status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,        &
     &                          var_name(i), var_id(i),                 &
     &                          0, gtype, Vsize,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, LatMin(ng), LatMax(ng),           &
     &                          GRID(ng) % rmask,                       &
     &                          GRID(ng) % latr)
              IF (FoundError(status, nf90_noerr, 1655, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'latitude of RHO-points lat_rho',   &
     &                              ng, TRIM(ncname),                   &
     &                              LatMin(ng), LatMax(ng)
                END IF
              END IF
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            .FALSE., .FALSE.,                     &
     &                            GRID(ng) % latr)
            END IF
!
!  Read in longitude at U-points.
!
          CASE ('lon_u')
            IF (spherical) THEN
              gtype=u2dvar
              status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,        &
     &                          var_name(i), var_id(i),                 &
     &                          0, gtype, Vsize,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
     &                          GRID(ng) % umask,                       &
     &                          GRID(ng) % lonu)
              IF (FoundError(status, nf90_noerr, 1707, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'longitude of U-points: lon_u',     &
     &                              ng, TRIM(ncname), Fmin, Fmax
                END IF
              END IF
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            .FALSE., .FALSE.,                     &
     &                            GRID(ng) % lonu)
            END IF
!
!  Read in latitude at U-points.
!
          CASE ('lat_u')
            IF (spherical) THEN
              gtype=u2dvar
              status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,        &
     &                          var_name(i), var_id(i),                 &
     &                          0, gtype, Vsize,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
     &                          GRID(ng) % umask,                       &
     &                          GRID(ng) % latu)
              IF (FoundError(status, nf90_noerr, 1758, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'latitude of U-points: lat_u',      &
     &                              ng, TRIM(ncname), Fmin, Fmax
                END IF
              END IF
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            .FALSE., .FALSE.,                     &
     &                            GRID(ng) % latu)
            END IF
!
!  Read in longitude at V-points.
!
          CASE ('lon_v')
            IF (spherical) THEN
              gtype=v2dvar
              status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,        &
     &                          var_name(i), var_id(i),                 &
     &                          0, gtype, Vsize,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
     &                          GRID(ng) % vmask,                       &
     &                          GRID(ng) % lonv)
              IF (FoundError(status, nf90_noerr, 1809, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'longitude of V-points: lon_v',     &
     &                              ng, TRIM(ncname), Fmin, Fmax
                END IF
              END IF
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            .FALSE., .FALSE.,                     &
     &                            GRID(ng) % lonv)
            END IF
!
!  Read in latitude at V-points.
!
          CASE ('lat_v')
            IF (spherical) THEN
              gtype=v2dvar
              status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,        &
     &                          var_name(i), var_id(i),                 &
     &                          0, gtype, Vsize,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
     &                          GRID(ng) % vmask,                       &
     &                          GRID(ng) % latv)
              IF (FoundError(status, nf90_noerr, 1860, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'latitude of V-points: lat_v',      &
     &                              ng, TRIM(ncname), Fmin, Fmax
                END IF
              END IF
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            .FALSE., .FALSE.,                     &
     &                            GRID(ng) % latv)
            END IF
!
!  Read in angle (radians) between XI-axis and EAST at RHO-points.
!
          CASE ('angle')
            gtype=r2dvar
            status=nf_fread2d(ng, model, ncname, GRD(ng)%ncid,          &
     &                        var_name(i), var_id(i),                   &
     &                        0, gtype, Vsize,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % rmask,                         &
     &                        GRID(ng) % angler)
            IF (FoundError(status, nf90_noerr, 1910, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'angle between XI-axis and EAST: '//  &
     &                            'angler',                             &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % angler)
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % angler)
        END SELECT
      END DO
      IF (FoundError(exit_flag, NoError, 2091, MyFile)) THEN
        IF (Master) WRITE (stdout,40) TRIM(var_name(i)), TRIM(ncname)
        RETURN
      END IF
!
! Close GRID NetCDF file.
!
      CALL netcdf_close (ng, model, GRD(ng)%ncid, ncname, .FALSE.)
      IF (FoundError(exit_flag, NoError, 2278, MyFile)) RETURN
!
  10  FORMAT (/,' GET_GRID_NF90 - unable to open grid NetCDF file: ',a)
  20  FORMAT (/,' GET_GRID_NF90 - unable to find grid variable: ',a,    &
     &        /,12x,'in grid NetCDF file: ',a)
  30  FORMAT (2x,'GET_GRID_NF90    - ',a,/,22x,                         &
     &        '(Grid = ',i2.2,', File: ',a,')',/,22x,                   &
     &        '(Min = ', 1p,e15.8,0p,' Max = ',1p,e15.8,0p,')')
  40  FORMAT (/,' GET_GRID_NF90 - error while reading variable: ',a,    &
     &        /,12x,'in grid NetCDF file: ',a)
  50  FORMAT (/,2x,'GET_GRID_NF90    - Reading adjoint sensitivity',    &
     &        ' scope arrays from file:',/22x,a,/)
!
      RETURN
      END SUBROUTINE get_grid_nf90
!
!***********************************************************************
      SUBROUTINE get_grid_pio (ng, tile, model,                         &
     &                         LBi, UBi, LBj, UBj)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: cr, i, status, vindex
      integer :: Vsize(4)
!
      real(dp), parameter :: Fscl = 1.0_dp
      real(r8) :: Fmax, Fmin
      real(r8) :: wtype(LBi:UBi,LBj:UBj)
!
      character (len=256) :: ncname
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/get_grid.F"//", get_grid_pio"
!
      TYPE (IO_desc_t), pointer :: ioDesc
      TYPE (My_VarDesc) :: pioVar
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Inquire about the contents of grid NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, 2567, MyFile)) RETURN
      ncname=GRD(ng)%name
!
!  Open grid NetCDF file for reading.
!
      IF (GRD(ng)%pioFile%fh.eq.-1) THEN
        CALL pio_netcdf_open (ng, model, ncname, 0, GRD(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, 2574, MyFile)) THEN
          WRITE (stdout,10) TRIM(ncname)
          RETURN
        END IF
      END IF
!
!  Check grid file dimensions for consitency.
!
      CALL pio_netcdf_check_dim (ng, model, ncname,                     &
     &                           pioFile = GRD(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, 2584, MyFile)) RETURN
!
!  Inquire about the variables.
!
      CALL pio_netcdf_inq_var (ng, model, ncname,                       &
     &                         pioFile = GRD(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, 2590, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Check if required variables are available.
!-----------------------------------------------------------------------
!
      IF (.not.find_string(var_name,n_var,'xl',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'xl', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'el',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'el', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'spherical',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'spherical', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'h',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'h', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'f',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'f', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'pm',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'pm', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'pn',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'pn', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'dndx',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'dndx', TRIM(ncname)
!       exit_flag=2
!       RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'dmde',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'dmde', TRIM(ncname)
!       exit_flag=2
!       RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'angle',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'angle', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'mask_rho',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'mask_rho', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'mask_u',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'mask_u', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'mask_v',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'mask_v', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'mask_psi',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'mask_psi', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.find_string(var_name,n_var,'wtype_grid',vindex)) THEN
        IF (Master) WRITE (stdout,20) 'wtype_grid', TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (LuvSponge(ng)) THEN
        IF (.not.find_string(var_name,n_var,'visc_factor',vindex)) THEN
          IF (Master) WRITE (stdout,20) 'visc_factor', TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
      END IF
      IF (ANY(LtracerSponge(:,ng))) THEN
        IF (.not.find_string(var_name,n_var,'diff_factor',vindex)) THEN
          IF (Master) WRITE (stdout,20) 'diff_factor', TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
      END IF
!
!  Read in logical switch for spherical grid configuration.
!
      spherical=.FALSE.
      IF (find_string(var_name,n_var,'spherical',vindex)) THEN
        CALL pio_netcdf_get_lvar (ng, model, ncname,                    &
     &                            'spherical', spherical,               &
     &                            pioFile = GRD(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, 2765, MyFile)) RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Read in grid variables.
!-----------------------------------------------------------------------
!
!  Set Vsize to zero to deativate interpolation of input data to model
!  grid in "nf_fread2d".
!
      DO i=1,4
        Vsize(i)=0
      END DO
!
!  Scan the variable list and read in needed variables.
!
      IF (Master) WRITE (stdout,'(1x)')
!
      DO i=1,n_var
        SELECT CASE (TRIM(ADJUSTL(var_name(i))))
!
!  Read in basin X-length.
!
          CASE ('xl')
            CALL pio_netcdf_get_fvar (ng, model, ncname,                &
     &                                'xl', xl(ng),                     &
     &                                pioFile = GRD(ng)%pioFile)
            IF (FoundError(exit_flag, NoError, 2793, MyFile)) EXIT
!
!  Read in basin Y-length.
!
          CASE ('el')
            CALL pio_netcdf_get_fvar (ng, model, ncname,                &
     &                                'el', el(ng),                     &
     &                                pioFile = GRD(ng)%pioFile)
            IF (FoundError(exit_flag, NoError, 2801, MyFile)) EXIT
!
!  Read in bathymetry.
!
          CASE ('h')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(GRID(ng)%h).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % rmask,                         &
     &                        GRID(ng) % h)
            IF (FoundError(status, nf90_noerr, 2830, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              Hmin(ng)=REAL(Fmin,dp)
              Hmax(ng)=REAL(Fmax,dp)
              IF (Master) THEN
                WRITE (stdout,30) 'bathymetry at RHO-points: h',        &
     &                            ng, TRIM(ncname), hmin(ng), hmax(ng)
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % h)
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % h)
!
!  Read in Land/Sea masking at RHO-points.
!
          CASE ('mask_rho')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(GRID(ng)%rmask).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % rmask,                         &
     &                        GRID(ng) % rmask)
            IF (FoundError(status, nf90_noerr, 2898, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'mask on RHO-points: mask_rho',       &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % rmask)
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % rmask)
!
!  Read in Land/Sea masking at U-points.
!
          CASE ('mask_u')
            pioVar%vd=var_desc(i)
            pioVar%gtype=u2dvar
            IF (KIND(GRID(ng)%umask).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_u2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_u2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % umask,                         &
     &                        GRID(ng) % umask)
            IF (FoundError(status, nf90_noerr, 2958, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'mask on U-points: mask_u',           &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_u2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % umask)
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % umask)
!
!  Read in Land/Sea masking at V-points.
!
          CASE ('mask_v')
            pioVar%vd=var_desc(i)
            pioVar%gtype=v2dvar
            IF (KIND(GRID(ng)%vmask).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_v2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_v2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % vmask,                         &
     &                        GRID(ng) % vmask)
            IF (FoundError(status, nf90_noerr, 3018, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'mask on V-points: mask_v',           &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_v2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % vmask)
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % vmask)
!
!  Read in Land/Sea masking at PSI-points.
!
          CASE ('mask_psi')
            pioVar%vd=var_desc(i)
            pioVar%gtype=p2dvar
            IF (KIND(GRID(ng)%pmask).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_p2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_p2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % pmask,                         &
     &                        GRID(ng) % pmask)
            IF (FoundError(status, nf90_noerr, 3078, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'mask on PSI-points: mask_psi',       &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_p2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % pmask)
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % pmask)
!
!  Read in Jerlov water type.
!
          CASE ('wtype_grid')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(wtype).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % rmask,                         &
     &                        wtype)
            IF (FoundError(status, nf90_noerr, 3190, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'Jerlov water type: wtype_grid',      &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                wtype)
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          wtype)
            MIXING(ng) % Jwtype = NINT(wtype)   ! nearest integer
!
!  Read in horizontal, spatially varying factor to increase/decrease
!  viscosity (nondimensional) in specific areas of the domain.
!
          CASE ('visc_factor')
            IF (LuvSponge(ng)) THEN
              pioVar%vd=var_desc(i)
              pioVar%gtype=r2dvar
              IF (KIND(MIXING(ng)%visc_factor).eq.8) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_r2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_r2dvar(ng)
              END IF
!
              status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,     &
     &                          var_name(i), pioVar,                    &
     &                          0, ioDesc, Vsize,                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
     &                          GRID(ng) % rmask,                       &
     &                          MIXING(ng) % visc_factor)
              IF (FoundError(status, nf90_noerr, 3248, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'horizontal viscosity sponge '//    &
     &                              'factor: visc_factor',              &
     &                              ng, TRIM(ncname), Fmin, Fmax
                END IF
              END IF
              IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
                CALL exchange_r2d_tile (ng, tile,                       &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  MIXING(ng) % visc_factor)
              END IF
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            MIXING(ng) % visc_factor)
            END IF
!
!  Read in horizontal, spatially varying factor to increase/decrease
!  diffusivity (nondimensional) in specific areas of the domain.
!
          CASE ('diff_factor')
            IF (ANY(LtracerSponge(:,ng))) THEN
              pioVar%vd=var_desc(i)
              pioVar%gtype=r2dvar
              IF (KIND(MIXING(ng)%diff_factor).eq.8) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_r2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_r2dvar(ng)
              END IF
!
              status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,     &
     &                          var_name(i), pioVar,                    &
     &                          0, ioDesc, Vsize,                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
     &                          GRID(ng) % rmask,                       &
     &                          MIXING(ng) % diff_factor)
              IF (FoundError(status, nf90_noerr, 3307, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'horizontal diffusivity sponge '//  &
     &                              'factor: diff_factor',              &
     &                              ng, TRIM(ncname), Fmin, Fmax
                END IF
              END IF
              IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
                CALL exchange_r2d_tile (ng, tile,                       &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  MIXING(ng) % diff_factor)
              END IF
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            MIXING(ng) % diff_factor)
            END IF
!
!  Read in Coriolis parameter.
!
          CASE ('f')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(GRID(ng)%f).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % rmask,                         &
     &                        GRID(ng) % f)
            IF (FoundError(status, nf90_noerr, 3364, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'Coriolis parameter at RHO-points: f',&
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % f)
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % f)
!
!  Read in coordinate transfomation metrics (m) associated with the
!  differential distances in XI.
!
          CASE ('pm')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(GRID(ng)%pm).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % rmask,                         &
     &                        GRID(ng) % pm)
            IF (FoundError(status, nf90_noerr, 3427, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'reciprocal XI-grid spacing: pm',     &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % pm)
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % pm)
!
!  Read in coordinate transfomation metrics (n) associated with the
!  differential distances in ETA.
!
          CASE ('pn')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(GRID(ng)%pn).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % rmask,                         &
     &                        GRID(ng) % pn)
            IF (FoundError(status, nf90_noerr, 3490, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'reciprocal ETA-grid spacing: pn',    &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % pn)
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % pn)
!
!  Read in derivatives of inverse metrics factors: d(m)/d(eta).
!
          CASE ('dmde')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(GRID(ng)%dmde).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % rmask,                         &
     &                        GRID(ng) % dmde)
            IF (FoundError(status, nf90_noerr, 3553, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'ETA-derivative of inverse metric '// &
     &                            'factor pm: dmde',                    &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % dmde)
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % dmde)
!
!  Read in derivatives of inverse metrics factors: d(n)/d(xi).
!
          CASE ('dndx')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(GRID(ng)%dndx).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % rmask,                         &
     &                        GRID(ng) % dndx)
            IF (FoundError(status, nf90_noerr, 3616, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'XI-derivative of inverse metric '//  &
     &                            'factor pn: dndx',                    &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % dndx)
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % dndx)
!
!  Read in X-coordinates at PSI-points.
!
          CASE ('x_psi')
            pioVar%vd=var_desc(i)
            pioVar%gtype=p2dvar
            IF (KIND(GRID(ng)%xp).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_p2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_p2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % pmask,                         &
     &                        GRID(ng) % xp)
            IF (FoundError(status, nf90_noerr, 3680, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'x-location of PSI-points: x_psi',    &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          .FALSE., .FALSE.,                       &
     &                          GRID(ng) % xp)
!
!  Read in Y-coordinates at PSI-points.
!
          CASE ('y_psi')
            pioVar%vd=var_desc(i)
            pioVar%gtype=p2dvar
            IF (KIND(GRID(ng)%yp).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_p2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_p2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % pmask,                         &
     &                        GRID(ng) % yp)
            IF (FoundError(status, nf90_noerr, 3728, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'y-location of PSI-points: y-psi',    &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          .FALSE., .FALSE.,                       &
     &                          GRID(ng) % yp)
!
!  Read in X-coordinates at RHO-points.
!
          CASE ('x_rho')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(GRID(ng)%xr).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % rmask,                         &
     &                        GRID(ng) % xr)
            IF (FoundError(status, nf90_noerr, 3776, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'x-location of RHO-points: x-rho',    &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          .FALSE., .FALSE.,                       &
     &                          GRID(ng) % xr)
!
!  Read in Y-coordinates at RHO-points.
!
          CASE ('y_rho')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(GRID(ng)%yr).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % rmask,                         &
     &                        GRID(ng) % yr)
            IF (FoundError(status, nf90_noerr, 3836, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'y-location of RHO-points: y_rho',    &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          .FALSE., .FALSE.,                       &
     &                          GRID(ng) % yr)
!
!  Read in X-coordinates at U-points.
!
          CASE ('x_u')
            pioVar%vd=var_desc(i)
            pioVar%gtype=u2dvar
            IF (KIND(GRID(ng)%xu).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_u2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_u2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % umask,                         &
     &                        GRID(ng) % xu)
            IF (FoundError(status, nf90_noerr, 3896, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'x-location of U-points: x_u',        &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          .FALSE., .FALSE.,                       &
     &                          GRID(ng) % xu)
!
!  Read in Y-coordinates at U-points.
!
          CASE ('y_u')
            pioVar%vd=var_desc(i)
            pioVar%gtype=u2dvar
            IF (KIND(GRID(ng)%yu).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_u2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_u2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % umask,                         &
     &                        GRID(ng) % yu)
            IF (FoundError(status, nf90_noerr, 3956, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'y-location of U-points: y_u',        &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          .FALSE., .FALSE.,                       &
     &                          GRID(ng) % yu)
!
!  Read in X-coordinates at V-points.
!
          CASE ('x_v')
            pioVar%vd=var_desc(i)
            pioVar%gtype=v2dvar
            IF (KIND(GRID(ng)%xv).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_v2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_v2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % vmask,                         &
     &                        GRID(ng) % xv)
            IF (FoundError(status, nf90_noerr, 4016, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'x-location of V-points: x_v',        &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          .FALSE., .FALSE.,                       &
     &                          GRID(ng) % xv)
!
!  Read in Y-coordinates at V-points.
!
          CASE ('y_v')
            pioVar%vd=var_desc(i)
            pioVar%gtype=v2dvar
            IF (KIND(GRID(ng)%yv).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_v2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_v2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % vmask,                         &
     &                        GRID(ng) % yv)
            IF (FoundError(status, nf90_noerr, 4076, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'y-location of V-points: y_v',        &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          .FALSE., .FALSE.,                       &
     &                          GRID(ng) % yv)
!
!  Read in longitude at PSI-points.
!
          CASE ('lon_psi')
            IF (spherical) THEN
              pioVar%vd=var_desc(i)
              pioVar%gtype=p2dvar
              IF (KIND(GRID(ng)%lonp).eq.8) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_p2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_p2dvar(ng)
              END IF
!
              status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,     &
     &                          var_name(i), pioVar,                    &
     &                          0, ioDesc, Vsize,                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
     &                          GRID(ng) % pmask,                       &
     &                          GRID(ng) % lonp)
              IF (FoundError(status, nf90_noerr, 4137, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'longitude of PSI-points: lon_psi', &
     &                              ng, TRIM(ncname), Fmin, Fmax
                END IF
              END IF
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            .FALSE., .FALSE.,                     &
     &                            GRID(ng) % lonp)
            END IF
!
!  Read in latitude at PSI-points.
!
          CASE ('lat_psi')
            IF (spherical) THEN
              pioVar%vd=var_desc(i)
              pioVar%gtype=p2dvar
              IF (KIND(GRID(ng)%latp).eq.8) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_p2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_p2dvar(ng)
              END IF
!
              status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,     &
     &                          var_name(i), pioVar,                    &
     &                          0, ioDesc, Vsize,                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
     &                          GRID(ng) % pmask,                       &
     &                          GRID(ng) % latp)
              IF (FoundError(status, nf90_noerr, 4187, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'latitude of PSI-points lat_psi',   &
     &                              ng, TRIM(ncname), Fmin, Fmax
                END IF
              END IF
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            .FALSE., .FALSE.,                     &
     &                            GRID(ng) % latp)
            END IF
!
!  Read in longitude at RHO-points.
!
          CASE ('lon_rho')
            IF (spherical) THEN
              pioVar%vd=var_desc(i)
              pioVar%gtype=r2dvar
              IF (KIND(GRID(ng)%lonr).eq.8) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_r2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_r2dvar(ng)
              END IF
!
              status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,     &
     &                          var_name(i), pioVar,                    &
     &                          0, ioDesc, Vsize,                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, LonMin(ng), LonMax(ng),           &
     &                          GRID(ng) % rmask,                       &
     &                          GRID(ng) % lonr)
              IF (FoundError(status, nf90_noerr, 4237, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'longitude of RHO-points: lon_rho', &
     &                              ng, TRIM(ncname),                   &
     &                              LonMin(ng), LonMax(ng)
                END IF
              END IF
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            .FALSE., .FALSE.,                     &
     &                            GRID(ng) % lonr)
            END IF
!
!  Read in latitude at RHO-points.
!
          CASE ('lat_rho')
            IF (spherical) THEN
              pioVar%vd=var_desc(i)
              pioVar%gtype=r2dvar
              IF (KIND(GRID(ng)%latr).eq.8) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_r2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_r2dvar(ng)
              END IF
!
              status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,     &
     &                          var_name(i), pioVar,                    &
     &                          0, ioDesc, Vsize,                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, LatMin(ng), LatMax(ng),           &
     &                          GRID(ng) % rmask,                       &
     &                          GRID(ng) % latr)
              IF (FoundError(status, nf90_noerr, 4298, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'latitude of RHO-points lat_rho',   &
     &                              ng, TRIM(ncname),                   &
     &                              LatMin(ng), LatMax(ng)
                END IF
              END IF
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            .FALSE., .FALSE.,                     &
     &                            GRID(ng) % latr)
            END IF
!
!  Read in longitude at U-points.
!
          CASE ('lon_u')
            IF (spherical) THEN
              pioVar%vd=var_desc(i)
              pioVar%gtype=u2dvar
              IF (KIND(GRID(ng)%lonu).eq.8) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_u2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_u2dvar(ng)
              END IF
!
              status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,     &
     &                          var_name(i), pioVar,                    &
     &                          0, ioDesc, Vsize,                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
     &                          GRID(ng) % umask,                       &
     &                          GRID(ng) % lonu)
              IF (FoundError(status, nf90_noerr, 4359, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'longitude of U-points: lon_u',     &
     &                              ng, TRIM(ncname), Fmin, Fmax
                END IF
              END IF
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            .FALSE., .FALSE.,                     &
     &                            GRID(ng) % lonu)
            END IF
!
!  Read in latitude at U-points.
!
          CASE ('lat_u')
            IF (spherical) THEN
              pioVar%vd=var_desc(i)
              pioVar%gtype=u2dvar
              IF (KIND(GRID(ng)%latu).eq.8) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_u2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_u2dvar(ng)
              END IF
!
              status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,     &
     &                          var_name(i), pioVar,                    &
     &                          0, ioDesc, Vsize,                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
     &                          GRID(ng) % umask,                       &
     &                          GRID(ng) % latu)
              IF (FoundError(status, nf90_noerr, 4419, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'latitude of U-points: lat_u',      &
     &                              ng, TRIM(ncname), Fmin, Fmax
                END IF
              END IF
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            .FALSE., .FALSE.,                     &
     &                            GRID(ng) % latu)
            END IF
!
!  Read in longitude at V-points.
!
          CASE ('lon_v')
            IF (spherical) THEN
              pioVar%vd=var_desc(i)
              pioVar%gtype=v2dvar
              IF (KIND(GRID(ng)%lonv).eq.8) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_v2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_v2dvar(ng)
              END IF
!
              status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,     &
     &                          var_name(i), pioVar,                    &
     &                          0, ioDesc, Vsize,                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
     &                          GRID(ng) % vmask,                       &
     &                          GRID(ng) % lonv)
              IF (FoundError(status, nf90_noerr, 4479, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'longitude of V-points: lon_v',     &
     &                              ng, TRIM(ncname), Fmin, Fmax
                END IF
              END IF
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            .FALSE., .FALSE.,                     &
     &                            GRID(ng) % lonv)
            END IF
!
!  Read in latitude at V-points.
!
          CASE ('lat_v')
            IF (spherical) THEN
              pioVar%vd=var_desc(i)
              pioVar%gtype=v2dvar
              IF (KIND(GRID(ng)%latv).eq.8) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_v2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_v2dvar(ng)
              END IF
!
              status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,     &
     &                          var_name(i), pioVar,                    &
     &                          0, ioDesc, Vsize,                       &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Fscl, Fmin, Fmax,                       &
     &                          GRID(ng) % vmask,                       &
     &                          GRID(ng) % latv)
              IF (FoundError(status, nf90_noerr, 4539, MyFile)) THEN
                exit_flag=2
                ioerror=status
                EXIT
              ELSE
                IF (Master) THEN
                  WRITE (stdout,30) 'latitude of V-points: lat_v',      &
     &                              ng, TRIM(ncname), Fmin, Fmax
                END IF
              END IF
              CALL mp_exchange2d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            .FALSE., .FALSE.,                     &
     &                            GRID(ng) % latv)
            END IF
!
!  Read in angle (radians) between XI-axis and EAST at RHO-points.
!
          CASE ('angle')
            pioVar%vd=var_desc(i)
            pioVar%gtype=r2dvar
            IF (KIND(GRID(ng)%angler).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r2dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r2dvar(ng)
            END IF
!
            status=nf_fread2d(ng, model, ncname, GRD(ng)%pioFile,       &
     &                        var_name(i), pioVar,                      &
     &                        0, ioDesc, Vsize,                         &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
     &                        GRID(ng) % rmask,                         &
     &                        GRID(ng) % angler)
            IF (FoundError(status, nf90_noerr, 4598, MyFile)) THEN
              exit_flag=2
              ioerror=status
              EXIT
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) 'angle between XI-axis and EAST: '//  &
     &                            'angler',                             &
     &                            ng, TRIM(ncname), Fmin, Fmax
              END IF
            END IF
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                GRID(ng) % angler)
            END IF
            CALL mp_exchange2d (ng, tile, model, 1,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          GRID(ng) % angler)
        END SELECT
      END DO
      IF (FoundError(exit_flag, NoError, 4806, MyFile)) THEN
        IF (Master) WRITE (stdout,40) TRIM(var_name(i)), TRIM(ncname)
        RETURN
      END IF
!
! Close GRID NetCDF file.
!
      CALL pio_netcdf_close (ng, model, GRD(ng)%pioFile, ncname, .FALSE.)
      IF (FoundError(exit_flag, NoError, 4980, MyFile)) RETURN
!
  10  FORMAT (/,' GET_GRID_PIO - unable to open grid NetCDF file: ',a)
  20  FORMAT (/,' GET_GRID_PIO - unable to find grid variable: ',a,     &
     &        /,16x,'in grid NetCDF file: ',a)
  30  FORMAT (2x,'GET_GRID_PIO     - ',a,/,22x,                         &
     &        '(Grid = ',i2.2,', File: ',a,')',/,22x,                   &
     &        '(Min = ', 1p,e15.8,0p,' Max = ',1p,e15.8,0p,')')
  40  FORMAT (/,' GET_GRID_PIO - error while reading variable: ',a,     &
     &        /,12x,'in grid NetCDF file: ',a)
  50  FORMAT (/,2x,'GET_GRID_PIO     - Reading adjoint sensitivity',    &
     &        ' scope arrays from file:',/22x,a,/)
!
      RETURN
      END SUBROUTINE get_grid_pio
      END MODULE get_grid_mod
