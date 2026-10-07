      MODULE wrt_info_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This routine defines information variables in requested NetCDF      !
!  file.                                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      Use mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_sources
!
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
      USE strings_mod,     ONLY : FoundError, find_string
!
      implicit none
!
      INTERFACE wrt_info
        MODULE PROCEDURE wrt_info_nf90
        MODULE PROCEDURE wrt_info_pio
      END INTERFACE wrt_info
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wrt_info_nf90 (ng, model, ncid, ncname)
!***********************************************************************
!                                                                      !
!  This routine writes out information variables into requested        !
!  NetCDF file using the standard NetCDF-3 or NetCDF-4 library.        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncid         NetCDF file ID (integer)                            !
!     ncname       NetCDF filename (string)                            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!***********************************************************************
!
      USE distribute_mod,  ONLY : mp_bcasti
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid
!
      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      logical :: Cgrid = .TRUE.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: i, j, k, ibry, ilev, itrc, status, varid
      integer, dimension(2) :: ibuffer
      integer :: ifield = 0
!
      real(dp) :: scale
      real(r8), dimension(NT(ng)) :: nudg
      real(r8), dimension(NT(ng),4) :: Tobc
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/wrt_info.F"//", wrt_info_nf90"
!
      SourceFile=MyFile
!
      IF (ncid.ne.XTR(ng)%ncid) THEN
        LBi=LBOUND(GRID(ng)%h,DIM=1)
        UBi=UBOUND(GRID(ng)%h,DIM=1)
        LBj=LBOUND(GRID(ng)%h,DIM=2)
        UBj=UBOUND(GRID(ng)%h,DIM=2)
      END IF
!
!-----------------------------------------------------------------------
!  Write out running parameters.
!-----------------------------------------------------------------------
!
!  Inquire about the variables.
!
      CALL netcdf_inq_var (ng, model, ncname, ncid)
      IF (FoundError(exit_flag, NoError, 144, MyFile)) RETURN
!
!  Time stepping parameters.
!
      CALL netcdf_put_ivar (ng, model, ncname, 'ntimes',                &
     &                      ntimes(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 151, MyFile)) RETURN
      CALL netcdf_put_ivar (ng, model, ncname, 'ndtfast',               &
     &                      ndtfast(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 156, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'dt',                    &
     &                      dt(ng), (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 161, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'dtfast',                &
     &                      dtfast(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 166, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'dstart',                &
     &                      dstart, (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 171, MyFile)) RETURN
      CALL netcdf_put_ivar (ng, model, ncname, 'nHIS',                  &
     &                      nHIS(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 205, MyFile)) RETURN
      CALL netcdf_put_ivar (ng, model, ncname, 'ndefHIS',               &
     &                      ndefHIS(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 210, MyFile)) RETURN
      CALL netcdf_put_ivar (ng, model, ncname, 'nRST',                  &
     &                      nRST(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 232, MyFile)) RETURN
      CALL netcdf_put_ivar (ng, model, ncname, 'ntsAVG',                &
     &                      ntsAVG(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 241, MyFile)) RETURN
      CALL netcdf_put_ivar (ng, model, ncname, 'nAVG',                  &
     &                      nAVG(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 246, MyFile)) RETURN
      CALL netcdf_put_ivar (ng, model, ncname, 'ndefAVG',               &
     &                      ndefAVG(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 251, MyFile)) RETURN
!
!  Power-law shape filter parameters for time-averaging of barotropic
!  fields.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'Falpha',                &
     &                      Falpha, (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 373, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'Fbeta',                 &
     &                      Fbeta, (/0/), (/0/),                        &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 378, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'Fgamma',                &
     &                      Fgamma, (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 383, MyFile)) RETURN
!
!  Horizontal mixing coefficients.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'nl_tnu2',               &
     &                      nl_tnu2(:,ng), (/1/), (/NT(ng)/),           &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 392, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'nl_visc2',              &
     &                      nl_visc2(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 444, MyFile)) RETURN
      CALL netcdf_put_lvar (ng, model, ncname, 'LuvSponge',             &
     &                      LuvSponge(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 502, MyFile)) RETURN
      CALL netcdf_put_lvar (ng, model, ncname, 'LtracerSponge',         &
     &                      LtracerSponge(:,ng), (/1/), (/NT(ng)/),     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 509, MyFile)) RETURN
!
!  Background vertical mixing coefficients.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'Akt_bak',               &
     &                      Akt_bak(:,ng), (/1/), (/NT(ng)/),           &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 519, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'Akv_bak',               &
     &                      Akv_bak(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 524, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'Akk_bak',               &
     &                      Akk_bak(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 530, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'Akp_bak',               &
     &                      Akp_bak(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 535, MyFile)) RETURN
!
!  Drag coefficients.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'rdrg',                  &
     &                      rdrg(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 577, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'rdrg2',                 &
     &                      rdrg2(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 582, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'Zob',                   &
     &                      Zob(ng), (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 588, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'Zos',                   &
     &                      Zos(ng), (/0/), (/0/),                      &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 593, MyFile)) RETURN
!
!  Generic length-scale parameters.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'gls_p',                 &
     &                      gls_p(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 603, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'gls_m',                 &
     &                      gls_m(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 608, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'gls_n',                 &
     &                      gls_n(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 613, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'gls_cmu0',              &
     &                      gls_cmu0(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 618, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'gls_c1',                &
     &                      gls_c1(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 623, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'gls_c2',                &
     &                      gls_c2(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 628, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'gls_c3m',               &
     &                      gls_c3m(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 633, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'gls_c3p',               &
     &                      gls_c3p(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 638, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'gls_sigk',              &
     &                      gls_sigk(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 643, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'gls_sigp',              &
     &                      gls_sigp(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 648, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'gls_Kmin',              &
     &                      gls_Kmin(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 653, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'gls_Pmin',              &
     &                      gls_Pmin(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 658, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'Charnok_alpha',         &
     &                      charnok_alpha(ng), (/0/), (/0/),            &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 663, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'Zos_hsig_alpha',        &
     &                      zos_hsig_alpha(ng), (/0/), (/0/),           &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 668, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'sz_alpha',              &
     &                      sz_alpha(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 673, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'CrgBan_cw',             &
     &                      crgban_cw(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 678, MyFile)) RETURN
!
!  Nudging inverse time scales used in various tasks.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'Znudg',                 &
     &                      Znudg(ng)/sec2day, (/0/), (/0/),            &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 692, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'M2nudg',                &
     &                      M2nudg(ng)/sec2day, (/0/), (/0/),           &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 697, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'M3nudg',                &
     &                      M3nudg(ng)/sec2day, (/0/), (/0/),           &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 703, MyFile)) RETURN
      DO itrc=1,NT(ng)
        nudg(itrc)=Tnudg(itrc,ng)/sec2day
      END DO
      CALL netcdf_put_fvar (ng, model, ncname, 'Tnudg',                 &
     &                      nudg, (/1/), (/NT(ng)/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 711, MyFile)) RETURN
!
!  Open boundary nudging, inverse time scales.
!
      IF (NudgingCoeff(ng)) THEN
        CALL netcdf_put_fvar (ng, model, ncname, 'FSobc_in',            &
     &                        FSobc_in(ng,:), (/1/), (/4/),             &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, 722, MyFile)) RETURN
        CALL netcdf_put_fvar (ng, model, ncname, 'FSobc_out',           &
     &                        FSobc_out(ng,:), (/1/), (/4/),            &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, 727, MyFile)) RETURN
        CALL netcdf_put_fvar (ng, model, ncname, 'M2obc_in',            &
     &                        M2obc_in(ng,:), (/1/), (/4/),             &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, 732, MyFile)) RETURN
        CALL netcdf_put_fvar (ng, model, ncname, 'M2obc_out',           &
     &                        M2obc_out(ng,:), (/1/), (/4/),            &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, 737, MyFile)) RETURN
        DO ibry=1,4
          DO itrc=1,NT(ng)
            Tobc(itrc,ibry)=Tobc_in(itrc,ng,ibry)
          END DO
        END DO
        CALL netcdf_put_fvar (ng, model, ncname, 'Tobc_in',             &
     &                        Tobc, (/1,1/), (/NT(ng),4/),              &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, 748, MyFile)) RETURN
        DO ibry=1,4
          DO itrc=1,NT(ng)
            Tobc(itrc,ibry)=Tobc_out(itrc,ng,ibry)
          END DO
        END DO
        CALL netcdf_put_fvar (ng, model, ncname, 'Tobc_out',            &
     &                        Tobc, (/1,1/), (/NT(ng),4/),              &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, 758, MyFile)) RETURN
        CALL netcdf_put_fvar (ng, model, ncname, 'M3obc_in',            &
     &                        M3obc_in(ng,:), (/1/), (/4/),             &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, 763, MyFile)) RETURN
        CALL netcdf_put_fvar (ng, model, ncname, 'M3obc_out',           &
     &                        M3obc_out(ng,:), (/1/), (/4/),            &
     &                      ncid = ncid)
        IF (FoundError(exit_flag, NoError, 768, MyFile)) RETURN
      END IF
!
!  Equation of State parameters.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'rho0',                  &
     &                      rho0, (/0/), (/0/),                         &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 778, MyFile)) RETURN
!
!  Slipperiness parameters.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'gamma2',                &
     &                      gamma2(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 826, MyFile)) RETURN
!
! Logical switches to activate horizontal momentum transport
! point Sources/Sinks (like river runoff transport) and mass point
! Sources/Sinks (like volume vertical influx).
!
      CALL netcdf_put_lvar (ng, model, ncname, 'LuvSrc',                &
     &                      LuvSrc(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 835, MyFile)) RETURN
      CALL netcdf_put_lvar (ng, model, ncname, 'LwSrc',                 &
     &                      LwSrc(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 840, MyFile)) RETURN
!
!  Logical switches to activate tracer point Sources/Sinks.
!
      CALL netcdf_put_lvar (ng, model, ncname, 'LtracerSrc',            &
     &                      LtracerSrc(:,ng), (/1/), (/NT(ng)/),        &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 849, MyFile)) RETURN
!
!  Logical switches to process climatology fields.
!
      CALL netcdf_put_lvar (ng, model, ncname, 'LsshCLM',               &
     &                      LsshCLM(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 857, MyFile)) RETURN
      CALL netcdf_put_lvar (ng, model, ncname, 'Lm2CLM',                &
     &                      Lm2CLM(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 862, MyFile)) RETURN
      CALL netcdf_put_lvar (ng, model, ncname, 'Lm3CLM',                &
     &                      Lm3CLM(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 868, MyFile)) RETURN
      CALL netcdf_put_lvar (ng, model, ncname, 'LtracerCLM',            &
     &                      LtracerCLM(:,ng), (/1/), (/NT(ng)/),        &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 873, MyFile)) RETURN
!
!  Logical switches for nudging climatology fields.
!
      CALL netcdf_put_lvar (ng, model, ncname, 'LnudgeM2CLM',           &
     &                      LnudgeM2CLM(ng), (/0/), (/0/),              &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 881, MyFile)) RETURN
      CALL netcdf_put_lvar (ng, model, ncname, 'LnudgeM3CLM',           &
     &                      LnudgeM3CLM(ng), (/0/), (/0/),              &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 887, MyFile)) RETURN
      CALL netcdf_put_lvar (ng, model, ncname, 'LnudgeTCLM',            &
     &                      LnudgeTCLM(:,ng), (/1/), (/NT(ng)/),        &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 892, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Write out grid variables.
!-----------------------------------------------------------------------
!
!  Grid type switch. Writing characters in parallel I/O is extremely
!  inefficient.  It is better to write this as an integer switch:
!  0=Cartesian, 1=spherical.
!
      CALL netcdf_put_lvar (ng, model, ncname, 'spherical',             &
     &                      spherical, (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1433, MyFile)) RETURN
!
!  Domain Length.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'xl',                    &
     &                      xl(ng), (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1440, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'el',                    &
     &                      el(ng), (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1445, MyFile)) RETURN
!
!  S-coordinate parameters.
!
      CALL netcdf_put_ivar (ng, model, ncname, 'Vtransform',            &
     &                      Vtransform(ng), (/0/), (/0/),               &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1454, MyFile)) RETURN
      CALL netcdf_put_ivar (ng, model, ncname, 'Vstretching',           &
     &                      Vstretching(ng), (/0/), (/0/),              &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1459, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'theta_s',               &
     &                      theta_s(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1464, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'theta_b',               &
     &                      theta_b(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1469, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'Tcline',                &
     &                      Tcline(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1474, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'hc',                    &
     &                      hc(ng), (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1479, MyFile)) RETURN
!
!  SGRID conventions for staggered data on structured grids. The value
!  is arbitrary but is set to unity so it can be used as logical during
!  post-processing.
!
      CALL netcdf_put_ivar (ng, model, ncname, 'grid',                  &
     &                      (/1/), (/0/), (/0/),                        &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1488, MyFile)) RETURN
!
!  S-coordinate non-dimensional independent variables.
!
      CALL netcdf_put_fvar (ng, model, ncname, 's_rho',                 &
     &                      SCALARS(ng)%sc_r(:), (/1/), (/N(ng)/),      &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1495, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 's_w',                   &
     &                      SCALARS(ng)%sc_w(0:), (/1/), (/N(ng)+1/),   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1500, MyFile)) RETURN
!
!  S-coordinate non-dimensional stretching curves.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'Cs_r',                  &
     &                      SCALARS(ng)%Cs_r(:), (/1/), (/N(ng)/),      &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1507, MyFile)) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'Cs_w',                  &
     &                      SCALARS(ng)%Cs_w(0:), (/1/), (/N(ng)+1/),   &
     &                      ncid = ncid)
      IF (FoundError(exit_flag, NoError, 1512, MyFile)) RETURN
!
!  Depths of horizontal slices.
!
      IF (Nslice.gt.0) THEN
        CALL netcdf_put_fvar (ng, model, ncname, 'z_slice',             &
     &                        Zslice, (/1/), (/Nslice/),                &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, 1521, MyFile)) RETURN
      END IF
!
!  User generic parameters.
!
      IF (Nuser.gt.0) THEN
        CALL netcdf_put_fvar (ng, model, ncname, 'user',                &
     &                        user, (/1/), (/Nuser/),                   &
     &                        ncid = ncid)
        IF (FoundError(exit_flag, NoError, 1530, MyFile)) RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Write out grid tiled variables.
!-----------------------------------------------------------------------
!
      GRID_VARS : IF (ncid.ne.FLT(ng)%ncid) THEN
!
!  Bathymetry.
!
        IF (exit_flag.eq.NoError) THEN
          scale=1.0_dp
          IF ((ncid.ne.STA(ng)%ncid).and.                               &
     &        (ncid.ne.XTR(ng)%ncid)) THEN
            IF (find_string(var_name, n_var, TRIM(Vname(1,idtopo)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idtopo,               &
     &                           varid, 0, r2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % rmask,                      &
     &                           GRID(ng) % h,                          &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, 1577, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idtopo)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idtopo)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
!  Coriolis parameter.
!
        IF (exit_flag.eq.NoError) THEN
          IF ((ncid.ne.STA(ng)%ncid).and.                               &
     &        (ncid.ne.XTR(ng)%ncid)) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, TRIM(Vname(1,idfcor)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idfcor,               &
     &                           varid, 0, r2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % rmask,                      &
     &                           GRID(ng) % f,                          &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, 1647, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idfcor)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idfcor)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
!  Curvilinear transformation metrics.
!
        IF (exit_flag.eq.NoError) THEN
          IF ((ncid.ne.STA(ng)%ncid).and.                               &
     &        (ncid.ne.XTR(ng)%ncid)) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, TRIM(Vname(1,idpmdx)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idpmdx,               &
     &                           varid, 0, r2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % rmask,                      &
     &                           GRID(ng) % pm,                         &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, 1705, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idpmdx)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idpmdx)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
        IF (exit_flag.eq.NoError) THEN
          IF ((ncid.ne.STA(ng)%ncid).and.                               &
     &        (ncid.ne.XTR(ng)%ncid)) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, TRIM(Vname(1,idpndy)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idpndy,               &
     &                           varid, 0, r2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % rmask,                      &
     &                           GRID(ng) % pn,                         &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, 1761, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idpndy)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idpndy)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
!  Grid coordinates of RHO-points.
!
        IF (spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            scale=1.0_dp
            IF ((ncid.ne.STA(ng)%ncid).and.                             &
     &          (ncid.ne.XTR(ng)%ncid)) THEN
              IF (find_string(var_name, n_var, TRIM(Vname(1,idLonR)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idLonR,             &
     &                             varid, 0, r2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % rmask,                    &
     &                             GRID(ng) % lonr,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         1821, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLonR)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLonR)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
!
          IF (exit_flag.eq.NoError) THEN
            scale=1.0_dp
            IF ((ncid.ne.STA(ng)%ncid).and.                             &
     &          (ncid.ne.XTR(ng)%ncid)) THEN
              IF (find_string(var_name, n_var, TRIM(Vname(1,idLatR)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idLatR,             &
     &                             varid, 0, r2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % rmask,                    &
     &                             GRID(ng) % latr,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         1891, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLatR)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLatR)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
        IF (.not.spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            scale=1.0_dp
            IF ((ncid.ne.STA(ng)%ncid).and.                             &
     &          (ncid.ne.XTR(ng)%ncid)) THEN
              IF (find_string(var_name, n_var, TRIM(Vname(1,idXgrR)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idXgrR,             &
     &                             varid, 0, r2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % rmask,                    &
     &                             GRID(ng) % xr,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         1963, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idXgrR)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idXgrR)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
!
          IF (exit_flag.eq.NoError) THEN
            scale=1.0_dp
            IF ((ncid.ne.STA(ng)%ncid).and.                             &
     &          (ncid.ne.XTR(ng)%ncid)) THEN
              IF (find_string(var_name, n_var, TRIM(Vname(1,idYgrR)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idYgrR,             &
     &                             varid, 0, r2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % rmask,                    &
     &                             GRID(ng) % yr,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         2033, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idYgrR)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idYgrR)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
!  Grid coordinates of U-points.
!
        IF (spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF ((ncid.ne.STA(ng)%ncid).and.                             &
     &          (ncid.ne.XTR(ng)%ncid)) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idLonU)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idLonU,             &
     &                             varid, 0, u2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % umask,                    &
     &                             GRID(ng) % lonu,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         2107, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLonU)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLonU)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
!
          IF (exit_flag.eq.NoError) THEN
            IF ((ncid.ne.STA(ng)%ncid).and.                             &
     &          (ncid.ne.XTR(ng)%ncid)) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idLatU)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idLatU,             &
     &                             varid, 0, u2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % umask,                    &
     &                             GRID(ng) % latu,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         2165, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLatU)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLatU)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
        IF (.not.spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF ((ncid.ne.STA(ng)%ncid).and.                             &
     &          (ncid.ne.XTR(ng)%ncid)) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idXgrU)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idXgrU,             &
     &                             varid, 0, u2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % umask,                    &
     &                             GRID(ng) % xu,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         2225, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idXgrU)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idXgrU)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
!
          IF (exit_flag.eq.NoError) THEN
            IF ((ncid.ne.STA(ng)%ncid).and.                             &
     &          (ncid.ne.XTR(ng)%ncid)) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idYgrU)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idYgrU,             &
     &                             varid, 0, u2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % umask,                    &
     &                             GRID(ng) % yu,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         2283, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idYgrU)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idYgrU)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
!  Grid coordinates of V-points.
!
        IF (spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF ((ncid.ne.STA(ng)%ncid).and.                             &
     &          (ncid.ne.XTR(ng)%ncid)) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idLonV)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idLonV,             &
     &                             varid, 0, v2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % vmask,                    &
     &                             GRID(ng) % lonv,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         2345, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLonV)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLonV)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
!
          IF (exit_flag.eq.NoError) THEN
            IF ((ncid.ne.STA(ng)%ncid).and.                             &
     &          (ncid.ne.XTR(ng)%ncid)) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idLatV)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idLatV,             &
     &                             varid, 0, v2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % vmask,                    &
     &                             GRID(ng) % latv,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         2403, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLatV)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLatV)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
        IF (.not.spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF ((ncid.ne.STA(ng)%ncid).and.                             &
     &          (ncid.ne.XTR(ng)%ncid)) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idXgrV)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idXgrV,             &
     &                             varid, 0, v2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % vmask,                    &
     &                             GRID(ng) % xv,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         2463, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idXgrV)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idXgrV)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
!
          IF (exit_flag.eq.NoError) THEN
            IF ((ncid.ne.STA(ng)%ncid).and.                             &
     &          (ncid.ne.XTR(ng)%ncid)) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idYgrV)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idYgrV,             &
     &                             varid, 0, v2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % vmask,                    &
     &                             GRID(ng) % yv,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         2521, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idYgrV)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idYgrV)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
!  Grid coordinates of PSI-points.
!
        IF (spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF ((ncid.ne.STA(ng)%ncid).and.                             &
     &          (ncid.ne.XTR(ng)%ncid)) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idLonP)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idLonP,             &
     &                             varid, 0, p2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % pmask,                    &
     &                             GRID(ng) % lonp,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         2583, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLonP)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLonP)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
          IF (exit_flag.eq.NoError) THEN
            IF ((ncid.ne.STA(ng)%ncid).and.                             &
     &          (ncid.ne.XTR(ng)%ncid)) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idLatP)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idLatP,             &
     &                             varid, 0, p2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % pmask,                    &
     &                             GRID(ng) % latp,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         2641, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLatP)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLatP)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
        IF (.not.spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF ((ncid.ne.STA(ng)%ncid).and.                             &
     &          (ncid.ne.XTR(ng)%ncid)) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idXgrP)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idXgrP,             &
     &                             varid, 0, p2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % pmask,                    &
     &                             GRID(ng) % xp,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         2701, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idXgrP)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idXgrP)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
!
          IF (exit_flag.eq.NoError) THEN
            IF ((ncid.ne.STA(ng)%ncid).and.                             &
     &          (ncid.ne.XTR(ng)%ncid)) THEN
              scale=1.0_dp
              IF (find_string(var_name, n_var, TRIM(Vname(1,idYgrP)),   &
     &                        varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, idYgrP,             &
     &                             varid, 0, p2dvar,                    &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % pmask,                    &
     &                             GRID(ng) % yp,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, nf90_noerr,                      &
     &                         2759, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idYgrP)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idYgrP)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
!  Angle between XI-axis and EAST at RHO-points.
!
        IF (exit_flag.eq.NoError) THEN
          scale=1.0_dp
          IF ((ncid.ne.STA(ng)%ncid).and.                               &
     &        (ncid.ne.XTR(ng)%ncid)) THEN
            IF (find_string(var_name, n_var, TRIM(Vname(1,idangR)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idangR,               &
     &                           varid, 0, r2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % rmask,                      &
     &                           GRID(ng) % angler,                     &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, 2821, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idangR)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idangR)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
!  Masking fields at RHO-, U-, V-points, and PSI-points.
!
        IF (exit_flag.eq.NoError) THEN
          IF ((ncid.ne.STA(ng)%ncid).and.                               &
     &        (ncid.ne.XTR(ng)%ncid)) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, TRIM(Vname(1,idmskR)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idmskR,               &
     &                           varid, 0, r2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % rmask,                      &
     &                           GRID(ng) % rmask,                      &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, 2891, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idmskR)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idmskR)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
        IF (exit_flag.eq.NoError) THEN
          IF ((ncid.ne.STA(ng)%ncid).and.                               &
     &        (ncid.ne.XTR(ng)%ncid)) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, TRIM(Vname(1,idmskU)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idmskU,               &
     &                           varid, 0, u2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % umask,                      &
     &                           GRID(ng) % umask,                      &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, 2943, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idmskU)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idmskU)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
        IF (exit_flag.eq.NoError) THEN
          IF ((ncid.ne.STA(ng)%ncid).and.                               &
     &        (ncid.ne.XTR(ng)%ncid)) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, TRIM(Vname(1,idmskV)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idmskV,               &
     &                           varid, 0, v2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % vmask,                      &
     &                           GRID(ng) % vmask,                      &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, 2995, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idmskV)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idmskV)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
        IF (exit_flag.eq.NoError) THEN
          IF ((ncid.ne.STA(ng)%ncid).and.                               &
     &        (ncid.ne.XTR(ng)%ncid)) THEN
            scale=1.0_dp
            IF (find_string(var_name, n_var, TRIM(Vname(1,idmskP)),     &
     &                      varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, idmskP,               &
     &                           varid, 0, p2dvar,                      &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % pmask,                      &
     &                           GRID(ng) % pmask,                      &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, nf90_noerr, 3047, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idmskP)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idmskP)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
      END IF GRID_VARS
!
!-----------------------------------------------------------------------
!  Synchronize NetCDF file to disk to allow other processes to access
!  data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, model, ncname, ncid)
      IF (FoundError(exit_flag, NoError, 3443, MyFile)) RETURN
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
!
  10  FORMAT (/,' WRT_INFO_NF90 - error while writing variable: ',a,/,  &
     &        17x,'into file: ',a)
  20  FORMAT (/,' WRT_INFO_NF90 - error while inquiring ID for',        &
     &        ' variable: ',a,/,17x,'in file: ',a)
  30  FORMAT (/,' WRT_INFO_NF90 - unable to synchronize to disk file:', &
     &        /,17x,a)
!
      RETURN
      END SUBROUTINE wrt_info_nf90
!
!***********************************************************************
      SUBROUTINE wrt_info_pio (ng, model, pioFile, ncname)
!***********************************************************************
!                                                                      !
!  This routine writes out information variables into requested        !
!  NetCDF file using the standard NetCDF-3 or NetCDF-4 library.        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     pioFile      PIO file descriptor structure, TYPE(File_desc_t)    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     ncname       PIO filename (string)                               !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
      character (len=*), intent(in) :: ncname
!
      TYPE (File_desc_t), intent(inout) :: pioFile
!
!  Local variable declarations.
!
      logical :: Cgrid = .TRUE.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: i, j, k, ibry, ilev, itrc, status
      integer :: ival
      integer :: FileH, MY_FOUT
      integer :: ifield = 0
!
      real(dp) :: scale
      real(r8), dimension(NT(ng)) :: nudg
      real(r8), dimension(NT(ng),4) :: Tobc
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/wrt_info.F"//", wrt_info_pio"
!
      TYPE (IO_desc_t), pointer :: ioDesc
      TYPE (My_VarDesc) :: pioVar
!
      SourceFile=MyFile
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out running parameters.
!-----------------------------------------------------------------------
!
!  Get NetCDF file handler from descriptor.
!
      FileH=ABS(pioFile%fh)
!
!  Time stepping parameters.
!
      CALL pio_netcdf_put_ivar (ng, model, ncname, 'ntimes',            &
     &                          ntimes(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 3550, MyFile)) RETURN
      CALL pio_netcdf_put_ivar (ng, model, ncname, 'ndtfast',           &
     &                          ndtfast(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 3555, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'dt',                &
     &                          dt(ng), (/0/), (/0/),                   &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 3560, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'dtfast',            &
     &                          dtfast(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 3565, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'dstart',            &
     &                          dstart, (/0/), (/0/),                   &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 3570, MyFile)) RETURN
      CALL pio_netcdf_put_ivar (ng, model, ncname, 'nHIS',              &
     &                          nHIS(ng), (/0/), (/0/),                 &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 3604, MyFile)) RETURN
      CALL pio_netcdf_put_ivar (ng, model, ncname, 'ndefHIS',           &
     &                          ndefHIS(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 3609, MyFile)) RETURN
      CALL pio_netcdf_put_ivar (ng, model, ncname, 'nRST',              &
     &                          nRST(ng), (/0/), (/0/),                 &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 3631, MyFile)) RETURN
      CALL pio_netcdf_put_ivar (ng, model, ncname, 'ntsAVG',            &
     &                          ntsAVG(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 3640, MyFile)) RETURN
      CALL pio_netcdf_put_ivar (ng, model, ncname, 'nAVG',              &
     &                          nAVG(ng), (/0/), (/0/),                 &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 3645, MyFile)) RETURN
      CALL pio_netcdf_put_ivar (ng, model, ncname, 'ndefAVG',           &
     &                          ndefAVG(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 3650, MyFile)) RETURN
!
!  Power-law shape filter parameters for time-averaging of barotropic
!  fields.
!
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Falpha',            &
     &                          Falpha, (/0/), (/0/),                   &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 3772, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Fbeta',             &
     &                          Fbeta, (/0/), (/0/),                    &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 3777, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Fgamma',            &
     &                          Fgamma, (/0/), (/0/),                   &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 3782, MyFile)) RETURN
!
!  Horizontal mixing coefficients.
!
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'nl_tnu2',           &
     &                          nl_tnu2(:,ng), (/1/), (/NT(ng)/),       &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 3791, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'nl_visc2',          &
     &                          nl_visc2(ng), (/0/), (/0/),             &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 3843, MyFile)) RETURN
      CALL pio_netcdf_put_lvar (ng, model, ncname, 'LuvSponge',         &
     &                          LuvSponge(ng), (/0/), (/0/),            &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 3901, MyFile)) RETURN
      CALL pio_netcdf_put_lvar (ng, model, ncname, 'LtracerSponge',     &
     &                          LtracerSponge(:,ng), (/1/), (/NT(ng)/), &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 3908, MyFile)) RETURN
!
!  Background vertical mixing coefficients.
!
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Akt_bak',           &
     &                          Akt_bak(:,ng), (/1/), (/NT(ng)/),       &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 3918, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Akv_bak',           &
     &                          Akv_bak(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 3923, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Akk_bak',           &
     &                          Akk_bak(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 3929, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Akp_bak',           &
     &                          Akp_bak(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 3934, MyFile)) RETURN
!
!  Drag coefficients.
!
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'rdrg',              &
     &                          rdrg(ng), (/0/), (/0/),                 &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 3976, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'rdrg2',             &
     &                          rdrg2(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 3981, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Zob',               &
     &                          Zob(ng), (/0/), (/0/),                  &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 3987, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Zos',               &
     &                          Zos(ng), (/0/), (/0/),                  &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 3992, MyFile)) RETURN
!
!  Generic length-scale parameters.
!
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'gls_p',             &
     &                          gls_p(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4002, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'gls_m',             &
     &                          gls_m(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4007, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'gls_n',             &
     &                          gls_n(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4012, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'gls_cmu0',          &
     &                          gls_cmu0(ng), (/0/), (/0/),             &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4017, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'gls_c1',            &
     &                          gls_c1(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4022, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'gls_c2',            &
     &                          gls_c2(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4027, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'gls_c3m',           &
     &                          gls_c3m(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4032, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'gls_c3p',           &
     &                          gls_c3p(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4037, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'gls_sigk',          &
     &                          gls_sigk(ng), (/0/), (/0/),             &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4042, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'gls_sigp',          &
     &                          gls_sigp(ng), (/0/), (/0/),             &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4047, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'gls_Kmin',          &
     &                          gls_Kmin(ng), (/0/), (/0/),             &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4052, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'gls_Pmin',          &
     &                          gls_Pmin(ng), (/0/), (/0/),             &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4057, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Charnok_alpha',     &
     &                          charnok_alpha(ng), (/0/), (/0/),        &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4062, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Zos_hsig_alpha',    &
     &                          zos_hsig_alpha(ng), (/0/), (/0/),       &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4067, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'sz_alpha',          &
     &                          sz_alpha(ng), (/0/), (/0/),             &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4072, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'CrgBan_cw',         &
     &                          crgban_cw(ng), (/0/), (/0/),            &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4077, MyFile)) RETURN
!
!  Nudging inverse time scales used in various tasks.
!
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Znudg',             &
     &                          Znudg(ng)/sec2day, (/0/), (/0/),        &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4091, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'M2nudg',            &
     &                          M2nudg(ng)/sec2day, (/0/), (/0/),       &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4096, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'M3nudg',            &
     &                          M3nudg(ng)/sec2day, (/0/), (/0/),       &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4102, MyFile)) RETURN
      DO itrc=1,NT(ng)
        nudg(itrc)=Tnudg(itrc,ng)/sec2day
      END DO
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Tnudg',             &
     &                          nudg, (/1/), (/NT(ng)/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4110, MyFile)) RETURN
!
!  Open boundary nudging, inverse time scales.
!
      IF (NudgingCoeff(ng)) THEN
        CALL pio_netcdf_put_fvar (ng, model, ncname, 'FSobc_in',        &
     &                            FSobc_in(ng,:), (/1/), (/4/),         &
     &                            pioFile = pioFile)
        IF (FoundError(exit_flag, NoError, 4121, MyFile)) RETURN
        CALL pio_netcdf_put_fvar (ng, model, ncname, 'FSobc_out',       &
     &                            FSobc_out(ng,:), (/1/), (/4/),        &
     &                            pioFile = pioFile)
        IF (FoundError(exit_flag, NoError, 4126, MyFile)) RETURN
        CALL pio_netcdf_put_fvar (ng, model, ncname, 'M2obc_in',        &
     &                            M2obc_in(ng,:), (/1/), (/4/),         &
     &                            pioFile = pioFile)
        IF (FoundError(exit_flag, NoError, 4131, MyFile)) RETURN
        CALL pio_netcdf_put_fvar (ng, model, ncname, 'M2obc_out',       &
     &                            M2obc_out(ng,:), (/1/), (/4/),        &
     &                            pioFile = pioFile)
        IF (FoundError(exit_flag, NoError, 4136, MyFile)) RETURN
        DO ibry=1,4
          DO itrc=1,NT(ng)
            Tobc(itrc,ibry)=Tobc_in(itrc,ng,ibry)
          END DO
        END DO
        CALL pio_netcdf_put_fvar (ng, model, ncname, 'Tobc_in',         &
     &                            Tobc, (/1,1/), (/NT(ng),4/),          &
     &                            pioFile = pioFile)
        IF (FoundError(exit_flag, NoError, 4147, MyFile)) RETURN
        DO ibry=1,4
          DO itrc=1,NT(ng)
            Tobc(itrc,ibry)=Tobc_out(itrc,ng,ibry)
          END DO
        END DO
        CALL pio_netcdf_put_fvar (ng, model, ncname, 'Tobc_out',        &
     &                            Tobc, (/1,1/), (/NT(ng),4/),          &
     &                            pioFile = pioFile)
        IF (FoundError(exit_flag, NoError, 4157, MyFile)) RETURN
        CALL pio_netcdf_put_fvar (ng, model, ncname, 'M3obc_in',        &
     &                            M3obc_in(ng,:), (/1/), (/4/),         &
     &                            pioFile = pioFile)
        IF (FoundError(exit_flag, NoError, 4162, MyFile)) RETURN
        CALL pio_netcdf_put_fvar (ng, model, ncname, 'M3obc_out',       &
     &                            M3obc_out(ng,:), (/1/), (/4/),        &
     &                            pioFile = pioFile)
        IF (FoundError(exit_flag, NoError, 4167, MyFile)) RETURN
      END IF
!
!  Equation of State parameters.
!
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'rho0',              &
     &                          rho0, (/0/), (/0/),                     &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4177, MyFile)) RETURN
!
!  Slipperiness parameters.
!
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'gamma2',            &
     &                          gamma2(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4225, MyFile)) RETURN
!
! Logical switches to activate horizontal momentum transport
! point Sources/Sinks (like river runoff transport) and mass point
! Sources/Sinks (like volume vertical influx).
!
      CALL pio_netcdf_put_lvar (ng, model, ncname, 'LuvSrc',            &
     &                          LuvSrc(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4234, MyFile)) RETURN
      CALL pio_netcdf_put_lvar (ng, model, ncname, 'LwSrc',             &
     &                          LwSrc(ng), (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4239, MyFile)) RETURN
!
!  Logical switches to activate tracer point Sources/Sinks.
!
      CALL pio_netcdf_put_lvar (ng, model, ncname, 'LtracerSrc',        &
     &                          LtracerSrc(:,ng), (/1/), (/NT(ng)/),    &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4248, MyFile)) RETURN
!
!  Logical switches to process climatology fields.
!
      CALL pio_netcdf_put_lvar (ng, model, ncname, 'LsshCLM',           &
     &                          LsshCLM(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4256, MyFile)) RETURN
      CALL pio_netcdf_put_lvar (ng, model, ncname, 'Lm2CLM',            &
     &                          Lm2CLM(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4261, MyFile)) RETURN
      CALL pio_netcdf_put_lvar (ng, model, ncname, 'Lm3CLM',            &
     &                          Lm3CLM(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4267, MyFile)) RETURN
      CALL pio_netcdf_put_lvar (ng, model, ncname, 'LtracerCLM',        &
     &                          LtracerCLM(:,ng), (/1/), (/NT(ng)/),    &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4272, MyFile)) RETURN
!
!  Logical switches for nudging climatology fields.
!
      CALL pio_netcdf_put_lvar (ng, model, ncname, 'LnudgeM2CLM',       &
     &                          LnudgeM2CLM(ng), (/0/), (/0/),          &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4280, MyFile)) RETURN
      CALL pio_netcdf_put_lvar (ng, model, ncname, 'LnudgeM3CLM',       &
     &                          LnudgeM3CLM(ng), (/0/), (/0/),          &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4286, MyFile)) RETURN
      CALL pio_netcdf_put_lvar (ng, model, ncname, 'LnudgeTCLM',        &
     &                          LnudgeTCLM(:,ng), (/1/), (/NT(ng)/),    &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4291, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Write out grid variables.
!-----------------------------------------------------------------------
!
!  Grid type switch. Writing characters in parallel I/O is extremely
!  inefficient.  It is better to write this as an integer switch:
!  0=Cartesian, 1=spherical.
!
      CALL pio_netcdf_put_lvar (ng, model, ncname, 'spherical',         &
     &                          spherical, (/0/), (/0/),                &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4843, MyFile)) RETURN
!
!  Domain Length.
!
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'xl',                &
     &                          xl(ng), (/0/), (/0/),                   &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4850, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'el',                &
     &                          el(ng), (/0/), (/0/),                   &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4855, MyFile)) RETURN
!
!  S-coordinate parameters.
!
      CALL pio_netcdf_put_ivar (ng, model, ncname, 'Vtransform',        &
     &                          Vtransform(ng), (/0/), (/0/),           &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4864, MyFile)) RETURN
      CALL pio_netcdf_put_ivar (ng, model, ncname, 'Vstretching',       &
     &                          Vstretching(ng), (/0/), (/0/),          &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4869, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'theta_s',           &
     &                          theta_s(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4874, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'theta_b',           &
     &                          theta_b(ng), (/0/), (/0/),              &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4879, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Tcline',            &
     &                          Tcline(ng), (/0/), (/0/),               &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4884, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'hc',                &
     &                          hc(ng), (/0/), (/0/),                   &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4889, MyFile)) RETURN
!
!  SGRID conventions for staggered data on structured grids. The value
!  is arbitrary but is set to unity so it can be used as logical during
!  post-processing.
!
      ival=1
      CALL pio_netcdf_put_ivar (ng, model, ncname, 'grid',              &
     &                          ival, (/0/), (/0/),                     &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4899, MyFile)) RETURN
!
!  S-coordinate non-dimensional independent variables.
!
      CALL pio_netcdf_put_fvar (ng, model, ncname, 's_rho',             &
     &                          SCALARS(ng)%sc_r(:), (/1/), (/N(ng)/),  &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4906, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 's_w',               &
     &                          SCALARS(ng)%sc_w(0:),                   &
     &                          (/1/), (/N(ng)+1/),                     &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4912, MyFile)) RETURN
!
!  S-coordinate non-dimensional stretching curves.
!
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Cs_r',              &
     &                          SCALARS(ng)%Cs_r(:), (/1/), (/N(ng)/),  &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4919, MyFile)) RETURN
      CALL pio_netcdf_put_fvar (ng, model, ncname, 'Cs_w',              &
     &                          SCALARS(ng)%Cs_w(0:),                   &
     &                          (/1/), (/N(ng)+1/),                     &
     &                          pioFile = pioFile)
      IF (FoundError(exit_flag, NoError, 4925, MyFile)) RETURN
!
!  Depths of horizontal slices.
!
      IF (Nslice.gt.0) THEN
        CALL pio_netcdf_put_fvar (ng, model, ncname, 'z_slice',         &
     &                            Zslice, (/1/), (/Nslice/),            &
     &                            pioFile = pioFile)
        IF (FoundError(exit_flag, NoError, 4934, MyFile)) RETURN
      END IF
!
!  User generic parameters.
!
      IF (Nuser.gt.0) THEN
        CALL pio_netcdf_put_fvar (ng, model, ncname, 'user',            &
     &                            user, (/1/), (/Nuser/),               &
     &                            pioFile = pioFile)
        IF (FoundError(exit_flag, NoError, 4943, MyFile)) RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Write out grid tiled variables.
!-----------------------------------------------------------------------
!
      GRID_VARS : IF (FileH.ne.ABS(FLT(ng)%pioFile%fh)) THEN
!
!  Bathymetry.
!
        IF (exit_flag.eq.NoError) THEN
          scale=1.0_dp
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            IF (pio_netcdf_find_var(ng, model, pioFile,                 &
     &                              TRIM(Vname(1,idtopo)),              &
     &                              pioVar%vd)) THEN
              pioVar%gtype=r2dvar
              IF (PIO_TYPE.eq.PIO_double) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_r2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_r2dvar(ng)
              END IF
              status=nf_fwrite2d(ng, model, pioFile, idtopo,            &
     &                           pioVar, 0, ioDesc,                     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % rmask,                      &
     &                           GRID(ng) % h,                          &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, PIO_noerr, 4998, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idtopo)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idtopo)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
!  Coriolis parameter.
!
        IF (exit_flag.eq.NoError) THEN
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            scale=1.0_dp
            IF (pio_netcdf_find_var(ng, model, pioFile,                 &
     &                              TRIM(Vname(1,idfcor)),              &
     &                              pioVar%vd)) THEN
              pioVar%gtype=r2dvar
              IF (PIO_TYPE.eq.PIO_double) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_r2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_r2dvar(ng)
              END IF
              status=nf_fwrite2d(ng, model, pioFile, idfcor,            &
     &                           pioVar, 0, ioDesc,                     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % rmask,                      &
     &                           GRID(ng) % f,                          &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, PIO_noerr, 5051, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idfcor)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idfcor)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
!  Curvilinear transformation metrics.
!
        IF (exit_flag.eq.NoError) THEN
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            scale=1.0_dp
            IF (pio_netcdf_find_var(ng, model, pioFile,                 &
     &                              TRIM(Vname(1,idpmdx)),              &
     &                              pioVar%vd)) THEN
              pioVar%gtype=r2dvar
              IF (PIO_TYPE.eq.PIO_double) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_r2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_r2dvar(ng)
              END IF
              status=nf_fwrite2d(ng, model, pioFile, idpmdx,            &
     &                           pioVar, 0, ioDesc,                     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % rmask,                      &
     &                           GRID(ng) % pm,                         &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, PIO_noerr, 5090, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idpmdx)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idpmdx)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
        IF (exit_flag.eq.NoError) THEN
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            scale=1.0_dp
            IF (pio_netcdf_find_var(ng, model, pioFile,                 &
     &                              TRIM(Vname(1,idpndy)),              &
     &                              pioVar%vd)) THEN
              pioVar%gtype=r2dvar
              IF (PIO_TYPE.eq.PIO_double) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_r2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_r2dvar(ng)
              END IF
              status=nf_fwrite2d(ng, model, pioFile, idpndy,            &
     &                           pioVar, 0, ioDesc,                     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % rmask,                      &
     &                           GRID(ng) % pn,                         &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, PIO_noerr, 5127, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idpndy)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idpndy)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
!  Grid coordinates of RHO-points.
!
        IF (spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            scale=1.0_dp
            IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
              IF (pio_netcdf_find_var(ng, model, pioFile,               &
     &                                TRIM(Vname(1,idLonR)),            &
     &                                pioVar%vd)) THEN
                pioVar%gtype=r2dvar
                IF (PIO_TYPE.eq.PIO_double) THEN
                  pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_r2dvar(ng)
                ELSE
                  pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_r2dvar(ng)
                END IF
                status=nf_fwrite2d(ng, model, pioFile, idLonR,          &
     &                             pioVar, 0, ioDesc,                   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % rmask,                    &
     &                             GRID(ng) % lonr,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, PIO_noerr,                       &
     &                         5168, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLonR)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLonR)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
!
          IF (exit_flag.eq.NoError) THEN
            scale=1.0_dp
            IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
              IF (pio_netcdf_find_var(ng, model, pioFile,               &
     &                                TRIM(Vname(1,idLatR)),            &
     &                                pioVar%vd)) THEN
                pioVar%gtype=r2dvar
                IF (PIO_TYPE.eq.PIO_double) THEN
                  pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_r2dvar(ng)
                ELSE
                  pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_r2dvar(ng)
                END IF
                status=nf_fwrite2d(ng, model, pioFile, idLatR,          &
     &                             pioVar, 0, ioDesc,                   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % rmask,                    &
     &                             GRID(ng) % latr,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, PIO_noerr,                       &
     &                         5220, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLatR)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLatR)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
        IF (.not.spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            scale=1.0_dp
            IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
              IF (pio_netcdf_find_var(ng, model, pioFile,               &
     &                                TRIM(Vname(1,idXgrR)),            &
     &                                pioVar%vd)) THEN
                pioVar%gtype=r2dvar
                IF (PIO_TYPE.eq.PIO_double) THEN
                  pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_r2dvar(ng)
                ELSE
                  pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_r2dvar(ng)
                END IF
                status=nf_fwrite2d(ng, model, pioFile, idXgrR,          &
     &                             pioVar, 0, ioDesc,                   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % rmask,                    &
     &                             GRID(ng) % xr,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, PIO_noerr,                       &
     &                         5274, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idXgrR)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idXgrR)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
!
          IF (exit_flag.eq.NoError) THEN
            scale=1.0_dp
            IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
              IF (pio_netcdf_find_var(ng, model, pioFile,               &
     &                                TRIM(Vname(1,idYgrR)),            &
     &                                pioVar%vd)) THEN
                pioVar%gtype=r2dvar
                IF (PIO_TYPE.eq.PIO_double) THEN
                  pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_r2dvar(ng)
                ELSE
                  pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_r2dvar(ng)
                END IF
                status=nf_fwrite2d(ng, model, pioFile, idYgrR,          &
     &                             pioVar, 0, ioDesc,                   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % rmask,                    &
     &                             GRID(ng) % yr,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, PIO_noerr,                       &
     &                         5326, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idYgrR)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idYgrR)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
!  Grid coordinates of U-points.
!
        IF (spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
              scale=1.0_dp
              IF (pio_netcdf_find_var(ng, model, pioFile,               &
     &                                TRIM(Vname(1,idLonU)),            &
     &                                pioVar%vd)) THEN
                pioVar%gtype=u2dvar
                IF (PIO_TYPE.eq.PIO_double) THEN
                  pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_u2dvar(ng)
                ELSE
                  pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_u2dvar(ng)
                END IF
                status=nf_fwrite2d(ng, model, pioFile, idLonU,          &
     &                             pioVar, 0, ioDesc,                   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % umask,                    &
     &                             GRID(ng) % lonu,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, PIO_noerr,                       &
     &                         5382, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLonU)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLonU)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
!
          IF (exit_flag.eq.NoError) THEN
            IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
              scale=1.0_dp
              IF (pio_netcdf_find_var(ng, model, pioFile,               &
     &                                TRIM(Vname(1,idLatU)),            &
     &                                pioVar%vd)) THEN
                pioVar%gtype=u2dvar
                IF (PIO_TYPE.eq.PIO_double) THEN
                  pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_u2dvar(ng)
                ELSE
                  pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_u2dvar(ng)
                END IF
                status=nf_fwrite2d(ng, model, pioFile, idLatU,          &
     &                             pioVar, 0, ioDesc,                   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % umask,                    &
     &                             GRID(ng) % latu,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, PIO_noerr,                       &
     &                         5420, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLatU)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLatU)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
        IF (.not.spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
              scale=1.0_dp
              IF (pio_netcdf_find_var(ng, model, pioFile,               &
     &                                TRIM(Vname(1,idXgrU)),            &
     &                                pioVar%vd)) THEN
                pioVar%gtype=u2dvar
                IF (PIO_TYPE.eq.PIO_double) THEN
                  pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_u2dvar(ng)
                ELSE
                  pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_u2dvar(ng)
                END IF
                status=nf_fwrite2d(ng, model, pioFile, idXgrU,          &
     &                             pioVar, 0, ioDesc,                   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % umask,                    &
     &                             GRID(ng) % xu,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, PIO_noerr,                       &
     &                         5460, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idXgrU)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idXgrU)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
!
          IF (exit_flag.eq.NoError) THEN
            IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
              scale=1.0_dp
              IF (pio_netcdf_find_var(ng, model, pioFile,               &
     &                                TRIM(Vname(1,idYgrU)),            &
     &                                pioVar%vd)) THEN
                pioVar%gtype=u2dvar
                IF (PIO_TYPE.eq.PIO_double) THEN
                  pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_u2dvar(ng)
                ELSE
                  pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_u2dvar(ng)
                END IF
                status=nf_fwrite2d(ng, model, pioFile, idYgrU,          &
     &                             pioVar, 0, ioDesc,                   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % umask,                    &
     &                             GRID(ng) % yu,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, PIO_noerr,                       &
     &                         5498, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idXgrU)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idXgrU)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
!  Grid coordinates of V-points.
!
        IF (spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
              scale=1.0_dp
              IF (pio_netcdf_find_var(ng, model, pioFile,               &
     &                                TRIM(Vname(1,idLonV)),            &
     &                                pioVar%vd)) THEN
                pioVar%gtype=v2dvar
                IF (PIO_TYPE.eq.PIO_double) THEN
                  pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_v2dvar(ng)
                ELSE
                  pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_v2dvar(ng)
                END IF
                status=nf_fwrite2d(ng, model, pioFile, idLonV,          &
     &                             pioVar, 0, ioDesc,                   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % vmask,                    &
     &                             GRID(ng) % lonv,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, PIO_noerr,                       &
     &                         5540, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLonV)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLonV)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
!
          IF (exit_flag.eq.NoError) THEN
            IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
              scale=1.0_dp
              IF (pio_netcdf_find_var(ng, model, pioFile,               &
     &                                TRIM(Vname(1,idLatV)),            &
     &                                pioVar%vd)) THEN
                pioVar%gtype=v2dvar
                IF (PIO_TYPE.eq.PIO_double) THEN
                  pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_v2dvar(ng)
                ELSE
                  pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_v2dvar(ng)
                END IF
                status=nf_fwrite2d(ng, model, pioFile, idLatV,          &
     &                             pioVar, 0, ioDesc,                   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % vmask,                    &
     &                             GRID(ng) % latv,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, PIO_noerr,                       &
     &                         5578, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLatV)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLatV)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
        IF (.not.spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
              scale=1.0_dp
              IF (pio_netcdf_find_var(ng, model, pioFile,               &
     &                                TRIM(Vname(1,idXgrV)),            &
     &                                pioVar%vd)) THEN
                pioVar%gtype=v2dvar
                IF (PIO_TYPE.eq.PIO_double) THEN
                  pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_v2dvar(ng)
                ELSE
                  pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_v2dvar(ng)
                END IF
                status=nf_fwrite2d(ng, model, pioFile, idXgrV,          &
     &                             pioVar, 0, ioDesc,                   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % vmask,                    &
     &                             GRID(ng) % xv,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, PIO_noerr,                       &
     &                         5618, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idXgrV)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idXgrV)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
!
          IF (exit_flag.eq.NoError) THEN
            IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
              scale=1.0_dp
              IF (pio_netcdf_find_var(ng, model, pioFile,               &
     &                                TRIM(Vname(1,idYgrV)),            &
     &                                pioVar%vd)) THEN
                pioVar%gtype=v2dvar
                IF (PIO_TYPE.eq.PIO_double) THEN
                  pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_v2dvar(ng)
                ELSE
                  pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_v2dvar(ng)
                END IF
                status=nf_fwrite2d(ng, model, pioFile, idYgrV,          &
     &                             pioVar, 0, ioDesc,                   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % vmask,                    &
     &                             GRID(ng) % yv,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, PIO_noerr,                       &
     &                         5656, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idYgrV)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idYgrV)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
!  Grid coordinates of PSI-points.
!
        IF (spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
              scale=1.0_dp
              IF (pio_netcdf_find_var(ng, model, pioFile,               &
     &                                TRIM(Vname(1,idLonP)),            &
     &                                pioVar%vd)) THEN
                pioVar%gtype=p2dvar
                IF (PIO_TYPE.eq.PIO_double) THEN
                  pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_p2dvar(ng)
                ELSE
                  pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_p2dvar(ng)
                END IF
                status=nf_fwrite2d(ng, model, pioFile, idLonP,          &
     &                             pioVar, 0, ioDesc,                   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % pmask,                    &
     &                             GRID(ng) % lonp,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, PIO_noerr,                       &
     &                         5698, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLonP)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLonP)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
!
          IF (exit_flag.eq.NoError) THEN
            IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
              scale=1.0_dp
              IF (pio_netcdf_find_var(ng, model, pioFile,               &
     &                                TRIM(Vname(1,idLatP)),            &
     &                                pioVar%vd)) THEN
                pioVar%gtype=p2dvar
                IF (PIO_TYPE.eq.PIO_double) THEN
                  pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_p2dvar(ng)
                ELSE
                  pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_p2dvar(ng)
                END IF
                status=nf_fwrite2d(ng, model, pioFile, idLatP,          &
     &                             pioVar, 0, ioDesc,                   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % pmask,                    &
     &                             GRID(ng) % latp,                     &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, PIO_noerr,                       &
     &                         5736, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idLatP)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idLatP)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
        IF (.not.spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
              scale=1.0_dp
              IF (pio_netcdf_find_var(ng, model, pioFile,               &
     &                                TRIM(Vname(1,idXgrP)),            &
     &                                pioVar%vd)) THEN
                pioVar%gtype=p2dvar
                IF (PIO_TYPE.eq.PIO_double) THEN
                  pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_p2dvar(ng)
                ELSE
                  pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_p2dvar(ng)
                END IF
                status=nf_fwrite2d(ng, model, pioFile, idXgrP,          &
     &                             pioVar, 0, ioDesc,                   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % pmask,                    &
     &                             GRID(ng) % xp,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, PIO_noerr,                       &
     &                         5776, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idXgrP)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'x_psi', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
          IF (exit_flag.eq.NoError) THEN
            IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
              scale=1.0_dp
              IF (pio_netcdf_find_var(ng, model, pioFile,               &
     &                                TRIM(Vname(1,idYgrP)),            &
     &                                pioVar%vd)) THEN
                pioVar%gtype=p2dvar
                IF (PIO_TYPE.eq.PIO_double) THEN
                  pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_p2dvar(ng)
                ELSE
                  pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_p2dvar(ng)
                END IF
                status=nf_fwrite2d(ng, model, pioFile, idYgrP,          &
     &                             pioVar, 0, ioDesc,                   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % pmask,                    &
     &                             GRID(ng) % yp,                       &
     &                             SetFillVal = .FALSE.)
                IF (FoundError(status, PIO_noerr,                       &
     &                         5813, MyFile)) THEN
                  IF (Master) WRITE (stdout,10) TRIM(Vname(1,idYgrP)),  &
     &                                          TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) TRIM(Vname(1,idYgrP)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
!  Angle between XI-axis and EAST at RHO-points.
!
        IF (exit_flag.eq.NoError) THEN
          scale=1.0_dp
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            IF (pio_netcdf_find_var(ng, model, pioFile,                 &
     &                              TRIM(Vname(1,idangR)),              &
     &                              pioVar%vd)) THEN
              pioVar%gtype=r2dvar
              IF (PIO_TYPE.eq.PIO_double) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_r2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_r2dvar(ng)
              END IF
              status=nf_fwrite2d(ng, model, pioFile, idangR,            &
     &                           pioVar, 0, ioDesc,                     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % rmask,                      &
     &                           GRID(ng) % angler,                     &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, PIO_noerr, 5855, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idangR)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idangR)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
!  Masking fields at RHO-, U-, V-points, and PSI-points.
!
        IF (exit_flag.eq.NoError) THEN
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            scale=1.0_dp
            IF (pio_netcdf_find_var(ng, model, pioFile,                 &
     &                              TRIM(Vname(1,idmskR)),              &
     &                              pioVar%vd)) THEN
              pioVar%gtype=r2dvar
              IF (PIO_TYPE.eq.PIO_double) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_r2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_r2dvar(ng)
              END IF
              status=nf_fwrite2d(ng, model, pioFile, idmskR,            &
     &                           pioVar, 0, ioDesc,                     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % rmask,                      &
     &                           GRID(ng) % rmask,                      &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, PIO_noerr, 5907, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idmskR)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idmskR)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
        IF (exit_flag.eq.NoError) THEN
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            scale=1.0_dp
            IF (pio_netcdf_find_var(ng, model, pioFile,                 &
     &                              TRIM(Vname(1,idmskU)),              &
     &                              pioVar%vd)) THEN
              pioVar%gtype=u2dvar
              IF (PIO_TYPE.eq.PIO_double) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_u2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_u2dvar(ng)
              END IF
              status=nf_fwrite2d(ng, model, pioFile, idmskU,            &
     &                           pioVar, 0, ioDesc,                     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % umask,                      &
     &                           GRID(ng) % umask,                      &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, PIO_noerr, 5942, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idmskU)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idmskU)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
        IF (exit_flag.eq.NoError) THEN
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            scale=1.0_dp
            IF (pio_netcdf_find_var(ng, model, pioFile,                 &
     &                              TRIM(Vname(1,idmskV)),              &
     &                              pioVar%vd)) THEN
              pioVar%gtype=v2dvar
              IF (PIO_TYPE.eq.PIO_double) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_v2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_v2dvar(ng)
              END IF
              status=nf_fwrite2d(ng, model, pioFile, idmskV,            &
     &                           pioVar, 0, ioDesc,                     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % vmask,                      &
     &                           GRID(ng) % vmask,                      &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, PIO_noerr, 5977, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idmskV)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idmskV)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
        IF (exit_flag.eq.NoError) THEN
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            scale=1.0_dp
            IF (pio_netcdf_find_var(ng, model, pioFile,                 &
     &                              TRIM(Vname(1,idmskP)),              &
     &                              pioVar%vd)) THEN
              pioVar%gtype=p2dvar
              IF (PIO_TYPE.eq.PIO_double) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_p2dvar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_p2dvar(ng)
              END IF
              status=nf_fwrite2d(ng, model, pioFile, idmskP,            &
     &                           pioVar, 0, ioDesc,                     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % pmask,                      &
     &                           GRID(ng) % pmask,                      &
     &                           SetFillVal = .FALSE.)
              IF (FoundError(status, PIO_noerr, 6012, MyFile)) THEN
                IF (Master) WRITE (stdout,10) TRIM(Vname(1,idmskP)),    &
     &                                        TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) TRIM(Vname(1,idmskP)),      &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
      END IF GRID_VARS
!
!-----------------------------------------------------------------------
!  Synchronize NetCDF file to disk to allow other processes to access
!  data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL pio_netcdf_sync (ng, model, ncname, pioFile)
      IF (FoundError(exit_flag, NoError, 6261, MyFile)) RETURN
!
  10  FORMAT (/,' WRT_INFO_PIO - error while writing variable: ',a,/,   &
     &        16x,'into file: ',a)
  20  FORMAT (/,' WRT_INFO_PIO - error while inquiring ID for',         &
     &        ' variable: ',a,/,16x,'in file: ',a)
  30  FORMAT (/,' WRT_INFO_PIO - unable to synchronize to disk file:',  &
     &        /,16x,a)
!
      RETURN
      END SUBROUTINE wrt_info_pio
      END MODULE wrt_info_mod
