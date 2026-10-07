      MODULE def_info_mod
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
!  On input the NetCDF dimesions IDs is an interger vector as follows: !
!                                                                      !
!    DimIDs( 1) => XI-dimension at RHO-points                          !
!    DimIDs( 2) => XI-dimension at U-points                            !
!    DimIDs( 3) => XI-dimension at V-points                            !
!    DimIDs( 4) => XI-dimension at PSI-points                          !
!    DimIDs( 5) => ETA-dimension at RHO-points                         !
!    DimIDs( 6) => ETA-dimension at U-points                           !
!    DimIDs( 7) => ETA-dimension at V-points                           !
!    DimIDs( 8) => ETA-dimension at PSI-points                         !
!    DimIDs( 9) => S-dimension at RHO-points                           !
!    DimIDs(10) => S-dimension at W-points                             !
!    DimIDs(11) => Number of tracers dimension                         !
!    DimIDs(12) => Unlimited time record dimension                     !
!    DimIDs(13) => Number of stations dimension                        !
!    DimIDs(14) => Boundary dimension                                  !
!    DimIDs(15) => Number of floats dimension                          !
!    DimIDs(16) => Number sediment bed layers dimension                !
!    DimIDs(17) => Dimension 2D water RHO-points                       !
!    DimIDs(18) => Dimension 2D water U-points                         !
!    DimIDs(19) => Dimension 2D water V-points                         !
!    DimIDs(20) => Dimension 3D water RHO-points                       !
!    DimIDs(21) => Dimension 3D water U-points                         !
!    DimIDs(23) => Dimension 3D water W-points                         !
!    DimIDs(24) => Dimension sediment bed water points                 !
!    DimIDs(25) => Number of EcoSim phytoplankton groups               !
!    DimIDs(26) => Number of EcoSim bacteria groups                    !
!    DimIDs(27) => Number of EcoSim DOM groups                         !
!    DimIDs(28) => Number of EcoSim fecal groups                       !
!    DimIDs(29) => Number of state variables                           !
!    DimIDs(30) => Number of 3D variables time levels (2)              !
!    DimIDs(31) => Number of 2D variables time levels (3)              !
!    DimIDs(32) => Number of sediment tracers                          !
!    DimIDs(33) => Number of light spectral bands                      !
!    DimIDs(34) => Number constant depth slices in 'extract_slice'     !
!    DimIDs(35) => Number of submerged acuatic vegetation types        !                     
!    DimIDs(36) => Number of correlation functions to be combined      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
      USE mod_strings
!
      USE def_dim_mod, ONLY : def_dim
      USE def_var_mod, ONLY : def_var
      USE lbc_mod,     ONLY : lbc_putatt
      USE strings_mod, ONLY : FoundError, join_string
      USE tadv_mod,    ONLY : tadv_putatt
!
      implicit none
!
      INTERFACE def_info
        MODULE PROCEDURE def_info_nf90
        MODULE PROCEDURE def_info_pio
      END INTERFACE def_info
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE def_info_nf90 (ng, model, ncid, ncname, DimIDs)
!***********************************************************************
!                                                                      !
!  This routine defines information variables for the requested NetCDF !
!  file using the standard NetCDF-3 or NetCDF-4 library.               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng       Nested grid number (integer)                            !
!     model    Calling model identifier (integer)                      !
!     ncid     NetCDF file ID (integer)                                !
!     ncname   NetCDF filename (character)                             !
!     DimIDs   NetCDF dimensions IDs (integer vector of size nDimID)   !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!***********************************************************************
!
      USE mod_netcdf
!
      USE distribute_mod, ONLY : mp_bcasti
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid
      integer, intent(in) :: DimIDs(nDimID)
!
      character (*), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer :: brydim, i, ie, is, j, lstr, varid
      integer :: slicedim, srdim, stadim, status, swdim, trcdim, usrdim
      integer :: ibuffer(2)
      integer :: p2dgrd(2), tbrydim(2)
      integer :: t2dgrd(3), u2dgrd(3), v2dgrd(3)
!
      real(r8) :: Aval(6)
!
      character (len=11 )    :: bryatt, clmatt, frcatt
      character (len=50 )    :: tiling
      character (len=80 )    :: type
      character (len=4096)   :: string
      character (len=MaxLen) :: Vinfo(Natt)
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/def_info.F"//", def_info_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set dimension variables.
!-----------------------------------------------------------------------
!
      p2dgrd(1)=DimIDs(4)
      p2dgrd(2)=DimIDs(8)
      t2dgrd(1)=DimIDs(1)
      t2dgrd(2)=DimIDs(5)
      u2dgrd(1)=DimIDs(2)
      u2dgrd(2)=DimIDs(6)
      v2dgrd(1)=DimIDs(3)
      v2dgrd(2)=DimIDs(7)
      srdim=DimIDs(9)
      swdim=DimIDs(10)
      trcdim=DimIDs(11)
      slicedim=DimIDs(34)
      stadim=DimIDs(13)
      brydim=DimIDs(14)
      tbrydim(1)=DimIDs(11)
      tbrydim(2)=DimIDs(14)
!
!  Set dimension for generic user parameters.
!
      IF ((Nuser.gt.0).and.(ncid.ne.GST(ng)%ncid)) THEN
        status=def_dim(ng, model, ncid, ncname, 'Nuser',                &
     &                 Nuser, usrdim)
        IF (FoundError(exit_flag, NoError, 206, MyFile)) RETURN
      END IF
!
!  Initialize local information variable arrays.
!
      DO i=1,Natt
        DO j=1,LEN(Vinfo(1))
          Vinfo(i)(j:j)=' '
        END DO
      END DO
      DO i=1,6
        Aval(i)=0.0_r8
      END DO
!
!-----------------------------------------------------------------------
!  Define global attributes.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
!
!  Define history global attribute.
!
        IF (LEN_TRIM(date_str).gt.0) THEN
          WRITE (history,'(a,1x,a,", ",a)') 'ROMS, Version',            &
     &                                      TRIM( version),             &
     &                                      TRIM(date_str)
        ELSE
          WRITE (history,'(a,1x,a)') 'ROMS, Version',                   &
     &                               TRIM(version)
        END IF
!
!  Set tile decomposition global attribute.
!
        WRITE (tiling,10) NtileI(ng), NtileJ(ng)
!
!  Define file name global attribute.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'file',                &
     &                        TRIM(ncname))
          IF (FoundError(status, nf90_noerr, 246, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Define NetCDF format type.
!
        IF (exit_flag.eq.NoError) THEN
          IF (CMODE.eq.nf90_netcdf4) THEN
            status=nf90_put_att(ncid, nf90_global, 'format',            &
     &                          'netCDF-4/HDF5 file')
          ELSE IF (CMODE.eq.nf90_64bit_offset) THEN
            status=nf90_put_att(ncid, nf90_global, 'format',            &
     &                          'netCDF-3 64bit offset file')
          END IF
          IF (FoundError(status, nf90_noerr, 265, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'format', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Define file climate and forecast metadata convention global
!  attribute.
!
        type='CF-1.4, SGRID-0.3'
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'Conventions',         &
     &                        TRIM(type))
          IF (FoundError(status, nf90_noerr, 280, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'Conventions', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Define file type global attribute.
!
        IF (ncid.eq.ADM(ng)%ncid) THEN
          type='ROMS adjoint history file'
        ELSE IF (ncid.eq.AVG(ng)%ncid) THEN
          type='ROMS nonlinear model averages file'
        ELSE IF (ncid.eq.DIA(ng)%ncid) THEN
          type='ROMS diagnostics file'
        ELSE IF (ncid.eq.FLT(ng)%ncid) THEN
          type='ROMS floats file'
        ELSE IF (ncid.eq.ERR(ng)%ncid) THEN
          type='ROMS posterior analysis error covariance matrix'
        ELSE IF (ncid.eq.GST(ng)%ncid) THEN
          type='ROMS GST check pointing restart file'
        ELSE IF (ncid.eq.HAR(ng)%ncid) THEN
          type='ROMS Least-squared Detiding Harmonics file'
        ELSE IF (ncid.eq.HSS(ng)%ncid) THEN
          type='ROMS 4D-Var Hessian eigenvectors file'
        ELSE IF (ncid.eq.HIS(ng)%ncid) THEN
          type='ROMS history file'
        ELSE IF (ncid.eq.ITL(ng)%ncid) THEN
          type='ROMS tangent linear model initial file'
        ELSE IF (ncid.eq.LCZ(ng)%ncid) THEN
          type='ROMS 4D-Var Lanczos vectors file'
        ELSE IF (ncid.eq.LZE(ng)%ncid) THEN
          type='ROMS 4D-Var Evolved Lanczos vectors file'
        ELSE IF (ncid.eq.NRM(1,ng)%ncid) THEN
          type='ROMS initial conditions error covariance norm file'
        ELSE IF (ncid.eq.NRM(2,ng)%ncid) THEN
          type='ROMS model error covariance norm file'
        ELSE IF (ncid.eq.NRM(3,ng)%ncid) THEN
         type='ROMS boundary conditions error covariance norm file'
        ELSE IF (ncid.eq.NRM(4,ng)%ncid) THEN
          type='ROMS surface forcing error covariance norm file'
        ELSE IF (ncid.eq.QCK(ng)%ncid) THEN
          type='ROMS quicksave file'
        ELSE IF (ncid.eq.RST(ng)%ncid) THEN
          type='ROMS restart file'
        ELSE IF (ncid.eq.STA(ng)%ncid) THEN
          type='ROMS station file'
        ELSE IF (ncid.eq.TLF(ng)%ncid) THEN
          type='ROMS tangent linear impulse forcing file'
        ELSE IF (ncid.eq.TLM(ng)%ncid) THEN
          type='ROMS tangent linear history file'
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'type',                &
     &                        TRIM(type))
          IF (FoundError(status, nf90_noerr, 353, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'type', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Define other global attributes to NetCDF file.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'title',               &
     &                        TRIM(title))
          IF (FoundError(status, nf90_noerr, 379, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'title', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'var_info',            &
     &                        TRIM(varname))
          IF (FoundError(status, nf90_noerr, 389, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'var_info', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'rst_file',            &
     &                        TRIM(RST(ng)%name))
          IF (FoundError(status, nf90_noerr, 424, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'rst_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          IF (LdefHIS(ng)) THEN
            IF (ndefHIS(ng).gt.0) THEN
              status=nf90_put_att(ncid, nf90_global, 'his_base',        &
     &                            TRIM(HIS(ng)%base))
            ELSE
              status=nf90_put_att(ncid, nf90_global, 'his_file',        &
     &                            TRIM(HIS(ng)%name))
            END IF
            IF (FoundError(status, nf90_noerr, 440, MyFile)) THEN
              IF (Master) WRITE (stdout,20) 'his_file', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          IF (ndefAVG(ng).gt.0) THEN
            status=nf90_put_att(ncid, nf90_global, 'avg_base',          &
     &                          TRIM(AVG(ng)%base))
          ELSE
            status=nf90_put_att(ncid, nf90_global, 'avg_file',          &
     &                          TRIM(AVG(ng)%name))
          END IF
          IF (FoundError(status, nf90_noerr, 489, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'avg_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'grd_file',            &
     &                        TRIM(GRD(ng)%name))
          IF (FoundError(status, nf90_noerr, 568, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'grd_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'ini_file',            &
     &                        TRIM(INI(ng)%name))
          IF (FoundError(status, nf90_noerr, 581, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'ini_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          IF (LuvSrc(ng).or.LwSrc(ng).or.(ANY(LtracerSrc(:,ng)))) THEN
            status=nf90_put_att(ncid, nf90_global, 'river_file',        &
     &                          TRIM(SSF(ng)%name))
            IF (FoundError(status, nf90_noerr, 748, MyFile)) THEN
              IF (Master) WRITE (stdout,20) 'river_file', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          DO i=1,nFfiles(ng)
            CALL join_string (FRC(i,ng)%files, FRC(i,ng)%Nfiles,        &
     &                        string, lstr)
            WRITE (frcatt,30) 'frc_file_', i
            status=nf90_put_att(ncid, nf90_global, frcatt,              &
     &                          string(1:lstr))
            IF (FoundError(status, nf90_noerr, 779, MyFile)) THEN
              IF (Master) WRITE (stdout,20) TRIM(frcatt), TRIM(ncname)
              exit_flag=3
              ioerror=status
              EXIT
            END IF
          END DO
        END IF
        IF (ObcData(ng)) THEN
          DO i=1,nBCfiles(ng)
            IF (exit_flag.eq.NoError) THEN
              CALL join_string (BRY(i,ng)%files, BRY(i,ng)%Nfiles,      &
     &                          string, lstr)
              WRITE (bryatt,30) 'bry_file_', i
              status=nf90_put_att(ncid, nf90_global, bryatt,            &
     &                            string(1:lstr))
              IF (FoundError(status, nf90_noerr, 797, MyFile)) THEN
                IF (Master) WRITE (stdout,20) TRIM(bryatt), TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            END IF
          END DO
        END IF
        IF (Lclimatology(ng)) THEN
          DO i=1,nCLMfiles(ng)
            IF (exit_flag.eq.NoError) THEN
              CALL join_string (CLM(i,ng)%files, CLM(i,ng)%Nfiles,      &
     &                          string, lstr)
              WRITE (clmatt,30) 'clm_file_', i
              status=nf90_put_att(ncid, nf90_global, clmatt,            &
     &                            string(1:lstr))
              IF (FoundError(status, nf90_noerr, 816, MyFile)) THEN
                IF (Master) WRITE (stdout,20) TRIM(clmatt), TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            END IF
          END DO
        END IF
        IF (Lnudging(ng)) THEN
          IF (exit_flag.eq.NoError) THEN
            status=nf90_put_att(ncid, nf90_global, 'nud_file',          &
     &                        TRIM(NUD(ng)%name))
            IF (FoundError(status, nf90_noerr, 831, MyFile)) THEN
              IF (Master) WRITE (stdout,20) 'nud_file', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'script_file',         &
     &                        TRIM(Iname))
          IF (FoundError(status, nf90_noerr, 870, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'script_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  NLM tracer advection scheme.
!
        IF (exit_flag.eq.NoError) THEN
          CALL tadv_putatt (ng, ncid, ncname, 'NLM_TADV',               &
     &                      Hadvection, Vadvection, status)
          IF (FoundError(status, nf90_noerr, 934, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'NLM_TADV', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  NLM Lateral boundary conditions.
!
        IF (exit_flag.eq.NoError) THEN
          CALL lbc_putatt (ng, ncid, ncname, 'NLM_LBC', LBC, status)
          IF (FoundError(status, nf90_noerr, 958, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'NLM_LBC', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  GIT repository information.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'git_url',             &
     &                        TRIM(git_url))
          IF (FoundError(status, nf90_noerr, 986, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'git_url', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'git_rev',             &
     &                        TRIM(git_rev))
          IF (FoundError(status, nf90_noerr, 997, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'git_rev', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  SVN repository information.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'svn_url',             &
     &                        TRIM(svn_url))
          IF (FoundError(status, nf90_noerr, 1010, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'svn_url', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Local root directory, cpp header directory and file, and analytical
!  directory
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'code_dir',            &
     &                        TRIM(Rdir))
          IF (FoundError(status, nf90_noerr, 1037, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'code_dir', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'header_dir',          &
     &                        TRIM(Hdir))
          IF (FoundError(status, nf90_noerr, 1049, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'header_dir', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'header_file',         &
     &                        TRIM(Hfile))
          IF (FoundError(status, nf90_noerr, 1061, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'header_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Attributes describing platform and compiler
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'os',                  &
     &                        TRIM(my_os))
          IF (FoundError(status, nf90_noerr, 1076, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'os', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'cpu',                 &
     &                        TRIM(my_cpu))
          IF (FoundError(status, nf90_noerr, 1086, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'cpu', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'compiler_system',     &
     &                        TRIM(my_fort))
          IF (FoundError(status, nf90_noerr, 1096, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'compiler_system',            &
     &                                    TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'compiler_command',    &
     &                        TRIM(my_fc))
          IF (FoundError(status, nf90_noerr, 1107, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'compiler_command',           &
     &                                    TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          lstr=INDEX(my_fflags, 'free')-2
          IF (lstr.le.0) lstr=LEN_TRIM(my_fflags)
          status=nf90_put_att(ncid, nf90_global, 'compiler_flags',      &
     &                        my_fflags(1:lstr))
          IF (FoundError(status, nf90_noerr, 1120, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'compiler_flags', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Tiling and history attributes.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'tiling',              &
     &                        TRIM(tiling))
          IF (FoundError(status, nf90_noerr, 1132, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'tiling', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'history',             &
     &                        TRIM(history))
          IF (FoundError(status, nf90_noerr, 1142, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'history', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Analytical header files used.
!
        IF (exit_flag.eq.NoError) THEN
          CALL join_string (ANANAME, SIZE(ANANAME), string, lstr)
          IF (lstr.gt.0) THEN
            status=nf90_put_att(ncid, nf90_global, 'ana_file',          &
     &                          string(1:lstr))
            IF (FoundError(status, nf90_noerr, 1156, MyFile)) THEN
              IF (Master) WRITE (stdout,20) 'ana_file', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
!
!  Activated CPP options.
!
        IF (exit_flag.eq.NoError) THEN
          lstr=LEN_TRIM(Coptions)-1
          status=nf90_put_att(ncid, nf90_global, 'CPP_options',         &
     &                        TRIM(Coptions(1:lstr)))
          IF (FoundError(status, nf90_noerr, 1191, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'CPP_options', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (FoundError(exit_flag, NoError, 1207, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Define running parameters.
!-----------------------------------------------------------------------
!
!  Time stepping parameters.
!
      Vinfo( 1)='ntimes'
      Vinfo( 2)='number of long time-steps'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1228, MyFile)) RETURN
      Vinfo( 1)='ndtfast'
      Vinfo( 2)='number of short time-steps'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1235, MyFile)) RETURN
      Vinfo( 1)='dt'
      Vinfo( 2)='size of long time-steps'
      Vinfo( 3)='second'
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1243, MyFile)) RETURN
      Vinfo( 1)='dtfast'
      Vinfo( 2)='size of short time-steps'
      Vinfo( 3)='second'
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1251, MyFile)) RETURN
      Vinfo( 1)='dstart'
      Vinfo( 2)='time stamp assigned to model initilization'
      WRITE (Vinfo( 3),'(a,a)') 'days since ', TRIM(Rclock%string)
      Vinfo( 4)=TRIM(Rclock%calendar)
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1260, MyFile)) RETURN
      Vinfo( 1)='nHIS'
      Vinfo( 2)='number of time-steps between history records'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1308, MyFile)) RETURN
      Vinfo( 1)='ndefHIS'
      Vinfo( 2)=                                                        &
     &    'number of time-steps between the creation of history files'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1316, MyFile)) RETURN
      Vinfo( 1)='nRST'
      Vinfo( 2)='number of time-steps between restart records'
      IF (LcycleRST(ng)) THEN
        Vinfo(13)='only latest two records are maintained'
      END IF
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1350, MyFile)) RETURN
      Vinfo( 1)='ntsAVG'
      Vinfo( 2)=                                                        &
     &   'starting time-step for accumulation of time-averaged fields'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1362, MyFile)) RETURN
      Vinfo( 1)='nAVG'
      Vinfo( 2)='number of time-steps between time-averaged records'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1369, MyFile)) RETURN
      Vinfo( 1)='ndefAVG'
      Vinfo( 2)=                                                        &
     &    'number of time-steps between the creation of average files'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1377, MyFile)) RETURN
!
!  Power-law shape filter parameters for time-averaging of barotropic
!  fields.
!
      Vinfo( 1)='Falpha'
      Vinfo( 2)='Power-law shape barotropic filter parameter'
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1549, MyFile)) RETURN
      Vinfo( 1)='Fbeta'
      Vinfo( 2)='Power-law shape barotropic filter parameter'
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1556, MyFile)) RETURN
      Vinfo( 1)='Fgamma'
      Vinfo( 2)='Power-law shape barotropic filter parameter'
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1563, MyFile)) RETURN
!
!  Horizontal mixing coefficients.
!
      Vinfo( 1)='nl_tnu2'
      Vinfo( 2)='nonlinear model Laplacian mixing coefficient '//       &
     &          'for tracers'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1576, MyFile)) RETURN
      Vinfo( 1)='nl_visc2'
      Vinfo( 2)='nonlinear model Laplacian mixing coefficient '//       &
     &          'for momentum'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1642, MyFile)) RETURN
      Vinfo( 1)='LuvSponge'
      Vinfo( 2)='horizontal viscosity sponge activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1728, MyFile)) RETURN
      Vinfo( 1)='LtracerSponge'
      Vinfo( 2)='horizontal diffusivity sponge activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1738, MyFile)) RETURN
!
!  Background vertical mixing coefficients.
!
      Vinfo( 1)='Akt_bak'
      Vinfo( 2)='background vertical mixing coefficient for tracers'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1750, MyFile)) RETURN
      Vinfo( 1)='Akv_bak'
      Vinfo( 2)='background vertical mixing coefficient for momentum'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1758, MyFile)) RETURN
      Vinfo( 1)='Akk_bak'
      Vinfo( 2)=                                                        &
     &   'background vertical mixing coefficient for turbulent energy'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1768, MyFile)) RETURN
      Vinfo( 1)='Akp_bak'
      Vinfo( 2)=                                                        &
     &   'background vertical mixing coefficient for length scale'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1777, MyFile)) RETURN
!
!  Drag coefficients.
!
      Vinfo( 1)='rdrg'
      Vinfo( 2)='linear drag coefficient'
      Vinfo( 3)='meter second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1834, MyFile)) RETURN
      Vinfo( 1)='rdrg2'
      Vinfo( 2)='quadratic drag coefficient'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo ,ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1841, MyFile)) RETURN
      Vinfo( 1)='Zob'
      Vinfo( 2)='bottom roughness'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1850, MyFile)) RETURN
      Vinfo( 1)='Zos'
      Vinfo( 2)='surface roughness'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1858, MyFile)) RETURN
!
!  Generic length-scale parameters.
!
      Vinfo( 1)='gls_p'
      Vinfo( 2)='stability exponent'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1869, MyFile)) RETURN
      Vinfo( 1)='gls_m'
      Vinfo( 2)='turbulent kinetic energy exponent'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1876, MyFile)) RETURN
      Vinfo( 1)='gls_n'
      Vinfo( 2)='turbulent length scale exponent'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1883, MyFile)) RETURN
      Vinfo( 1)='gls_cmu0'
      Vinfo( 2)='stability coefficient'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1890, MyFile)) RETURN
      Vinfo( 1)='gls_c1'
      Vinfo( 2)='shear production coefficient'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1897, MyFile)) RETURN
      Vinfo( 1)='gls_c2'
      Vinfo( 2)='dissipation coefficient'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1904, MyFile)) RETURN
      Vinfo( 1)='gls_c3m'
      Vinfo( 2)='buoyancy production coefficient (minus)'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1911, MyFile)) RETURN
      Vinfo( 1)='gls_c3p'
      Vinfo( 2)='buoyancy production coefficient (plus)'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1918, MyFile)) RETURN
      Vinfo( 1)='gls_sigk'
      Vinfo( 2)='constant Schmidt number for TKE'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1925, MyFile)) RETURN
      Vinfo( 1)='gls_sigp'
      Vinfo( 2)='constant Schmidt number for PSI'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1932, MyFile)) RETURN
      Vinfo( 1)='gls_Kmin'
      Vinfo( 2)='minimum value of specific turbulent kinetic energy'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1939, MyFile)) RETURN
      Vinfo( 1)='gls_Pmin'
      Vinfo( 2)='minimum Value of dissipation'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1946, MyFile)) RETURN
      Vinfo( 1)='Charnok_alpha'
      Vinfo( 2)='Charnock factor for surface roughness'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1953, MyFile)) RETURN
      Vinfo( 1)='Zos_hsig_alpha'
      Vinfo( 2)='wave amplitude factor for surface roughness'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1960, MyFile)) RETURN
      Vinfo( 1)='sz_alpha'
      Vinfo( 2)='surface flux from wave dissipation'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1967, MyFile)) RETURN
      Vinfo( 1)='CrgBan_cw'
      Vinfo( 2)='surface flux due to Craig and Banner wave breaking'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1974, MyFile)) RETURN
!
!  Nudging inverse time scales used in various tasks.
!
      Vinfo( 1)='Znudg'
      Vinfo( 2)='free-surface nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 1995, MyFile)) RETURN
      Vinfo( 1)='M2nudg'
      Vinfo( 2)='2D momentum nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2003, MyFile)) RETURN
      Vinfo( 1)='M3nudg'
      Vinfo( 2)='3D momentum nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2012, MyFile)) RETURN
      Vinfo( 1)='Tnudg'
      Vinfo( 2)='Tracers nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2020, MyFile)) RETURN
!
!  Open boundary nudging, inverse time scales.
!
      IF (NudgingCoeff(ng)) THEN
        Vinfo( 1)='FSobc_in'
        Vinfo( 2)='free-surface inflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TOUT,                 &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 2034, MyFile)) RETURN
        Vinfo( 1)='FSobc_out'
        Vinfo( 2)='free-surface outflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TOUT,                 &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 2042, MyFile)) RETURN
        Vinfo( 1)='M2obc_in'
        Vinfo( 2)='2D momentum inflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TOUT,                 &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 2050, MyFile)) RETURN
        Vinfo( 1)='M2obc_out'
        Vinfo( 2)='2D momentum outflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TOUT,                 &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 2058, MyFile)) RETURN
        Vinfo( 1)='Tobc_in'
        Vinfo( 2)='tracers inflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TOUT,                 &
     &                 2, tbrydim, Aval, Vinfo, ncname,                 &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 2067, MyFile)) RETURN
        Vinfo( 1)='Tobc_out'
        Vinfo( 2)='tracers outflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TOUT,                 &
     &                 2, tbrydim, Aval, Vinfo, ncname,                 &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 2075, MyFile)) RETURN
        Vinfo( 1)='M3obc_in'
        Vinfo( 2)='3D momentum inflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TOUT,                 &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 2083, MyFile)) RETURN
        Vinfo( 1)='M3obc_out'
        Vinfo( 2)='3D momentum outflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TOUT,                 &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 2091, MyFile)) RETURN
      END IF
!
!  Equation of State parameters.
!
      Vinfo( 1)='rho0'
      Vinfo( 2)='mean density used in Boussinesq approximation'
      Vinfo( 3)='kilogram meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2104, MyFile)) RETURN
!
!  Various parameters.
!
!
!  Slipperiness parameters.
!
      Vinfo( 1)='gamma2'
      Vinfo( 2)='slipperiness parameter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2169, MyFile)) RETURN
!
! Logical switches to activate horizontal momentum transport
! point Sources/Sinks (like river runoff transport) and mass point
! Sources/Sinks (like volume vertical influx).
!
      Vinfo( 1)='LuvSrc'
      Vinfo( 2)='momentum point sources and sink activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2182, MyFile)) RETURN
      Vinfo( 1)='LwSrc'
      Vinfo( 2)='mass point sources and sink activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2191, MyFile)) RETURN
!
!  Logical switches indicating which tracer variables are processed
!  during point Sources/Sinks.
!
      Vinfo( 1)='LtracerSrc'
      Vinfo( 2)='tracer point sources and sink activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2205, MyFile)) RETURN
!
!  Logical switches to process climatology fields.
!
      Vinfo( 1)='LsshCLM'
      Vinfo( 2)='sea surface height climatology processing switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2217, MyFile)) RETURN
      Vinfo( 1)='Lm2CLM'
      Vinfo( 2)='2D momentum climatology processing switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2226, MyFile)) RETURN
      Vinfo( 1)='Lm3CLM'
      Vinfo( 2)='3D momentum climatology processing switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2236, MyFile)) RETURN
      Vinfo( 1)='LtracerCLM'
      Vinfo( 2)='tracer climatology processing switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2245, MyFile)) RETURN
!
!  Logical switches for nudging of climatology fields.
!
      Vinfo( 1)='LnudgeM2CLM'
      Vinfo( 2)='2D momentum climatology nudging activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2257, MyFile)) RETURN
!
      Vinfo( 1)='LnudgeM3CLM'
      Vinfo( 2)='3D momentum climatology nudging activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2267, MyFile)) RETURN
!
      Vinfo( 1)='LnudgeTCLM'
      Vinfo( 2)='tracer climatology nudging activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 2276, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Define grid variables.
!-----------------------------------------------------------------------
!
!  Grid type switch: Spherical or Cartesian. Writing characters in
!  parallel I/O is extremely inefficient.  It is better to write
!  this as an integer switch: 0=Cartesian, 1=spherical.
!
      Vinfo( 1)='spherical'
      Vinfo( 2)='grid type logical switch'
      Vinfo( 9)='Cartesian'
      Vinfo(10)='spherical'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 3073, MyFile)) RETURN
!
!  Domain Length.
!
      Vinfo( 1)='xl'
      Vinfo( 2)='domain length in the XI-direction'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 3083, MyFile)) RETURN
      Vinfo( 1)='el'
      Vinfo( 2)='domain length in the ETA-direction'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 3091, MyFile)) RETURN
!
!  S-coordinate parameters.
!
      Vinfo( 1)='Vtransform'
      Vinfo( 2)='vertical terrain-following transformation equation'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 3101, MyFile)) RETURN
      Vinfo( 1)='Vstretching'
      Vinfo( 2)='vertical terrain-following stretching function'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 3108, MyFile)) RETURN
      Vinfo( 1)='theta_s'
      Vinfo( 2)='S-coordinate surface control parameter'
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 3115, MyFile)) RETURN
      Vinfo( 1)='theta_b'
      Vinfo( 2)='S-coordinate bottom control parameter'
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 3122, MyFile)) RETURN
      Vinfo( 1)='Tcline'
      Vinfo( 2)='S-coordinate surface/bottom layer width'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 3130, MyFile)) RETURN
      Vinfo( 1)='hc'
      Vinfo( 2)='S-coordinate parameter, critical depth'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 3138, MyFile)) RETURN
!
!  SGRID conventions for staggered data on structured grids.
!
      Vinfo( 1)='grid'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 3146, MyFile)) RETURN
!
!  S-coordinate non-dimensional independent variable at RHO-points.
!
      Vinfo( 1)='s_rho'
      Vinfo( 2)='S-coordinate at RHO-points'
      Vinfo( 5)='valid_min'
      Vinfo( 6)='valid_max'
      IF (Vtransform(ng).eq.1) THEN
        Vinfo(21)='ocean_s_coordinate_g1'
      ELSE IF (Vtransform(ng).eq.2) THEN
        Vinfo(21)='ocean_s_coordinate_g2'
      END IF
      Vinfo(23)='s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc'
      vinfo(25)='up'
      Aval(2)=-1.0_r8
      Aval(3)=0.0_r8
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/srdim/), Aval, Vinfo, ncname,                 &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 3170, MyFile)) RETURN
!
!  S-coordinate non-dimensional independent variable at W-points.
!
      Vinfo( 1)='s_w'
      Vinfo( 2)='S-coordinate at W-points'
      Vinfo( 5)='valid_min'
      Vinfo( 6)='valid_max'
      Vinfo(21)='ocean_s_coordinate'
      IF (Vtransform(ng).eq.1) THEN
        Vinfo(21)='ocean_s_coordinate_g1'
      ELSE IF (Vtransform(ng).eq.2) THEN
        Vinfo(21)='ocean_s_coordinate_g2'
      END IF
      Vinfo(23)='s: s_w C: Cs_w eta: zeta depth: h depth_c: hc'
      vinfo(25)='up'
      Aval(2)=-1.0_r8
      Aval(3)=0.0_r8
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/swdim/), Aval, Vinfo, ncname,                 &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 3195, MyFile)) RETURN
!
!  S-coordinate non-dimensional stretching curves at RHO-points.
!
      Vinfo( 1)='Cs_r'
      Vinfo( 2)='S-coordinate stretching curves at RHO-points'
      Vinfo( 5)='valid_min'
      Vinfo( 6)='valid_max'
      Aval(2)=-1.0_r8
      Aval(3)=0.0_r8
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/srdim/), Aval, Vinfo, ncname,                 &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 3208, MyFile)) RETURN
!
!  S-coordinate non-dimensional stretching curves at W-points.
!
      Vinfo( 1)='Cs_w'
      Vinfo( 2)='S-coordinate stretching curves at W-points'
      Vinfo( 5)='valid_min'
      Vinfo( 6)='valid_max'
      Aval(2)=-1.0_r8
      Aval(3)=0.0_r8
      status=def_var(ng, model, ncid, varid, NF_TOUT,                   &
     &               1, (/swdim/), Aval, Vinfo, ncname,                 &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 3221, MyFile)) RETURN
!
!  Depth of horizontal slices.
!
      IF (Nslice.gt.0) THEN
        Vinfo( 1)='z_slice'
        Vinfo( 2)='constant depth of output fields horizontal slices'
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/slicedim/), Aval, Vinfo, ncname,            &
     &               SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 3234, MyFile)) RETURN
      END IF
!
!  User generic parameters.
!
      IF (Nuser.gt.0) THEN
        Vinfo( 1)='user'
        Vinfo( 2)='user generic parameters'
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/usrdim/), Aval, Vinfo, ncname,              &
     &               SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 3247, MyFile)) RETURN
      END IF
      IF (ncid.ne.FLT(ng)%ncid) THEN
!
!  Bathymetry.
!
        Vinfo( 1)=Vname(1,idtopo)
        Vinfo( 2)=Vname(2,idtopo)
        Vinfo( 3)=Vname(3,idtopo)
        Vinfo(14)=Vname(4,idtopo)
        Vinfo(21)=Vname(6,idtopo)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idtopo,ng),r8)
        IF (ncid.eq.STA(ng)%ncid) THEN
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   1, (/stadim/), Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3288, MyFile)) RETURN
        ELSE
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3292, MyFile)) RETURN
        END IF
!
!  Coriolis Parameter.
!
        IF (ncid.ne.STA(ng)%ncid) THEN
          Vinfo( 1)=Vname(1,idfcor)
          Vinfo( 2)=Vname(2,idfcor)
          Vinfo( 3)=Vname(3,idfcor)
          Vinfo(14)=Vname(4,idfcor)
          Vinfo(21)=Vname(6,idfcor)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idfcor,ng),r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3308, MyFile)) RETURN
        END IF
!
!  Curvilinear coordinate metrics.
!
        IF (ncid.ne.STA(ng)%ncid) THEN
          Vinfo( 1)=Vname(1,idpmdx)
          Vinfo( 2)=Vname(2,idpmdx)
          Vinfo( 3)=Vname(3,idpmdx)
          Vinfo(14)=Vname(4,idpmdx)
          Vinfo(21)=Vname(6,idpmdx)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idpmdx,ng),r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3323, MyFile)) RETURN
!
          Vinfo( 1)=Vname(1,idpndy)
          Vinfo( 2)=Vname(2,idpndy)
          Vinfo( 3)=Vname(3,idpndy)
          Vinfo(14)=Vname(4,idpndy)
          Vinfo(21)=Vname(6,idpndy)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idpndy,ng),r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3334, MyFile)) RETURN
        END IF
!
!  Grid coordinates of RHO-points.
!
        IF (spherical) THEN
          Vinfo( 1)=Vname(1,idLonR)
          Vinfo( 2)=Vname(2,idLonR)
          Vinfo( 3)=Vname(3,idLonR)
          Vinfo(14)=Vname(4,idLonR)
          Vinfo(21)=Vname(6,idLonR)
          IF (ncid.eq.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     1, (/stadim/), Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3348, MyFile)) RETURN
          ELSE
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, t2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3352, MyFile)) RETURN
          END IF
!
          Vinfo( 1)=Vname(1,idLatR)
          Vinfo( 2)=Vname(2,idLatR)
          Vinfo( 3)=Vname(3,idLatR)
          Vinfo(14)=Vname(4,idLatR)
          Vinfo(21)=Vname(6,idLatR)
          IF (ncid.eq.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     1, (/stadim/), Aval, Vinfo,  ncname)
            IF (FoundError(exit_flag, NoError, 3363, MyFile)) RETURN
          ELSE
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, t2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3367, MyFile)) RETURN
          END IF
        ELSE
          Vinfo( 1)=Vname(1,idXgrR)
          Vinfo( 2)=Vname(2,idXgrR)
          Vinfo( 3)=Vname(3,idXgrR)
          Vinfo(14)=Vname(4,idXgrR)
          Vinfo(21)=Vname(6,idXgrR)
          IF (ncid.eq.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     1, (/stadim/), Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3378, MyFile)) RETURN
          ELSE
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, t2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3382, MyFile)) RETURN
          END IF
!
          Vinfo( 1)=Vname(1,idYgrR)
          Vinfo( 2)=Vname(2,idYgrR)
          Vinfo( 3)=Vname(3,idYgrR)
          Vinfo(14)=Vname(4,idYgrR)
          Vinfo(21)=Vname(6,idYgrR)
          IF (ncid.eq.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     1, (/stadim/), Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3393, MyFile)) RETURN
          ELSE
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, t2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3397, MyFile)) RETURN
          END IF
        END IF
!
!  Grid coordinates of U-points.
!
        IF (spherical) THEN
          Vinfo( 1)=Vname(1,idLonU)
          Vinfo( 2)=Vname(2,idLonU)
          Vinfo( 3)=Vname(3,idLonU)
          Vinfo(14)=Vname(4,idLonU)
          Vinfo(21)=Vname(6,idLonU)
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, u2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3412, MyFile)) RETURN
          END IF
!
          Vinfo( 1)=Vname(1,idLatU)
          Vinfo( 2)=Vname(2,idLatU)
          Vinfo( 3)=Vname(3,idLatU)
          Vinfo(14)=Vname(4,idLatU)
          Vinfo(21)=Vname(6,idLatU)
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, u2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3423, MyFile)) RETURN
          END IF
        ELSE
          Vinfo( 1)=Vname(1,idXgrU)
          Vinfo( 2)=Vname(2,idXgrU)
          Vinfo( 3)=Vname(3,idXgrU)
          Vinfo(14)=Vname(4,idXgrU)
          Vinfo(21)=Vname(6,idXgrU)
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, u2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3434, MyFile)) RETURN
          END IF
!
          Vinfo( 1)=Vname(1,idYgrU)
          Vinfo( 2)=Vname(2,idYgrU)
          Vinfo( 3)=Vname(3,idYgrU)
          Vinfo(14)=Vname(4,idYgrU)
          Vinfo(21)=Vname(6,idYgrU)
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, u2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3445, MyFile)) RETURN
          END IF
        END IF
!
!  Grid coordinates of V-points.
!
        IF (spherical) THEN
          Vinfo( 1)=Vname(1,idLonV)
          Vinfo( 2)=Vname(2,idLonV)
          Vinfo( 3)=Vname(3,idLonV)
          Vinfo(14)=Vname(4,idLonV)
          Vinfo(21)=Vname(6,idLonV)
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, v2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3460, MyFile)) RETURN
          END IF
!
          Vinfo( 1)=Vname(1,idLatV)
          Vinfo( 2)=Vname(2,idLatV)
          Vinfo( 3)=Vname(3,idLatV)
          Vinfo(14)=Vname(4,idLatV)
          Vinfo(21)=Vname(6,idLatV)
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, v2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3471, MyFile)) RETURN
          END IF
        ELSE
          Vinfo( 1)=Vname(1,idXgrV)
          Vinfo( 2)=Vname(2,idXgrV)
          Vinfo( 3)=Vname(3,idXgrV)
          Vinfo(14)=Vname(4,idXgrV)
          Vinfo(21)=Vname(6,idXgrV)
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, v2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3482, MyFile)) RETURN
          END IF
!
          Vinfo( 1)=Vname(1,idYgrV)
          Vinfo( 2)=Vname(2,idYgrV)
          Vinfo( 3)=Vname(3,idYgrV)
          Vinfo(14)=Vname(4,idYgrV)
          Vinfo(21)=Vname(6,idYgrV)
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, v2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3493, MyFile)) RETURN
          END IF
        END IF
!
!  Grid coordinates of PSI-points.
!
        IF (spherical) THEN
          Vinfo( 1)=Vname(1,idLonP)
          Vinfo( 2)=Vname(2,idLonP)
          Vinfo( 3)=Vname(3,idLonP)
          Vinfo(14)=Vname(4,idLonP)
          Vinfo(21)=Vname(6,idLonP)
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, p2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3508, MyFile)) RETURN
          END IF
!
          Vinfo( 1)=Vname(1,idLatP)
          Vinfo( 2)=Vname(2,idLatP)
          Vinfo( 3)=Vname(3,idLatP)
          Vinfo(14)=Vname(4,idLatP)
          Vinfo(21)=Vname(6,idLatP)
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, p2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3519, MyFile)) RETURN
          END IF
        ELSE
          Vinfo( 1)=Vname(1,idXgrP)
          Vinfo( 2)=Vname(2,idXgrP)
          Vinfo( 3)=Vname(3,idXgrP)
          Vinfo(14)=Vname(4,idXgrP)
          Vinfo(21)=Vname(6,idXgrP)
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, p2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3530, MyFile)) RETURN
          END IF
!
          Vinfo( 1)=Vname(1,idYgrP)
          Vinfo( 2)=Vname(2,idYgrP)
          Vinfo( 3)=Vname(3,idYgrP)
          Vinfo(14)=Vname(4,idYgrP)
          Vinfo(21)=Vname(6,idYgrP)
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, p2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3541, MyFile)) RETURN
          END IF
        END IF
!
!  Angle between XI-axis and EAST at RHO-points.
!
        Vinfo( 1)=Vname(1,idangR)
        Vinfo( 2)=Vname(2,idangR)
        Vinfo( 3)=Vname(3,idangR)
        Vinfo(14)=Vname(4,idangR)
        Vinfo(21)=Vname(6,idangR)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idangR,ng),r8)
        IF (ncid.eq.STA(ng)%ncid) THEN
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   1, (/stadim/), Aval, Vinfo,  ncname)
          IF (FoundError(exit_flag, NoError, 3559, MyFile)) RETURN
        ELSE
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3563, MyFile)) RETURN
        END IF
!
!  Masking fields at RHO-, U-, V-points, and PSI-points.
!
        IF (ncid.ne.STA(ng)%ncid) THEN
          Vinfo( 1)=Vname(1,idmskR)
          Vinfo( 2)=Vname(2,idmskR)
          Vinfo( 9)='land'
          Vinfo(10)='water'
          Vinfo(21)=Vname(6,idmskR)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idmskR,ng),r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3580, MyFile)) RETURN
!
          Vinfo( 1)=Vname(1,idmskU)
          Vinfo( 2)=Vname(2,idmskU)
          Vinfo( 9)='land'
          Vinfo(10)='water'
          Vinfo(21)=Vname(6,idmskU)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idmskU,ng),r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3591, MyFile)) RETURN
!
          Vinfo( 1)=Vname(1,idmskV)
          Vinfo( 2)=Vname(2,idmskV)
          Vinfo( 9)='land'
          Vinfo(10)='water'
          Vinfo(21)=Vname(6,idmskV)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idmskV,ng),r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3602, MyFile)) RETURN
!
          Vinfo( 1)=Vname(1,idmskP)
          Vinfo( 2)=Vname(2,idmskP)
          Vinfo( 9)='land'
          Vinfo(10)='water'
          Vinfo(21)=Vname(6,idmskP)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idmskP,ng),r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, p2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3613, MyFile)) RETURN
        END IF
      END IF
!
  10  FORMAT (i0,'x',i0)
  20  FORMAT (/,' DEF_INFO_NF90 - error while creating global',         &
     &        ' attribute: ',a,/,17x,a)
  30  FORMAT (a,i2.2)
!
      RETURN
      END SUBROUTINE def_info_nf90
!
!***********************************************************************
      SUBROUTINE def_info_pio (ng, model, pioFile, ncname, DimIDs)
!***********************************************************************
!                                                                      !
!  This routine defines information variables for the requested NetCDF !
!  file using the NCAR Parallel-IO library.                            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng       Nested grid number (integer)                            !
!     model    Calling model identifier (integer)                      !
!     pioFile  PIO file descriptor structure, TYPE(file_desc_t)        !
!                pioFile%fh         file handler                       !
!                pioFile%iosystem   IO system descriptor (struct)      !
!     ncname   PIO filename (character)                                !
!     DimIDs   PIO dimensions IDs (integer vector of size nDimID)      !
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
      integer, intent(in) :: DimIDs(nDimID)
!
      character (*), intent(in) :: ncname
!
      TYPE (file_desc_t), intent(in) :: pioFile
!
!  Local variable declarations.
!
      integer :: brydim, FileH, i, ie, is, j, lstr
      integer :: slicedim, srdim, stadim, status, swdim, trcdim, usrdim
      integer :: ibuffer(2)
      integer :: p2dgrd(2), tbrydim(2)
      integer :: t2dgrd(3), u2dgrd(3), v2dgrd(3)
!
      real(r8) :: Aval(6)
!
      character (len=11 )    :: bryatt, clmatt, frcatt
      character (len=50 )    :: tiling
      character (len=80 )    :: type
      character (len=4096)   :: string
      character (len=MaxLen) :: Vinfo(Natt)
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/def_info.F"//", def_info_pio"
!
      TYPE (Var_desc_t) :: pioVar
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set dimension variables.
!-----------------------------------------------------------------------
!
      p2dgrd(1)=DimIDs(4)
      p2dgrd(2)=DimIDs(8)
      t2dgrd(1)=DimIDs(1)
      t2dgrd(2)=DimIDs(5)
      u2dgrd(1)=DimIDs(2)
      u2dgrd(2)=DimIDs(6)
      v2dgrd(1)=DimIDs(3)
      v2dgrd(2)=DimIDs(7)
      srdim=DimIDs(9)
      swdim=DimIDs(10)
      trcdim=DimIDs(11)
      slicedim=DimIDs(34)
      stadim=DimIDs(13)
      brydim=DimIDs(14)
      tbrydim(1)=DimIDs(11)
      tbrydim(2)=DimIDs(14)
!
!  Get NetCDF file handle from descriptor.
!
     FileH=ABS(pioFile%fh)
!
!  Set dimension for generic user parameters.
!
      IF ((Nuser.gt.0).and.(FileH.ne.ABS(GST(ng)%pioFile%fh))) THEN
        status=def_dim(ng, model, pioFile, ncname, 'Nuser',             &
     &                 Nuser, usrdim)
        IF (FoundError(exit_flag, NoError, 3840, MyFile)) RETURN
      END IF
!
!  Initialize local information variable arrays.
!
      DO i=1,Natt
        DO j=1,LEN(Vinfo(1))
          Vinfo(i)(j:j)=' '
        END DO
      END DO
      DO i=1,6
        Aval(i)=0.0_r8
      END DO
!
!-----------------------------------------------------------------------
!  Define global attributes.
!-----------------------------------------------------------------------
!
!  Define history global attribute.
!
      IF (LEN_TRIM(date_str).gt.0) THEN
        WRITE (history,'(a,1x,a,", ",a)') 'ROMS, Version',              &
     &                                    TRIM( version),               &
     &                                    TRIM(date_str)
      ELSE
        WRITE (history,'(a,1x,a)') 'ROMS, Version',                     &
     &                             TRIM(version)
      END IF
!
!  Set tile decomposition global attribute.
!
      WRITE (tiling,10) NtileI(ng), NtileJ(ng)
!
!  Define file name global attribute.
!
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_att(pioFile, PIO_global, 'file',                 &
     &                     TRIM(ncname))
        IF (FoundError(status, PIO_noerr, 3878, MyFile)) THEN
          IF (Master) WRITE (stdout,20) 'file', TRIM(ncname)
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
!  Define NetCDF format type.
!
      IF (exit_flag.eq.NoError) THEN
        SELECT CASE (pio_method)
          CASE (1, 2)
            type='netCDF-3 64bit offset file'
          CASE (3, 4)
            type='netCDF-4/HDF5 file'
        END SELECT
        status=PIO_put_att(pioFile, PIO_global, 'format',               &
     &                     TRIM(type))
        IF (FoundError(status, PIO_noerr, 3898, MyFile)) THEN
          IF (Master) WRITE (stdout,20) 'format', TRIM(ncname)
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
!  Define file climate and forecast metadata convention global
!  attribute.
!
      type='CF-1.4, SGRID-0.3'
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_att(pioFile, PIO_global, 'Conventions',          &
     &                     TRIM(type))
        IF (FoundError(status, PIO_noerr, 3913, MyFile)) THEN
          IF (Master) WRITE (stdout,20) 'Conventions', TRIM(ncname)
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
!  Define file type global attribute.
!
      IF (FileH.eq.ABS(ADM(ng)%pioFile%fh)) THEN
        type='ROMS adjoint history file'
      ELSE IF (FileH.eq.ABS(AVG(ng)%pioFile%fh)) THEN
        type='ROMS nonlinear model averages file'
      ELSE IF (FileH.eq.ABS(DIA(ng)%pioFile%fh)) THEN
        type='ROMS diagnostics file'
      ELSE IF (FileH.eq.ABS(FLT(ng)%pioFile%fh)) THEN
        type='ROMS floats file'
      ELSE IF (FileH.eq.ABS(ERR(ng)%pioFile%fh)) THEN
        type='ROMS posterior analysis error covariance matrix'
      ELSE IF (FileH.eq.ABS(GST(ng)%pioFile%fh)) THEN
        type='ROMS GST check pointing restart file'
      ELSE IF (FileH.eq.ABS(HAR(ng)%pioFile%fh)) THEN
        type='ROMS Least-squared Detiding Harmonics file'
      ELSE IF (FileH.eq.ABS(HSS(ng)%pioFile%fh)) THEN
        type='ROMS 4D-Var Hessian eigenvectors file'
      ELSE IF (FileH.eq.ABS(HIS(ng)%pioFile%fh)) THEN
        type='ROMS history file'
      ELSE IF (FileH.eq.ABS(ITL(ng)%pioFile%fh)) THEN
        type='ROMS tangent linear model initial file'
      ELSE IF (FileH.eq.ABS(LCZ(ng)%pioFile%fh)) THEN
        type='ROMS 4D-Var Lanczos vectors file'
      ELSE IF (FileH.eq.ABS(LZE(ng)%pioFile%fh)) THEN
        type='ROMS 4D-Var Evolved Lanczos vectors file'
      ELSE IF (FileH.eq.ABS(NRM(1,ng)%pioFile%fh)) THEN
        type='ROMS initial conditions error covariance norm file'
      ELSE IF (FileH.eq.ABS(NRM(2,ng)%pioFile%fh)) THEN
        type='ROMS model error covariance norm file'
      ELSE IF (FileH.eq.ABS(NRM(3,ng)%pioFile%fh)) THEN
        type='ROMS boundary conditions error covariance norm file'
      ELSE IF (FileH.eq.ABS(NRM(4,ng)%pioFile%fh)) THEN
        type='ROMS surface forcing error covariance norm file'
      ELSE IF (FileH.eq.ABS(QCK(ng)%pioFile%fh)) THEN
        type='ROMS quicksave file'
      ELSE IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
        type='ROMS restart file'
      ELSE IF (FileH.eq.ABS(STA(ng)%pioFile%fh)) THEN
        type='ROMS station file'
      ELSE IF (FileH.eq.ABS(TLF(ng)%pioFile%fh)) THEN
        type='ROMS tangent linear impulse forcing file'
      ELSE IF (FileH.eq.ABS(TLM(ng)%pioFile%fh)) THEN
        type='ROMS tangent linear history file'
      END IF
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_att(pioFile, PIO_global, 'type',                 &
     &                     TRIM(type))
        IF (FoundError(status, PIO_noerr, 3986, MyFile)) THEN
          IF (Master) WRITE (stdout,20) 'type', TRIM(ncname)
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
!  Define other global attributes to NetCDF file.
!
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_att(pioFile, PIO_global, 'title',                &
     &                     TRIM(title))
        IF (FoundError(status, PIO_noerr, 4012, MyFile)) THEN
          IF (Master) WRITE (stdout,20) 'title', TRIM(ncname)
          exit_flag=3
          ioerror=status
        END IF
      END IF
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_att(pioFile, PIO_global, 'var_info',             &
     &                     TRIM(varname))
        IF (FoundError(status, PIO_noerr, 4022, MyFile)) THEN
          IF (Master) WRITE (stdout,20) 'var_info', TRIM(ncname)
          exit_flag=3
          ioerror=status
        END IF
      END IF
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_att(pioFile, PIO_global, 'rst_file',             &
     &                     TRIM(RST(ng)%name))
        IF (FoundError(status, PIO_noerr, 4057, MyFile)) THEN
          IF (Master) WRITE (stdout,20) 'rst_file', TRIM(ncname)
          exit_flag=3
          ioerror=status
        END IF
      END IF
      IF (exit_flag.eq.NoError) THEN
        IF (LdefHIS(ng)) THEN
          IF (ndefHIS(ng).gt.0) THEN
            status=PIO_put_att(pioFile, PIO_global, 'his_base',         &
     &                         TRIM(HIS(ng)%base))
          ELSE
            status=PIO_put_att(pioFile, PIO_global, 'his_file',         &
     &                         TRIM(HIS(ng)%name))
          END IF
          IF (FoundError(status, PIO_noerr, 4073, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'his_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF
      IF (exit_flag.eq.NoError) THEN
        IF (ndefAVG(ng).gt.0) THEN
          status=PIO_put_att(pioFile, PIO_global, 'avg_base',           &
     &                       TRIM(AVG(ng)%base))
        ELSE
          status=PIO_put_att(pioFile, PIO_global, 'avg_file',           &
     &                       TRIM(AVG(ng)%name))
        END IF
        IF (FoundError(status, PIO_noerr, 4122, MyFile)) THEN
          IF (Master) WRITE (stdout,20) 'avg_file', TRIM(ncname)
          exit_flag=3
          ioerror=status
        END IF
      END IF
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_att(pioFile, PIO_global, 'grd_file',             &
     &                     TRIM(GRD(ng)%name))
        IF (FoundError(status, PIO_noerr, 4201, MyFile)) THEN
          IF (Master) WRITE (stdout,20) 'grd_file', TRIM(ncname)
          exit_flag=3
          ioerror=status
        END IF
      END IF
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_att(pioFile, PIO_global, 'ini_file',             &
     &                     TRIM(INI(ng)%name))
        IF (FoundError(status, PIO_noerr, 4214, MyFile)) THEN
          IF (Master) WRITE (stdout,20) 'ini_file', TRIM(ncname)
          exit_flag=3
          ioerror=status
        END IF
      END IF
      IF (exit_flag.eq.NoError) THEN
        IF (LuvSrc(ng).or.LwSrc(ng).or.(ANY(LtracerSrc(:,ng)))) THEN
          status=PIO_put_att(pioFile, PIO_global, 'river_file',         &
     &                       TRIM(SSF(ng)%name))
          IF (FoundError(status, PIO_noerr, 4381, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'river_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF
      IF (exit_flag.eq.NoError) THEN
        DO i=1,nFfiles(ng)
          CALL join_string (FRC(i,ng)%files, FRC(i,ng)%Nfiles,          &
     &                      string, lstr)
          WRITE (frcatt,30) 'frc_file_', i
          status=PIO_put_att(pioFile, PIO_global, frcatt,               &
     &                       string(1:lstr))
          IF (FoundError(status, PIO_noerr, 4412, MyFile)) THEN
            IF (Master) WRITE (stdout,20) TRIM(frcatt), TRIM(ncname)
            exit_flag=3
            ioerror=status
            EXIT
          END IF
        END DO
      END IF
      IF (ObcData(ng)) THEN
        DO i=1,nBCfiles(ng)
          IF (exit_flag.eq.NoError) THEN
            CALL join_string (BRY(i,ng)%files, BRY(i,ng)%Nfiles,        &
     &                        string, lstr)
            WRITE (bryatt,30) 'bry_file_', i
            status=PIO_put_att(pioFile, PIO_global, bryatt,             &
     &                          string(1:lstr))
            IF (FoundError(status, PIO_noerr, 4430, MyFile)) THEN
              IF (Master) WRITE (stdout,20) TRIM(bryatt), TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END DO
      END IF
      IF (Lclimatology(ng)) THEN
        DO i=1,nCLMfiles(ng)
          IF (exit_flag.eq.NoError) THEN
            CALL join_string (CLM(i,ng)%files, CLM(i,ng)%Nfiles,        &
     &                        string, lstr)
            WRITE (clmatt,30) 'clm_file_', i
            status=PIO_put_att(pioFile, PIO_global, clmatt,             &
     &                         string(1:lstr))
            IF (FoundError(status, PIO_noerr, 4449, MyFile)) THEN
              IF (Master) WRITE (stdout,20) TRIM(clmatt), TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END DO
      END IF
      IF (Lnudging(ng)) THEN
        IF (exit_flag.eq.NoError) THEN
          status=PIO_put_att(pioFile, PIO_global, 'nud_file',           &
     &                       TRIM(NUD(ng)%name))
          IF (FoundError(status, PIO_noerr, 4464, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'nud_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_att(pioFile, PIO_global, 'script_file',          &
     &                     TRIM(Iname))
        IF (FoundError(status, PIO_noerr, 4503, MyFile)) THEN
          IF (Master) WRITE (stdout,20) 'script_file', TRIM(ncname)
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
!  NLM tracer advection scheme.
!
      IF (exit_flag.eq.NoError) THEN
        CALL tadv_putatt (ng, pioFile, ncname, 'NLM_TADV',              &
     &                    Hadvection, Vadvection, status)
        IF (FoundError(status, PIO_noerr, 4567, MyFile)) THEN
          IF (Master) WRITE (stdout,20) 'NLM_TADV', TRIM(ncname)
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
!  NLM Lateral boundary conditions.
!
      IF (exit_flag.eq.NoError) THEN
        CALL lbc_putatt (ng, pioFile, ncname, 'NLM_LBC', LBC, status)
        IF (FoundError(status, PIO_noerr, 4591, MyFile)) THEN
          IF (Master) WRITE (stdout,20) 'NLM_LBC', TRIM(ncname)
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
!  GIT repository information.
!
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_att(pioFile, PIO_global, 'git_url',              &
     &                     TRIM(git_url))
        IF (FoundError(status, PIO_noerr, 4619, MyFile)) THEN
          IF (Master) WRITE (stdout,20) 'git_url', TRIM(ncname)
          exit_flag=3
          ioerror=status
        END IF
      END IF
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_att(pioFile, PIO_global, 'git_rev',              &
     &                     TRIM(git_rev))
        IF (FoundError(status, PIO_noerr, 4630, MyFile)) THEN
          IF (Master) WRITE (stdout,20) 'git_rev', TRIM(ncname)
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
!  SVN repository information.
!
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_att(pioFile, PIO_global, 'svn_url',              &
     &                     TRIM(svn_url))
        IF (FoundError(status, PIO_noerr, 4643, MyFile)) THEN
          IF (Master) WRITE (stdout,20) 'svn_url', TRIM(ncname)
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
!  Local root directory, cpp header directory and file, and analytical
!  directory
!
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_att(pioFile, PIO_global, 'code_dir',             &
     &                     TRIM(Rdir))
        IF (FoundError(status, PIO_noerr, 4670, MyFile)) THEN
          IF (Master) WRITE (stdout,20) 'code_dir', TRIM(ncname)
          exit_flag=3
          ioerror=status
        END IF
      END IF
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_att(pioFile, PIO_global, 'header_dir',           &
     &                     TRIM(Hdir))
        IF (FoundError(status, PIO_noerr, 4682, MyFile)) THEN
          IF (Master) WRITE (stdout,20) 'header_dir', TRIM(ncname)
          exit_flag=3
          ioerror=status
        END IF
      END IF
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_att(pioFile, PIO_global, 'header_file',          &
     &                     TRIM(Hfile))
        IF (FoundError(status, PIO_noerr, 4694, MyFile)) THEN
          IF (Master) WRITE (stdout,20) 'header_file', TRIM(ncname)
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
!  Attributes describing platform and compiler
!
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_att(pioFile, PIO_global, 'os',                   &
     &                     TRIM(my_os))
        IF (FoundError(status, PIO_noerr, 4709, MyFile)) THEN
          IF (Master) WRITE (stdout,20) 'os', TRIM(ncname)
          exit_flag=3
          ioerror=status
        END IF
      END IF
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_att(pioFile, PIO_global, 'cpu',                  &
     &                     TRIM(my_cpu))
        IF (FoundError(status, PIO_noerr, 4719, MyFile)) THEN
          IF (Master) WRITE (stdout,20) 'cpu', TRIM(ncname)
          exit_flag=3
          ioerror=status
        END IF
      END IF
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_att(pioFile, PIO_global, 'compiler_system',      &
     &                     TRIM(my_fort))
        IF (FoundError(status, PIO_noerr, 4729, MyFile)) THEN
          IF (Master) WRITE (stdout,20) 'compiler_system',              &
     &                                  TRIM(ncname)
          exit_flag=3
          ioerror=status
        END IF
      END IF
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_att(pioFile, PIO_global, 'compiler_command',     &
     &                     TRIM(my_fc))
        IF (FoundError(status, PIO_noerr, 4740, MyFile)) THEN
          IF (Master) WRITE (stdout,20) 'compiler_command',             &
     &                                    TRIM(ncname)
          exit_flag=3
          ioerror=status
        END IF
      END IF
      IF (exit_flag.eq.NoError) THEN
        lstr=INDEX(my_fflags, 'free')-2
        IF (lstr.le.0) lstr=LEN_TRIM(my_fflags)
        status=PIO_put_att(pioFile, PIO_global, 'compiler_flags',       &
     &                     my_fflags(1:lstr))
        IF (FoundError(status, PIO_noerr, 4753, MyFile)) THEN
          IF (Master) WRITE (stdout,20) 'compiler_flags', TRIM(ncname)
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
!  Tiling and history attributes.
!
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_att(pioFile, PIO_global, 'tiling',               &
     &                     TRIM(tiling))
        IF (FoundError(status, PIO_noerr, 4765, MyFile)) THEN
          IF (Master) WRITE (stdout,20) 'tiling', TRIM(ncname)
          exit_flag=3
          ioerror=status
        END IF
      END IF
      IF (exit_flag.eq.NoError) THEN
        status=PIO_put_att(pioFile, PIO_global, 'history',              &
     &                     TRIM(history))
        IF (FoundError(status, PIO_noerr, 4775, MyFile)) THEN
          IF (Master) WRITE (stdout,20) 'history', TRIM(ncname)
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
!  Analytical header files used.
!
      IF (exit_flag.eq.NoError) THEN
        CALL join_string (ANANAME, SIZE(ANANAME), string, lstr)
        IF (lstr.gt.0) THEN
          status=PIO_put_att(pioFile, PIO_global, 'ana_file',           &
     &                       string(1:lstr))
          IF (FoundError(status, PIO_noerr, 4789, MyFile)) THEN
            IF (Master) WRITE (stdout,20) 'ana_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF
!
!  Activated CPP options.
!
      IF (exit_flag.eq.NoError) THEN
        lstr=LEN_TRIM(Coptions)-1
        status=PIO_put_att(pioFile, PIO_global, 'CPP_options',          &
     &                     TRIM(Coptions(1:lstr)))
        IF (FoundError(status, PIO_noerr, 4824, MyFile)) THEN
          IF (Master) WRITE (stdout,20) 'CPP_options', TRIM(ncname)
          exit_flag=3
          ioerror=status
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Define running parameters.
!-----------------------------------------------------------------------
!
!  Time stepping parameters.
!
      Vinfo( 1)='ntimes'
      Vinfo( 2)='number of long time-steps'
      status=def_var(ng, model, pioFile, pioVar, PIO_int,               &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 4850, MyFile)) RETURN
      Vinfo( 1)='ndtfast'
      Vinfo( 2)='number of short time-steps'
      status=def_var(ng, model, pioFile, pioVar, PIO_int,               &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 4857, MyFile)) RETURN
      Vinfo( 1)='dt'
      Vinfo( 2)='size of long time-steps'
      Vinfo( 3)='second'
      status=def_var(ng, model, pioFile, pioVar, PIO_TOUT,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 4865, MyFile)) RETURN
      Vinfo( 1)='dtfast'
      Vinfo( 2)='size of short time-steps'
      Vinfo( 3)='second'
      status=def_var(ng, model, pioFile, pioVar, PIO_TOUT,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 4873, MyFile)) RETURN
      Vinfo( 1)='dstart'
      Vinfo( 2)='time stamp assigned to model initilization'
      WRITE (Vinfo( 3),'(a,a)') 'days since ', TRIM(Rclock%string)
      Vinfo( 4)=TRIM(Rclock%calendar)
      status=def_var(ng, model, pioFile, pioVar, PIO_TOUT,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 4882, MyFile)) RETURN
      Vinfo( 1)='nHIS'
      Vinfo( 2)='number of time-steps between history records'
      status=def_var(ng, model, pioFile, pioVar, PIO_int,               &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 4930, MyFile)) RETURN
      Vinfo( 1)='ndefHIS'
      Vinfo( 2)=                                                        &
     &    'number of time-steps between the creation of history files'
      status=def_var(ng, model, pioFile, pioVar, PIO_int,               &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 4938, MyFile)) RETURN
      Vinfo( 1)='nRST'
      Vinfo( 2)='number of time-steps between restart records'
      IF (LcycleRST(ng)) THEN
        Vinfo(13)='only latest two records are maintained'
      END IF
      status=def_var(ng, model, pioFile, pioVar, PIO_int,               &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 4972, MyFile)) RETURN
      Vinfo( 1)='ntsAVG'
      Vinfo( 2)=                                                        &
     &   'starting time-step for accumulation of time-averaged fields'
      status=def_var(ng, model, pioFile, pioVar, PIO_int,               &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 4984, MyFile)) RETURN
      Vinfo( 1)='nAVG'
      Vinfo( 2)='number of time-steps between time-averaged records'
      status=def_var(ng, model, pioFile, pioVar, PIO_int,               &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 4991, MyFile)) RETURN
      Vinfo( 1)='ndefAVG'
      Vinfo( 2)=                                                        &
     &    'number of time-steps between the creation of average files'
      status=def_var(ng, model, pioFile, pioVar, PIO_int,               &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 4999, MyFile)) RETURN
!
!  Power-law shape filter parameters for time-averaging of barotropic
!  fields.
!
      Vinfo( 1)='Falpha'
      Vinfo( 2)='Power-law shape barotropic filter parameter'
      status=def_var(ng, model, pioFile, pioVar, PIO_TOUT,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5171, MyFile)) RETURN
      Vinfo( 1)='Fbeta'
      Vinfo( 2)='Power-law shape barotropic filter parameter'
      status=def_var(ng, model, pioFile, pioVar, PIO_TOUT,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5178, MyFile)) RETURN
      Vinfo( 1)='Fgamma'
      Vinfo( 2)='Power-law shape barotropic filter parameter'
      status=def_var(ng, model, pioFile, pioVar, PIO_TOUT,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5185, MyFile)) RETURN
!
!  Horizontal mixing coefficients.
!
      Vinfo( 1)='nl_tnu2'
      Vinfo( 2)='nonlinear model Laplacian mixing coefficient '//       &
     &          'for tracers'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5198, MyFile)) RETURN
      Vinfo( 1)='nl_visc2'
      Vinfo( 2)='nonlinear model Laplacian mixing coefficient '//       &
     &          'for momentum'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5264, MyFile)) RETURN
      Vinfo( 1)='LuvSponge'
      Vinfo( 2)='horizontal viscosity sponge activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, pioFile, pioVar, PIO_int,               &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5350, MyFile)) RETURN
      Vinfo( 1)='LtracerSponge'
      Vinfo( 2)='horizontal diffusivity sponge activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, pioFile, pioVar, PIO_int,               &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5360, MyFile)) RETURN
!
!  Background vertical mixing coefficients.
!
      Vinfo( 1)='Akt_bak'
      Vinfo( 2)='background vertical mixing coefficient for tracers'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5372, MyFile)) RETURN
      Vinfo( 1)='Akv_bak'
      Vinfo( 2)='background vertical mixing coefficient for momentum'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5380, MyFile)) RETURN
      Vinfo( 1)='Akk_bak'
      Vinfo( 2)=                                                        &
     &   'background vertical mixing coefficient for turbulent energy'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5390, MyFile)) RETURN
      Vinfo( 1)='Akp_bak'
      Vinfo( 2)=                                                        &
     &   'background vertical mixing coefficient for length scale'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5399, MyFile)) RETURN
!
!  Drag coefficients.
!
      Vinfo( 1)='rdrg'
      Vinfo( 2)='linear drag coefficient'
      Vinfo( 3)='meter second-1'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5456, MyFile)) RETURN
      Vinfo( 1)='rdrg2'
      Vinfo( 2)='quadratic drag coefficient'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo ,ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5463, MyFile)) RETURN
      Vinfo( 1)='Zob'
      Vinfo( 2)='bottom roughness'
      Vinfo( 3)='meter'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5472, MyFile)) RETURN
      Vinfo( 1)='Zos'
      Vinfo( 2)='surface roughness'
      Vinfo( 3)='meter'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5480, MyFile)) RETURN
!
!  Generic length-scale parameters.
!
      Vinfo( 1)='gls_p'
      Vinfo( 2)='stability exponent'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5491, MyFile)) RETURN
      Vinfo( 1)='gls_m'
      Vinfo( 2)='turbulent kinetic energy exponent'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5498, MyFile)) RETURN
      Vinfo( 1)='gls_n'
      Vinfo( 2)='turbulent length scale exponent'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5505, MyFile)) RETURN
      Vinfo( 1)='gls_cmu0'
      Vinfo( 2)='stability coefficient'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5512, MyFile)) RETURN
      Vinfo( 1)='gls_c1'
      Vinfo( 2)='shear production coefficient'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5519, MyFile)) RETURN
      Vinfo( 1)='gls_c2'
      Vinfo( 2)='dissipation coefficient'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5526, MyFile)) RETURN
      Vinfo( 1)='gls_c3m'
      Vinfo( 2)='buoyancy production coefficient (minus)'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5533, MyFile)) RETURN
      Vinfo( 1)='gls_c3p'
      Vinfo( 2)='buoyancy production coefficient (plus)'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5540, MyFile)) RETURN
      Vinfo( 1)='gls_sigk'
      Vinfo( 2)='constant Schmidt number for TKE'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5547, MyFile)) RETURN
      Vinfo( 1)='gls_sigp'
      Vinfo( 2)='constant Schmidt number for PSI'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5554, MyFile)) RETURN
      Vinfo( 1)='gls_Kmin'
      Vinfo( 2)='minimum value of specific turbulent kinetic energy'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5561, MyFile)) RETURN
      Vinfo( 1)='gls_Pmin'
      Vinfo( 2)='minimum Value of dissipation'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5568, MyFile)) RETURN
      Vinfo( 1)='Charnok_alpha'
      Vinfo( 2)='Charnock factor for surface roughness'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5575, MyFile)) RETURN
      Vinfo( 1)='Zos_hsig_alpha'
      Vinfo( 2)='wave amplitude factor for surface roughness'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5582, MyFile)) RETURN
      Vinfo( 1)='sz_alpha'
      Vinfo( 2)='surface flux from wave dissipation'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5589, MyFile)) RETURN
      Vinfo( 1)='CrgBan_cw'
      Vinfo( 2)='surface flux due to Craig and Banner wave breaking'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5596, MyFile)) RETURN
!
!  Nudging inverse time scales used in various tasks.
!
      Vinfo( 1)='Znudg'
      Vinfo( 2)='free-surface nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, pioFile, pioVar, PIO_TOUT,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5617, MyFile)) RETURN
      Vinfo( 1)='M2nudg'
      Vinfo( 2)='2D momentum nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, pioFile, pioVar, PIO_TOUT,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5625, MyFile)) RETURN
      Vinfo( 1)='M3nudg'
      Vinfo( 2)='3D momentum nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, pioFile, pioVar, PIO_TOUT,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5634, MyFile)) RETURN
      Vinfo( 1)='Tnudg'
      Vinfo( 2)='Tracers nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, pioFile, pioVar, PIO_TOUT,              &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5642, MyFile)) RETURN
!
!  Open boundary nudging, inverse time scales.
!
      IF (NudgingCoeff(ng)) THEN
        Vinfo( 1)='FSobc_in'
        Vinfo( 2)='free-surface inflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, pioFile, pioVar, PIO_TOUT,            &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 5656, MyFile)) RETURN
        Vinfo( 1)='FSobc_out'
        Vinfo( 2)='free-surface outflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, pioFile, pioVar, PIO_TOUT,            &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 5664, MyFile)) RETURN
        Vinfo( 1)='M2obc_in'
        Vinfo( 2)='2D momentum inflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, pioFile, pioVar, PIO_TOUT,            &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 5672, MyFile)) RETURN
        Vinfo( 1)='M2obc_out'
        Vinfo( 2)='2D momentum outflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, pioFile, pioVar, PIO_TOUT,            &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 5680, MyFile)) RETURN
        Vinfo( 1)='Tobc_in'
        Vinfo( 2)='tracers inflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, pioFile, pioVar, PIO_TOUT,            &
     &                 2, tbrydim, Aval, Vinfo, ncname,                 &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 5689, MyFile)) RETURN
        Vinfo( 1)='Tobc_out'
        Vinfo( 2)='tracers outflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, pioFile, pioVar, PIO_TOUT,            &
     &                 2, tbrydim, Aval, Vinfo, ncname,                 &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 5697, MyFile)) RETURN
        Vinfo( 1)='M3obc_in'
        Vinfo( 2)='3D momentum inflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, pioFile, pioVar, PIO_TOUT,            &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 5705, MyFile)) RETURN
        Vinfo( 1)='M3obc_out'
        Vinfo( 2)='3D momentum outflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, pioFile, pioVar, PIO_TOUT,            &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 5713, MyFile)) RETURN
      END IF
!
!  Equation of State parameters.
!
      Vinfo( 1)='rho0'
      Vinfo( 2)='mean density used in Boussinesq approximation'
      Vinfo( 3)='kilogram meter-3'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5726, MyFile)) RETURN
!
!  Various parameters.
!
!
!  Slipperiness parameters.
!
      Vinfo( 1)='gamma2'
      Vinfo( 2)='slipperiness parameter'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5791, MyFile)) RETURN
!
! Logical switches to activate horizontal momentum transport
! point Sources/Sinks (like river runoff transport) and mass point
! Sources/Sinks (like volume vertical influx).
!
      Vinfo( 1)='LuvSrc'
      Vinfo( 2)='momentum point sources and sink activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, pioFile, pioVar, PIO_int,               &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5804, MyFile)) RETURN
      Vinfo( 1)='LwSrc'
      Vinfo( 2)='mass point sources and sink activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, pioFile, pioVar, PIO_int,               &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5813, MyFile)) RETURN
!
!  Logical switches indicating which tracer variables are processed
!  during point Sources/Sinks.
!
      Vinfo( 1)='LtracerSrc'
      Vinfo( 2)='tracer point sources and sink activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, pioFile, pioVar, PIO_int,               &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5827, MyFile)) RETURN
!
!  Logical switches to process climatology fields.
!
      Vinfo( 1)='LsshCLM'
      Vinfo( 2)='sea surface height climatology processing switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, pioFile, pioVar, PIO_int,               &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5839, MyFile)) RETURN
      Vinfo( 1)='Lm2CLM'
      Vinfo( 2)='2D momentum climatology processing switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, pioFile, pioVar, PIO_int,               &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5848, MyFile)) RETURN
      Vinfo( 1)='Lm3CLM'
      Vinfo( 2)='3D momentum climatology processing switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, pioFile, pioVar, PIO_int,               &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5858, MyFile)) RETURN
      Vinfo( 1)='LtracerCLM'
      Vinfo( 2)='tracer climatology processing switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, pioFile, pioVar, PIO_int,               &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5867, MyFile)) RETURN
!
!  Logical switches for nudging of climatology fields.
!
      Vinfo( 1)='LnudgeM2CLM'
      Vinfo( 2)='2D momentum climatology nudging activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, pioFile, pioVar, PIO_int,               &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5879, MyFile)) RETURN
!
      Vinfo( 1)='LnudgeM3CLM'
      Vinfo( 2)='3D momentum climatology nudging activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, pioFile, pioVar, PIO_int,               &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5889, MyFile)) RETURN
!
      Vinfo( 1)='LnudgeTCLM'
      Vinfo( 2)='tracer climatology nudging activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, pioFile, pioVar, PIO_int,               &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 5898, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Define grid variables.
!-----------------------------------------------------------------------
!
!  Grid type switch: Spherical or Cartesian. Writing characters in
!  parallel I/O is extremely inefficient.  It is better to write
!  this as an integer switch: 0=Cartesian, 1=spherical.
!
      Vinfo( 1)='spherical'
      Vinfo( 2)='grid type logical switch'
      Vinfo( 9)='Cartesian'
      Vinfo(10)='spherical'
      status=def_var(ng, model, pioFile, pioVar, PIO_int,               &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 6693, MyFile)) RETURN
!
!  Domain Length.
!
      Vinfo( 1)='xl'
      Vinfo( 2)='domain length in the XI-direction'
      Vinfo( 3)='meter'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 6703, MyFile)) RETURN
      Vinfo( 1)='el'
      Vinfo( 2)='domain length in the ETA-direction'
      Vinfo( 3)='meter'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 6711, MyFile)) RETURN
!
!  S-coordinate parameters.
!
      Vinfo( 1)='Vtransform'
      Vinfo( 2)='vertical terrain-following transformation equation'
      status=def_var(ng, model, pioFile, pioVar, PIO_int,               &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 6721, MyFile)) RETURN
      Vinfo( 1)='Vstretching'
      Vinfo( 2)='vertical terrain-following stretching function'
      status=def_var(ng, model, pioFile, pioVar, PIO_int,               &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 6728, MyFile)) RETURN
      Vinfo( 1)='theta_s'
      Vinfo( 2)='S-coordinate surface control parameter'
      status=def_var(ng, model, pioFile, pioVar, PIO_TOUT,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 6735, MyFile)) RETURN
      Vinfo( 1)='theta_b'
      Vinfo( 2)='S-coordinate bottom control parameter'
      status=def_var(ng, model, pioFile, pioVar, PIO_TOUT,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 6742, MyFile)) RETURN
      Vinfo( 1)='Tcline'
      Vinfo( 2)='S-coordinate surface/bottom layer width'
      Vinfo( 3)='meter'
      status=def_var(ng, model, pioFile, pioVar, PIO_TOUT,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 6750, MyFile)) RETURN
      Vinfo( 1)='hc'
      Vinfo( 2)='S-coordinate parameter, critical depth'
      Vinfo( 3)='meter'
      status=def_var(ng, model, pioFile, pioVar, PIO_TOUT,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 6758, MyFile)) RETURN
!
!  SGRID conventions for staggered data on structured grids.
!
      Vinfo( 1)='grid'
      status=def_var(ng, model, pioFile, pioVar, PIO_int,               &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 6766, MyFile)) RETURN
!
!  S-coordinate non-dimensional independent variable at RHO-points.
!
      Vinfo( 1)='s_rho'
      Vinfo( 2)='S-coordinate at RHO-points'
      Vinfo( 5)='valid_min'
      Vinfo( 6)='valid_max'
      IF (Vtransform(ng).eq.1) THEN
        Vinfo(21)='ocean_s_coordinate_g1'
      ELSE IF (Vtransform(ng).eq.2) THEN
        Vinfo(21)='ocean_s_coordinate_g2'
      END IF
      Vinfo(23)='s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc'
      vinfo(25)='up'
      Aval(2)=-1.0_r8
      Aval(3)=0.0_r8
      status=def_var(ng, model, pioFile, pioVar, PIO_TOUT,              &
     &               1, (/srdim/), Aval, Vinfo, ncname,                 &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 6790, MyFile)) RETURN
!
!  S-coordinate non-dimensional independent variable at W-points.
!
      Vinfo( 1)='s_w'
      Vinfo( 2)='S-coordinate at W-points'
      Vinfo( 5)='valid_min'
      Vinfo( 6)='valid_max'
      Vinfo(21)='ocean_s_coordinate'
      IF (Vtransform(ng).eq.1) THEN
        Vinfo(21)='ocean_s_coordinate_g1'
      ELSE IF (Vtransform(ng).eq.2) THEN
        Vinfo(21)='ocean_s_coordinate_g2'
      END IF
      Vinfo(23)='s: s_w C: Cs_w eta: zeta depth: h depth_c: hc'
      vinfo(25)='up'
      Aval(2)=-1.0_r8
      Aval(3)=0.0_r8
      status=def_var(ng, model, pioFile, pioVar, PIO_TOUT,              &
     &               1, (/swdim/), Aval, Vinfo, ncname,                 &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 6815, MyFile)) RETURN
!
!  S-coordinate non-dimensional stretching curves at RHO-points.
!
      Vinfo( 1)='Cs_r'
      Vinfo( 2)='S-coordinate stretching curves at RHO-points'
      Vinfo( 5)='valid_min'
      Vinfo( 6)='valid_max'
      Aval(2)=-1.0_r8
      Aval(3)=0.0_r8
      status=def_var(ng, model, pioFile, pioVar, PIO_TOUT,              &
     &               1, (/srdim/), Aval, Vinfo, ncname,                 &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 6828, MyFile)) RETURN
!
!  S-coordinate non-dimensional stretching curves at W-points.
!
      Vinfo( 1)='Cs_w'
      Vinfo( 2)='S-coordinate stretching curves at W-points'
      Vinfo( 5)='valid_min'
      Vinfo( 6)='valid_max'
      Aval(2)=-1.0_r8
      Aval(3)=0.0_r8
      status=def_var(ng, model, pioFile, pioVar, PIO_TOUT,              &
     &               1, (/swdim/), Aval, Vinfo, ncname,                 &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, 6841, MyFile)) RETURN
!
!  Depth of horizontal slices.
!
      IF (Nslice.gt.0) THEN
        Vinfo( 1)='z_slice'
        Vinfo( 2)='constant depth of output fields horizontal slices'
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,            &
     &                 1, (/slicedim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 6854, MyFile)) RETURN
      END IF
!
!  User generic parameters.
!
      IF (Nuser.gt.0) THEN
        Vinfo( 1)='user'
        Vinfo( 2)='user generic parameters'
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,            &
     &                 1, (/usrdim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, 6867, MyFile)) RETURN
      END IF
      IF (FileH.ne.ABS(FLT(ng)%pioFile%fh)) THEN
!
!  Bathymetry.
!
        Vinfo( 1)=Vname(1,idtopo)
        Vinfo( 2)=Vname(2,idtopo)
        Vinfo( 3)=Vname(3,idtopo)
        Vinfo(14)=Vname(4,idtopo)
        Vinfo(21)=Vname(6,idtopo)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idtopo,ng),r8)
        IF (FileH.eq.ABS(STA(ng)%pioFile%fh)) THEN
          status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,          &
     &                   1, (/stadim/), Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 6908, MyFile)) RETURN
        ELSE
          status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,          &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 6912, MyFile)) RETURN
        END IF
!
!  Coriolis Parameter.
!
        IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
          Vinfo( 1)=Vname(1,idfcor)
          Vinfo( 2)=Vname(2,idfcor)
          Vinfo( 3)=Vname(3,idfcor)
          Vinfo(14)=Vname(4,idfcor)
          Vinfo(21)=Vname(6,idfcor)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idfcor,ng),r8)
          status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,          &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 6928, MyFile)) RETURN
        END IF
!
!  Curvilinear coordinate metrics.
!
        IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
          Vinfo( 1)=Vname(1,idpmdx)
          Vinfo( 2)=Vname(2,idpmdx)
          Vinfo( 3)=Vname(3,idpmdx)
          Vinfo(14)=Vname(4,idpmdx)
          Vinfo(21)=Vname(6,idpmdx)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idpmdx,ng),r8)
          status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,          &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 6943, MyFile)) RETURN
!
          Vinfo( 1)=Vname(1,idpndy)
          Vinfo( 2)=Vname(2,idpndy)
          Vinfo( 3)=Vname(3,idpndy)
          Vinfo(14)=Vname(4,idpndy)
          Vinfo(21)=Vname(6,idpndy)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idpndy,ng),r8)
          status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,          &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 6954, MyFile)) RETURN
        END IF
!
!  Grid coordinates of RHO-points.
!
        IF (spherical) THEN
          Vinfo( 1)=Vname(1,idLonR)
          Vinfo( 2)=Vname(2,idLonR)
          Vinfo( 3)=Vname(3,idLonR)
          Vinfo(14)=Vname(4,idLonR)
          Vinfo(21)=Vname(6,idLonR)
          IF (FileH.eq.ABS(STA(ng)%pioFile%fh)) THEN
            status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,        &
     &                     1, (/stadim/), Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 6968, MyFile)) RETURN
          ELSE
            status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,        &
     &                     2, t2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 6972, MyFile)) RETURN
          END IF
!
          Vinfo( 1)=Vname(1,idLatR)
          Vinfo( 2)=Vname(2,idLatR)
          Vinfo( 3)=Vname(3,idLatR)
          Vinfo(14)=Vname(4,idLatR)
          Vinfo(21)=Vname(6,idLatR)
          IF (FileH.eq.ABS(STA(ng)%pioFile%fh)) THEN
            status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,        &
     &                     1, (/stadim/), Aval, Vinfo,  ncname)
            IF (FoundError(exit_flag, NoError, 6983, MyFile)) RETURN
          ELSE
            status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,        &
     &                     2, t2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 6987, MyFile)) RETURN
          END IF
        ELSE
          Vinfo( 1)=Vname(1,idXgrR)
          Vinfo( 2)=Vname(2,idXgrR)
          Vinfo( 3)=Vname(3,idXgrR)
          Vinfo(14)=Vname(4,idXgrR)
          Vinfo(21)=Vname(6,idXgrR)
          IF (FileH.eq.ABS(STA(ng)%pioFile%fh)) THEN
            status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,        &
     &                     1, (/stadim/), Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 6998, MyFile)) RETURN
          ELSE
            status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,        &
     &                     2, t2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 7002, MyFile)) RETURN
          END IF
!
          Vinfo( 1)=Vname(1,idYgrR)
          Vinfo( 2)=Vname(2,idYgrR)
          Vinfo( 3)=Vname(3,idYgrR)
          Vinfo(14)=Vname(4,idYgrR)
          Vinfo(21)=Vname(6,idYgrR)
          IF (FileH.eq.ABS(STA(ng)%pioFile%fh)) THEN
            status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,        &
     &                     1, (/stadim/), Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 7013, MyFile)) RETURN
          ELSE
            status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,        &
     &                     2, t2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 7017, MyFile)) RETURN
          END IF
        END IF
!
!  Grid coordinates of U-points.
!
        IF (spherical) THEN
          Vinfo( 1)=Vname(1,idLonU)
          Vinfo( 2)=Vname(2,idLonU)
          Vinfo( 3)=Vname(3,idLonU)
          Vinfo(14)=Vname(4,idLonU)
          Vinfo(21)=Vname(6,idLonU)
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,        &
     &                     2, u2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 7032, MyFile)) RETURN
          END IF
!
          Vinfo( 1)=Vname(1,idLatU)
          Vinfo( 2)=Vname(2,idLatU)
          Vinfo( 3)=Vname(3,idLatU)
          Vinfo(14)=Vname(4,idLatU)
          Vinfo(21)=Vname(6,idLatU)
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,        &
     &                     2, u2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 7043, MyFile)) RETURN
          END IF
        ELSE
          Vinfo( 1)=Vname(1,idXgrU)
          Vinfo( 2)=Vname(2,idXgrU)
          Vinfo( 3)=Vname(3,idXgrU)
          Vinfo(14)=Vname(4,idXgrU)
          Vinfo(21)=Vname(6,idXgrU)
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,        &
     &                     2, u2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 7054, MyFile)) RETURN
          END IF
!
          Vinfo( 1)=Vname(1,idYgrU)
          Vinfo( 2)=Vname(2,idYgrU)
          Vinfo( 3)=Vname(3,idYgrU)
          Vinfo(14)=Vname(4,idYgrU)
          Vinfo(21)=Vname(6,idYgrU)
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,        &
     &                     2, u2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 7065, MyFile)) RETURN
          END IF
        END IF
!
!  Grid coordinates of V-points.
!
        IF (spherical) THEN
          Vinfo( 1)=Vname(1,idLonV)
          Vinfo( 2)=Vname(2,idLonV)
          Vinfo( 3)=Vname(3,idLonV)
          Vinfo(14)=Vname(4,idLonV)
          Vinfo(21)=Vname(6,idLonV)
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,        &
     &                     2, v2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 7080, MyFile)) RETURN
          END IF
!
          Vinfo( 1)=Vname(1,idLatV)
          Vinfo( 2)=Vname(2,idLatV)
          Vinfo( 3)=Vname(3,idLatV)
          Vinfo(14)=Vname(4,idLatV)
          Vinfo(21)=Vname(6,idLatV)
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,        &
     &                     2, v2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 7091, MyFile)) RETURN
          END IF
        ELSE
          Vinfo( 1)=Vname(1,idXgrV)
          Vinfo( 2)=Vname(2,idXgrV)
          Vinfo( 3)=Vname(3,idXgrV)
          Vinfo(14)=Vname(4,idXgrV)
          Vinfo(21)=Vname(6,idXgrV)
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,        &
     &                     2, v2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 7102, MyFile)) RETURN
          END IF
!
          Vinfo( 1)=Vname(1,idYgrV)
          Vinfo( 2)=Vname(2,idYgrV)
          Vinfo( 3)=Vname(3,idYgrV)
          Vinfo(14)=Vname(4,idYgrV)
          Vinfo(21)=Vname(6,idYgrV)
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,        &
     &                     2, v2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 7113, MyFile)) RETURN
          END IF
        END IF
!
!  Grid coordinates of PSI-points.
!
        IF (spherical) THEN
          Vinfo( 1)=Vname(1,idLonP)
          Vinfo( 2)=Vname(2,idLonP)
          Vinfo( 3)=Vname(3,idLonP)
          Vinfo(14)=Vname(4,idLonP)
          Vinfo(21)=Vname(6,idLonP)
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,        &
     &                     2, p2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 7128, MyFile)) RETURN
          END IF
!
          Vinfo( 1)=Vname(1,idLatP)
          Vinfo( 2)=Vname(2,idLatP)
          Vinfo( 3)=Vname(3,idLatP)
          Vinfo(14)=Vname(4,idLatP)
          Vinfo(21)=Vname(6,idLatP)
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,        &
     &                     2, p2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 7139, MyFile)) RETURN
          END IF
        ELSE
          Vinfo( 1)=Vname(1,idXgrP)
          Vinfo( 2)=Vname(2,idXgrP)
          Vinfo( 3)=Vname(3,idXgrP)
          Vinfo(14)=Vname(4,idXgrP)
          Vinfo(21)=Vname(6,idXgrP)
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,        &
     &                     2, p2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 7150, MyFile)) RETURN
          END IF
!
          Vinfo( 1)=Vname(1,idYgrP)
          Vinfo( 2)=Vname(2,idYgrP)
          Vinfo( 3)=Vname(3,idYgrP)
          Vinfo(14)=Vname(4,idYgrP)
          Vinfo(21)=Vname(6,idYgrP)
          IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
            status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,        &
     &                     2, p2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 7161, MyFile)) RETURN
          END IF
        END IF
!
!  Angle between XI-axis and EAST at RHO-points.
!
        Vinfo( 1)=Vname(1,idangR)
        Vinfo( 2)=Vname(2,idangR)
        Vinfo( 3)=Vname(3,idangR)
        Vinfo(14)=Vname(4,idangR)
        Vinfo(21)=Vname(6,idangR)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idangR,ng),r8)
        IF (FileH.eq.ABS(STA(ng)%pioFile%fh)) THEN
          status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,          &
     &                   1, (/stadim/), Aval, Vinfo,  ncname)
          IF (FoundError(exit_flag, NoError, 7179, MyFile)) RETURN
        ELSE
          status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,          &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 7183, MyFile)) RETURN
        END IF
!
!  Masking fields at RHO-, U-, V-points, and PSI-points.
!
        IF (FileH.ne.ABS(STA(ng)%pioFile%fh)) THEN
          Vinfo( 1)=Vname(1,idmskR)
          Vinfo( 2)=Vname(2,idmskR)
          Vinfo( 9)='land'
          Vinfo(10)='water'
          Vinfo(21)=Vname(6,idmskR)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idmskR,ng),r8)
          status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,          &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 7200, MyFile)) RETURN
!
          Vinfo( 1)=Vname(1,idmskU)
          Vinfo( 2)=Vname(2,idmskU)
          Vinfo( 9)='land'
          Vinfo(10)='water'
          Vinfo(21)=Vname(6,idmskU)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idmskU,ng),r8)
          status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,          &
     &                   2, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 7211, MyFile)) RETURN
!
          Vinfo( 1)=Vname(1,idmskV)
          Vinfo( 2)=Vname(2,idmskV)
          Vinfo( 9)='land'
          Vinfo(10)='water'
          Vinfo(21)=Vname(6,idmskV)
          status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,          &
     &                   2, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 7220, MyFile)) RETURN
!
          Vinfo( 1)=Vname(1,idmskP)
          Vinfo( 2)=Vname(2,idmskP)
          Vinfo( 9)='land'
          Vinfo(10)='water'
          Vinfo(21)=Vname(6,idmskP)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idmskP,ng),r8)
          status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,          &
     &                   2, p2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 7231, MyFile)) RETURN
        END IF
      END IF
!
  10  FORMAT (i0,'x',i0)
  20  FORMAT (/,' DEF_INFO_PIO - error while creating global',          &
     &        ' attribute: ',a,/,16x,a)
  30  FORMAT (a,i2.2)
!
      RETURN
      END SUBROUTINE def_info_pio
      END MODULE def_info_mod
