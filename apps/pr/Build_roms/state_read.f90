      MODULE state_read_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This routine reads the ROMS ocean state from requested NetCDF file. !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     tile         Domaing partition (integer)                         !
!     model        Calling model identifier (integer)                  !
!     IOtype       I/O NetCDF library to use (integer)                 !
!     LBi          I-dimension Lower bound (integer)                   !
!     UBi          I-dimension Upper bound (integer)                   !
!     LBj          J-dimension Lower bound (integer)                   !
!     UBj          J-dimension Upper bound (integer)                   !
!     LBij         IJ-dimension Lower bound (integer)                  !
!     UBij         IJ-dimension Upper bound (integer)                  !
!     Lout         State array index to process (integer)              !
!     rec          NetCDF file time record to read (integer)           !
!     nopen        Flag to open NetCDF file (integer)                  !
!                    nopen = 0   => file already opened for access     !
!                    nopen > 0   => open at entry and close at exit    !
!     ncid         NetCDF file ID, if already opened (integer)         !
!     pioFile      PIO file descriptor, TYPE(File_desc_t)              !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     ncname       NetCDF filename (string)                            !
!     rmask        State land/sea mask on RHO-points (real 2D array)   !
!     umask        State land/sea mask on U-points   (real 2D array)   !
!     vmask        State land/sea mask on V-points   (real 2D array)   !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     s_t          State tracers               (:,:,N,Lout,NT)         !
!     s_u          State 3D U-velocity         (:,:,N,Lout)            !
!     s_v          State 3D V-velocity         (:,:,N,Lout)            !
!     s_zeta       State free-surface          (:,:,Lout)              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      USE dateclock_mod,      ONLY : time_string
      USE distribute_mod,     ONLY : mp_bcasti
      USE nf_fread2d_mod,     ONLY : nf_fread2d
      USE nf_fread3d_mod,     ONLY : nf_fread3d
      USE strings_mod,        ONLY : FoundError
!
      implicit none
!
      PUBLIC  :: state_read
      PRIVATE :: state_read_nf90
      PRIVATE :: state_read_pio
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE state_read (ng, tile, model, IOtype,                   &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       Lout, rec,                                 &
     &                       nopen, ncid,                               &
     &                       pioFile,                                   &
     &                       ncname,                                    &
     &                       rmask, umask, vmask,                       &
     &                       s_t, s_u, s_v,                             &
     &                       s_zeta)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, IOtype
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: Lout, rec, nopen
      integer, intent(inout) :: ncid
!
      TYPE (File_desc_t), intent(inout) :: pioFile
!
      character (len=*), intent(in) :: ncname
!
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
      real(r8), intent(inout) :: s_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: s_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s_v(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s_zeta(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      logical :: Lreport = .FALSE.
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/state_read.F"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Read in requested state vector
!-----------------------------------------------------------------------
!
      SELECT CASE (IOtype)
        CASE (io_nf90)
          CALL state_read_nf90 (ng, tile, model, Lreport,               &
     &                          LBi, UBi, LBj, UBj, LBij, UBij,         &
     &                          Lout, rec,                              &
     &                          nopen, ncid, ncname,                    &
     &                          rmask, umask, vmask,                    &
     &                          s_t, s_u, s_v,                          &
     &                          s_zeta)
        CASE (io_pio)
          CALL state_read_pio (ng, tile, model, Lreport,                &
     &                         LBi, UBi, LBj, UBj, LBij, UBij,          &
     &                         Lout, rec,                               &
     &                         nopen, pioFile, ncname,                  &
     &                         rmask, umask, vmask,                     &
     &                         s_t, s_u, s_v,                           &
     &                         s_zeta)
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) IOtype
          exit_flag=2
      END SELECT
      IF (FoundError(exit_flag, NoError, 299, MyFile)) RETURN
!
  10  FORMAT (' STATE_READ - Illegal inpput file type, io_type = ',i0,  &
     &        /,14x,'Check KeyWord ''INP_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE state_read
!
!***********************************************************************
      SUBROUTINE state_read_nf90 (ng, tile, model, Lreport,             &
     &                            LBi, UBi, LBj, UBj, LBij, UBij,       &
     &                            Lout, rec,                            &
     &                            nopen, ncid, ncname,                  &
     &                            rmask, umask, vmask,                  &
     &                            s_t, s_u, s_v,                        &
     &                            s_zeta)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: Lreport
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: Lout, rec, nopen
      integer, intent(inout) :: ncid
!
      character (len=*), intent(in) :: ncname
!
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
      real(r8), intent(inout) :: s_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: s_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s_v(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s_zeta(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      integer :: Sstr, Send
      integer :: i, j, k
      integer :: ifield, itrc
      integer :: gtype, my_ncid, status, varid
      integer, dimension(4) :: Vsize
!
      real(r8) :: Fmin, Fmax
      real(dp) :: stime, scale
!
      character (len=15) :: Tstring
      character (len=22) :: t_code
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/state_read.F"//", state_read_nf90"
!
!-----------------------------------------------------------------------
!  Set lower and upper tile bounds and staggered variables bounds for
!  this horizontal domain partition.  Notice that if tile=-1, it will
!  set the values for the global grid.
!-----------------------------------------------------------------------
!
      integer :: Istr, IstrB, IstrP, IstrR, IstrT, IstrM, IstrU
      integer :: Iend, IendB, IendP, IendR, IendT
      integer :: Jstr, JstrB, JstrP, JstrR, JstrT, JstrM, JstrV
      integer :: Jend, JendB, JendP, JendR, JendT
      integer :: Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1
      integer :: Iendp1, Iendp2, Iendp2i, Iendp3
      integer :: Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1
      integer :: Jendp1, Jendp2, Jendp2i, Jendp3
!
      Istr   =BOUNDS(ng) % Istr   (tile)
      IstrB  =BOUNDS(ng) % IstrB  (tile)
      IstrM  =BOUNDS(ng) % IstrM  (tile)
      IstrP  =BOUNDS(ng) % IstrP  (tile)
      IstrR  =BOUNDS(ng) % IstrR  (tile)
      IstrT  =BOUNDS(ng) % IstrT  (tile)
      IstrU  =BOUNDS(ng) % IstrU  (tile)
      Iend   =BOUNDS(ng) % Iend   (tile)
      IendB  =BOUNDS(ng) % IendB  (tile)
      IendP  =BOUNDS(ng) % IendP  (tile)
      IendR  =BOUNDS(ng) % IendR  (tile)
      IendT  =BOUNDS(ng) % IendT  (tile)
      Jstr   =BOUNDS(ng) % Jstr   (tile)
      JstrB  =BOUNDS(ng) % JstrB  (tile)
      JstrM  =BOUNDS(ng) % JstrM  (tile)
      JstrP  =BOUNDS(ng) % JstrP  (tile)
      JstrR  =BOUNDS(ng) % JstrR  (tile)
      JstrT  =BOUNDS(ng) % JstrT  (tile)
      JstrV  =BOUNDS(ng) % JstrV  (tile)
      Jend   =BOUNDS(ng) % Jend   (tile)
      JendB  =BOUNDS(ng) % JendB  (tile)
      JendP  =BOUNDS(ng) % JendP  (tile)
      JendR  =BOUNDS(ng) % JendR  (tile)
      JendT  =BOUNDS(ng) % JendT  (tile)
!
      Istrm3 =BOUNDS(ng) % Istrm3 (tile)            ! Istr-3
      Istrm2 =BOUNDS(ng) % Istrm2 (tile)            ! Istr-2
      Istrm1 =BOUNDS(ng) % Istrm1 (tile)            ! Istr-1
      IstrUm2=BOUNDS(ng) % IstrUm2(tile)            ! IstrU-2
      IstrUm1=BOUNDS(ng) % IstrUm1(tile)            ! IstrU-1
      Iendp1 =BOUNDS(ng) % Iendp1 (tile)            ! Iend+1
      Iendp2 =BOUNDS(ng) % Iendp2 (tile)            ! Iend+2
      Iendp2i=BOUNDS(ng) % Iendp2i(tile)            ! Iend+2 interior
      Iendp3 =BOUNDS(ng) % Iendp3 (tile)            ! Iend+3
      Jstrm3 =BOUNDS(ng) % Jstrm3 (tile)            ! Jstr-3
      Jstrm2 =BOUNDS(ng) % Jstrm2 (tile)            ! Jstr-2
      Jstrm1 =BOUNDS(ng) % Jstrm1 (tile)            ! Jstr-1
      JstrVm2=BOUNDS(ng) % JstrVm2(tile)            ! JstrV-2
      JstrVm1=BOUNDS(ng) % JstrVm1(tile)            ! JstrV-1
      Jendp1 =BOUNDS(ng) % Jendp1 (tile)            ! Jend+1
      Jendp2 =BOUNDS(ng) % Jendp2 (tile)            ! Jend+2
      Jendp2i=BOUNDS(ng) % Jendp2i(tile)            ! Jend+2 interior
      Jendp3 =BOUNDS(ng) % Jendp3 (tile)            ! Jend+3
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Read in requested model state record. Load data into state array
!  index Lout.
!-----------------------------------------------------------------------
!
!  Turn on time wall clock.
!
      CALL wclock_on (ng, model, 80, 450, MyFile)
!
!  Determine file and variables ids.
!
      IF ((nopen.gt.0).or.(ncid.eq.-1)) THEN
        CALL netcdf_open (ng, model, ncname, 0, my_ncid)
        IF (FoundError(exit_flag, NoError, 457, MyFile)) RETURN
        ncid=my_ncid
      ELSE
        my_ncid=ncid
      END IF
!
      DO i=1,4
        Vsize(i)=0
      END DO
!
!  Read in time.
!
      CALL netcdf_get_time (ng, model, ncname, Vname(1,idtime),         &
     &                      Rclock%DateNumber, stime,                   &
     &                      my_ncid, (/rec/), (/1/))
      IF (FoundError(exit_flag, NoError, 483, MyFile)) RETURN
!
!  Report information.
!
      IF (Master) THEN
        CALL time_string (stime, t_code)
        Sstr=SCAN(CalledFrom,'/',BACK=.TRUE.)+1
        Send=LEN_TRIM(CalledFrom)
        WRITE (Tstring,'(f15.4)') stime*sec2day
        WRITE (stdout,10) 'Reading state fields,', t_code,              &
     &                    ng, Tstring, TRIM(ncname), rec, Lout,         &
     &                    CalledFrom(Sstr:Send)
      END IF
!
!  Read in free-surface.
!
      gtype=r2dvar
      scale=1.0_dp
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idFsur),        &
     &                       my_ncid, varid)
      IF (FoundError(exit_flag, NoError, 503, MyFile)) RETURN
      status=nf_fread2d(ng, model, ncname, my_ncid,                     &
     &                  Vname(1,idFsur), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, Fmin, Fmax,                              &
     &                  rmask,                                          &
     &                  s_zeta(:,:,Lout))
      IF (FoundError(status, nf90_noerr, 514, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(1,idFsur)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master.and.Lreport) THEN
          WRITE (stdout,40) TRIM(Vname(2,idFsur)), Fmin, Fmax
        END IF
      END IF
!
!  Read in 3D U-momentum component.
!
      gtype=u3dvar
      scale=1.0_dp
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idUvel),        &
     &                       my_ncid, varid)
      IF (FoundError(exit_flag, NoError, 755, MyFile)) RETURN
      status=nf_fread3d(ng, model, ncname, my_ncid,                     &
     &                  Vname(1,idUvel), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  scale, Fmin, Fmax,                              &
     &                  umask,                                          &
     &                  s_u(:,:,:,Lout))
      IF (FoundError(status, nf90_noerr, 766, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(1,idUvel)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master.and.Lreport) THEN
          WRITE (stdout,40) TRIM(Vname(2,idUvel)), Fmin, Fmax
        END IF
      END IF
!
!  Read in 3D V-momentum component.
!
      gtype=v3dvar
      scale=1.0_dp
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idVvel),        &
     &                       my_ncid, varid)
      IF (FoundError(exit_flag, NoError, 818, MyFile)) RETURN
      status=nf_fread3d(ng, model, ncname, my_ncid,                     &
     &                  Vname(1,idVvel), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  scale, Fmin, Fmax,                              &
     &                  vmask,                                          &
     &                  s_v(:,:,:,Lout))
      IF (FoundError(status, nf90_noerr, 829, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(1,idVvel)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master.and.Lreport) THEN
          WRITE (stdout,40) TRIM(Vname(2,idVvel)), Fmin, Fmax
        END IF
      END IF
!
!  Read in tracers.
!
      gtype=r3dvar
      scale=1.0_dp
      DO itrc=1,NT(ng)
        CALL netcdf_inq_varid (ng, model, ncname,                       &
     &                         Vname(1,idTvar(itrc)), my_ncid, varid)
        IF (FoundError(exit_flag, NoError, 882, MyFile)) RETURN
        status=nf_fread3d(ng, model, ncname, my_ncid,                   &
     &                    Vname(1,idTvar(itrc)), varid,                 &
     &                    rec, gtype, Vsize,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    scale, Fmin, Fmax,                            &
     &                    rmask,                                        &
     &                    s_t(:,:,:,Lout,itrc))
        IF (FoundError(status, nf90_noerr, 893, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idTvar(itrc))), rec,         &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        ELSE
          IF (Master.and.Lreport) THEN
            WRITE (stdout,40) TRIM(Vname(2,idTvar(itrc))), Fmin, Fmax
          END IF
        END IF
      END DO
!
!  Close current file.
!
      IF (nopen.gt.0) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
        IF (FoundError(exit_flag, NoError, 986, MyFile)) RETURN
      END IF
!
!  Turn off time wall clock.
!
      CALL wclock_off (ng, model, 80, 993, MyFile)
!
  10  FORMAT ('STATE_READ_NF90 - ',a,t75,a,                             &
     &        /,19x,'(Grid ',i2.2,', t = ',a,', File: ',a,              &
     &        ', Rec=',i4.4,', Index=',i1,')',                          &
     &        /,19x,'Called from ''',a,'''')
  20  FORMAT (' STATE_READ_NF90 - unable to open NetCDF file: ',a)
  30  FORMAT (' STATE_READ_NF90 - error while reading variable: ',a,2x, &
     &        'at time record = ',i3,/,14x,'in NetCDF file: ',a)
  40  FORMAT (16x,'- ',a,/,19x,'(Min = ',1p,e15.8,                      &
     &        ' Max = ',1p,e15.8,')')
!
      RETURN
      END SUBROUTINE state_read_nf90
!
!***********************************************************************
      SUBROUTINE state_read_pio (ng, tile, model, Lreport,              &
     &                           LBi, UBi, LBj, UBj, LBij, UBij,        &
     &                           Lout, rec,                             &
     &                           nopen, pioFile, ncname,                &
     &                           rmask, umask, vmask,                   &
     &                           s_t, s_u, s_v,                         &
     &                           s_zeta)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: Lreport
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: Lout, rec, nopen
!
      TYPE (File_desc_t), intent(inout) :: pioFile
!
      character (len=*), intent(in) :: ncname
!
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
      real(r8), intent(inout) :: s_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: s_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s_v(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s_zeta(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      integer :: Sstr, Send
      integer :: i, j, k
      integer :: ifield, itrc
      integer :: status
      integer, dimension(4) :: Vsize
!
      real(r8) :: Fmin, Fmax
      real(dp) :: stime, scale
!
      character (len=15) :: Tstring
      character (len=22) :: t_code
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/state_read.F"//", state_read_pio"
!
      TYPE (IO_desc_t), pointer :: my_ioDesc
      TYPE (File_desc_t) :: my_pioFile
      TYPE (My_VarDesc)  :: my_pioVar
!
!-----------------------------------------------------------------------
!  Set lower and upper tile bounds and staggered variables bounds for
!  this horizontal domain partition.  Notice that if tile=-1, it will
!  set the values for the global grid.
!-----------------------------------------------------------------------
!
      integer :: Istr, IstrB, IstrP, IstrR, IstrT, IstrM, IstrU
      integer :: Iend, IendB, IendP, IendR, IendT
      integer :: Jstr, JstrB, JstrP, JstrR, JstrT, JstrM, JstrV
      integer :: Jend, JendB, JendP, JendR, JendT
      integer :: Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1
      integer :: Iendp1, Iendp2, Iendp2i, Iendp3
      integer :: Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1
      integer :: Jendp1, Jendp2, Jendp2i, Jendp3
!
      Istr   =BOUNDS(ng) % Istr   (tile)
      IstrB  =BOUNDS(ng) % IstrB  (tile)
      IstrM  =BOUNDS(ng) % IstrM  (tile)
      IstrP  =BOUNDS(ng) % IstrP  (tile)
      IstrR  =BOUNDS(ng) % IstrR  (tile)
      IstrT  =BOUNDS(ng) % IstrT  (tile)
      IstrU  =BOUNDS(ng) % IstrU  (tile)
      Iend   =BOUNDS(ng) % Iend   (tile)
      IendB  =BOUNDS(ng) % IendB  (tile)
      IendP  =BOUNDS(ng) % IendP  (tile)
      IendR  =BOUNDS(ng) % IendR  (tile)
      IendT  =BOUNDS(ng) % IendT  (tile)
      Jstr   =BOUNDS(ng) % Jstr   (tile)
      JstrB  =BOUNDS(ng) % JstrB  (tile)
      JstrM  =BOUNDS(ng) % JstrM  (tile)
      JstrP  =BOUNDS(ng) % JstrP  (tile)
      JstrR  =BOUNDS(ng) % JstrR  (tile)
      JstrT  =BOUNDS(ng) % JstrT  (tile)
      JstrV  =BOUNDS(ng) % JstrV  (tile)
      Jend   =BOUNDS(ng) % Jend   (tile)
      JendB  =BOUNDS(ng) % JendB  (tile)
      JendP  =BOUNDS(ng) % JendP  (tile)
      JendR  =BOUNDS(ng) % JendR  (tile)
      JendT  =BOUNDS(ng) % JendT  (tile)
!
      Istrm3 =BOUNDS(ng) % Istrm3 (tile)            ! Istr-3
      Istrm2 =BOUNDS(ng) % Istrm2 (tile)            ! Istr-2
      Istrm1 =BOUNDS(ng) % Istrm1 (tile)            ! Istr-1
      IstrUm2=BOUNDS(ng) % IstrUm2(tile)            ! IstrU-2
      IstrUm1=BOUNDS(ng) % IstrUm1(tile)            ! IstrU-1
      Iendp1 =BOUNDS(ng) % Iendp1 (tile)            ! Iend+1
      Iendp2 =BOUNDS(ng) % Iendp2 (tile)            ! Iend+2
      Iendp2i=BOUNDS(ng) % Iendp2i(tile)            ! Iend+2 interior
      Iendp3 =BOUNDS(ng) % Iendp3 (tile)            ! Iend+3
      Jstrm3 =BOUNDS(ng) % Jstrm3 (tile)            ! Jstr-3
      Jstrm2 =BOUNDS(ng) % Jstrm2 (tile)            ! Jstr-2
      Jstrm1 =BOUNDS(ng) % Jstrm1 (tile)            ! Jstr-1
      JstrVm2=BOUNDS(ng) % JstrVm2(tile)            ! JstrV-2
      JstrVm1=BOUNDS(ng) % JstrVm1(tile)            ! JstrV-1
      Jendp1 =BOUNDS(ng) % Jendp1 (tile)            ! Jend+1
      Jendp2 =BOUNDS(ng) % Jendp2 (tile)            ! Jend+2
      Jendp2i=BOUNDS(ng) % Jendp2i(tile)            ! Jend+2 interior
      Jendp3 =BOUNDS(ng) % Jendp3 (tile)            ! Jend+3
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Read in requested model state record. Load data into state array
!  index Lout.
!-----------------------------------------------------------------------
!
!  Turn on time wall clock.
!
      CALL wclock_on (ng, model, 80, 1158, MyFile)
!
!  Determine file and variables ids.
!
      IF ((nopen.gt.0).or.(pioFile%fh.eq.-1)) THEN
        CALL pio_netcdf_open (ng, model, ncname, 0, my_pioFile)
        IF (FoundError(exit_flag, NoError, 1165, MyFile)) RETURN
        pioFile=my_pioFile
      ELSE
        my_pioFile=pioFile
      END IF
!
      DO i=1,4
        Vsize(i)=0
      END DO
!
!  Read in time.
!
      CALL pio_netcdf_get_time (ng, model, ncname, Vname(1,idtime),     &
     &                          Rclock%DateNumber, stime,               &
     &                          my_pioFile, (/rec/), (/1/))
      IF (FoundError(exit_flag, NoError, 1191, MyFile)) RETURN
!
!  Report information.
!
      IF (Master) THEN
        CALL time_string (stime, t_code)
        Sstr=SCAN(CalledFrom,'/',BACK=.TRUE.)+1
        Send=LEN_TRIM(CalledFrom)
        WRITE (Tstring,'(f15.4)') stime*sec2day
        WRITE (stdout,10) 'Reading state fields,', t_code,              &
     &                    ng, Tstring, TRIM(ncname), rec, Lout,         &
     &                    CalledFrom(Sstr:Send)
      END IF
!
!  Read in free-surface.
!
      scale=1.0_dp
      CALL pio_netcdf_inq_varid (ng, model, ncname, Vname(1,idFsur),    &
     &                           my_pioFile, my_pioVar%vd)
      IF (FoundError(exit_flag, NoError, 1210, MyFile)) RETURN
!
      my_pioVar%gtype=r2dvar
      IF (KIND(s_zeta).eq.8) THEN
        my_pioVar%dkind=PIO_double
        my_ioDesc => ioDesc_dp_r2dvar(ng)
      ELSE
        my_pioVar%dkind=PIO_real
        my_ioDesc => ioDesc_sp_r2dvar(ng)
      END IF
!
      status=nf_fread2d(ng, model, ncname, my_pioFile,                  &
     &                  Vname(1,idFsur), my_pioVar,                     &
     &                  rec, my_ioDesc, Vsize,                          &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, Fmin, Fmax,                              &
     &                  rmask,                                          &
     &                  s_zeta(:,:,Lout))
      IF (FoundError(status, PIO_noerr, 1230, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(1,idFsur)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master.and.Lreport) THEN
          WRITE (stdout,40) TRIM(Vname(2,idFsur)), Fmin, Fmax
        END IF
      END IF
!
!  Read in 3D U-momentum component.
!
      scale=1.0_dp
      CALL pio_netcdf_inq_varid (ng, model, ncname, Vname(1,idUvel),    &
     &                           my_pioFile, my_pioVar%vd)
      IF (FoundError(exit_flag, NoError, 1520, MyFile)) RETURN
!
      my_pioVar%gtype=u3dvar
      IF (KIND(s_u).eq.8) THEN
        my_pioVar%dkind=PIO_double
        my_ioDesc => ioDesc_dp_u3dvar(ng)
      ELSE
        my_pioVar%dkind=PIO_real
        my_ioDesc => ioDesc_sp_u3dvar(ng)
      END IF
!
      status=nf_fread3d(ng, model, ncname, my_pioFile,                  &
     &                  Vname(1,idUvel), my_pioVar,                     &
     &                  rec, my_ioDesc, Vsize,                          &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  scale, Fmin, Fmax,                              &
     &                  umask,                                          &
     &                  s_u(:,:,:,Lout))
      IF (FoundError(status, PIO_noerr, 1540, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(1,idUvel)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master.and.Lreport) THEN
          WRITE (stdout,40) TRIM(Vname(2,idUvel)), Fmin, Fmax
        END IF
      END IF
!
!  Read in 3D V-momentum component.
!
      scale=1.0_dp
      CALL pio_netcdf_inq_varid (ng, model, ncname, Vname(1,idVvel),    &
     &                           my_pioFile, my_pioVar%vd)
      IF (FoundError(exit_flag, NoError, 1597, MyFile)) RETURN
!
      my_pioVar%gtype=v3dvar
      IF (KIND(s_v).eq.8) THEN
        my_pioVar%dkind=PIO_double
        my_ioDesc => ioDesc_dp_v3dvar(ng)
      ELSE
        my_pioVar%dkind=PIO_real
        my_ioDesc => ioDesc_sp_v3dvar(ng)
      END IF
!
      status=nf_fread3d(ng, model, ncname, my_pioFile,                  &
     &                  Vname(1,idVvel), my_pioVar,                     &
     &                  rec, my_ioDesc, Vsize,                          &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  scale, Fmin, Fmax,                              &
     &                  vmask,                                          &
     &                  s_v(:,:,:,Lout))
      IF (FoundError(status, PIO_noerr, 1617, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(1,idVvel)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master.and.Lreport) THEN
          WRITE (stdout,40) TRIM(Vname(2,idVvel)), Fmin, Fmax
        END IF
      END IF
!
!  Read in tracers.
!
      scale=1.0_dp
      DO itrc=1,NT(ng)
        CALL pio_netcdf_inq_varid (ng, model, ncname,                   &
     &                             Vname(1,idTvar(itrc)),               &
     &                             my_pioFile, my_pioVar%vd)
        IF (FoundError(exit_flag, NoError, 1676, MyFile)) RETURN
!
        my_pioVar%gtype=r3dvar
        IF (KIND(s_t).eq.8) THEN
          my_pioVar%dkind=PIO_double
          my_ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          my_pioVar%dkind=PIO_real
          my_ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
!
        status=nf_fread3d(ng, model, ncname, my_pioFile,                &
     &                    Vname(1,idTvar(itrc)), my_pioVar,             &
     &                    rec, my_ioDesc, Vsize,                        &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    scale, Fmin, Fmax,                            &
     &                    rmask,                                        &
     &                    s_t(:,:,:,Lout,itrc))
        IF (FoundError(status, PIO_noerr, 1696, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idTvar(itrc))), rec,         &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        ELSE
          IF (Master.and.Lreport) THEN
            WRITE (stdout,40) TRIM(Vname(2,ifield)), Fmin, Fmax
          END IF
        END IF
      END DO
!
!  Close current file.
!
      IF (nopen.gt.0) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
        IF (FoundError(exit_flag, NoError, 1805, MyFile)) RETURN
      END IF
!
!  Turn off time wall clock.
!
      CALL wclock_off (ng, model, 80, 1812, MyFile)
!
  10  FORMAT ('STATE_READ_PIO  - ',a,t75,a,                             &
     &        /,19x,'(Grid ',i2.2,', t = ',a,', File: ',a,              &
     &        ', Rec=',i4.4,', Index=',i1,')',                          &
     &        /,19x,'Called from ''',a,'''')
  20  FORMAT (' STATE_READ_PIO  - unable to open NetCDF file: ',a)
  30  FORMAT (' STATE_READ_PIO  - error while reading variable: ',a,2x, &
     &        'at time record = ',i3,/,14x,'in NetCDF file: ',a)
  40  FORMAT (16x,'- ',a,/,19x,'(Min = ',1p,e15.8,                      &
     &        ' Max = ',1p,e15.8,')')
!
      RETURN
      END SUBROUTINE state_read_pio
      END MODULE state_read_mod
