      MODULE def_his_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module creates output HISTORY file using either the standard   !
!  NetCDF library or the Parallel-IO (PIO) library.  It defines its    !
!  dimensions, attributes, and variables.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      USE def_dim_mod,           ONLY : def_dim
      USE def_info_mod,          ONLY : def_info
      USE def_var_mod,           ONLY : def_var
      USE strings_mod,           ONLY : FoundError
      USE wrt_info_mod,        ONLY : wrt_info
!
      implicit none
!
      PUBLIC  :: def_his
      PRIVATE :: def_his_nf90
      PRIVATE :: def_his_pio
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE def_his (ng, ldef)
!***********************************************************************
!
!  Imported variable declarations.
!
      logical, intent(in) :: ldef
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/def_his.F"
!
!-----------------------------------------------------------------------
!  Create a new history file according to IO type.
!-----------------------------------------------------------------------
!
      SELECT CASE (HIS(ng)%IOtype)
        CASE (io_nf90)
          CALL def_his_nf90 (ng, iNLM, ldef)
        CASE (io_pio)
          CALL def_his_pio (ng, iNLM, ldef)
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) HIS(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, 112, MyFile)) RETURN
!
  10  FORMAT (' DEF_HIS - Illegal output file type, io_type = ',i0,     &
     &        /,11x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE def_his
!
!***********************************************************************
      SUBROUTINE def_his_nf90 (ng, model, ldef)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: ldef
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      logical :: got_var(NV)
!
      integer :: i, j, ifield, itrc, nvd3, nvd4, varid
      integer :: recdim, status
      integer :: DimIDs(nDimID)
      integer :: p2dgrd(3), t2dgrd(3), u2dgrd(3), v2dgrd(3)
      integer :: p3dgrd(4),  t3dgrd(4),  u3dgrd(4),  v3dgrd(4)
      integer :: p3dzgrd(4), t3dzgrd(4), u3dzgrd(4), v3dzgrd(4)
      integer :: w3dgrd(4)
!
      real(r8) :: Aval(6)
!
      character (len=256)    :: ncname
      character (len=MaxLen) :: Vinfo(Natt)
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/def_his.F"//", def_his_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, 175, MyFile)) RETURN
      ncname=HIS(ng)%name
!
      IF (Master) THEN
        IF (ldef) THEN
          WRITE (stdout,10) ng, TRIM(ncname)
        ELSE
          WRITE (stdout,20) ng, TRIM(ncname)
        END IF
      END IF
!
!=======================================================================
!  Create a new history file.
!=======================================================================
!
      DEFINE : IF (ldef) THEN
        CALL netcdf_create (ng, model, TRIM(ncname), HIS(ng)%ncid)
        IF (FoundError(exit_flag, NoError, 192, MyFile)) THEN
          IF (Master) WRITE (stdout,30) TRIM(ncname)
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Define file dimensions.
!-----------------------------------------------------------------------
!
        DimIDs=0
!
        status=def_dim(ng, model, HIS(ng)%ncid, ncname, 'xi_rho',       &
     &                 IOBOUNDS(ng)%xi_rho, DimIDs( 1))
        IF (FoundError(exit_flag, NoError, 205, MyFile)) RETURN
        status=def_dim(ng, model, HIS(ng)%ncid, ncname, 'xi_u',         &
     &                 IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (FoundError(exit_flag, NoError, 209, MyFile)) RETURN
        status=def_dim(ng, model, HIS(ng)%ncid, ncname, 'xi_v',         &
     &                 IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (FoundError(exit_flag, NoError, 213, MyFile)) RETURN
        status=def_dim(ng, model, HIS(ng)%ncid, ncname, 'xi_psi',       &
     &                 IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (FoundError(exit_flag, NoError, 217, MyFile)) RETURN
        status=def_dim(ng, model, HIS(ng)%ncid, ncname, 'eta_rho',      &
     &                 IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (FoundError(exit_flag, NoError, 221, MyFile)) RETURN
        status=def_dim(ng, model, HIS(ng)%ncid, ncname, 'eta_u',        &
     &                 IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (FoundError(exit_flag, NoError, 225, MyFile)) RETURN
        status=def_dim(ng, model, HIS(ng)%ncid, ncname, 'eta_v',        &
     &                 IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (FoundError(exit_flag, NoError, 229, MyFile)) RETURN
        status=def_dim(ng, model, HIS(ng)%ncid, ncname, 'eta_psi',      &
     &                 IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (FoundError(exit_flag, NoError, 233, MyFile)) RETURN
        status=def_dim(ng, model, HIS(ng)%ncid, ncname, 'N',            &
     &                 N(ng), DimIDs( 9))
        IF (FoundError(exit_flag, NoError, 276, MyFile)) RETURN
        status=def_dim(ng, model, HIS(ng)%ncid, ncname, 's_rho',        &
     &                 N(ng), DimIDs( 9))
        IF (FoundError(exit_flag, NoError, 280, MyFile)) RETURN
        status=def_dim(ng, model, HIS(ng)%ncid, ncname, 's_w',          &
     &                 N(ng)+1, DimIDs(10))
        IF (FoundError(exit_flag, NoError, 284, MyFile)) RETURN
        IF (Nslice.gt.0) THEN
          status=def_dim(ng, model, HIS(ng)%ncid, ncname, 'z_slice',    &
     &                   Nslice, DimIDs(34))
          IF (FoundError(exit_flag, NoError, 289, MyFile)) RETURN
        END IF
        status=def_dim(ng, model, HIS(ng)%ncid, ncname, 'tracer',       &
     &                 NT(ng), DimIDs(11))
        IF (FoundError(exit_flag, NoError, 294, MyFile)) RETURN
        status=def_dim(ng, model, HIS(ng)%ncid, ncname, 'boundary',     &
     &                 4, DimIDs(14))
        IF (FoundError(exit_flag, NoError, 343, MyFile)) RETURN
        status=def_dim(ng, model, HIS(ng)%ncid, ncname,                 &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 nf90_unlimited, DimIDs(12))
        IF (FoundError(exit_flag, NoError, 366, MyFile)) RETURN
        recdim=DimIDs(12)
!
!  Set number of dimensions for output variables.
!
        nvd3=3
        nvd4=4
!
!  Define dimension vectors for staggered tracer type variables.
!
        t2dgrd(1)=DimIDs( 1)
        t2dgrd(2)=DimIDs( 5)
        t2dgrd(3)=DimIDs(12)
        t3dgrd(1)=DimIDs( 1)
        t3dgrd(2)=DimIDs( 5)
        t3dgrd(3)=DimIDs( 9)
        t3dgrd(4)=DimIDs(12)
!
        t3dzgrd(1)=DimIDs( 1)
        t3dzgrd(2)=DimIDs( 5)
        t3dzgrd(3)=DimIDs(34)          ! selected constant depth slices
        t3dzgrd(4)=DimIDs(12)
!
!  Define dimension vectors for staggered type variables at PSI-points.
!
        p2dgrd(1)=DimIDs( 4)
        p2dgrd(2)=DimIDs( 8)
        p2dgrd(3)=DimIDs(12)
        p3dgrd(1)=DimIDs( 4)
        p3dgrd(2)=DimIDs( 8)
        p3dgrd(3)=DimIDs( 9)
        p3dgrd(4)=DimIDs(12)
!
        p3dzgrd(1)=DimIDs( 4)
        p3dzgrd(2)=DimIDs( 8)
        p3dzgrd(3)=DimIDs(34)          ! selected constant depth slices
        p3dzgrd(4)=DimIDs(12)
!
!  Define dimension vectors for staggered u-momentum type variables.
!
        u2dgrd(1)=DimIDs( 2)
        u2dgrd(2)=DimIDs( 6)
        u2dgrd(3)=DimIDs(12)
        u3dgrd(1)=DimIDs( 2)
        u3dgrd(2)=DimIDs( 6)
        u3dgrd(3)=DimIDs( 9)
        u3dgrd(4)=DimIDs(12)
!
        u3dzgrd(1)=DimIDs( 2)
        u3dzgrd(2)=DimIDs( 6)
        u3dzgrd(3)=DimIDs(34)          ! selected constant depth slices
        u3dzgrd(4)=DimIDs(12)
!
!  Define dimension vectors for staggered v-momentum type variables.
!
        v2dgrd(1)=DimIDs( 3)
        v2dgrd(2)=DimIDs( 7)
        v2dgrd(3)=DimIDs(12)
        v3dgrd(1)=DimIDs( 3)
        v3dgrd(2)=DimIDs( 7)
        v3dgrd(3)=DimIDs( 9)
        v3dgrd(4)=DimIDs(12)
!
        v3dzgrd(1)=DimIDs( 3)
        v3dzgrd(2)=DimIDs( 7)
        v3dzgrd(3)=DimIDs(34)          ! selected constant depth slices
        v3dzgrd(4)=DimIDs(12)
!
!  Define dimension vector for staggered w-momentum type variables.
!
        w3dgrd(1)=DimIDs( 1)
        w3dgrd(2)=DimIDs( 5)
        w3dgrd(3)=DimIDs(10)
        w3dgrd(4)=DimIDs(12)
!
!  Initialize unlimited time record dimension.
!
        HIS(ng)%Rindex=0
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
!  Define time-recordless information variables.
!-----------------------------------------------------------------------
!
        CALL def_info (ng, model, HIS(ng)%ncid, ncname, DimIDs)
        IF (FoundError(exit_flag, NoError, 548, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Define time-varying variables.
!-----------------------------------------------------------------------
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        Vinfo( 2)=Vname(2,idtime)
        WRITE (Vinfo( 3),'(a,a)') 'seconds since ', TRIM(Rclock%string)
        Vinfo( 4)=TRIM(Rclock%calendar)
        Vinfo(14)=Vname(4,idtime)
        Vinfo(21)=Vname(6,idtime)
        status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idtime),    &
     &                 NF_TOUT, 1, (/recdim/), Aval, Vinfo, ncname,     &
     &                 SetParAccess = .TRUE.)
        IF (FoundError(exit_flag, NoError, 565, MyFile)) RETURN
!
!  Define time-varying depth of RHO-points.
!
        IF (Hout(idpthR,ng)) THEN
          Vinfo( 1)=Vname(1,idpthR)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idpthR))
          Vinfo( 3)=Vname(3,idpthR)
          Vinfo(14)=Vname(4,idpthR)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idpthR)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idpthR,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idpthR),  &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, 675, MyFile)) RETURN
        END IF
!
!  Define time-varying depth of U-points.
!
        IF (Hout(idpthU,ng)) THEN
          Vinfo( 1)=Vname(1,idpthU)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idpthU))
          Vinfo( 3)=Vname(3,idpthU)
          Vinfo(14)=Vname(4,idpthU)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idpthU)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idpthU,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idpthU),  &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, 695, MyFile)) RETURN
        END IF
!
!  Define time-varying depth of V-points.
!
        IF (Hout(idpthV,ng)) THEN
          Vinfo( 1)=Vname(1,idpthV)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idpthV))
          Vinfo( 3)=Vname(3,idpthV)
          Vinfo(14)=Vname(4,idpthV)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idpthV)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idpthV,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idpthV),  &
     &                   NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, 715, MyFile)) RETURN
        END IF
!
!  Define time-varying depth of W-points.
!
        IF (Hout(idpthW,ng)) THEN
          Vinfo( 1)=Vname(1,idpthW)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idpthW))
          Vinfo( 3)=Vname(3,idpthW)
          Vinfo(14)=Vname(4,idpthW)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idpthW)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idpthW,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idpthW),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, 735, MyFile)) RETURN
        END IF
!
!  Define free-surface.
!
        IF (Hout(idFsur,ng)) THEN
          Vinfo( 1)=Vname(1,idFsur)
          Vinfo( 2)=Vname(2,idFsur)
          Vinfo( 3)=Vname(3,idFsur)
          Vinfo(14)=Vname(4,idFsur)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idFsur)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idFsur),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 760, MyFile)) RETURN
        END IF
!
!  Define 2D U-momentum component.
!
        IF (Hout(idUbar,ng)) THEN
          Vinfo( 1)=Vname(1,idUbar)
          Vinfo( 2)=Vname(2,idUbar)
          Vinfo( 3)=Vname(3,idUbar)
          Vinfo(14)=Vname(4,idUbar)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUbar)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbar,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idUbar),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 817, MyFile)) RETURN
        END IF
!
!  Define 2D V-momentum component.
!
        IF (Hout(idVbar,ng)) THEN
          Vinfo( 1)=Vname(1,idVbar)
          Vinfo( 2)=Vname(2,idVbar)
          Vinfo( 3)=Vname(3,idVbar)
          Vinfo(14)=Vname(4,idVbar)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVbar)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbar,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idVbar),  &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 925, MyFile)) RETURN
        END IF
!
!  Define 2D Eastward momentum component at RHO-points.
!
        IF (Hout(idu2dE,ng)) THEN
          Vinfo( 1)=Vname(1,idu2dE)
          Vinfo( 2)=Vname(2,idu2dE)
          Vinfo( 3)=Vname(3,idu2dE)
          Vinfo(14)=Vname(4,idu2dE)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idu2dE)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu2dE,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idu2dE),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1033, MyFile)) RETURN
        END IF
!
!  Define 2D Northward momentum component at RHO-points.
!
        IF (Hout(idv2dN,ng)) THEN
          Vinfo( 1)=Vname(1,idv2dN)
          Vinfo( 2)=Vname(2,idv2dN)
          Vinfo( 3)=Vname(3,idv2dN)
          Vinfo(14)=Vname(4,idv2dN)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idv2dN)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv2dN,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idv2dN),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1052, MyFile)) RETURN
        END IF
!
!  Define 3D U-momentum component.
!
        IF (Hout(idUvel,ng)) THEN
          Vinfo( 1)=Vname(1,idUvel)
          Vinfo( 2)=Vname(2,idUvel)
          Vinfo( 3)=Vname(3,idUvel)
          Vinfo(14)=Vname(4,idUvel)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUvel,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idUvel),  &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1073, MyFile)) RETURN
        END IF
!
!  Define 3D U-momentum component at specified constant depth slices.
!
        IF (Hout(idUzsl,ng).and.(Nslice.gt.0)) THEN
          Vinfo( 1)=Vname(1,idUzsl)
          Vinfo( 2)=Vname(2,idUzsl)
          Vinfo( 3)=Vname(3,idUzsl)
          Vinfo(14)=Vname(4,idUzsl)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUzsl)
          Vinfo(22)='coordinates, z_slice'
          Aval(5)=REAL(Iinfo(1,idUzsl,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idUzsl),  &
     &                   NF_FOUT, nvd4, u3dzgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1110, MyFile)) RETURN
        END IF
!
!  Define 3D V-momentum component.
!
        IF (Hout(idVvel,ng)) THEN
          Vinfo( 1)=Vname(1,idVvel)
          Vinfo( 2)=Vname(2,idVvel)
          Vinfo( 3)=Vname(3,idVvel)
          Vinfo(14)=Vname(4,idVvel)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvel,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idVvel),  &
     &                   NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1149, MyFile)) RETURN
        END IF
!
!  Define 3D V-momentum component at specified constant depth slices.
!
        IF (Hout(idVzsl,ng).and.(Nslice.gt.0)) THEN
          Vinfo( 1)=Vname(1,idVzsl)
          Vinfo( 2)=Vname(2,idVzsl)
          Vinfo( 3)=Vname(3,idVzsl)
          Vinfo(14)=Vname(4,idVzsl)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVzsl)
          Vinfo(22)='coordinates, z_slice'
          Aval(5)=REAL(Iinfo(1,idVzsl,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idVzsl),  &
     &                   NF_FOUT, nvd4, v3dzgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1186, MyFile)) RETURN
        END IF
!
!  Define 3D Eastward momentum at RHO-points, A-grid.
!
        IF (Hout(idu3dE,ng)) THEN
          Vinfo( 1)=Vname(1,idu3dE)
          Vinfo( 2)=Vname(2,idu3dE)
          Vinfo( 3)=Vname(3,idu3dE)
          Vinfo(14)=Vname(4,idu3dE)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idu3dE)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu3dE,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idu3dE),  &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1225, MyFile)) RETURN
        END IF
!
!  Define 3D Eastward momentum component at RHO-point for specified
!  constant depth slices.
!
        IF (Hout(idUzsE,ng).and.(Nslice.gt.0)) THEN
          Vinfo( 1)=Vname(1,idUzsE)
          Vinfo( 2)=Vname(2,idUzsE)
          Vinfo( 3)=Vname(3,idUzsE)
          Vinfo(14)=Vname(4,idUzsE)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUzsE)
          Vinfo(22)='coordinates, z_slice'
          Aval(5)=REAL(Iinfo(1,idUzsE,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idUzsE),  &
     &                   NF_FOUT, nvd4, t3dzgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1245, MyFile)) RETURN
        END IF
!
!  Define 3D Northward momentum at RHO-points, A-grid.
!
        IF (Hout(idv3dN,ng)) THEN
          Vinfo( 1)=Vname(1,idv3dN)
          Vinfo( 2)=Vname(2,idv3dN)
          Vinfo( 3)=Vname(3,idv3dN)
          Vinfo(14)=Vname(4,idv3dN)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idv3dN)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv3dN,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idv3dN),  &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1264, MyFile)) RETURN
        END IF
!
!  Define 3D Northward momentum component at RHO-points for specified
!  constant depth slices.
!
        IF (Hout(idVzsN,ng).and.(Nslice.gt.0)) THEN
          Vinfo( 1)=Vname(1,idVzsN)
          Vinfo( 2)=Vname(2,idVzsN)
          Vinfo( 3)=Vname(3,idVzsN)
          Vinfo(14)=Vname(4,idVzsN)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVzsN)
          Vinfo(22)='coordinates, zslice'
          Aval(5)=REAL(Iinfo(1,idVzsN,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idVzsN),  &
     &                   NF_FOUT, nvd4, t3dzgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1284, MyFile)) RETURN
        END IF
!
!  Define model potential vorticity at PSI-points.
!
        IF (Hout(id3dPV,ng)) THEN
          Vinfo( 1)=Vname(1,id3dPV)
          Vinfo( 2)=Vname(2,id3dPV)
          Vinfo( 3)=Vname(3,id3dPV)
          Vinfo(14)=Vname(4,id3dPV)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,id3dPV)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,id3dPV,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(id3dPV),  &
     &                   NF_FOUT, nvd4, p3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1303, MyFile)) RETURN
        END IF
!
!  Define model relative vorticity at PSI-points.
!
        IF (Hout(id3dRV,ng)) THEN
          Vinfo( 1)=Vname(1,id3dRV)
          Vinfo( 2)=Vname(2,id3dRV)
          Vinfo( 3)=Vname(3,id3dRV)
          Vinfo(14)=Vname(4,id3dRV)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,id3dRV)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,id3dRV,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(id3dRV),  &
     &                   NF_FOUT, nvd4, p3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1322, MyFile)) RETURN
        END IF
!
!  Define model potential vorticity at PSI-points for specified constant
!  depth slices.
!
        IF (Hout(idPVzs,ng)) THEN
          Vinfo( 1)=Vname(1,idPVzs)
          Vinfo( 2)=Vname(2,idPVzs)
          Vinfo( 3)=Vname(3,idPVzs)
          Vinfo(14)=Vname(4,idPVzs)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idPVzs)
          Vinfo(22)='coordinates, z_slice'
          Aval(5)=REAL(Iinfo(1,idPVzs,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idPVzs),  &
     &                   NF_FOUT, nvd4, p3dzgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1342, MyFile)) RETURN
        END IF
!
!  Define model relative vorticity at PSI-points for specified constant
!  depth slices.
!
        IF (Hout(idRVzs,ng)) THEN
          Vinfo( 1)=Vname(1,idRVzs)
          Vinfo( 2)=Vname(2,idRVzs)
          Vinfo( 3)=Vname(3,idRVzs)
          Vinfo(14)=Vname(4,idRVzs)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idRVzs)
          Vinfo(22)='coordinates, z_slice'
          Aval(5)=REAL(Iinfo(1,id3dRV,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idRVzs),  &
     &                   NF_FOUT, nvd4, p3dzgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1362, MyFile)) RETURN
        END IF
!
!  Define 3D momentum component in the Z-direction.
!
        IF (Hout(idWvel,ng)) THEN
          Vinfo( 1)=Vname(1,idWvel)
          Vinfo( 2)=Vname(2,idWvel)
          Vinfo( 3)=Vname(3,idWvel)
          Vinfo(14)=Vname(4,idWvel)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idWvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWvel,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idWvel),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1381, MyFile)) RETURN
        END IF
!
!  Define S-coordinate vertical "omega" momentum component.
!
        IF (Hout(idOvel,ng)) THEN
          Vinfo( 1)=Vname(1,idOvel)
          Vinfo( 2)=Vname(2,idOvel)
          Vinfo( 3)='meter second-1'
          Vinfo(14)=Vname(4,idOvel)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idOvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idOvel,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idOvel),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1400, MyFile)) RETURN
        END IF
!
!  Define S-coordinate implicit vertical "omega" momentum component.
!
        IF (Hout(idOvil,ng)) THEN
          Vinfo( 1)=Vname(1,idOvil)
          Vinfo( 2)=Vname(2,idOvil)
          Vinfo( 3)='meter second-1'
          Vinfo(14)=Vname(4,idOvil)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idOvil)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idOvil,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idOvil),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1421, MyFile)) RETURN
        END IF
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (Hout(idTvar(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTvar(itrc))
            Vinfo( 2)=Vname(2,idTvar(itrc))
            Vinfo( 3)=Vname(3,idTvar(itrc))
            Vinfo(14)=Vname(4,idTvar(itrc))
            Vinfo(16)=Vname(1,idtime)
            Vinfo(21)=Vname(6,idTvar(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idTvar(itrc),ng),r8)
            status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Tid(itrc),  &
     &                     NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1449, MyFile)) RETURN
          END IF
        END DO
!
!  Define tracer type variables at specified depth slices.
!
        DO itrc=1,NT(ng)
          IF (Hout(idzslT(itrc),ng).and.(Nslice.gt.0)) THEN
            Vinfo( 1)=Vname(1,idzslT(itrc))
            Vinfo( 2)=Vname(2,idzslT(itrc))
            Vinfo( 3)=Vname(3,idzslT(itrc))
            Vinfo(14)=Vname(4,idzslT(itrc))
            Vinfo(16)=Vname(1,idtime)
            Vinfo(21)=Vname(6,idzslT(itrc))
            Vinfo(22)='coordinates, z_slice'
            Aval(5)=REAL(r3dvar,r8)
            status=def_var(ng, model, HIS(ng)%ncid,                     &
     &                     HIS(ng)%Vid(idzslT(itrc)),                   &
     &                     NF_FOUT, nvd4, t3dzgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1478, MyFile)) RETURN
          END IF
        END DO
!
!  Define density anomaly.
!
        IF (Hout(idDano,ng)) THEN
          Vinfo( 1)=Vname(1,idDano)
          Vinfo( 2)=Vname(2,idDano)
          Vinfo( 3)=Vname(3,idDano)
          Vinfo(14)=Vname(4,idDano)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idDano)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idDano,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idDano),  &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1527, MyFile)) RETURN
        END IF
!
!  Define vertical viscosity coefficient.
!
        IF (Hout(idVvis,ng)) THEN
          Vinfo( 1)=Vname(1,idVvis)
          Vinfo( 2)=Vname(2,idVvis)
          Vinfo( 3)=Vname(3,idVvis)
          Vinfo(14)=Vname(4,idVvis)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVvis)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvis,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idVvis),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, 1614, MyFile)) RETURN
        END IF
!
!  Define vertical diffusion coefficient for potential temperature.
!
        IF (Hout(idTdif,ng)) THEN
          Vinfo( 1)=Vname(1,idTdif)
          Vinfo( 2)=Vname(2,idTdif)
          Vinfo( 3)=Vname(3,idTdif)
          Vinfo(14)=Vname(4,idTdif)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idTdif)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTdif,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idTdif),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, 1634, MyFile)) RETURN
        END IF
!
!  Define vertical diffusion coefficient for salinity.
!
        IF (Hout(idSdif,ng)) THEN
          Vinfo( 1)=Vname(1,idSdif)
          Vinfo( 2)=Vname(2,idSdif)
          Vinfo( 3)=Vname(3,idSdif)
          Vinfo(14)=Vname(4,idSdif)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idSdif)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSdif,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idSdif),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, 1656, MyFile)) RETURN
        END IF
!
!  Define turbulent kinetic energy.
!
        IF (Hout(idMtke,ng)) THEN
          Vinfo( 1)=Vname(1,idMtke)
          Vinfo( 2)=Vname(2,idMtke)
          Vinfo( 3)=Vname(3,idMtke)
          Vinfo(14)=Vname(4,idMtke)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idMtke)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idMtke,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idMtke),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, 1678, MyFile)) RETURN
        END IF
!
!  Define turbulent kinetic energy time length scale.
!
        IF (Hout(idMtls,ng)) THEN
          Vinfo( 1)=Vname(1,idMtls)
          Vinfo( 2)=Vname(2,idMtls)
          Vinfo( 3)=Vname(3,idMtls)
          Vinfo(14)=Vname(4,idMtls)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idMtls)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idMtls,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idMtls),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, 1716, MyFile)) RETURN
        END IF
!
!  Define surface air pressure.
!
        IF (Hout(idPair,ng)) THEN
          Vinfo( 1)=Vname(1,idPair)
          Vinfo( 2)=Vname(2,idPair)
          Vinfo( 3)=Vname(3,idPair)
          Vinfo(14)=Vname(4,idPair)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idPair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idPair,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idPair),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1772, MyFile)) RETURN
        END IF
!
!  Define surface winds.
!
        IF (Hout(idUair,ng)) THEN
          Vinfo( 1)=Vname(1,idUair)
          Vinfo( 2)=Vname(2,idUair)
          Vinfo( 3)=Vname(3,idUair)
          Vinfo(14)=Vname(4,idUair)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUair,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idUair),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1793, MyFile)) RETURN
        END IF
!
        IF (Hout(idVair,ng)) THEN
          Vinfo( 1)=Vname(1,idVair)
          Vinfo( 2)=Vname(2,idVair)
          Vinfo( 3)=Vname(3,idVair)
          Vinfo(14)=Vname(4,idVair)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVair,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idVair),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1810, MyFile)) RETURN
        END IF
!
!  Define Eastward/Northward surface winds at RHO-points.
!
        IF (Hout(idUaiE,ng)) THEN
          Vinfo( 1)=Vname(1,idUaiE)
          Vinfo( 2)=Vname(2,idUaiE)
          Vinfo( 3)=Vname(3,idUaiE)
          Vinfo(14)=Vname(4,idUaiE)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUaiE)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUaiE,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idUaiE),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1829, MyFile)) RETURN
        END IF
!
        IF (Hout(idVaiN,ng)) THEN
          Vinfo( 1)=Vname(1,idVaiN)
          Vinfo( 2)=Vname(2,idVaiN)
          Vinfo( 3)=Vname(3,idVaiN)
          Vinfo(14)=Vname(4,idVaiN)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVaiN)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVaiN,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idVaiN),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1846, MyFile)) RETURN
        END IF
!
!  Define surface active tracer fluxes.
!
        DO itrc=1,NAT
          IF (Hout(idTsur(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTsur(itrc))
            Vinfo( 2)=Vname(2,idTsur(itrc))
            Vinfo( 3)=Vname(3,idTsur(itrc))
            IF (itrc.eq.itemp) THEN
              Vinfo(11)='upward flux, cooling'
              Vinfo(12)='downward flux, heating'
            ELSE IF (itrc.eq.isalt) THEN
              Vinfo(11)='upward flux, freshening (net precipitation)'
              Vinfo(12)='downward flux, salting (net evaporation)'
            END IF
            Vinfo(14)=Vname(4,idTsur(itrc))
            Vinfo(16)=Vname(1,idtime)
            Vinfo(21)=Vname(6,idTsur(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idTsur(itrc),ng),r8)
            status=def_var(ng, model, HIS(ng)%ncid,                     &
     &                     HIS(ng)%Vid(idTsur(itrc)), NF_FOUT,          &
     &                     nvd3, t2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 1875, MyFile)) RETURN
          END IF
        END DO
!
!  Define latent heat flux.
!
        IF (Hout(idLhea,ng)) THEN
          Vinfo( 1)=Vname(1,idLhea)
          Vinfo( 2)=Vname(2,idLhea)
          Vinfo( 3)=Vname(3,idLhea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLhea)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idLhea)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idLhea,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idLhea),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1899, MyFile)) RETURN
        END IF
!
!  Define sensible heat flux.
!
        IF (Hout(idShea,ng)) THEN
          Vinfo( 1)=Vname(1,idShea)
          Vinfo( 2)=Vname(2,idShea)
          Vinfo( 3)=Vname(3,idShea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idShea)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idShea)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idShea,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idShea),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1920, MyFile)) RETURN
        END IF
!
!  Define net longwave radiation flux.
!
        IF (Hout(idLrad,ng)) THEN
          Vinfo( 1)=Vname(1,idLrad)
          Vinfo( 2)=Vname(2,idLrad)
          Vinfo( 3)=Vname(3,idLrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLrad)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idLrad)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idLrad,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idLrad),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1941, MyFile)) RETURN
        END IF
!
!  Define atmospheric air temperature.
!
        IF (Hout(idTair,ng)) THEN
          Vinfo( 1)=Vname(1,idTair)
          Vinfo( 2)=Vname(2,idTair)
          Vinfo( 3)=Vname(3,idTair)
          Vinfo(14)=Vname(4,idTair)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idTair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTair,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idTair),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1963, MyFile)) RETURN
        END IF
!
!  Define evaporation rate.
!
        IF (Hout(idevap,ng)) THEN
          Vinfo( 1)=Vname(1,idevap)
          Vinfo( 2)=Vname(2,idevap)
          Vinfo( 3)=Vname(3,idevap)
          Vinfo(11)='downward flux, freshening (condensation)'
          Vinfo(12)='upward flux, salting (evaporation)'
          Vinfo(14)=Vname(4,idevap)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idevap)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idevap,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idevap),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1986, MyFile)) RETURN
        END IF
!
!  Define precipitation rate.
!
        IF (Hout(idrain,ng)) THEN
          Vinfo( 1)=Vname(1,idrain)
          Vinfo( 2)=Vname(2,idrain)
          Vinfo( 3)=Vname(3,idrain)
          Vinfo(11)='upward flux, salting (NOT POSSIBLE)'
          Vinfo(12)='downward flux, freshening (precipitation)'
          Vinfo(14)=Vname(4,idrain)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idrain)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idrain,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idrain),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 2007, MyFile)) RETURN
        END IF
!
!  Define E-P flux.
!
        IF (Hout(idEmPf,ng)) THEN
          Vinfo( 1)=Vname(1,idEmPf)
          Vinfo( 2)=Vname(2,idEmPf)
          Vinfo( 3)=Vname(3,idEmPf)
          Vinfo(11)='upward flux, freshening (net precipitation)'
          Vinfo(12)='downward flux, salting (net evaporation)'
          Vinfo(14)=Vname(4,idEmPf)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idEmPf)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idEmPf,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idEmPf),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 2030, MyFile)) RETURN
        END IF
!
!  Define net shortwave radiation flux.
!
        IF (Hout(idSrad,ng)) THEN
          Vinfo( 1)=Vname(1,idSrad)
          Vinfo( 2)=Vname(2,idSrad)
          Vinfo( 3)=Vname(3,idSrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idSrad)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idSrad)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSrad,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idSrad),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 2053, MyFile)) RETURN
        END IF
!
!  Define surface U-momentum stress.
!
        IF (Hout(idUsms,ng)) THEN
          Vinfo( 1)=Vname(1,idUsms)
          Vinfo( 2)=Vname(2,idUsms)
          Vinfo( 3)=Vname(3,idUsms)
          Vinfo(14)=Vname(4,idUsms)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUsms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUsms,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idUsms),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 2074, MyFile)) RETURN
        END IF
!
!  Define surface V-momentum stress.
!
        IF (Hout(idVsms,ng)) THEN
          Vinfo( 1)=Vname(1,idVsms)
          Vinfo( 2)=Vname(2,idVsms)
          Vinfo( 3)=Vname(3,idVsms)
          Vinfo(14)=Vname(4,idVsms)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVsms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVsms,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idVsms),  &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 2093, MyFile)) RETURN
        END IF
!
!  Define bottom U-momentum stress.
!
        IF (Hout(idUbms,ng)) THEN
          Vinfo( 1)=Vname(1,idUbms)
          Vinfo( 2)=Vname(2,idUbms)
          Vinfo( 3)=Vname(3,idUbms)
          Vinfo(14)=Vname(4,idUbms)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUbms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbms,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idUbms),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 2112, MyFile)) RETURN
        END IF
!
!  Define bottom V-momentum stress.
!
        IF (Hout(idVbms,ng)) THEN
          Vinfo( 1)=Vname(1,idVbms)
          Vinfo( 2)=Vname(2,idVbms)
          Vinfo( 3)=Vname(3,idVbms)
          Vinfo(14)=Vname(4,idVbms)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVbms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbms,ng),r8)
          status=def_var(ng, model, HIS(ng)%ncid, HIS(ng)%Vid(idVbms),  &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 2131, MyFile)) RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, model, ncname, HIS(ng)%ncid)
        IF (FoundError(exit_flag, NoError, 2196, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, model, HIS(ng)%ncid, ncname)
        IF (FoundError(exit_flag, NoError, 2203, MyFile)) RETURN
      END IF DEFINE
!
!=======================================================================
!  Open an existing history file, check its contents, and prepare for
!  appending data.
!=======================================================================
!
      QUERY : IF (.not.ldef) THEN
        ncname=HIS(ng)%name
!
!  Open history file for read/write.
!
        CALL netcdf_open (ng, model, ncname, 1, HIS(ng)%ncid)
        IF (FoundError(exit_flag, NoError, 2218, MyFile)) THEN
          WRITE (stdout,60) TRIM(ncname)
          RETURN
        END IF
!
!  Inquire about the dimensions and check for consistency.
!
        CALL netcdf_check_dim (ng, model, ncname,                       &
     &                         ncid = HIS(ng)%ncid)
        IF (FoundError(exit_flag, NoError, 2227, MyFile)) RETURN
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, model, ncname,                         &
     &                       ncid = HIS(ng)%ncid)
        IF (FoundError(exit_flag, NoError, 2233, MyFile)) RETURN
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  history variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            HIS(ng)%Vid(idtime)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idpthR))) THEN
            got_var(idpthR)=.TRUE.
            HIS(ng)%Vid(idpthR)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idpthU))) THEN
            got_var(idpthU)=.TRUE.
            HIS(ng)%Vid(idpthU)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idpthV))) THEN
            got_var(idpthV)=.TRUE.
            HIS(ng)%Vid(idpthV)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idpthW))) THEN
            got_var(idpthW)=.TRUE.
            HIS(ng)%Vid(idpthW)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            HIS(ng)%Vid(idFsur)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            HIS(ng)%Vid(idUbar)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            HIS(ng)%Vid(idVbar)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu2dE))) THEN
            got_var(idu2dE)=.TRUE.
            HIS(ng)%Vid(idu2dE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv2dN))) THEN
            got_var(idv2dN)=.TRUE.
            HIS(ng)%Vid(idv2dN)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            HIS(ng)%Vid(idUvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            HIS(ng)%Vid(idVvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUzsl))) THEN
            got_var(idUzsl)=.TRUE.
            HIS(ng)%Vid(idUzsl)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVzsl))) THEN
            got_var(idVzsl)=.TRUE.
            HIS(ng)%Vid(idVzsl)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu3dE))) THEN
            got_var(idu3dE)=.TRUE.
            HIS(ng)%Vid(idu3dE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv3dN))) THEN
            got_var(idv3dN)=.TRUE.
            HIS(ng)%Vid(idv3dN)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUzsE))) THEN
            got_var(idUzsE)=.TRUE.
            HIS(ng)%Vid(idUzsE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVzsN))) THEN
            got_var(idVzsN)=.TRUE.
            HIS(ng)%Vid(idVzsN)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,id3dPV))) THEN
            got_var(id3dPV)=.TRUE.
            HIS(ng)%Vid(id3dPV)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,id3dRV))) THEN
            got_var(id3dRV)=.TRUE.
            HIS(ng)%Vid(id3dRV)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPVzs))) THEN
            got_var(idPVzs)=.TRUE.
            HIS(ng)%Vid(idRVzs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRVzs))) THEN
            got_var(idRVzs)=.TRUE.
            HIS(ng)%Vid(idRVzs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWvel))) THEN
            got_var(idWvel)=.TRUE.
            HIS(ng)%Vid(idWvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOvel))) THEN
            got_var(idOvel)=.TRUE.
            HIS(ng)%Vid(idOvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOvil))) THEN
            got_var(idOvil)=.TRUE.
            HIS(ng)%Vid(idOvil)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            HIS(ng)%Vid(idDano)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvis))) THEN
            got_var(idVvis)=.TRUE.
            HIS(ng)%Vid(idVvis)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTdif))) THEN
            got_var(idTdif)=.TRUE.
            HIS(ng)%Vid(idTdif)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSdif))) THEN
            got_var(idSdif)=.TRUE.
            HIS(ng)%Vid(idSdif)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idMtke))) THEN
            got_var(idMtke)=.TRUE.
            HIS(ng)%Vid(idMtke)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idMtls))) THEN
            got_var(idMtls)=.TRUE.
            HIS(ng)%Vid(idMtls)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPair))) THEN
            got_var(idPair)=.TRUE.
            HIS(ng)%Vid(idPair)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUair))) THEN
            got_var(idUair)=.TRUE.
            HIS(ng)%Vid(idUair)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVair))) THEN
            got_var(idVair)=.TRUE.
            HIS(ng)%Vid(idVair)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUaiE))) THEN
            got_var(idUaiE)=.TRUE.
            HIS(ng)%Vid(idUaiE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVaiN))) THEN
            got_var(idVaiN)=.TRUE.
            HIS(ng)%Vid(idVaiN)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLhea))) THEN
            got_var(idLhea)=.TRUE.
            HIS(ng)%Vid(idLhea)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idShea))) THEN
            got_var(idShea)=.TRUE.
            HIS(ng)%Vid(idShea)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLrad))) THEN
            got_var(idLrad)=.TRUE.
            HIS(ng)%Vid(idLrad)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTair))) THEN
            got_var(idTair)=.TRUE.
            HIS(ng)%Vid(idTair)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idevap))) THEN
            got_var(idevap)=.TRUE.
            HIS(ng)%Vid(idevap)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idrain))) THEN
            got_var(idrain)=.TRUE.
            HIS(ng)%Vid(idrain)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idEmPf))) THEN
            got_var(idEmPf)=.TRUE.
            HIS(ng)%Vid(idEmPf)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSrad))) THEN
            got_var(idSrad)=.TRUE.
            HIS(ng)%Vid(idSrad)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            HIS(ng)%Vid(idUsms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            HIS(ng)%Vid(idVsms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbms))) THEN
            got_var(idUbms)=.TRUE.
            HIS(ng)%Vid(idUbms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbms))) THEN
            got_var(idVbms)=.TRUE.
            HIS(ng)%Vid(idVbms)=var_id(i)
          END IF
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              HIS(ng)%Tid(itrc)=var_id(i)
            END IF
          END DO
          DO itrc=1,NAT
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idzslT(itrc)))) THEN
              got_var(idzslT(itrc))=.TRUE.
              HIS(ng)%Vid(idzslT(itrc))=var_id(i)
            END IF
          END DO
          DO itrc=1,NAT
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTsur(itrc)))) THEN
              got_var(idTsur(itrc))=.TRUE.
              HIS(ng)%Vid(idTsur(itrc))=var_id(i)
            END IF
          END DO
        END DO
!
!  Check if history variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idpthR).and.Hout(idpthR,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idpthR)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idpthU).and.Hout(idpthU,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idpthU)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idpthV).and.Hout(idpthV,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idpthV)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idpthW).and.Hout(idpthW,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idpthW)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsur).and.Hout(idFsur,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar).and.Hout(idUbar,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar).and.Hout(idVbar,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idu2dE).and.Hout(idu2dE,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idu2dE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv2dN).and.Hout(idv2dN,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idv2dN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUvel).and.Hout(idUvel,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel).and.Hout(idVvel,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUzsl).and.Hout(idUzsl,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUzsl)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVzsl).and.Hout(idVzsl,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVzsl)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idu3dE).and.Hout(idu3dE,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idu3dE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv3dN).and.Hout(idv3dN,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idv3dN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUzsE).and.Hout(idUzsE,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUzsE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVzsN).and.Hout(idVzsN,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVzsN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(id3dPV).and.Hout(id3dPV,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,id3dPV)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(id3dRV).and.Hout(id3dRV,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,id3dRV)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idPVzs).and.Hout(idPVzs,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idPVzs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idRVzs).and.Hout(idRVzs,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idRVzs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWvel).and.Hout(idWvel,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idWvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idOvel).and.Hout(idOvel,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idOvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idOvil).and.Hout(idOvil,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idOvil)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idDano).and.Hout(idDano,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idDano)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvis).and.Hout(idVvis,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVvis)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTdif).and.Hout(idTdif,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSdif).and.Hout(idSdif,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idMtke).and.Hout(idMtke,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idMtke)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idMtls).and.Hout(idMtls,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idMtls)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idPair).and.Hout(idPair,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idPair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUair).and.Hout(idUair,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVair).and.Hout(idVair,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUaiE).and.Hout(idUaiE,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUaiE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVaiN).and.Hout(idVaiN,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVaiN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLhea).and.Hout(idLhea,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idLhea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idShea).and.Hout(idShea,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idShea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLrad).and.Hout(idLrad,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idLrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTair).and.Hout(idTair,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idevap).and.Hout(idevap,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idevap)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idrain).and.Hout(idrain,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idrain)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idEmPf).and.Hout(idEmPf,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idEmPf)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSrad).and.Hout(idSrad,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUsms).and.Hout(idUsms,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms).and.Hout(idVsms,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbms).and.Hout(idUbms,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbms).and.Hout(idVbms,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc)).and.Hout(idTvar(itrc),ng)) THEN
            IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO itrc=1,NAT
          IF (.not.got_var(idzslT(itrc)).and.Hout(idzslT(itrc),ng)) THEN
            IF (Master) WRITE (stdout,70) TRIM(Vname(1,idzslT(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO itrc=1,NAT
          IF (.not.got_var(idTsur(itrc)).and.Hout(idTsur(itrc),ng)) THEN
            IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTsur(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
!
!  Set unlimited time record dimension to the appropriate value.
!
        IF (ndefHIS(ng).gt.0) THEN
          HIS(ng)%Rindex=((ntstart(ng)-1)-                              &
     &                    ndefHIS(ng)*((ntstart(ng)-1)/ndefHIS(ng)))/   &
     &                   nHIS(ng)
        ELSE
          HIS(ng)%Rindex=(ntstart(ng)-1)/nHIS(ng)
        END IF
        HIS(ng)%Rindex=MIN(HIS(ng)%Rindex,rec_size)
      END IF QUERY
!
  10  FORMAT (2x,'DEF_HIS_NF90     - creating history file,',t56,       &
     &        'Grid ',i2.2,': ',a)
  20  FORMAT (2x,'DEF_HIS_NF90     - inquiring history file,',t56,      &
     &        'Grid ',i2.2,': ',a)
  30  FORMAT (/,' DEF_HIS_NF90 - unable to create history NetCDF',      &
     &        ' file: ',a)
  40  FORMAT ('time dependent',1x,a)
  50  FORMAT (1pe11.4,1x,'millimeter')
  60  FORMAT (/,' DEF_HIS_NF90 - unable to open history NetCDF',        &
     &        ' file: ',a)
  70  FORMAT (/,' DEF_HIS_NF90 - unable to find variable: ',a,2x,       &
     &        ' in history NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE def_his_nf90
!
!***********************************************************************
      SUBROUTINE def_his_pio (ng, model, ldef)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: ldef
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      logical :: got_var(NV)
!
      integer :: i, j, ifield, itrc, nvd3, nvd4
      integer :: recdim, status
      integer :: DimIDs(nDimID)
      integer :: p2dgrd(3), t2dgrd(3), u2dgrd(3), v2dgrd(3)
      integer :: p3dgrd(4),  t3dgrd(4),  u3dgrd(4),  v3dgrd(4)
      integer :: p3dzgrd(4), t3dzgrd(4), u3dzgrd(4), v3dzgrd(4)
      integer :: w3dgrd(4)
!
      real(r8) :: Aval(6)
!
      character (len=256)    :: ncname
      character (len=MaxLen) :: Vinfo(Natt)
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/def_his.F"//", def_his_pio"
!
      TYPE (Var_desc_t) :: varDesc
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, 3185, MyFile)) RETURN
      ncname=HIS(ng)%name
!
      IF (Master) THEN
        IF (ldef) THEN
          WRITE (stdout,10) ng, TRIM(ncname)
        ELSE
          WRITE (stdout,20) ng, TRIM(ncname)
        END IF
      END IF
!
!=======================================================================
!  Create a new history file.
!=======================================================================
!
      DEFINE : IF (ldef) THEN
        CALL pio_netcdf_create (ng, model, TRIM(ncname), HIS(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, 3202, MyFile)) THEN
          IF (Master) WRITE (stdout,30) TRIM(ncname)
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Define file dimensions.
!-----------------------------------------------------------------------
!
        DimIDs=0
!
        status=def_dim(ng, model, HIS(ng)%pioFile, ncname, 'xi_rho',    &
     &                 IOBOUNDS(ng)%xi_rho, DimIDs( 1))
        IF (FoundError(exit_flag, NoError, 3215, MyFile)) RETURN
        status=def_dim(ng, model, HIS(ng)%pioFile, ncname, 'xi_u',      &
     &                 IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (FoundError(exit_flag, NoError, 3219, MyFile)) RETURN
        status=def_dim(ng, model, HIS(ng)%pioFile, ncname, 'xi_v',      &
     &                 IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (FoundError(exit_flag, NoError, 3223, MyFile)) RETURN
        status=def_dim(ng, model, HIS(ng)%pioFile, ncname, 'xi_psi',    &
     &                 IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (FoundError(exit_flag, NoError, 3227, MyFile)) RETURN
        status=def_dim(ng, model, HIS(ng)%pioFile, ncname, 'eta_rho',   &
     &                 IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (FoundError(exit_flag, NoError, 3231, MyFile)) RETURN
        status=def_dim(ng, model, HIS(ng)%pioFile, ncname, 'eta_u',     &
     &                 IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (FoundError(exit_flag, NoError, 3235, MyFile)) RETURN
        status=def_dim(ng, model, HIS(ng)%pioFile, ncname, 'eta_v',     &
     &                 IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (FoundError(exit_flag, NoError, 3239, MyFile)) RETURN
        status=def_dim(ng, model, HIS(ng)%pioFile, ncname, 'eta_psi',   &
     &                 IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (FoundError(exit_flag, NoError, 3243, MyFile)) RETURN
        status=def_dim(ng, model, HIS(ng)%pioFile, ncname, 'N',         &
     &                 N(ng), DimIDs( 9))
        IF (FoundError(exit_flag, NoError, 3286, MyFile)) RETURN
        status=def_dim(ng, model, HIS(ng)%pioFile, ncname, 's_rho',     &
     &                 N(ng), DimIDs( 9))
        IF (FoundError(exit_flag, NoError, 3290, MyFile)) RETURN
        status=def_dim(ng, model, HIS(ng)%pioFile, ncname, 's_w',       &
     &                 N(ng)+1, DimIDs(10))
        IF (FoundError(exit_flag, NoError, 3294, MyFile)) RETURN
        IF (Nslice.gt.0) THEN
          status=def_dim(ng, model, HIS(ng)%pioFile, ncname, 'z_slice', &
     &                   Nslice, DimIDs(34))
          IF (FoundError(exit_flag, NoError, 3299, MyFile)) RETURN
        END IF
        status=def_dim(ng, model, HIS(ng)%pioFile, ncname, 'tracer',    &
     &                 NT(ng), DimIDs(11))
        IF (FoundError(exit_flag, NoError, 3304, MyFile)) RETURN
        status=def_dim(ng, model, HIS(ng)%pioFile, ncname, 'boundary',  &
     &                 4, DimIDs(14))
        IF (FoundError(exit_flag, NoError, 3353, MyFile)) RETURN
        status=def_dim(ng, model, HIS(ng)%pioFile, ncname,              &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 PIO_unlimited, DimIDs(12))
        IF (FoundError(exit_flag, NoError, 3376, MyFile)) RETURN
        recdim=DimIDs(12)
!
!  Set number of dimensions for output variables.
!
        nvd3=3
        nvd4=4
!
!  Define dimension vectors for staggered tracer type variables.
!
        t2dgrd(1)=DimIDs( 1)
        t2dgrd(2)=DimIDs( 5)
        t2dgrd(3)=DimIDs(12)
        t3dgrd(1)=DimIDs( 1)
        t3dgrd(2)=DimIDs( 5)
        t3dgrd(3)=DimIDs( 9)
        t3dgrd(4)=DimIDs(12)
!
        t3dzgrd(1)=DimIDs( 1)
        t3dzgrd(2)=DimIDs( 5)
        t3dzgrd(3)=DimIDs(34)          ! selected constant depth slices
        t3dzgrd(4)=DimIDs(12)
!
!  Define dimension vectors for staggered type variables at PSI-points.
!
        p2dgrd(1)=DimIDs( 4)
        p2dgrd(2)=DimIDs( 8)
        p2dgrd(3)=DimIDs(12)
        p3dgrd(1)=DimIDs( 4)
        p3dgrd(2)=DimIDs( 8)
        p3dgrd(3)=DimIDs( 9)
        p3dgrd(4)=DimIDs(12)
!
        p3dzgrd(1)=DimIDs( 4)
        p3dzgrd(2)=DimIDs( 8)
        p3dzgrd(3)=DimIDs(34)          ! selected constant depth slices
        p3dzgrd(4)=DimIDs(12)
!
!  Define dimension vectors for staggered u-momentum type variables.
!
        u2dgrd(1)=DimIDs( 2)
        u2dgrd(2)=DimIDs( 6)
        u2dgrd(3)=DimIDs(12)
        u3dgrd(1)=DimIDs( 2)
        u3dgrd(2)=DimIDs( 6)
        u3dgrd(3)=DimIDs( 9)
        u3dgrd(4)=DimIDs(12)
!
        u3dzgrd(1)=DimIDs( 1)
        u3dzgrd(2)=DimIDs( 5)
        u3dzgrd(3)=DimIDs(34)          ! selected constant depth slices
        u3dzgrd(4)=DimIDs(12)
!
!  Define dimension vectors for staggered v-momentum type variables.
!
        v2dgrd(1)=DimIDs( 3)
        v2dgrd(2)=DimIDs( 7)
        v2dgrd(3)=DimIDs(12)
        v3dgrd(1)=DimIDs( 3)
        v3dgrd(2)=DimIDs( 7)
        v3dgrd(3)=DimIDs( 9)
        v3dgrd(4)=DimIDs(12)
!
        v3dzgrd(1)=DimIDs( 1)
        v3dzgrd(2)=DimIDs( 5)
        v3dzgrd(3)=DimIDs(34)          ! selected constant depth slices
        v3dzgrd(4)=DimIDs(12)
!
!  Define dimension vector for staggered w-momentum type variables.
!
        w3dgrd(1)=DimIDs( 1)
        w3dgrd(2)=DimIDs( 5)
        w3dgrd(3)=DimIDs(10)
        w3dgrd(4)=DimIDs(12)
!
!  Initialize unlimited time record dimension.
!
        HIS(ng)%Rindex=0
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
!  Define time-recordless information variables.
!-----------------------------------------------------------------------
!
        CALL def_info (ng, model, HIS(ng)%pioFile, ncname, DimIDs)
        IF (FoundError(exit_flag, NoError, 3558, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Define time-varying variables.
!-----------------------------------------------------------------------
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        Vinfo( 2)=Vname(2,idtime)
        WRITE (Vinfo( 3),'(a,a)') 'seconds since ', TRIM(Rclock%string)
        Vinfo( 4)=TRIM(Rclock%calendar)
        Vinfo(14)=Vname(4,idtime)
        Vinfo(21)=Vname(6,idtime)
        HIS(ng)%pioVar(idtime)%dkind=PIO_TOUT
        HIS(ng)%pioVar(idtime)%gtype=0
!
        status=def_var(ng, model, HIS(ng)%pioFile,                      &
     &                 HIS(ng)%pioVar(idtime)%vd,                       &
     &                 PIO_TOUT, 1, (/recdim/), Aval, Vinfo, ncname,    &
     &                 SetParAccess = .TRUE.)
        IF (FoundError(exit_flag, NoError, 3579, MyFile)) RETURN
!
!  Define time-varying depth of RHO-points.
!
        IF (Hout(idpthR,ng)) THEN
          Vinfo( 1)=Vname(1,idpthR)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idpthR))
          Vinfo( 3)=Vname(3,idpthR)
          Vinfo(14)=Vname(4,idpthR)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idpthR)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idpthR,ng),r8)
          HIS(ng)%pioVar(idpthR)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idpthR)%gtype=r3dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idpthR)%vd,                     &
     &                   PIO_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, 3708, MyFile)) RETURN
        END IF
!
!  Define time-varying depth of U-points.
!
        IF (Hout(idpthU,ng)) THEN
          Vinfo( 1)=Vname(1,idpthU)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idpthU))
          Vinfo( 3)=Vname(3,idpthU)
          Vinfo(14)=Vname(4,idpthU)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idpthU)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idpthU,ng),r8)
          HIS(ng)%pioVar(idpthU)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idpthU)%gtype=u3dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idpthU)%vd,                     &
     &                   PIO_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, 3732, MyFile)) RETURN
        END IF
!
!  Define time-varying depth of V-points.
!
        IF (Hout(idpthV,ng)) THEN
          Vinfo( 1)=Vname(1,idpthV)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idpthV))
          Vinfo( 3)=Vname(3,idpthV)
          Vinfo(14)=Vname(4,idpthV)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idpthV)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idpthV,ng),r8)
          HIS(ng)%pioVar(idpthV)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idpthV)%gtype=v3dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idpthV)%vd,                     &
     &                   PIO_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, 3756, MyFile)) RETURN
        END IF
!
!  Define time-varying depth of W-points.
!
        IF (Hout(idpthW,ng)) THEN
          Vinfo( 1)=Vname(1,idpthW)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idpthW))
          Vinfo( 3)=Vname(3,idpthW)
          Vinfo(14)=Vname(4,idpthW)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idpthW)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idpthW,ng),r8)
          HIS(ng)%pioVar(idpthW)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idpthW)%gtype=w3dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idpthW)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, 3780, MyFile)) RETURN
        END IF
!
!  Define free-surface.
!
        IF (Hout(idFsur,ng)) THEN
          Vinfo( 1)=Vname(1,idFsur)
          Vinfo( 2)=Vname(2,idFsur)
          Vinfo( 3)=Vname(3,idFsur)
          Vinfo(14)=Vname(4,idFsur)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idFsur)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
          HIS(ng)%pioVar(idFsur)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idFsur)%gtype=r2dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idFsur)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3809, MyFile)) RETURN
        END IF
!
!  Define 2D U-momentum component.
!
        IF (Hout(idUbar,ng)) THEN
          Vinfo( 1)=Vname(1,idUbar)
          Vinfo( 2)=Vname(2,idUbar)
          Vinfo( 3)=Vname(3,idUbar)
          Vinfo(14)=Vname(4,idUbar)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUbar)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbar,ng),r8)
          HIS(ng)%pioVar(idUbar)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idUbar)%gtype=u2dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idUbar)%vd,                     &
     &                   PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3878, MyFile)) RETURN
        END IF
!
!  Define 2D V-momentum component.
!
        IF (Hout(idVbar,ng)) THEN
          Vinfo( 1)=Vname(1,idVbar)
          Vinfo( 2)=Vname(2,idVbar)
          Vinfo( 3)=Vname(3,idVbar)
          Vinfo(14)=Vname(4,idVbar)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVbar)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbar,ng),r8)
          HIS(ng)%pioVar(idVbar)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idVbar)%gtype=v2dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idVbar)%vd,                     &
     &                   PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4010, MyFile)) RETURN
        END IF
!
!  Define 2D Eastward momentum component at RHO-points.
!
        IF (Hout(idu2dE,ng)) THEN
          Vinfo( 1)=Vname(1,idu2dE)
          Vinfo( 2)=Vname(2,idu2dE)
          Vinfo( 3)=Vname(3,idu2dE)
          Vinfo(14)=Vname(4,idu2dE)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idu2dE)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu2dE,ng),r8)
          HIS(ng)%pioVar(idu2dE)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idu2dE)%gtype=r2dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idu2dE)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4142, MyFile)) RETURN
        END IF
!
!  Define 2D Northward momentum component at RHO-points.
!
        IF (Hout(idv2dN,ng)) THEN
          Vinfo( 1)=Vname(1,idv2dN)
          Vinfo( 2)=Vname(2,idv2dN)
          Vinfo( 3)=Vname(3,idv2dN)
          Vinfo(14)=Vname(4,idv2dN)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idv2dN)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv2dN,ng),r8)
          HIS(ng)%pioVar(idv2dN)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idv2dN)%gtype=r2dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idv2dN)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4165, MyFile)) RETURN
        END IF
!
!  Define 3D U-momentum component.
!
        IF (Hout(idUvel,ng)) THEN
          Vinfo( 1)=Vname(1,idUvel)
          Vinfo( 2)=Vname(2,idUvel)
          Vinfo( 3)=Vname(3,idUvel)
          Vinfo(14)=Vname(4,idUvel)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUvel,ng),r8)
          HIS(ng)%pioVar(idUvel)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idUvel)%gtype=u3dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idUvel)%vd,                     &
     &                   PIO_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4190, MyFile)) RETURN
        END IF
!
!  Define 3D U-momentum component at specified constant depth slices.
!
        IF (Hout(idUzsl,ng).and.(Nslice.gt.0)) THEN
          Vinfo( 1)=Vname(1,idUzsl)
          Vinfo( 2)=Vname(2,idUzsl)
          Vinfo( 3)=Vname(3,idUzsl)
          Vinfo(14)=Vname(4,idUzsl)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUzsl)
          Vinfo(22)='coordinates, z_slice'
          Aval(5)=REAL(Iinfo(1,idUzsl,ng),r8)
          HIS(ng)%pioVar(idUzsl)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idUzsl)%gtype=u3dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idUzsl)%vd,                     &
     &                   PIO_FOUT, nvd4, u3dzgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4235, MyFile)) RETURN
        END IF
!
!  Define 3D V-momentum component.
!
        IF (Hout(idVvel,ng)) THEN
          Vinfo( 1)=Vname(1,idVvel)
          Vinfo( 2)=Vname(2,idVvel)
          Vinfo( 3)=Vname(3,idVvel)
          Vinfo(14)=Vname(4,idVvel)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvel,ng),r8)
          HIS(ng)%pioVar(idVvel)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idVvel)%gtype=v3dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idVvel)%vd,                     &
     &                   PIO_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4282, MyFile)) RETURN
        END IF
!
!  Define 3D V-momentum component at specified constant depth slices.
!
        IF (Hout(idVzsl,ng).and.(Nslice.gt.0)) THEN
          Vinfo( 1)=Vname(1,idVzsl)
          Vinfo( 2)=Vname(2,idVzsl)
          Vinfo( 3)=Vname(3,idVzsl)
          Vinfo(14)=Vname(4,idVzsl)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVzsl)
          Vinfo(22)='coordinates, z_slice'
          Aval(5)=REAL(Iinfo(1,idVzsl,ng),r8)
          HIS(ng)%pioVar(idVzsl)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idVzsl)%gtype=v3dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idVzsl)%vd,                     &
     &                   PIO_FOUT, nvd4, v3dzgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4327, MyFile)) RETURN
        END IF
!
!  Define 3D Eastward momentum at RHO-points, A-grid.
!
        IF (Hout(idu3dE,ng)) THEN
          Vinfo( 1)=Vname(1,idu3dE)
          Vinfo( 2)=Vname(2,idu3dE)
          Vinfo( 3)=Vname(3,idu3dE)
          Vinfo(14)=Vname(4,idu3dE)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idu3dE)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu3dE,ng),r8)
          HIS(ng)%pioVar(idu3dE)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idu3dE)%gtype=r3dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idu3dE)%vd,                     &
     &                   PIO_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4374, MyFile)) RETURN
        END IF
!
!  Define 3D Eastward momentum component at RHO-points for selected
!  constant depth slices.
!
        IF (Hout(idUzsE,ng).and.(Nslice.gt.0)) THEN
          Vinfo( 1)=Vname(1,idUzsE)
          Vinfo( 2)=Vname(2,idUzsE)
          Vinfo( 3)=Vname(3,idUzsE)
          Vinfo(14)=Vname(4,idUzsE)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUzsE)
          Vinfo(22)='coordinates, z_slice'
          Aval(5)=REAL(Iinfo(1,idUzsE,ng),r8)
          HIS(ng)%pioVar(idUzsE)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idUzsE)%gtype=r3dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idUzsE)%vd,                     &
     &                   PIO_FOUT, nvd4, t3dzgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4398, MyFile)) RETURN
        END IF
!
!  Define 3D Northward momentum at RHO-points, A-grid.
!
        IF (Hout(idv3dN,ng)) THEN
          Vinfo( 1)=Vname(1,idv3dN)
          Vinfo( 2)=Vname(2,idv3dN)
          Vinfo( 3)=Vname(3,idv3dN)
          Vinfo(14)=Vname(4,idv3dN)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idv3dN)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv3dN,ng),r8)
          HIS(ng)%pioVar(idv3dN)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idv3dN)%gtype=r3dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idv3dN)%vd,                     &
     &                   PIO_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4421, MyFile)) RETURN
        END IF
!
!  Define 3D Northward momentum component at RHO-points for selected
!  constant depth slices.
!
        IF (Hout(idVzsN,ng).and.(Nslice.gt.0)) THEN
          Vinfo( 1)=Vname(1,idVzsN)
          Vinfo( 2)=Vname(2,idVzsN)
          Vinfo( 3)=Vname(3,idVzsN)
          Vinfo(14)=Vname(4,idVzsN)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVzsN)
          Vinfo(22)='coordinates, z_slice'
          Aval(5)=REAL(Iinfo(1,idVzsN,ng),r8)
          HIS(ng)%pioVar(idVzsN)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idVzsN)%gtype=r3dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idVzsN)%vd,                     &
     &                   PIO_FOUT, nvd4, t3dzgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4445, MyFile)) RETURN
        END IF
!
!  Define model potential vorticity at PSI-points.
!
        IF (Hout(id3dPV,ng)) THEN
          Vinfo( 1)=Vname(1,id3dPV)
          Vinfo( 2)=Vname(2,id3dPV)
          Vinfo( 3)=Vname(3,id3dPV)
          Vinfo(14)=Vname(4,id3dPV)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,id3dPV)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,id3dPV,ng),r8)
          HIS(ng)%pioVar(id3dPV)%dkind=PIO_FOUT
          HIS(ng)%pioVar(id3dPV)%gtype=p3dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(id3dPV)%vd,                     &
     &                   PIO_FOUT, nvd4, p3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4468, MyFile)) RETURN
        END IF
!
!  Define model relative vorticity at PSI-points.
!
        IF (Hout(id3dRV,ng)) THEN
          Vinfo( 1)=Vname(1,id3dRV)
          Vinfo( 2)=Vname(2,id3dRV)
          Vinfo( 3)=Vname(3,id3dRV)
          Vinfo(14)=Vname(4,id3dRV)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,id3dRV)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,id3dRV,ng),r8)
          HIS(ng)%pioVar(id3dRV)%dkind=PIO_FOUT
          HIS(ng)%pioVar(id3dRV)%gtype=p3dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(id3dRV)%vd,                     &
     &                   PIO_FOUT, nvd4, p3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4491, MyFile)) RETURN
        END IF
!
!  Define model potential vorticity at PSI-points for selected constant
!  depth slices.
!
        IF (Hout(idPVzs,ng)) THEN
          Vinfo( 1)=Vname(1,idPVzs)
          Vinfo( 2)=Vname(2,idPVzs)
          Vinfo( 3)=Vname(3,idPVzs)
          Vinfo(14)=Vname(4,idPVzs)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idPVzs)
          Vinfo(22)='coordinates, z_slice'
          Aval(5)=REAL(Iinfo(1,idPVzs,ng),r8)
          HIS(ng)%pioVar(idPVzs)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idPVzs)%gtype=p3dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idPVzs)%vd,                     &
     &                   PIO_FOUT, nvd4, p3dzgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4515, MyFile)) RETURN
        END IF
!
!  Define model relative vorticity at PSI-points for selected constant
!  depth slices.
!
        IF (Hout(idRVzs,ng)) THEN
          Vinfo( 1)=Vname(1,idRVzs)
          Vinfo( 2)=Vname(2,idRVzs)
          Vinfo( 3)=Vname(3,idRVzs)
          Vinfo(14)=Vname(4,idRVzs)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idRVzs)
          Vinfo(22)='coordinates, z_slice'
          Aval(5)=REAL(Iinfo(1,id3dRV,ng),r8)
          HIS(ng)%pioVar(idRVzs)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idRVzs)%gtype=p3dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idRVzs)%vd,                     &
     &                   PIO_FOUT, nvd4, p3dzgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4539, MyFile)) RETURN
        END IF
!
!  Define 3D momentum component in the Z-direction.
!
        IF (Hout(idWvel,ng)) THEN
          Vinfo( 1)=Vname(1,idWvel)
          Vinfo( 2)=Vname(2,idWvel)
          Vinfo( 3)=Vname(3,idWvel)
          Vinfo(14)=Vname(4,idWvel)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idWvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWvel,ng),r8)
          HIS(ng)%pioVar(idWvel)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idWvel)%gtype=w3dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idWvel)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4562, MyFile)) RETURN
        END IF
!
!  Define S-coordinate vertical "omega" momentum component.
!
        IF (Hout(idOvel,ng)) THEN
          Vinfo( 1)=Vname(1,idOvel)
          Vinfo( 2)=Vname(2,idOvel)
          Vinfo( 3)='meter second-1'
          Vinfo(14)=Vname(4,idOvel)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idOvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idOvel,ng),r8)
          HIS(ng)%pioVar(idOvel)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idOvel)%gtype=w3dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idOvel)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4585, MyFile)) RETURN
        END IF
!
!  Define S-coordinate implicit vertical "omega" momentum component.
!
        IF (Hout(idOvil,ng)) THEN
          Vinfo( 1)=Vname(1,idOvil)
          Vinfo( 2)=Vname(2,idOvil)
          Vinfo( 3)='meter second-1'
          Vinfo(14)=Vname(4,idOvil)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idOvil)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idOvil,ng),r8)
          HIS(ng)%pioVar(idOvil)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idOvil)%gtype=w3dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idOvil)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4610, MyFile)) RETURN
        END IF
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (Hout(idTvar(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTvar(itrc))
            Vinfo( 2)=Vname(2,idTvar(itrc))
            Vinfo( 3)=Vname(3,idTvar(itrc))
            Vinfo(14)=Vname(4,idTvar(itrc))
            Vinfo(16)=Vname(1,idtime)
            Vinfo(21)=Vname(6,idTvar(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idTvar(itrc),ng),r8)
            HIS(ng)%pioTrc(itrc)%dkind=PIO_FOUT
            HIS(ng)%pioTrc(itrc)%gtype=r3dvar
!
            status=def_var(ng, model, HIS(ng)%pioFile,                  &
     &                     HIS(ng)%pioTrc(itrc)%vd,                     &
     &                     PIO_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4642, MyFile)) RETURN
          END IF
        END DO
!
!  Define surface tracer type variables at specified depth slices.
!
        DO itrc=1,NT(ng)
          IF (Hout(idzslT(itrc),ng).and.(Nslice.gt.0)) THEN
            Vinfo( 1)=Vname(1,idzslT(itrc))
            Vinfo( 2)=Vname(2,idzslT(itrc))
            Vinfo( 3)=Vname(3,idzslT(itrc))
            Vinfo(14)=Vname(4,idzslT(itrc))
            Vinfo(16)=Vname(1,idtime)
            Vinfo(21)=Vname(6,idzslT(itrc))
            Vinfo(22)='coordinates, z_slice'
            Aval(5)=REAL(r3dvar,r8)
            HIS(ng)%pioVar(idzslT(itrc))%dkind=PIO_FOUT
            HIS(ng)%pioVar(idzslT(itrc))%gtype=r2dvar
!
            status=def_var(ng, model, HIS(ng)%pioFile,                  &
     &                     HIS(ng)%pioVar(idzslT(itrc))%vd,             &
     &                     PIO_FOUT, nvd4, t3dzgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4674, MyFile)) RETURN
          END IF
        END DO
!
!  Define density anomaly.
!
        IF (Hout(idDano,ng)) THEN
          Vinfo( 1)=Vname(1,idDano)
          Vinfo( 2)=Vname(2,idDano)
          Vinfo( 3)=Vname(3,idDano)
          Vinfo(14)=Vname(4,idDano)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idDano)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idDano,ng),r8)
          HIS(ng)%pioVar(idDano)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idDano)%gtype=r3dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idDano)%vd,                     &
     &                   PIO_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4731, MyFile)) RETURN
        END IF
!
!  Define vertical viscosity coefficient.
!
        IF (Hout(idVvis,ng)) THEN
          Vinfo( 1)=Vname(1,idVvis)
          Vinfo( 2)=Vname(2,idVvis)
          Vinfo( 3)=Vname(3,idVvis)
          Vinfo(14)=Vname(4,idVvis)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVvis)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvis,ng),r8)
          HIS(ng)%pioVar(idVvis)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idVvis)%gtype=w3dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idVvis)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, 4833, MyFile)) RETURN
        END IF
!
!  Define vertical diffusion coefficient for potential temperature.
!
        IF (Hout(idTdif,ng)) THEN
          Vinfo( 1)=Vname(1,idTdif)
          Vinfo( 2)=Vname(2,idTdif)
          Vinfo( 3)=Vname(3,idTdif)
          Vinfo(14)=Vname(4,idTdif)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idTdif)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTdif,ng),r8)
          HIS(ng)%pioVar(idTdif)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idTdif)%gtype=w3dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idTdif)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, 4857, MyFile)) RETURN
        END IF
!
!  Define vertical diffusion coefficient for salinity.
!
        IF (Hout(idSdif,ng)) THEN
          Vinfo( 1)=Vname(1,idSdif)
          Vinfo( 2)=Vname(2,idSdif)
          Vinfo( 3)=Vname(3,idSdif)
          Vinfo(14)=Vname(4,idSdif)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idSdif)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSdif,ng),r8)
          HIS(ng)%pioVar(idSdif)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idSdif)%gtype=w3dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idSdif)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, 4883, MyFile)) RETURN
        END IF
!
!  Define turbulent kinetic energy.
!
        IF (Hout(idMtke,ng)) THEN
          Vinfo( 1)=Vname(1,idMtke)
          Vinfo( 2)=Vname(2,idMtke)
          Vinfo( 3)=Vname(3,idMtke)
          Vinfo(14)=Vname(4,idMtke)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idMtke)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idMtke,ng),r8)
          HIS(ng)%pioVar(idMtke)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idMtke)%gtype=w3dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idMtke)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, 4909, MyFile)) RETURN
        END IF
!
!  Define turbulent kinetic energy time length scale.
!
        IF (Hout(idMtls,ng)) THEN
          Vinfo( 1)=Vname(1,idMtls)
          Vinfo( 2)=Vname(2,idMtls)
          Vinfo( 3)=Vname(3,idMtls)
          Vinfo(14)=Vname(4,idMtls)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idMtls)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idMtls,ng),r8)
          HIS(ng)%pioVar(idMtls)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idMtls)%gtype=w3dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idMtls)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, 4955, MyFile)) RETURN
        END IF
!
!  Define surface air pressure.
!
        IF (Hout(idPair,ng)) THEN
          Vinfo( 1)=Vname(1,idPair)
          Vinfo( 2)=Vname(2,idPair)
          Vinfo( 3)=Vname(3,idPair)
          Vinfo(14)=Vname(4,idPair)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idPair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idPair,ng),r8)
          HIS(ng)%pioVar(idPair)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idPair)%gtype=r2dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idPair)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 5023, MyFile)) RETURN
        END IF
!
!  Define surface winds.
!
        IF (Hout(idUair,ng)) THEN
          Vinfo( 1)=Vname(1,idUair)
          Vinfo( 2)=Vname(2,idUair)
          Vinfo( 3)=Vname(3,idUair)
          Vinfo(14)=Vname(4,idUair)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUair,ng),r8)
          HIS(ng)%pioVar(idUair)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idUair)%gtype=r2dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idUair)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 5048, MyFile)) RETURN
        END IF
!
        IF (Hout(idVair,ng)) THEN
          Vinfo( 1)=Vname(1,idVair)
          Vinfo( 2)=Vname(2,idVair)
          Vinfo( 3)=Vname(3,idVair)
          Vinfo(14)=Vname(4,idVair)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVair,ng),r8)
          HIS(ng)%pioVar(idVair)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idVair)%gtype=r2dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idVair)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 5069, MyFile)) RETURN
        END IF
!
!  Define Eastward/Northward surface winds at RHO-points.
!
        IF (Hout(idUaiE,ng)) THEN
          Vinfo( 1)=Vname(1,idUaiE)
          Vinfo( 2)=Vname(2,idUaiE)
          Vinfo( 3)=Vname(3,idUaiE)
          Vinfo(14)=Vname(4,idUaiE)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUaiE)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUaiE,ng),r8)
          HIS(ng)%pioVar(idUaiE)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idUaiE)%gtype=r2dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idUaiE)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 5092, MyFile)) RETURN
        END IF
!
        IF (Hout(idVaiN,ng)) THEN
          Vinfo( 1)=Vname(1,idVaiN)
          Vinfo( 2)=Vname(2,idVaiN)
          Vinfo( 3)=Vname(3,idVaiN)
          Vinfo(14)=Vname(4,idVaiN)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVaiN)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVaiN,ng),r8)
          HIS(ng)%pioVar(idVaiN)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idVaiN)%gtype=r2dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idVaiN)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 5113, MyFile)) RETURN
        END IF
!
!  Define surface active tracer fluxes.
!
        DO itrc=1,NAT
          IF (Hout(idTsur(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTsur(itrc))
            Vinfo( 2)=Vname(2,idTsur(itrc))
            Vinfo( 3)=Vname(3,idTsur(itrc))
            IF (itrc.eq.itemp) THEN
              Vinfo(11)='upward flux, cooling'
              Vinfo(12)='downward flux, heating'
            ELSE IF (itrc.eq.isalt) THEN
              Vinfo(11)='upward flux, freshening (net precipitation)'
              Vinfo(12)='downward flux, salting (net evaporation)'
            END IF
            Vinfo(14)=Vname(4,idTsur(itrc))
            Vinfo(16)=Vname(1,idtime)
            Vinfo(21)=Vname(6,idTsur(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idTsur(itrc),ng),r8)
            HIS(ng)%pioVar(idTsur(itrc))%dkind=PIO_FOUT
            HIS(ng)%pioVar(idTsur(itrc))%gtype=r2dvar
!
            status=def_var(ng, model, HIS(ng)%pioFile,                  &
     &                     HIS(ng)%pioVar(idTsur(itrc))%vd,             &
     &                     PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 5145, MyFile)) RETURN
          END IF
        END DO
!
!  Define latent heat flux.
!
        IF (Hout(idLhea,ng)) THEN
          Vinfo( 1)=Vname(1,idLhea)
          Vinfo( 2)=Vname(2,idLhea)
          Vinfo( 3)=Vname(3,idLhea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLhea)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idLhea)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idLhea,ng),r8)
          HIS(ng)%pioVar(idLhea)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idLhea)%gtype=r2dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idLhea)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 5173, MyFile)) RETURN
        END IF
!
!  Define sensible heat flux.
!
        IF (Hout(idShea,ng)) THEN
          Vinfo( 1)=Vname(1,idShea)
          Vinfo( 2)=Vname(2,idShea)
          Vinfo( 3)=Vname(3,idShea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idShea)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idShea)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idShea,ng),r8)
          HIS(ng)%pioVar(idShea)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idShea)%gtype=r2dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idShea)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 5198, MyFile)) RETURN
        END IF
!
!  Define net longwave radiation flux.
!
        IF (Hout(idLrad,ng)) THEN
          Vinfo( 1)=Vname(1,idLrad)
          Vinfo( 2)=Vname(2,idLrad)
          Vinfo( 3)=Vname(3,idLrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLrad)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idLrad)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idLrad,ng),r8)
          HIS(ng)%pioVar(idLrad)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idLrad)%gtype=r2dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idLrad)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 5223, MyFile)) RETURN
        END IF
!
!  Define atmospheric air temperature.
!
        IF (Hout(idTair,ng)) THEN
          Vinfo( 1)=Vname(1,idTair)
          Vinfo( 2)=Vname(2,idTair)
          Vinfo( 3)=Vname(3,idTair)
          Vinfo(14)=Vname(4,idTair)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idTair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTair,ng),r8)
          HIS(ng)%pioVar(idTair)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idTair)%gtype=r2dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idTair)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 5249, MyFile)) RETURN
        END IF
!
!  Define evaporation rate.
!
        IF (Hout(idevap,ng)) THEN
          Vinfo( 1)=Vname(1,idevap)
          Vinfo( 2)=Vname(2,idevap)
          Vinfo( 3)=Vname(3,idevap)
          Vinfo(11)='downward flux, freshening (condensation)'
          Vinfo(12)='upward flux, salting (evaporation)'
          Vinfo(14)=Vname(4,idevap)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idevap)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idevap,ng),r8)
          HIS(ng)%pioVar(idevap)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idevap)%gtype=r2dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idevap)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 5276, MyFile)) RETURN
        END IF
!
!  Define precipitation rate.
!
        IF (Hout(idrain,ng)) THEN
          Vinfo( 1)=Vname(1,idrain)
          Vinfo( 2)=Vname(2,idrain)
          Vinfo( 3)=Vname(3,idrain)
          Vinfo(11)='upward flux, salting (NOT POSSIBLE)'
          Vinfo(12)='downward flux, freshening (precipitation)'
          Vinfo(14)=Vname(4,idrain)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idrain)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idrain,ng),r8)
          HIS(ng)%pioVar(idrain)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idrain)%gtype=r2dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idrain)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 5301, MyFile)) RETURN
        END IF
!
!  Define E-P flux.
!
        IF (Hout(idEmPf,ng)) THEN
          Vinfo( 1)=Vname(1,idEmPf)
          Vinfo( 2)=Vname(2,idEmPf)
          Vinfo( 3)=Vname(3,idEmPf)
          Vinfo(11)='upward flux, freshening (net precipitation)'
          Vinfo(12)='downward flux, salting (net evaporation)'
          Vinfo(14)=Vname(4,idEmPf)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idEmPf)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idEmPf,ng),r8)
          HIS(ng)%pioVar(idEmPf)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idEmPf)%gtype=r2dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idEmPf)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 5328, MyFile)) RETURN
        END IF
!
!  Define net shortwave radiation flux.
!
        IF (Hout(idSrad,ng)) THEN
          Vinfo( 1)=Vname(1,idSrad)
          Vinfo( 2)=Vname(2,idSrad)
          Vinfo( 3)=Vname(3,idSrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idSrad)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idSrad)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSrad,ng),r8)
          HIS(ng)%pioVar(idSrad)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idSrad)%gtype=r2dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idSrad)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 5355, MyFile)) RETURN
        END IF
!
!  Define surface U-momentum stress.
!
        IF (Hout(idUsms,ng)) THEN
          Vinfo( 1)=Vname(1,idUsms)
          Vinfo( 2)=Vname(2,idUsms)
          Vinfo( 3)=Vname(3,idUsms)
          Vinfo(14)=Vname(4,idUsms)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUsms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUsms,ng),r8)
          HIS(ng)%pioVar(idUsms)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idUsms)%gtype=u2dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idUsms)%vd,                     &
     &                   PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 5380, MyFile)) RETURN
        END IF
!
!  Define surface V-momentum stress.
!
        IF (Hout(idVsms,ng)) THEN
          Vinfo( 1)=Vname(1,idVsms)
          Vinfo( 2)=Vname(2,idVsms)
          Vinfo( 3)=Vname(3,idVsms)
          Vinfo(14)=Vname(4,idVsms)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVsms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVsms,ng),r8)
          HIS(ng)%pioVar(idVsms)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idVsms)%gtype=v2dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idVsms)%vd,                     &
     &                   PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 5403, MyFile)) RETURN
        END IF
!
!  Define bottom U-momentum stress.
!
        IF (Hout(idUbms,ng)) THEN
          Vinfo( 1)=Vname(1,idUbms)
          Vinfo( 2)=Vname(2,idUbms)
          Vinfo( 3)=Vname(3,idUbms)
          Vinfo(14)=Vname(4,idUbms)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUbms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbms,ng),r8)
          HIS(ng)%pioVar(idUbms)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idUbms)%gtype=u2dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idUbms)%vd,                     &
     &                   PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 5426, MyFile)) RETURN
        END IF
!
!  Define bottom V-momentum stress.
!
        IF (Hout(idVbms,ng)) THEN
          Vinfo( 1)=Vname(1,idVbms)
          Vinfo( 2)=Vname(2,idVbms)
          Vinfo( 3)=Vname(3,idVbms)
          Vinfo(14)=Vname(4,idVbms)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVbms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbms,ng),r8)
          HIS(ng)%pioVar(idVbms)%dkind=PIO_FOUT
          HIS(ng)%pioVar(idVbms)%gtype=v2dvar
!
          status=def_var(ng, model, HIS(ng)%pioFile,                    &
     &                   HIS(ng)%pioVar(idVbms)%vd,                     &
     &                   PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 5449, MyFile)) RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL pio_netcdf_enddef (ng, model, ncname, HIS(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, 5514, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, model, HIS(ng)%pioFile, ncname)
        IF (FoundError(exit_flag, NoError, 5521, MyFile)) RETURN
      END IF DEFINE
!
!=======================================================================
!  Open an existing history file, check its contents, and prepare for
!  appending data.
!=======================================================================
!
      QUERY : IF (.not.ldef) THEN
        ncname=HIS(ng)%name
!
!  Open history file for read/write.
!
        CALL pio_netcdf_open (ng, model, ncname, 1, HIS(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, 5535, MyFile)) THEN
          WRITE (stdout,60) TRIM(ncname)
          RETURN
        END IF
!
!  Inquire about the dimensions and check for consistency.
!
        CALL pio_netcdf_check_dim (ng, model, ncname,                   &
     &                             pioFile = HIS(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, 5544, MyFile)) RETURN
!
!  Inquire about the variables.
!
        CALL pio_netcdf_inq_var (ng, model, ncname,                     &
     &                           pioFile = HIS(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, 5550, MyFile)) RETURN
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  history variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            HIS(ng)%pioVar(idtime)%vd=var_desc(i)
            HIS(ng)%pioVar(idtime)%dkind=PIO_TOUT
            HIS(ng)%pioVar(idtime)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idpthR))) THEN
            got_var(idpthR)=.TRUE.
            HIS(ng)%pioVar(idpthR)%vd=var_desc(i)
            HIS(ng)%pioVar(idpthR)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idpthR)%gtype=r3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idpthU))) THEN
            got_var(idpthU)=.TRUE.
            HIS(ng)%pioVar(idpthU)%vd=var_desc(i)
            HIS(ng)%pioVar(idpthU)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idpthU)%gtype=u3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idpthV))) THEN
            got_var(idpthV)=.TRUE.
            HIS(ng)%pioVar(idpthV)%vd=var_desc(i)
            HIS(ng)%pioVar(idpthV)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idpthV)%gtype=v3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idpthW))) THEN
            got_var(idpthW)=.TRUE.
            HIS(ng)%pioVar(idpthW)%vd=var_desc(i)
            HIS(ng)%pioVar(idpthW)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idpthW)%gtype=w3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            HIS(ng)%pioVar(idFsur)%vd=var_desc(i)
            HIS(ng)%pioVar(idFsur)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idFsur)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            HIS(ng)%pioVar(idUbar)%vd=var_desc(i)
            HIS(ng)%pioVar(idUbar)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idUbar)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            HIS(ng)%pioVar(idVbar)%vd=var_desc(i)
            HIS(ng)%pioVar(idVbar)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idVbar)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu2dE))) THEN
            got_var(idu2dE)=.TRUE.
            HIS(ng)%pioVar(idu2dE)%vd=var_desc(i)
            HIS(ng)%pioVar(idu2dE)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idu2dE)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv2dN))) THEN
            got_var(idv2dN)=.TRUE.
            HIS(ng)%pioVar(idv2dN)%vd=var_desc(i)
            HIS(ng)%pioVar(idv2dN)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idv2dN)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            HIS(ng)%pioVar(idUvel)%vd=var_desc(i)
            HIS(ng)%pioVar(idUvel)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idUvel)%gtype=u3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            HIS(ng)%pioVar(idVvel)%vd=var_desc(i)
            HIS(ng)%pioVar(idVvel)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idVvel)%gtype=v3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUzsl))) THEN
            got_var(idUzsl)=.TRUE.
            HIS(ng)%pioVar(idUzsl)%vd=var_desc(i)
            HIS(ng)%pioVar(idUzsl)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idUzsl)%gtype=u3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVzsl))) THEN
            got_var(idVzsl)=.TRUE.
            HIS(ng)%pioVar(idVzsl)%vd=var_desc(i)
            HIS(ng)%pioVar(idVzsl)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idVzsl)%gtype=v3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu3dE))) THEN
            got_var(idu3dE)=.TRUE.
            HIS(ng)%pioVar(idu3dE)%vd=var_desc(i)
            HIS(ng)%pioVar(idu3dE)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idu3dE)%gtype=r3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv3dN))) THEN
            got_var(idv3dN)=.TRUE.
            HIS(ng)%pioVar(idv3dN)%vd=var_desc(i)
            HIS(ng)%pioVar(idv3dN)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idv3dN)%gtype=r3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUzsE))) THEN
            got_var(idUzsE)=.TRUE.
            HIS(ng)%pioVar(idUzsE)%vd=var_desc(i)
            HIS(ng)%pioVar(idUzsE)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idUzsE)%gtype=r3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVzsN))) THEN
            got_var(idVzsN)=.TRUE.
            HIS(ng)%pioVar(idVzsN)%vd=var_desc(i)
            HIS(ng)%pioVar(idVzsN)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idVzsN)%gtype=r3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,id3dPV))) THEN
            got_var(id3dPV)=.TRUE.
            HIS(ng)%pioVar(id3dPV)%vd=var_desc(i)
            HIS(ng)%pioVar(id3dPV)%dkind=PIO_FOUT
            HIS(ng)%pioVar(id3dPV)%gtype=p3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,id3dRV))) THEN
            got_var(id3dRV)=.TRUE.
            HIS(ng)%pioVar(id3dRV)%vd=var_desc(i)
            HIS(ng)%pioVar(id3dRV)%dkind=PIO_FOUT
            HIS(ng)%pioVar(id3dRV)%gtype=p3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPVzs))) THEN
            got_var(idPVzs)=.TRUE.
            HIS(ng)%pioVar(idPVzs)%vd=var_desc(i)
            HIS(ng)%pioVar(idPVzs)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idPVzs)%gtype=p3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRVzs))) THEN
            got_var(idRVzs)=.TRUE.
            HIS(ng)%pioVar(idRVzs)%vd=var_desc(i)
            HIS(ng)%pioVar(idRVzs)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idRVzs)%gtype=p3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWvel))) THEN
            got_var(idWvel)=.TRUE.
            HIS(ng)%pioVar(idWvel)%vd=var_desc(i)
            HIS(ng)%pioVar(idWvel)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idWvel)%gtype=w3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOvel))) THEN
            got_var(idOvel)=.TRUE.
            HIS(ng)%pioVar(idOvel)%vd=var_desc(i)
            HIS(ng)%pioVar(idOvel)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idOvel)%gtype=w3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOvil))) THEN
            got_var(idOvil)=.TRUE.
            HIS(ng)%pioVar(idOvil)%vd=var_desc(i)
            HIS(ng)%pioVar(idOvil)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idOvil)%gtype=w3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            HIS(ng)%pioVar(idDano)%vd=var_desc(i)
            HIS(ng)%pioVar(idDano)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idDano)%gtype=r3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvis))) THEN
            got_var(idVvis)=.TRUE.
            HIS(ng)%pioVar(idVvis)%vd=var_desc(i)
            HIS(ng)%pioVar(idVvis)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idVvis)%gtype=w3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTdif))) THEN
            got_var(idTdif)=.TRUE.
            HIS(ng)%pioVar(idTdif)%vd=var_desc(i)
            HIS(ng)%pioVar(idTdif)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idTdif)%gtype=w3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSdif))) THEN
            got_var(idSdif)=.TRUE.
            HIS(ng)%pioVar(idSdif)%vd=var_desc(i)
            HIS(ng)%pioVar(idSdif)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idSdif)%gtype=w3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idMtke))) THEN
            got_var(idMtke)=.TRUE.
            HIS(ng)%pioVar(idMtke)%vd=var_desc(i)
            HIS(ng)%pioVar(idMtke)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idMtke)%gtype=w3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idMtls))) THEN
            got_var(idMtls)=.TRUE.
            HIS(ng)%pioVar(idMtls)%vd=var_desc(i)
            HIS(ng)%pioVar(idMtls)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idMtls)%gtype=w3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPair))) THEN
            got_var(idPair)=.TRUE.
            HIS(ng)%pioVar(idPair)%vd=var_desc(i)
            HIS(ng)%pioVar(idPair)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idPair)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUair))) THEN
            got_var(idUair)=.TRUE.
            HIS(ng)%pioVar(idUair)%vd=var_desc(i)
            HIS(ng)%pioVar(idUair)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idUair)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVair))) THEN
            got_var(idVair)=.TRUE.
            HIS(ng)%pioVar(idVair)%vd=var_desc(i)
            HIS(ng)%pioVar(idVair)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idVair)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUaiE))) THEN
            got_var(idUair)=.TRUE.
            HIS(ng)%pioVar(idUaiE)%vd=var_desc(i)
            HIS(ng)%pioVar(idUaiE)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idUaiE)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVaiN))) THEN
            got_var(idVair)=.TRUE.
            HIS(ng)%pioVar(idVaiN)%vd=var_desc(i)
            HIS(ng)%pioVar(idVaiN)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idVaiN)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLhea))) THEN
            got_var(idLhea)=.TRUE.
            HIS(ng)%pioVar(idLhea)%vd=var_desc(i)
            HIS(ng)%pioVar(idLhea)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idLhea)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idShea))) THEN
            got_var(idShea)=.TRUE.
            HIS(ng)%pioVar(idShea)%vd=var_desc(i)
            HIS(ng)%pioVar(idShea)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idShea)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLrad))) THEN
            got_var(idLrad)=.TRUE.
            HIS(ng)%pioVar(idLrad)%vd=var_desc(i)
            HIS(ng)%pioVar(idLrad)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idLrad)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTair))) THEN
            got_var(idTair)=.TRUE.
            HIS(ng)%pioVar(idTair)%vd=var_desc(i)
            HIS(ng)%pioVar(idTair)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idTair)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idevap))) THEN
            got_var(idevap)=.TRUE.
            HIS(ng)%pioVar(idevap)%vd=var_desc(i)
            HIS(ng)%pioVar(idevap)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idevap)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idrain))) THEN
            got_var(idrain)=.TRUE.
            HIS(ng)%pioVar(idrain)%vd=var_desc(i)
            HIS(ng)%pioVar(idrain)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idrain)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idEmPf))) THEN
            got_var(idEmPf)=.TRUE.
            HIS(ng)%pioVar(idEmPf)%vd=var_desc(i)
            HIS(ng)%pioVar(idEmPf)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idEmPf)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSrad))) THEN
            got_var(idSrad)=.TRUE.
            HIS(ng)%pioVar(idSrad)%vd=var_desc(i)
            HIS(ng)%pioVar(idSrad)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idSrad)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            HIS(ng)%pioVar(idUsms)%vd=var_desc(i)
            HIS(ng)%pioVar(idUsms)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idUsms)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            HIS(ng)%pioVar(idVsms)%vd=var_desc(i)
            HIS(ng)%pioVar(idVsms)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idVsms)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbms))) THEN
            got_var(idUbms)=.TRUE.
            HIS(ng)%pioVar(idUbms)%vd=var_desc(i)
            HIS(ng)%pioVar(idUbms)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idUbms)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbms))) THEN
            got_var(idVbms)=.TRUE.
            HIS(ng)%pioVar(idVbms)%vd=var_desc(i)
            HIS(ng)%pioVar(idVbms)%dkind=PIO_FOUT
            HIS(ng)%pioVar(idVbms)%gtype=v2dvar
          END IF
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              HIS(ng)%pioTrc(itrc)%vd=var_desc(i)
              HIS(ng)%pioTrc(itrc)%dkind=PIO_FOUT
              HIS(ng)%pioTrc(itrc)%gtype=r3dvar
            END IF
          END DO
          DO itrc=1,NAT
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idzslT(itrc)))) THEN
              got_var(idzslT(itrc))=.TRUE.
              HIS(ng)%pioVar(idzslT(itrc))%vd=var_desc(i)
              HIS(ng)%pioVar(idzslT(itrc))%dkind=PIO_FOUT
              HIS(ng)%pioVar(idzslT(itrc))%gtype=r3dvar
            END IF
          END DO
          DO itrc=1,NAT
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTsur(itrc)))) THEN
              got_var(idTsur(itrc))=.TRUE.
              HIS(ng)%pioVar(idTsur(itrc))%vd=var_desc(i)
              HIS(ng)%pioVar(idTsur(itrc))%dkind=PIO_FOUT
              HIS(ng)%pioVar(idTsur(itrc))%gtype=r2dvar
            END IF
          END DO
        END DO
!
!  Check if history variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idpthR).and.Hout(idpthR,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idpthR)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idpthU).and.Hout(idpthU,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idpthU)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idpthV).and.Hout(idpthV,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idpthV)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idpthW).and.Hout(idpthW,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idpthW)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsur).and.Hout(idFsur,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar).and.Hout(idUbar,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar).and.Hout(idVbar,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idu2dE).and.Hout(idu2dE,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idu2dE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv2dN).and.Hout(idv2dN,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idv2dN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUvel).and.Hout(idUvel,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel).and.Hout(idVvel,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUzsl).and.Hout(idUzsl,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUzsl)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVzsl).and.Hout(idVzsl,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVzsl)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idu3dE).and.Hout(idu3dE,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idu3dE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv3dN).and.Hout(idv3dN,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idv3dN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUzsE).and.Hout(idUzsE,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUzsE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVzsN).and.Hout(idVzsN,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVzsN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(id3dPV).and.Hout(id3dPV,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,id3dPV)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(id3dRV).and.Hout(id3dRV,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,id3dRV)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idPVzs).and.Hout(idPVzs,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idPVzs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idRVzs).and.Hout(idRVzs,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idRVzs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWvel).and.Hout(idWvel,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idWvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idOvel).and.Hout(idOvel,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idOvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idOvil).and.Hout(idOvil,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idOvil)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idDano).and.Hout(idDano,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idDano)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvis).and.Hout(idVvis,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVvis)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTdif).and.Hout(idTdif,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSdif).and.Hout(idSdif,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idMtke).and.Hout(idMtke,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idMtke)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idMtls).and.Hout(idMtls,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idMtls)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idPair).and.Hout(idPair,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idPair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUair).and.Hout(idUair,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVair).and.Hout(idVair,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUaiE).and.Hout(idUaiE,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUaiE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVaiN).and.Hout(idVaiN,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVaiN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLhea).and.Hout(idLhea,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idLhea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idShea).and.Hout(idShea,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idShea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLrad).and.Hout(idLrad,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idLrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTair).and.Hout(idTair,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idevap).and.Hout(idevap,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idevap)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idrain).and.Hout(idrain,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idrain)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idEmPf).and.Hout(idEmPf,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idEmPf)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSrad).and.Hout(idSrad,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUsms).and.Hout(idUsms,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms).and.Hout(idVsms,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbms).and.Hout(idUbms,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbms).and.Hout(idVbms,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc)).and.Hout(idTvar(itrc),ng)) THEN
            IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO itrc=1,NAT
          IF (.not.got_var(idzslT(itrc)).and.Hout(idzslT(itrc),ng)) THEN
            IF (Master) WRITE (stdout,70) TRIM(Vname(1,idzslT(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO itrc=1,NAT
          IF (.not.got_var(idTsur(itrc)).and.Hout(idTsur(itrc),ng)) THEN
            IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTsur(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
!
!  Set unlimited time record dimension to the appropriate value.
!
        IF (ndefHIS(ng).gt.0) THEN
          HIS(ng)%Rindex=((ntstart(ng)-1)-                              &
     &                    ndefHIS(ng)*((ntstart(ng)-1)/ndefHIS(ng)))/   &
     &                   nHIS(ng)
        ELSE
          HIS(ng)%Rindex=(ntstart(ng)-1)/nHIS(ng)
        END IF
        HIS(ng)%Rindex=MIN(HIS(ng)%Rindex,rec_size)
      END IF QUERY
!
  10  FORMAT (2x,'DEF_HIS_PIO      - creating history file,',t56,       &
     &        'Grid ',i2.2,': ',a)
  20  FORMAT (2x,'DEF_HIS_PIO      - inquiring history file,',t56,      &
     &        'Grid ',i2.2,': ',a)
  30  FORMAT (/,' DEF_HIS_PIO - unable to create history NetCDF',       &
     &        ' file: ',a)
  40  FORMAT ('time dependent',1x,a)
  50  FORMAT (1pe11.4,1x,'millimeter')
  60  FORMAT (/,' DEF_HIS_PIO - unable to open history NetCDF file: ',a)
  70  FORMAT (/,' DEF_HIS_PIO - unable to find variable: ',a,2x,        &
     &        ' in history NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE def_his_pio
!
      END MODULE def_his_mod
