      MODULE def_avg_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module creates output TIME-AVERAGED file using either the      !
!  standard NetCDF library or the Parallel-IO (PIO) library.  It       !
!  defines its dimensions, attributes, and variables.                  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE def_dim_mod,         ONLY : def_dim
      USE def_info_mod,        ONLY : def_info
      USE def_var_mod,         ONLY : def_var
      USE strings_mod,         ONLY : FoundError
      USE wrt_info_mod,        ONLY : wrt_info
!
      implicit none
!
      PUBLIC  :: def_avg
      PRIVATE :: def_avg_nf90
      PRIVATE :: def_avg_pio
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE def_avg (ng, ldef)
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
      integer :: model
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/def_avg.F"
!
!-----------------------------------------------------------------------
!  Create a new time-averaged file according to IO type.
!-----------------------------------------------------------------------
!
      model=iNLM
!
      SELECT CASE (AVG(ng)%IOtype)
        CASE (io_nf90)
          CALL def_avg_nf90 (ng, model, ldef)
        CASE (io_pio)
          CALL def_avg_pio (ng, model, ldef)
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) AVG(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, 116, MyFile)) RETURN
!
  10  FORMAT (' DEF_AVG - Illegal output file type, io_type = ',i0,     &
     &        /,11x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE def_avg
!
!***********************************************************************
      SUBROUTINE def_avg_nf90 (ng, model, ldef)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
      logical, intent(in) :: ldef
!
!  Local variable declarations.
!
      logical :: got_var(NV)
!
      integer :: i, ifield, itrc, j, nvd3, nvd4
      integer :: recdim, status
      integer :: DimIDs(nDimID)
      integer :: p2dgrd(3), t2dgrd(3), u2dgrd(3), v2dgrd(3)
      integer :: p3dgrd(4), t3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)
!
      real(r8) :: Aval(6)
!
      character (len=13)     :: Prefix
      character (len=256)    :: ncname
      character (len=MaxLen) :: Vinfo(Natt)
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/def_avg.F"//", def_avg_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, 175, MyFile)) RETURN
      ncname=AVG(ng)%name
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
!  Create a new averages file.
!=======================================================================
!
      DEFINE : IF (ldef) THEN
        CALL netcdf_create (ng, model, TRIM(ncname), AVG(ng)%ncid)
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
        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'xi_rho',       &
     &                 IOBOUNDS(ng)%xi_rho, DimIDs( 1))
        IF (FoundError(exit_flag, NoError, 205, MyFile)) RETURN
        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'xi_u',         &
     &                 IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (FoundError(exit_flag, NoError, 209, MyFile)) RETURN
        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'xi_v',         &
     &                 IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (FoundError(exit_flag, NoError, 213, MyFile)) RETURN
        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'xi_psi',       &
     &                 IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (FoundError(exit_flag, NoError, 217, MyFile)) RETURN
        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'eta_rho',      &
     &                 IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (FoundError(exit_flag, NoError, 221, MyFile)) RETURN
        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'eta_u',        &
     &                 IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (FoundError(exit_flag, NoError, 225, MyFile)) RETURN
        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'eta_v',        &
     &                 IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (FoundError(exit_flag, NoError, 229, MyFile)) RETURN
        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'eta_psi',      &
     &                 IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (FoundError(exit_flag, NoError, 233, MyFile)) RETURN
        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 's_rho',        &
     &                 N(ng), DimIDs( 9))
        IF (FoundError(exit_flag, NoError, 278, MyFile)) RETURN
        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 's_w',          &
     &                 N(ng)+1, DimIDs(10))
        IF (FoundError(exit_flag, NoError, 282, MyFile)) RETURN
        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'tracer',       &
     &                 NT(ng), DimIDs(11))
        IF (FoundError(exit_flag, NoError, 286, MyFile)) RETURN
        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'boundary',     &
     &                 4, DimIDs(14))
        IF (FoundError(exit_flag, NoError, 329, MyFile)) RETURN
        status=def_dim(ng, model, AVG(ng)%ncid, ncname,                 &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 nf90_unlimited, DimIDs(12))
        IF (FoundError(exit_flag, NoError, 346, MyFile)) RETURN
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
!  Define dimension vectors for staggered variables at PSI-points.
!
        p2dgrd(1)=DimIDs( 4)
        p2dgrd(2)=DimIDs( 8)
        p2dgrd(3)=DimIDs(12)
        p3dgrd(1)=DimIDs( 4)
        p3dgrd(2)=DimIDs( 8)
        p3dgrd(3)=DimIDs( 9)
        p3dgrd(4)=DimIDs(12)
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
        AVG(ng)%Rindex=0
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
!  Set long name prefix string.
!
        Prefix=CHAR(32)                                   ! blank space
!
!-----------------------------------------------------------------------
!  Define time-recordless information variables.
!-----------------------------------------------------------------------
!
        CALL def_info (ng, model, AVG(ng)%ncid, ncname, DimIDs)
        IF (FoundError(exit_flag, NoError, 502, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Define time-varying variables.
!-----------------------------------------------------------------------
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        WRITE (Vinfo( 2),'(a,a)') 'averaged ', TRIM(Vname(2,idtime))
        WRITE (Vinfo( 3),'(a,a)') 'seconds since ', TRIM(Rclock%string)
        Vinfo( 4)=TRIM(Rclock%calendar)
        Vinfo(14)=Vname(4,idtime)
        Vinfo(21)=Vname(6,idtime)
        status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idtime),    &
     &                 NF_TOUT, 1, (/recdim/), Aval, Vinfo, ncname,     &
     &                 SetParAccess = .TRUE.)
        IF (FoundError(exit_flag, NoError, 519, MyFile)) RETURN
!
!  Define free-surface.
!
        IF (Aout(idFsur,ng)) THEN
          Vinfo( 1)=Vname(1,idFsur)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idFsur))
          Vinfo( 3)=Vname(3,idFsur)
          Vinfo(14)=Vname(4,idFsur)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idFsur)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idFsur),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 537, MyFile)) RETURN
        END IF
!
!  Define 2D momentum in the XI-direction.
!
        IF (Aout(idUbar,ng)) THEN
          Vinfo( 1)=Vname(1,idUbar)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUbar))
          Vinfo( 3)=Vname(3,idUbar)
          Vinfo(14)=Vname(4,idUbar)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUbar)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbar,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idUbar),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 578, MyFile)) RETURN
        END IF
!
!  Define 2D momentum in the ETA-direction.
!
        IF (Aout(idVbar,ng)) THEN
          Vinfo( 1)=Vname(1,idVbar)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVbar))
          Vinfo( 3)=Vname(3,idVbar)
          Vinfo(14)=Vname(4,idVbar)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVbar)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbar,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idVbar),  &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 619, MyFile)) RETURN
        END IF
!
!  Define 2D Eastward momentum component at RHO-points.
!
        IF (Aout(idu2dE,ng)) THEN
          Vinfo( 1)=Vname(1,idu2dE)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idu2dE))
          Vinfo( 3)=Vname(3,idu2dE)
          Vinfo(14)=Vname(4,idu2dE)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idu2dE)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu2dE,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idu2dE),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 660, MyFile)) RETURN
        END IF
!
!  Define 2D Northward momentum component at RHO-points.
!
        IF (Aout(idv2dN,ng)) THEN
          Vinfo( 1)=Vname(1,idv2dN)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idv2dN))
          Vinfo( 3)=Vname(3,idv2dN)
          Vinfo(14)=Vname(4,idv2dN)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idv2dN)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv2dN,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idv2dN),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 679, MyFile)) RETURN
        END IF
!
!  Define 3D momentum component in the XI-direction.
!
        IF (Aout(idUvel,ng)) THEN
          Vinfo( 1)=Vname(1,idUvel)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUvel))
          Vinfo( 3)=Vname(3,idUvel)
          Vinfo(14)=Vname(4,idUvel)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUvel,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idUvel),  &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 772, MyFile)) RETURN
        END IF
!
!  Define 3D momentum component in the ETA-direction.
!
        IF (Aout(idVvel,ng)) THEN
          Vinfo( 1)=Vname(1,idVvel)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVvel))
          Vinfo( 3)=Vname(3,idVvel)
          Vinfo(14)=Vname(4,idVvel)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvel,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idVvel),  &
     &                   NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 813, MyFile)) RETURN
        END IF
!
!  Define 3D Eastward momentum component at RHO-points.
!
        IF (Aout(idu3dE,ng)) THEN
          Vinfo( 1)=Vname(1,idu3dE)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idu3dE))
          Vinfo( 3)=Vname(3,idu3dE)
          Vinfo(14)=Vname(4,idu3dE)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idu3dE)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu3dE,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idu3dE),  &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 854, MyFile)) RETURN
        END IF
!
!  Define 3D Northward momentum component at RHO-points.
!
        IF (Aout(idv3dN,ng)) THEN
          Vinfo( 1)=Vname(1,idv3dN)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idv3dN))
          Vinfo( 3)=Vname(3,idv3dN)
          Vinfo(14)=Vname(4,idv3dN)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idv3dN)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv3dN,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idv3dN),  &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 873, MyFile)) RETURN
        END IF
!
!  Define S-coordinate vertical "omega" momentum component.
!
        IF (Aout(idOvel,ng)) THEN
          Vinfo( 1)=Vname(1,idOvel)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idOvel))
          Vinfo( 3)=Vname(3,idOvel)
          Vinfo(14)=Vname(4,idOvel)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idOvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idOvel,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idOvel),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 892, MyFile)) RETURN
        END IF
!
!  Define "true" vertical momentum component.
!
        IF (Aout(idWvel,ng)) THEN
          Vinfo( 1)=Vname(1,idWvel)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWvel))
          Vinfo( 3)=Vname(3,idWvel)
          Vinfo(14)=Vname(4,idWvel)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idWvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWvel,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idWvel),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 911, MyFile)) RETURN
        END IF
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (Aout(idTvar(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTvar(itrc))
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idTvar(itrc)))
            Vinfo( 3)=Vname(3,idTvar(itrc))
            Vinfo(14)=Vname(4,idTvar(itrc))
            Vinfo(16)=Vname(1,idtime)
            Vinfo(21)=Vname(6,idTvar(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Tid(itrc),  &
     &                     NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 939, MyFile)) RETURN
          END IF
        END DO
!
!  Define density anomaly.
!
        IF (Aout(idDano,ng)) THEN
          Vinfo( 1)=Vname(1,idDano)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idDano))
          Vinfo( 3)=Vname(3,idDano)
          Vinfo(14)=Vname(4,idDano)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idDano)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idDano,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idDano),  &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 985, MyFile)) RETURN
        END IF
!
!  Define 2D potential vorticity.
!
        IF (Aout(id2dPV,ng)) THEN
          Vinfo( 1)=Vname(1,id2dPV)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,id2dPV))
          Vinfo( 3)=Vname(3,id2dPV)
          Vinfo(14)=Vname(4,id2dPV)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,id2dPV)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(p2dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(id2dPV),  &
     &                   NF_FOUT, nvd3, p2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1048, MyFile)) RETURN
        END IF
!
!  Define 2D relative vorticity.
!
        IF (Aout(id2dRV,ng)) THEN
          Vinfo( 1)=Vname(1,id2dRV)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,id2dRV))
          Vinfo( 3)=Vname(3,id2dRV)
          Vinfo(14)=Vname(4,id2dRV)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,id2dRV)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(p2dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(id2dRV),  &
     &                   NF_FOUT, nvd3, p2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1067, MyFile)) RETURN
        END IF
!
!  Define 3D potential vorticity.
!
        IF (Aout(id3dPV,ng)) THEN
          Vinfo( 1)=Vname(1,id3dPV)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,id3dPV))
          Vinfo( 3)=Vname(3,id3dPV)
          Vinfo(14)=Vname(4,id3dPV)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,id3dPV)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(p3dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(id3dPV),  &
     &                   NF_FOUT, nvd4, p3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1088, MyFile)) RETURN
        END IF
!
!  Define 3D relative vorticity.
!
        IF (Aout(id3dRV,ng)) THEN
          Vinfo( 1)=Vname(1,id3dRV)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,id3dRV))
          Vinfo( 3)=Vname(3,id3dRV)
          Vinfo(14)=Vname(4,id3dRV)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,id3dRV)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(p3dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(id3dRV),  &
     &                   NF_FOUT, nvd4, p3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1107, MyFile)) RETURN
        END IF
!
!  Define quadratic <zeta*zeta> term.
!
        IF (Aout(idZZav,ng)) THEN
          Vinfo( 1)=Vname(1,idZZav)
          Vinfo( 2)=TRIM(Vname(2,idZZav))
          Vinfo( 3)=Vname(3,idZZav)
          Vinfo(14)=Vname(4,idZZav)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idZZav)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idZZav),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1127, MyFile)) RETURN
        END IF
!
!  Define quadratic <ubar*ubar> term.
!
        IF (Aout(idU2av,ng)) THEN
          Vinfo( 1)=Vname(1,idU2av)
          Vinfo( 2)=TRIM(Vname(2,idU2av))
          Vinfo( 3)=Vname(3,idU2av)
          Vinfo(14)=Vname(4,idU2av)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idU2av)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idU2av),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1146, MyFile)) RETURN
        END IF
!
!  Define quadratic <vbar*vbar> term.
!
        IF (Aout(idV2av,ng)) THEN
          Vinfo( 1)=Vname(1,idV2av)
          Vinfo( 2)=TRIM(Vname(2,idV2av))
          Vinfo( 3)=Vname(3,idV2av)
          Vinfo(14)=Vname(4,idV2av)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idV2av)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idV2av),  &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1165, MyFile)) RETURN
        END IF
!
!  Define u-volume flux.
!
        IF (Aout(idHUav,ng)) THEN
          Vinfo( 1)=Vname(1,idHUav)
          Vinfo( 2)=TRIM(Vname(2,idHUav))
          Vinfo( 3)=Vname(3,idHUav)
          Vinfo(14)=Vname(4,idHUav)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idHUav)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idHUav),  &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1186, MyFile)) RETURN
        END IF
!
!  Define v-volume flux.
!
        IF (Aout(idHVav,ng)) THEN
          Vinfo( 1)=Vname(1,idHVav)
          Vinfo( 2)=TRIM(Vname(2,idHVav))
          Vinfo( 3)=Vname(3,idHVav)
          Vinfo(14)=Vname(4,idHVav)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idHVav)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v3dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idHVav),  &
     &                   NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1205, MyFile)) RETURN
        END IF
!
!  Define quadratic <u*u> term.
!
        IF (Aout(idUUav,ng)) THEN
          Vinfo( 1)=Vname(1,idUUav)
          Vinfo( 2)=TRIM(Vname(2,idUUav))
          Vinfo( 3)=Vname(3,idUUav)
          Vinfo(14)=Vname(4,idUUav)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUUav)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idUUav),  &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1224, MyFile)) RETURN
        END IF
!
!  Define quadratic <u*v> term.
!
        IF (Aout(idUVav,ng)) THEN
          Vinfo( 1)=Vname(1,idUVav)
          Vinfo( 2)=TRIM(Vname(2,idUVav))
          Vinfo( 3)=Vname(3,idUVav)
          Vinfo(14)=Vname(4,idUVav)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUVav)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idUVav),  &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1243, MyFile)) RETURN
        END IF
!
!  Define quadratic <v*v> term.
!
        IF (Aout(idVVav,ng)) THEN
          Vinfo( 1)=Vname(1,idVVav)
          Vinfo( 2)=TRIM(Vname(2,idVVav))
          Vinfo( 3)=Vname(3,idVVav)
          Vinfo(14)=Vname(4,idVVav)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVVav)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v3dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idVVav),  &
     &                   NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1262, MyFile)) RETURN
        END IF
!
!  Define quadratic <t*t> terms.
!
        DO itrc=1,NT(ng)
          IF (Aout(idTTav(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTTav(itrc))
            Vinfo( 2)=TRIM(Vname(2,idTTav(itrc)))
            Vinfo( 3)=Vname(3,idTTav(itrc))
            Vinfo(14)=Vname(4,idTTav(itrc))
            Vinfo(16)=Vname(1,idtime)
            Vinfo(21)=Vname(6,idTTav(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            status=def_var(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idTTav(itrc)), NF_FOUT,          &
     &                     nvd4, t3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 1283, MyFile)) RETURN
          END IF
        END DO
!
!  Define active tracers volume fluxes.
!
        DO itrc=1,NT(ng)
          IF (Aout(iHUTav(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,iHUTav(itrc))
            Vinfo( 2)=TRIM(Vname(2,iHUTav(itrc)))
            Vinfo( 3)=Vname(3,iHUTav(itrc))
            Vinfo(14)=Vname(4,iHUTav(itrc))
            Vinfo(16)=Vname(1,idtime)
            Vinfo(21)=Vname(6,iHUTav(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(u3dvar,r8)
            status=def_var(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(iHUTav(itrc)), NF_FOUT,          &
     &                     nvd4, u3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 1305, MyFile)) RETURN
          END IF
!
          IF (Aout(iHVTav(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,iHVTav(itrc))
            Vinfo( 2)=TRIM(Vname(2,iHVTav(itrc)))
            Vinfo( 3)=Vname(3,iHVTav(itrc))
            Vinfo(14)=Vname(4,iHVTav(itrc))
            Vinfo(16)=Vname(1,idtime)
            Vinfo(21)=Vname(6,iHVTav(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(v3dvar,r8)
            status=def_var(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(iHVTav(itrc)), NF_FOUT,          &
     &                     nvd4, v3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 1323, MyFile)) RETURN
          END IF
        END DO
!
!  Define quadratic <u*t> and <v*t> terms.
!
        DO itrc=1,NT(ng)
          IF (Aout(idUTav(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idUTav(itrc))
            Vinfo( 2)=TRIM(Vname(2,idUTav(itrc)))
            Vinfo( 3)=Vname(3,idUTav(itrc))
            Vinfo(14)=Vname(4,idUTav(itrc))
            Vinfo(16)=Vname(1,idtime)
            Vinfo(21)=Vname(6,idUTav(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(u3dvar,r8)
            status=def_var(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idUTav(itrc)), NF_FOUT,          &
     &                     nvd4, u3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 1345, MyFile)) RETURN
          END IF
!
          IF (Aout(idVTav(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idVTav(itrc))
            Vinfo( 2)=TRIM(Vname(2,idVTav(itrc)))
            Vinfo( 3)=Vname(3,idVTav(itrc))
            Vinfo(14)=Vname(4,idVTav(itrc))
            Vinfo(16)=Vname(1,idtime)
            Vinfo(21)=Vname(6,idVTav(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(v3dvar,r8)
            status=def_var(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idVTav(itrc)), NF_FOUT,          &
     &                     nvd4, v3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 1363, MyFile)) RETURN
          END IF
        END DO
!
!  Define vertical viscosity coefficient.
!
        IF (Aout(idVvis,ng)) THEN
          Vinfo( 1)=Vname(1,idVvis)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVvis))
          Vinfo( 3)=Vname(3,idVvis)
          Vinfo(14)=Vname(4,idVvis)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVvis)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvis,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idVvis),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, 1383, MyFile)) RETURN
        END IF
!
!  Define vertical diffusion coefficient for potential temperature.
!
        IF (Aout(idTdif,ng)) THEN
          Vinfo( 1)=Vname(1,idTdif)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idTdif))
          Vinfo( 3)=Vname(3,idTdif)
          Vinfo(14)=Vname(4,idTdif)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idTdif)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTdif,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idTdif),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, 1400, MyFile)) RETURN
        END IF
!
!  Define vertical diffusion coefficient for salinity.
!
        IF (Aout(idSdif,ng)) THEN
          Vinfo( 1)=Vname(1,idSdif)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idSdif))
          Vinfo( 3)=Vname(3,idSdif)
          Vinfo(14)=Vname(4,idSdif)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idSdif)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSdif,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idSdif),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, 1419, MyFile)) RETURN
        END IF
!
!  Define surface air pressure.
!
        IF (Aout(idPair,ng)) THEN
          Vinfo( 1)=Vname(1,idPair)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idPair))
          Vinfo( 3)=Vname(3,idPair)
          Vinfo(14)=Vname(4,idPair)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idPair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idPair,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idPair),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1441, MyFile)) RETURN
        END IF
!
!  Define surface air temperature.
!
        IF (Aout(idTair,ng)) THEN
          Vinfo( 1)=Vname(1,idTair)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idTair))
          Vinfo( 3)=Vname(3,idTair)
          Vinfo(14)=Vname(4,idTair)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idTair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTair,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idTair),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1462, MyFile)) RETURN
        END IF
!
!  Define surface winds.
!
        IF (Aout(idUair,ng)) THEN
          Vinfo( 1)=Vname(1,idUair)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUair))
          Vinfo( 3)=Vname(3,idUair)
          Vinfo(14)=Vname(4,idUair)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUair,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idUair),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1483, MyFile)) RETURN
        END IF
!
        IF (Aout(idVair,ng)) THEN
          Vinfo( 1)=Vname(1,idVair)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVair))
          Vinfo( 3)=Vname(3,idVair)
          Vinfo(14)=Vname(4,idVair)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVair,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idVair),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1500, MyFile)) RETURN
        END IF
!
!  Define Eastward/Northward surface wind at RHO-points.
!
        IF (Aout(idUaiE,ng)) THEN
          Vinfo( 1)=Vname(1,idUaiE)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUaiE))
          Vinfo( 3)=Vname(3,idUaiE)
          Vinfo(14)=Vname(4,idUaiE)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUaiE)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUaiE,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idUaiE),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1519, MyFile)) RETURN
        END IF
!
        IF (Aout(idVaiN,ng)) THEN
          Vinfo( 1)=Vname(1,idVaiN)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVaiN))
          Vinfo( 3)=Vname(3,idVaiN)
          Vinfo(14)=Vname(4,idVaiN)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVaiN)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVaiN,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idVaiN),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1536, MyFile)) RETURN
        END IF
!
!  Define surface net heat flux.
!
        IF (Aout(idTsur(itemp),ng)) THEN
          Vinfo( 1)=Vname(1,idTsur(itemp))
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                          &
     &                                 TRIM(Vname(2,idTsur(itemp)))
          Vinfo( 3)=Vname(3,idTsur(itemp))
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idTsur(itemp))
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idTsur(itemp))
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTsur(itemp),ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idTsur(itemp)), NF_FOUT,           &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1560, MyFile)) RETURN
        END IF
!
!  Define surface net salt flux.
!
        IF (Aout(idTsur(isalt),ng)) THEN
          Vinfo( 1)=Vname(1,idTsur(isalt))
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                          &
     &                                 TRIM(Vname(2,idTsur(isalt)))
          Vinfo( 3)=Vname(3,idTsur(isalt))
          Vinfo(11)='upward flux, freshening (net precipitation)'
          Vinfo(12)='downward flux, salting (net evaporation)'
          Vinfo(14)=Vname(4,idTsur(isalt))
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idTsur(isalt))
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTsur(isalt),ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idTsur(isalt)), NF_FOUT,           &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1585, MyFile)) RETURN
        END IF
!
!  Define latent heat flux.
!
        IF (Aout(idLhea,ng)) THEN
          Vinfo( 1)=Vname(1,idLhea)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idLhea))
          Vinfo( 3)=Vname(3,idLhea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLhea)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idLhea)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idLhea,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idLhea),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1609, MyFile)) RETURN
        END IF
!
!  Define sensible heat flux.
!
        IF (Aout(idShea,ng)) THEN
          Vinfo( 1)=Vname(1,idShea)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idShea))
          Vinfo( 3)=Vname(3,idShea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idShea)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idShea)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idShea,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idShea),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1630, MyFile)) RETURN
        END IF
!
!  Define net longwave radiation flux.
!
        IF (Aout(idLrad,ng)) THEN
          Vinfo( 1)=Vname(1,idLrad)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idLrad))
          Vinfo( 3)=Vname(3,idLrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLrad)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idLrad)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idLrad,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idLrad),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1651, MyFile)) RETURN
        END IF
!
!  Define evaporation rate.
!
        IF (Aout(idevap,ng)) THEN
          Vinfo( 1)=Vname(1,idevap)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idevap))
          Vinfo( 3)=Vname(3,idevap)
          Vinfo(11)='downward flux, freshening (condensation)'
          Vinfo(12)='upward flux, salting (evaporation)'
          Vinfo(14)=Vname(4,idevap)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idevap)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idevap,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idevap),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1676, MyFile)) RETURN
        END IF
!
!  Define precipitation rate.
!
        IF (Aout(idrain,ng)) THEN
          Vinfo( 1)=Vname(1,idrain)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idrain))
          Vinfo( 3)=Vname(3,idrain)
          Vinfo(12)='downward flux, freshening (precipitation)'
          Vinfo(14)=Vname(4,idrain)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idrain)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idrain,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idrain),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1696, MyFile)) RETURN
        END IF
!
!  Define net shortwave radiation flux.
!
        IF (Aout(idSrad,ng)) THEN
          Vinfo( 1)=Vname(1,idSrad)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idSrad))
          Vinfo( 3)=Vname(3,idSrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idSrad)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idSrad)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSrad,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idSrad),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1720, MyFile)) RETURN
        END IF
!
!  Define surface u-momentum stress.
!
        IF (Aout(idUsms,ng)) THEN
          Vinfo( 1)=Vname(1,idUsms)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUsms))
          Vinfo( 3)=Vname(3,idUsms)
          Vinfo(14)=Vname(4,idUsms)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUsms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUsms,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idUsms),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1741, MyFile)) RETURN
        END IF
!
!  Define surface v-momentum stress.
!
        IF (Aout(idVsms,ng)) THEN
          Vinfo( 1)=Vname(1,idVsms)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVsms))
          Vinfo( 3)=Vname(3,idVsms)
          Vinfo(14)=Vname(4,idVsms)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVsms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVsms,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idVsms),  &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1760, MyFile)) RETURN
        END IF
!
!  Define bottom u-momentum stress.
!
        IF (Aout(idUbms,ng)) THEN
          Vinfo( 1)=Vname(1,idUbms)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUbms))
          Vinfo( 3)=Vname(3,idUbms)
          Vinfo(14)=Vname(4,idUbms)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUbms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbms,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idUbms),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1779, MyFile)) RETURN
        END IF
!
!  Define bottom v-momentum stress.
!
        IF (Aout(idVbms,ng)) THEN
          Vinfo( 1)=Vname(1,idVbms)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVbms))
          Vinfo( 3)=Vname(3,idVbms)
          Vinfo(14)=Vname(4,idVbms)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVbms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbms,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idVbms),  &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 1798, MyFile)) RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, model, ncname, AVG(ng)%ncid)
        IF (FoundError(exit_flag, NoError, 1852, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, model, AVG(ng)%ncid, ncname)
        IF (FoundError(exit_flag, NoError, 1859, MyFile)) RETURN
      END IF DEFINE
!
!=======================================================================
!  Open an existing averages file, check its contents, and prepare
!  for appending data.
!=======================================================================
!
      QUERY : IF (.not.ldef) THEN
        ncname=AVG(ng)%name
!
!  Open averages file for read/write.
!
        CALL netcdf_open (ng, model, ncname, 1, AVG(ng)%ncid)
        IF (FoundError(exit_flag, NoError, 1873, MyFile)) THEN
          WRITE (stdout,50) TRIM(ncname)
          RETURN
        END IF
!
!  Inquire about the dimensions and check for consistency.
!
        CALL netcdf_check_dim (ng, model, ncname,                       &
     &                         ncid = AVG(ng)%ncid)
        IF (FoundError(exit_flag, NoError, 1882, MyFile)) RETURN
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, model, ncname,                         &
     &                       ncid = AVG(ng)%ncid)
        IF (FoundError(exit_flag, NoError, 1888, MyFile)) RETURN
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  average variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            AVG(ng)%Vid(idtime)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            AVG(ng)%Vid(idFsur)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            AVG(ng)%Vid(idUbar)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            AVG(ng)%Vid(idVbar)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu2dE))) THEN
            got_var(idu2dE)=.TRUE.
            AVG(ng)%Vid(idu2dE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv2dN))) THEN
            got_var(idv2dN)=.TRUE.
            AVG(ng)%Vid(idv2dN)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            AVG(ng)%Vid(idUvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            AVG(ng)%Vid(idVvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu3dE))) THEN
            got_var(idu3dE)=.TRUE.
            AVG(ng)%Vid(idu3dE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv3dN))) THEN
            got_var(idv3dN)=.TRUE.
            AVG(ng)%Vid(idv3dN)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOvel))) THEN
            got_var(idOvel)=.TRUE.
            AVG(ng)%Vid(idOvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWvel))) THEN
            got_var(idWvel)=.TRUE.
            AVG(ng)%Vid(idWvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            AVG(ng)%Vid(idDano)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,id2dPV))) THEN
            got_var(id2dPV)=.TRUE.
            AVG(ng)%Vid(id2dPV)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,id2dRV))) THEN
            got_var(id2dRV)=.TRUE.
            AVG(ng)%Vid(id2dRV)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,id3dPV))) THEN
            got_var(id3dPV)=.TRUE.
            AVG(ng)%Vid(id3dPV)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,id3dRV))) THEN
            got_var(id3dRV)=.TRUE.
            AVG(ng)%Vid(id3dRV)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idZZav))) THEN
            got_var(idZZav)=.TRUE.
            AVG(ng)%Vid(idZZav)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU2av))) THEN
            got_var(idU2av)=.TRUE.
            AVG(ng)%Vid(idU2av)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV2av))) THEN
            got_var(idV2av)=.TRUE.
            AVG(ng)%Vid(idV2av)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHUav))) THEN
            got_var(idHUav)=.TRUE.
            AVG(ng)%Vid(idHUav)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHVav))) THEN
            got_var(idHVav)=.TRUE.
            AVG(ng)%Vid(idHVav)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUUav))) THEN
            got_var(idUUav)=.TRUE.
            AVG(ng)%Vid(idUUav)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUVav))) THEN
            got_var(idUVav)=.TRUE.
            AVG(ng)%Vid(idUVav)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVVav))) THEN
            got_var(idVVav)=.TRUE.
            AVG(ng)%Vid(idVVav)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvis))) THEN
            got_var(idVvis)=.TRUE.
            AVG(ng)%Vid(idVvis)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTdif))) THEN
            got_var(idTdif)=.TRUE.
            AVG(ng)%Vid(idTdif)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSdif))) THEN
            got_var(idSdif)=.TRUE.
            AVG(ng)%Vid(idSdif)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPair))) THEN
            got_var(idPair)=.TRUE.
            AVG(ng)%Vid(idPair)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUair))) THEN
            got_var(idUair)=.TRUE.
            AVG(ng)%Vid(idUair)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVair))) THEN
            got_var(idVair)=.TRUE.
            AVG(ng)%Vid(idVair)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUaiE))) THEN
            got_var(idUaiE)=.TRUE.
            AVG(ng)%Vid(idUaiE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVaiN))) THEN
            got_var(idVaiN)=.TRUE.
            AVG(ng)%Vid(idVaiN)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idTsur(itemp)))) THEN
            got_var(idTsur(itemp))=.TRUE.
            AVG(ng)%Vid(idTsur(itemp))=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idTsur(isalt)))) THEN
            got_var(idTsur(isalt))=.TRUE.
            AVG(ng)%Vid(idTsur(isalt))=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLhea))) THEN
            got_var(idLhea)=.TRUE.
            AVG(ng)%Vid(idLhea)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idShea))) THEN
            got_var(idShea)=.TRUE.
            AVG(ng)%Vid(idShea)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLrad))) THEN
            got_var(idLrad)=.TRUE.
            AVG(ng)%Vid(idLrad)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTair))) THEN
            got_var(idTair)=.TRUE.
            AVG(ng)%Vid(idTair)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idevap))) THEN
            got_var(idevap)=.TRUE.
            AVG(ng)%Vid(idevap)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idrain))) THEN
            got_var(idrain)=.TRUE.
            AVG(ng)%Vid(idrain)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSrad))) THEN
            got_var(idSrad)=.TRUE.
            AVG(ng)%Vid(idSrad)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            AVG(ng)%Vid(idUsms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            AVG(ng)%Vid(idVsms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbms))) THEN
            got_var(idUbms)=.TRUE.
            AVG(ng)%Vid(idUbms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbms))) THEN
            got_var(idVbms)=.TRUE.
            AVG(ng)%Vid(idVbms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbrs))) THEN
            got_var(idUbrs)=.TRUE.
            AVG(ng)%Vid(idUbrs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbrs))) THEN
            got_var(idVbrs)=.TRUE.
            AVG(ng)%Vid(idVbrs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbws))) THEN
            got_var(idUbws)=.TRUE.
            AVG(ng)%Vid(idUbws)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbws))) THEN
            got_var(idVbws)=.TRUE.
            AVG(ng)%Vid(idVbws)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbcs))) THEN
            got_var(idUbcs)=.TRUE.
            AVG(ng)%Vid(idUbcs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbcs))) THEN
            got_var(idVbcs)=.TRUE.
            AVG(ng)%Vid(idVbcs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUVwc))) THEN
            got_var(idUVwc)=.TRUE.
            AVG(ng)%Vid(idUVwc)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbot))) THEN
            got_var(idUbot)=.TRUE.
            AVG(ng)%Vid(idUbot)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbot))) THEN
            got_var(idVbot)=.TRUE.
            AVG(ng)%Vid(idVbot)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbur))) THEN
            got_var(idUbur)=.TRUE.
            AVG(ng)%Vid(idUbur)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbvr))) THEN
            got_var(idVbvr)=.TRUE.
            AVG(ng)%Vid(idVbvr)=var_id(i)
          END IF
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
             got_var(idTvar(itrc))=.TRUE.
             AVG(ng)%Tid(itrc)=var_id(i)
            END IF
          END DO
          DO itrc=1,NAT
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,iHUTav(itrc)))) THEN
              got_var(iHUTav(itrc))=.TRUE.
              AVG(ng)%Vid(iHUTav(itrc))=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,iHVTav(itrc)))) THEN
              got_var(iHVTav(itrc))=.TRUE.
              AVG(ng)%Vid(iHVTav(itrc))=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idUTav(itrc)))) THEN
              got_var(idUTav(itrc))=.TRUE.
              AVG(ng)%Vid(idUTav(itrc))=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idVTav(itrc)))) THEN
              got_var(idVTav(itrc))=.TRUE.
              AVG(ng)%Vid(idVTav(itrc))=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idTTav(itrc)))) THEN
              got_var(idTTav(itrc))=.TRUE.
              AVG(ng)%Vid(idTTav(itrc))=var_id(i)
            END IF
          END DO
        END DO
!
!  Check if averages variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsur).and.Aout(idFsur,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar).and.Aout(idUbar,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar).and.Aout(idVbar,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idu2dE).and.Aout(idu2dE,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idu2dE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv2dN).and.Aout(idv2dN,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idv2dN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUvel).and.Aout(idUvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel).and.Aout(idVvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idu3dE).and.Aout(idu3dE,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idu3dE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv3dN).and.Aout(idv3dN,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idv3dN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idOvel).and.Aout(idOvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idOvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWvel).and.Aout(idWvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idDano).and.Aout(idDano,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idDano)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(id2dPV).and.Aout(id2dPV,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,id2dPV)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(id2dRV).and.Aout(id2dRV,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,id2dRV)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(id3dPV).and.Aout(id3dPV,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,id3dPV)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(id3dRV).and.Aout(id3dRV,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,id3dRV)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idZZav).and.Aout(idZZav,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idZZav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idU2av).and.Aout(idU2av,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idU2av)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV2av).and.Aout(idV2av,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idV2av)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idHUav).and.Aout(idHUav,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHUav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idHVav).and.Aout(idHVav,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHVav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUUav).and.Aout(idUUav,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUUav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUVav).and.Aout(idUVav,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUVav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVVav).and.Aout(idVVav,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVVav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvis).and.Aout(idVvis,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVvis)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTdif).and.Aout(idTdif,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSdif).and.Aout(idSdif,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idPair).and.Aout(idPair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idPair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUair).and.Aout(idUair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVair).and.Aout(idVair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUaiE).and.Aout(idUaiE,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUaiE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVaiN).and.Aout(idVaiN,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVaiN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTsur(itemp)).and.Aout(idTsur(itemp),ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTsur(itemp))),   &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTsur(isalt)).and.Aout(idTsur(isalt),ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTsur(isalt))),   &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLhea).and.Aout(idLhea,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idLhea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idShea).and.Aout(idShea,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idShea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLrad).and.Aout(idLrad,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idLrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTair).and.Aout(idTair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idevap).and.Aout(idevap,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idevap)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idrain).and.Aout(idrain,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idrain)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSrad).and.Aout(idSrad,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUsms).and.Aout(idUsms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms).and.Aout(idVsms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbms).and.Aout(idUbms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbms).and.Aout(idVbms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc)).and.Aout(idTvar(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO itrc=1,NAT
          IF (.not.got_var(iHUTav(itrc)).and.Aout(iHUTav(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,iHUTav(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(iHVTav(itrc)).and.Aout(iHVTav(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,iHVTav(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idUTav(itrc)).and.Aout(idUTav(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUTav(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idVTav(itrc)).and.Aout(idVTav(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVTav(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idTTav(itrc)).and.Aout(idTTav(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTTav(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
!
!  Set unlimited time record dimension to the appropriate value.
!
        IF (nRST(ng).eq.nAVG(ng)) THEN
          IF (ndefAVG(ng).gt.0) THEN
            AVG(ng)%Rindex=((ntstart(ng)-1)-                            &
     &                      ndefAVG(ng)*((ntstart(ng)-1)/ndefAVG(ng)))/ &
     &                     nAVG(ng)
          ELSE
            AVG(ng)%Rindex=(ntstart(ng)-1)/nAVG(ng)
          END IF
        ELSE
          AVG(ng)%Rindex=rec_size
        END IF
      END IF QUERY
!
!  Set initial average time. Notice that the value is offset by half
!  nAVG*dt so there is not a special case when computing its value
!  in "set_avg".
!
      IF (ntsAVG(ng).eq.1) THEN
        AVGtime(ng)=time(ng)-0.5_r8*REAL(nAVG(ng),r8)*dt(ng)
      ELSE
        AVGtime(ng)=time(ng)+REAL(ntsAVG(ng),r8)*dt(ng)-                &
     &              0.5_r8*REAL(nAVG(ng),r8)*dt(ng)
      END IF
!
  10  FORMAT (2x,'DEF_AVG_NF90     - creating average file,',t56,       &
     &        'Grid ',i2.2,': ',a)
  20  FORMAT (2x,'DEF_AVG_NF90     - inquiring average file,',t56,      &
     &        'Grid ',i2.2,': ',a)
  30  FORMAT (/,' DEF_AVG_NF90 - unable to create averages NetCDF',     &
     &        ' file: ',a)
  40  FORMAT (1pe11.4,1x,'millimeter')
  50  FORMAT (/,' DEF_AVG_NF90 - unable to open averages NetCDF',       &
     &        ' file: ',a)
  60  FORMAT (/,' DEF_AVG_NF90 - unable to find variable: ',a,2x,       &
     &        ' in averages NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE def_avg_nf90
!
!***********************************************************************
      SUBROUTINE def_avg_pio (ng, model, ldef)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
      logical, intent(in) :: ldef
!
!  Local variable declarations.
!
      logical :: got_var(NV)
!
      integer :: i, ifield, itrc, j, nvd3, nvd4
      integer :: recdim, status
      integer :: DimIDs(nDimID)
      integer :: p2dgrd(3), t2dgrd(3), u2dgrd(3), v2dgrd(3)
      integer :: p3dgrd(4), t3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)
!
      real(r8) :: Aval(6)
!
      character (len=13)     :: Prefix
      character (len=256)    :: ncname
      character (len=MaxLen) :: Vinfo(Natt)
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/def_avg.F"//", def_avg_pio"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, 2754, MyFile)) RETURN
      ncname=AVG(ng)%name
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
!  Create a new averages file.
!=======================================================================
!
      DEFINE : IF (ldef) THEN
        CALL pio_netcdf_create (ng, model, TRIM(ncname),                &
     &                          AVG(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, 2772, MyFile)) THEN
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
        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'xi_rho',    &
     &                 IOBOUNDS(ng)%xi_rho, DimIDs( 1))
        IF (FoundError(exit_flag, NoError, 2785, MyFile)) RETURN
        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'xi_u',      &
     &                 IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (FoundError(exit_flag, NoError, 2789, MyFile)) RETURN
        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'xi_v',      &
     &                 IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (FoundError(exit_flag, NoError, 2793, MyFile)) RETURN
        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'xi_psi',    &
     &                 IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (FoundError(exit_flag, NoError, 2797, MyFile)) RETURN
        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'eta_rho',   &
     &                 IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (FoundError(exit_flag, NoError, 2801, MyFile)) RETURN
        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'eta_u',     &
     &                 IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (FoundError(exit_flag, NoError, 2805, MyFile)) RETURN
        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'eta_v',     &
     &                 IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (FoundError(exit_flag, NoError, 2809, MyFile)) RETURN
        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'eta_psi',   &
     &                 IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (FoundError(exit_flag, NoError, 2813, MyFile)) RETURN
        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 's_rho',     &
     &                 N(ng), DimIDs( 9))
        IF (FoundError(exit_flag, NoError, 2858, MyFile)) RETURN
        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 's_w',       &
     &                 N(ng)+1, DimIDs(10))
        IF (FoundError(exit_flag, NoError, 2862, MyFile)) RETURN
        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'tracer',    &
     &                 NT(ng), DimIDs(11))
        IF (FoundError(exit_flag, NoError, 2866, MyFile)) RETURN
        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'boundary',  &
     &                 4, DimIDs(14))
        IF (FoundError(exit_flag, NoError, 2909, MyFile)) RETURN
        status=def_dim(ng, model, AVG(ng)%pioFile, ncname,              &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 PIO_unlimited, DimIDs(12))
        IF (FoundError(exit_flag, NoError, 2926, MyFile)) RETURN
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
!  Define dimension vectors for staggered variables at PSI-points.
!
        p2dgrd(1)=DimIDs( 4)
        p2dgrd(2)=DimIDs( 8)
        p2dgrd(3)=DimIDs(12)
        p3dgrd(1)=DimIDs( 4)
        p3dgrd(2)=DimIDs( 8)
        p3dgrd(3)=DimIDs( 9)
        p3dgrd(4)=DimIDs(12)
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
        AVG(ng)%Rindex=0
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
!  Set long name prefix string.
!
        Prefix=CHAR(32)                                   ! blank space
!
!-----------------------------------------------------------------------
!  Define time-recordless information variables.
!-----------------------------------------------------------------------
!
        CALL def_info (ng, model, AVG(ng)%pioFile, ncname, DimIDs)
        IF (FoundError(exit_flag, NoError, 3082, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Define time-varying variables.
!-----------------------------------------------------------------------
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        WRITE (Vinfo( 2),'(a,a)') 'averaged ', TRIM(Vname(2,idtime))
        WRITE (Vinfo( 3),'(a,a)') 'seconds since ', TRIM(Rclock%string)
        Vinfo( 4)=TRIM(Rclock%calendar)
        Vinfo(14)=Vname(4,idtime)
        Vinfo(21)=Vname(6,idtime)
        AVG(ng)%pioVar(idtime)%dkind=PIO_TOUT
        AVG(ng)%pioVar(idtime)%gtype=0
!
        status=def_var(ng, model, AVG(ng)%pioFile,                      &
     &                 AVG(ng)%pioVar(idtime)%vd,                       &
     &                 PIO_TOUT, 1, (/recdim/), Aval, Vinfo, ncname,    &
     &                 SetParAccess = .TRUE.)
        IF (FoundError(exit_flag, NoError, 3103, MyFile)) RETURN
!
!  Define free-surface.
!
        IF (Aout(idFsur,ng)) THEN
          Vinfo( 1)=Vname(1,idFsur)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idFsur))
          Vinfo( 3)=Vname(3,idFsur)
          Vinfo(14)=Vname(4,idFsur)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idFsur)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
          AVG(ng)%pioVar(idFsur)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idFsur)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idFsur)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3125, MyFile)) RETURN
        END IF
!
!  Define 2D momentum in the XI-direction.
!
        IF (Aout(idUbar,ng)) THEN
          Vinfo( 1)=Vname(1,idUbar)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUbar))
          Vinfo( 3)=Vname(3,idUbar)
          Vinfo(14)=Vname(4,idUbar)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUbar)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbar,ng),r8)
          AVG(ng)%pioVar(idUbar)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idUbar)%gtype=u2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idUbar)%vd,                     &
     &                   PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3174, MyFile)) RETURN
        END IF
!
!  Define 2D momentum in the ETA-direction.
!
        IF (Aout(idVbar,ng)) THEN
          Vinfo( 1)=Vname(1,idVbar)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVbar))
          Vinfo( 3)=Vname(3,idVbar)
          Vinfo(14)=Vname(4,idVbar)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVbar)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbar,ng),r8)
          AVG(ng)%pioVar(idVbar)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idVbar)%gtype=v2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idVbar)%vd,                     &
     &                   PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3223, MyFile)) RETURN
        END IF
!
!  Define 2D Eastward momentum component at RHO-points.
!
        IF (Aout(idu2dE,ng)) THEN
          Vinfo( 1)=Vname(1,idu2dE)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idu2dE))
          Vinfo( 3)=Vname(3,idu2dE)
          Vinfo(14)=Vname(4,idu2dE)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idu2dE)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu2dE,ng),r8)
          AVG(ng)%pioVar(idu2dE)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idu2dE)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idu2dE)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3272, MyFile)) RETURN
        END IF
!
!  Define 2D Northward momentum component at RHO-points.
!
        IF (Aout(idv2dN,ng)) THEN
          Vinfo( 1)=Vname(1,idv2dN)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idv2dN))
          Vinfo( 3)=Vname(3,idv2dN)
          Vinfo(14)=Vname(4,idv2dN)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idv2dN)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv2dN,ng),r8)
          AVG(ng)%pioVar(idv2dN)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idv2dN)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idv2dN)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3295, MyFile)) RETURN
        END IF
!
!  Define 3D momentum component in the XI-direction.
!
        IF (Aout(idUvel,ng)) THEN
          Vinfo( 1)=Vname(1,idUvel)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUvel))
          Vinfo( 3)=Vname(3,idUvel)
          Vinfo(14)=Vname(4,idUvel)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUvel,ng),r8)
          AVG(ng)%pioVar(idUvel)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idUvel)%gtype=u3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idUvel)%vd,                     &
     &                   PIO_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3408, MyFile)) RETURN
        END IF
!
!  Define 3D momentum component in the ETA-direction.
!
        IF (Aout(idVvel,ng)) THEN
          Vinfo( 1)=Vname(1,idVvel)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVvel))
          Vinfo( 3)=Vname(3,idVvel)
          Vinfo(14)=Vname(4,idVvel)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvel,ng),r8)
          AVG(ng)%pioVar(idVvel)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idVvel)%gtype=v3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idVvel)%vd,                     &
     &                   PIO_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3457, MyFile)) RETURN
        END IF
!
!  Define 3D Eastward momentum component at RHO-points.
!
        IF (Aout(idu3dE,ng)) THEN
          Vinfo( 1)=Vname(1,idu3dE)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idu3dE))
          Vinfo( 3)=Vname(3,idu3dE)
          Vinfo(14)=Vname(4,idu3dE)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idu3dE)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu3dE,ng),r8)
          AVG(ng)%pioVar(idu3dE)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idu3dE)%gtype=r3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idu3dE)%vd,                     &
     &                   PIO_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3506, MyFile)) RETURN
        END IF
!
!  Define 3D Northward momentum component at RHO-points.
!
        IF (Aout(idv3dN,ng)) THEN
          Vinfo( 1)=Vname(1,idv3dN)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idv3dN))
          Vinfo( 3)=Vname(3,idv3dN)
          Vinfo(14)=Vname(4,idv3dN)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idv3dN)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv3dN,ng),r8)
          AVG(ng)%pioVar(idv3dN)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idv3dN)%gtype=r3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idv3dN)%vd,                     &
     &                   PIO_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3529, MyFile)) RETURN
        END IF
!
!  Define S-coordinate vertical "omega" momentum component.
!
        IF (Aout(idOvel,ng)) THEN
          Vinfo( 1)=Vname(1,idOvel)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idOvel))
          Vinfo( 3)=Vname(3,idOvel)
          Vinfo(14)=Vname(4,idOvel)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idOvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idOvel,ng),r8)
          AVG(ng)%pioVar(idOvel)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idOvel)%gtype=w3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idOvel)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3552, MyFile)) RETURN
        END IF
!
!  Define "true" vertical momentum component.
!
        IF (Aout(idWvel,ng)) THEN
          Vinfo( 1)=Vname(1,idWvel)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWvel))
          Vinfo( 3)=Vname(3,idWvel)
          Vinfo(14)=Vname(4,idWvel)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idWvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWvel,ng),r8)
          AVG(ng)%pioVar(idWvel)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idWvel)%gtype=w3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idWvel)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3575, MyFile)) RETURN
        END IF
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (Aout(idTvar(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTvar(itrc))
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idTvar(itrc)))
            Vinfo( 3)=Vname(3,idTvar(itrc))
            Vinfo(14)=Vname(4,idTvar(itrc))
            Vinfo(16)=Vname(1,idtime)
            Vinfo(21)=Vname(6,idTvar(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            AVG(ng)%pioTrc(itrc)%dkind=PIO_FOUT
            AVG(ng)%pioTrc(itrc)%gtype=r3dvar
!
            status=def_var(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioTrc(itrc)%vd,                     &
     &                     PIO_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 3607, MyFile)) RETURN
          END IF
        END DO
!
!  Define density anomaly.
!
        IF (Aout(idDano,ng)) THEN
          Vinfo( 1)=Vname(1,idDano)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idDano))
          Vinfo( 3)=Vname(3,idDano)
          Vinfo(14)=Vname(4,idDano)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idDano)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idDano,ng),r8)
          AVG(ng)%pioVar(idDano)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idDano)%gtype=r3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idDano)%vd,                     &
     &                   PIO_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3660, MyFile)) RETURN
        END IF
!
!  Define 2D potential vorticity.
!
        IF (Aout(id2dPV,ng)) THEN
          Vinfo( 1)=Vname(1,id2dPV)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,id2dPV))
          Vinfo( 3)=Vname(3,id2dPV)
          Vinfo(14)=Vname(4,id2dPV)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,id2dPV)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(p2dvar,r8)
          AVG(ng)%pioVar(id2dPV)%dkind=PIO_FOUT
          AVG(ng)%pioVar(id2dPV)%gtype=p2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(id2dPV)%vd,                     &
     &                   PIO_FOUT, nvd3, p2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3735, MyFile)) RETURN
        END IF
!
!  Define 2D relative vorticity.
!
        IF (Aout(id2dRV,ng)) THEN
          Vinfo( 1)=Vname(1,id2dRV)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,id2dRV))
          Vinfo( 3)=Vname(3,id2dRV)
          Vinfo(14)=Vname(4,id2dRV)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,id2dRV)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(p2dvar,r8)
          AVG(ng)%pioVar(id2dRV)%dkind=PIO_FOUT
          AVG(ng)%pioVar(id2dRV)%gtype=p2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(id2dRV)%vd,                     &
     &                   PIO_FOUT, nvd3, p2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3758, MyFile)) RETURN
        END IF
!
!  Define 3D potential vorticity.
!
        IF (Aout(id3dPV,ng)) THEN
          Vinfo( 1)=Vname(1,id3dPV)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,id3dPV))
          Vinfo( 3)=Vname(3,id3dPV)
          Vinfo(14)=Vname(4,id3dPV)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,id3dPV)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(p3dvar,r8)
          AVG(ng)%pioVar(id3dPV)%dkind=PIO_FOUT
          AVG(ng)%pioVar(id3dPV)%gtype=p3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(id3dPV)%vd,                     &
     &                   PIO_FOUT, nvd4, p3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3783, MyFile)) RETURN
        END IF
!
!  Define 3D relative vorticity.
!
        IF (Aout(id3dRV,ng)) THEN
          Vinfo( 1)=Vname(1,id3dRV)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,id3dRV))
          Vinfo( 3)=Vname(3,id3dRV)
          Vinfo(14)=Vname(4,id3dRV)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,id3dRV)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(p3dvar,r8)
          AVG(ng)%pioVar(id3dRV)%dkind=PIO_FOUT
          AVG(ng)%pioVar(id3dRV)%gtype=p3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(id3dRV)%vd,                     &
     &                   PIO_FOUT, nvd4, p3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3806, MyFile)) RETURN
        END IF
!
!  Define quadratic <zeta*zeta> term.
!
        IF (Aout(idZZav,ng)) THEN
          Vinfo( 1)=Vname(1,idZZav)
          Vinfo( 2)=TRIM(Vname(2,idZZav))
          Vinfo( 3)=Vname(3,idZZav)
          Vinfo(14)=Vname(4,idZZav)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idZZav)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          AVG(ng)%pioVar(idZZav)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idZZav)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idZZav)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3830, MyFile)) RETURN
        END IF
!
!  Define quadratic <ubar*ubar> term.
!
        IF (Aout(idU2av,ng)) THEN
          Vinfo( 1)=Vname(1,idU2av)
          Vinfo( 2)=TRIM(Vname(2,idU2av))
          Vinfo( 3)=Vname(3,idU2av)
          Vinfo(14)=Vname(4,idU2av)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idU2av)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          AVG(ng)%pioVar(idU2av)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idU2av)%gtype=u2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idU2av)%vd,                     &
     &                   PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3853, MyFile)) RETURN
        END IF
!
!  Define quadratic <vbar*vbar> term.
!
        IF (Aout(idV2av,ng)) THEN
          Vinfo( 1)=Vname(1,idV2av)
          Vinfo( 2)=TRIM(Vname(2,idV2av))
          Vinfo( 3)=Vname(3,idV2av)
          Vinfo(14)=Vname(4,idV2av)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idV2av)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          AVG(ng)%pioVar(idV2av)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idV2av)%gtype=v2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idV2av)%vd,                     &
     &                   PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3876, MyFile)) RETURN
        END IF
!
!  Define u-volume flux.
!
        IF (Aout(idHUav,ng)) THEN
          Vinfo( 1)=Vname(1,idHUav)
          Vinfo( 2)=TRIM(Vname(2,idHUav))
          Vinfo( 3)=Vname(3,idHUav)
          Vinfo(14)=Vname(4,idHUav)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idHUav)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          AVG(ng)%pioVar(idHUav)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idHUav)%gtype=u3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idHUav)%vd,                     &
     &                   PIO_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3901, MyFile)) RETURN
        END IF
!
!  Define v-volume flux.
!
        IF (Aout(idHVav,ng)) THEN
          Vinfo( 1)=Vname(1,idHVav)
          Vinfo( 2)=TRIM(Vname(2,idHVav))
          Vinfo( 3)=Vname(3,idHVav)
          Vinfo(14)=Vname(4,idHVav)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idHVav)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v3dvar,r8)
          AVG(ng)%pioVar(idHVav)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idHVav)%gtype=v3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idHVav)%vd,                     &
     &                   PIO_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3924, MyFile)) RETURN
        END IF
!
!  Define quadratic <u*u> term.
!
        IF (Aout(idUUav,ng)) THEN
          Vinfo( 1)=Vname(1,idUUav)
          Vinfo( 2)=TRIM(Vname(2,idUUav))
          Vinfo( 3)=Vname(3,idUUav)
          Vinfo(14)=Vname(4,idUUav)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUUav)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          AVG(ng)%pioVar(idUUav)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idUUav)%gtype=u3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idUUav)%vd,                     &
     &                   PIO_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3947, MyFile)) RETURN
        END IF
!
!  Define quadratic <u*v> term.
!
        IF (Aout(idUVav,ng)) THEN
          Vinfo( 1)=Vname(1,idUVav)
          Vinfo( 2)=TRIM(Vname(2,idUVav))
          Vinfo( 3)=Vname(3,idUVav)
          Vinfo(14)=Vname(4,idUVav)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUVav)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dvar,r8)
          AVG(ng)%pioVar(idUVav)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idUVav)%gtype=r3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idUVav)%vd,                     &
     &                   PIO_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3970, MyFile)) RETURN
        END IF
!
!  Define quadratic <v*v> term.
!
        IF (Aout(idVVav,ng)) THEN
          Vinfo( 1)=Vname(1,idVVav)
          Vinfo( 2)=TRIM(Vname(2,idVVav))
          Vinfo( 3)=Vname(3,idVVav)
          Vinfo(14)=Vname(4,idVVav)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVVav)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v3dvar,r8)
          AVG(ng)%pioVar(idVVav)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idVVav)%gtype=v3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idVVav)%vd,                     &
     &                   PIO_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 3993, MyFile)) RETURN
        END IF
!
!  Define quadratic <t*t> terms.
!
        DO itrc=1,NT(ng)
          IF (Aout(idTTav(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTTav(itrc))
            Vinfo( 2)=TRIM(Vname(2,idTTav(itrc)))
            Vinfo( 3)=Vname(3,idTTav(itrc))
            Vinfo(14)=Vname(4,idTTav(itrc))
            Vinfo(16)=Vname(1,idtime)
            Vinfo(21)=Vname(6,idTTav(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            AVG(ng)%pioVar(idTTav(itrc))%dkind=PIO_FOUT
            AVG(ng)%pioVar(idTTav(itrc))%gtype=r3dvar
!
            status=def_var(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idTTav(itrc))%vd,             &
     &                     PIO_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 4017, MyFile)) RETURN
          END IF
        END DO
!
!  Define active tracers volume fluxes.
!
        DO itrc=1,NT(ng)
          IF (Aout(iHUTav(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,iHUTav(itrc))
            Vinfo( 2)=TRIM(Vname(2,iHUTav(itrc)))
            Vinfo( 3)=Vname(3,iHUTav(itrc))
            Vinfo(14)=Vname(4,iHUTav(itrc))
            Vinfo(16)=Vname(1,idtime)
            Vinfo(21)=Vname(6,iHUTav(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(u3dvar,r8)
            AVG(ng)%pioVar(iHUTav(itrc))%dkind=PIO_FOUT
            AVG(ng)%pioVar(iHUTav(itrc))%gtype=u3dvar
!
            status=def_var(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(iHUTav(itrc))%vd,             &
     &                     PIO_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 4042, MyFile)) RETURN
          END IF
!
          IF (Aout(iHVTav(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,iHVTav(itrc))
            Vinfo( 2)=TRIM(Vname(2,iHVTav(itrc)))
            Vinfo( 3)=Vname(3,iHVTav(itrc))
            Vinfo(14)=Vname(4,iHVTav(itrc))
            Vinfo(16)=Vname(1,idtime)
            Vinfo(21)=Vname(6,iHVTav(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(v3dvar,r8)
            AVG(ng)%pioVar(iHVTav(itrc))%dkind=PIO_FOUT
            AVG(ng)%pioVar(iHVTav(itrc))%gtype=v3dvar
!
            status=def_var(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(iHVTav(itrc))%vd,             &
     &                     PIO_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 4063, MyFile)) RETURN
          END IF
        END DO
!
!  Define quadratic <u*t> and <v*t> terms.
!
        DO itrc=1,NT(ng)
          IF (Aout(idUTav(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idUTav(itrc))
            Vinfo( 2)=TRIM(Vname(2,idUTav(itrc)))
            Vinfo( 3)=Vname(3,idUTav(itrc))
            Vinfo(14)=Vname(4,idUTav(itrc))
            Vinfo(16)=Vname(1,idtime)
            Vinfo(21)=Vname(6,idUTav(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(u3dvar,r8)
            AVG(ng)%pioVar(idUTav(itrc))%dkind=PIO_FOUT
            AVG(ng)%pioVar(idUTav(itrc))%gtype=u3dvar
!
            status=def_var(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idUTav(itrc))%vd,             &
     &                     PIO_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 4088, MyFile)) RETURN
          END IF
!
          IF (Aout(idVTav(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idVTav(itrc))
            Vinfo( 2)=TRIM(Vname(2,idVTav(itrc)))
            Vinfo( 3)=Vname(3,idVTav(itrc))
            Vinfo(14)=Vname(4,idVTav(itrc))
            Vinfo(16)=Vname(1,idtime)
            Vinfo(21)=Vname(6,idVTav(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(v3dvar,r8)
            AVG(ng)%pioVar(idVTav(itrc))%dkind=PIO_FOUT
            AVG(ng)%pioVar(idVTav(itrc))%gtype=v3dvar
!
            status=def_var(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idVTav(itrc))%vd,             &
     &                     PIO_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, 4109, MyFile)) RETURN
          END IF
        END DO
!
!  Define vertical viscosity coefficient.
!
        IF (Aout(idVvis,ng)) THEN
          Vinfo( 1)=Vname(1,idVvis)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVvis))
          Vinfo( 3)=Vname(3,idVvis)
          Vinfo(14)=Vname(4,idVvis)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVvis)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvis,ng),r8)
          AVG(ng)%pioVar(idVvis)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idVvis)%gtype=w3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idVvis)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, 4133, MyFile)) RETURN
        END IF
!
!  Define vertical diffusion coefficient for potential temperature.
!
        IF (Aout(idTdif,ng)) THEN
          Vinfo( 1)=Vname(1,idTdif)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idTdif))
          Vinfo( 3)=Vname(3,idTdif)
          Vinfo(14)=Vname(4,idTdif)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idTdif)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTdif,ng),r8)
          AVG(ng)%pioVar(idTdif)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idTdif)%gtype=w3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idTdif)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, 4154, MyFile)) RETURN
        END IF
!
!  Define vertical diffusion coefficient for salinity.
!
        IF (Aout(idSdif,ng)) THEN
          Vinfo( 1)=Vname(1,idSdif)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idSdif))
          Vinfo( 3)=Vname(3,idSdif)
          Vinfo(14)=Vname(4,idSdif)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idSdif)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSdif,ng),r8)
          AVG(ng)%pioVar(idSdif)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idSdif)%gtype=w3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idSdif)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, 4177, MyFile)) RETURN
        END IF
!
!  Define surface air pressure.
!
        IF (Aout(idPair,ng)) THEN
          Vinfo( 1)=Vname(1,idPair)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idPair))
          Vinfo( 3)=Vname(3,idPair)
          Vinfo(14)=Vname(4,idPair)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idPair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idPair,ng),r8)
          AVG(ng)%pioVar(idPair)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idPair)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idPair)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4203, MyFile)) RETURN
        END IF
!
!  Define surface air temperature.
!
        IF (Aout(idTair,ng)) THEN
          Vinfo( 1)=Vname(1,idTair)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idTair))
          Vinfo( 3)=Vname(3,idTair)
          Vinfo(14)=Vname(4,idTair)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idTair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTair,ng),r8)
          AVG(ng)%pioVar(idTair)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idTair)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idTair)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4228, MyFile)) RETURN
        END IF
!
!  Define surface winds.
!
        IF (Aout(idUair,ng)) THEN
          Vinfo( 1)=Vname(1,idUair)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUair))
          Vinfo( 3)=Vname(3,idUair)
          Vinfo(14)=Vname(4,idUair)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUair,ng),r8)
          AVG(ng)%pioVar(idUair)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idUair)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idUair)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4253, MyFile)) RETURN
        END IF
!
        IF (Aout(idVair,ng)) THEN
          Vinfo( 1)=Vname(1,idVair)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVair))
          Vinfo( 3)=Vname(3,idVair)
          Vinfo(14)=Vname(4,idVair)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVair,ng),r8)
          AVG(ng)%pioVar(idVair)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idVair)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idVair)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4274, MyFile)) RETURN
        END IF
!
!  Define Eastward/Northward surface wind at RHO-points.
!
        IF (Aout(idUaiE,ng)) THEN
          Vinfo( 1)=Vname(1,idUaiE)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUaiE))
          Vinfo( 3)=Vname(3,idUaiE)
          Vinfo(14)=Vname(4,idUaiE)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUaiE)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUaiE,ng),r8)
          AVG(ng)%pioVar(idUaiE)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idUaiE)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idUaiE)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4297, MyFile)) RETURN
        END IF
!
        IF (Aout(idVaiN,ng)) THEN
          Vinfo( 1)=Vname(1,idVaiN)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVaiN))
          Vinfo( 3)=Vname(3,idVaiN)
          Vinfo(14)=Vname(4,idVaiN)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVaiN)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVaiN,ng),r8)
          AVG(ng)%pioVar(idVaiN)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idVaiN)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idVaiN)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4318, MyFile)) RETURN
        END IF
!
!  Define surface net heat flux.
!
        IF (Aout(idTsur(itemp),ng)) THEN
          Vinfo( 1)=Vname(1,idTsur(itemp))
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                          &
     &                                 TRIM(Vname(2,idTsur(itemp)))
          Vinfo( 3)=Vname(3,idTsur(itemp))
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idTsur(itemp))
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idTsur(itemp))
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTsur(itemp),ng),r8)
          AVG(ng)%pioVar(idTsur(itemp))%dkind=PIO_FOUT
          AVG(ng)%pioVar(idTsur(itemp))%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idTsur(itemp))%vd,              &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4345, MyFile)) RETURN
        END IF
!
!  Define surface net salt flux.
!
        IF (Aout(idTsur(isalt),ng)) THEN
          Vinfo( 1)=Vname(1,idTsur(isalt))
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                          &
     &                                 TRIM(Vname(2,idTsur(isalt)))
          Vinfo( 3)=Vname(3,idTsur(isalt))
          Vinfo(11)='upward flux, freshening (net precipitation)'
          Vinfo(12)='downward flux, salting (net evaporation)'
          Vinfo(14)=Vname(4,idTsur(isalt))
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idTsur(isalt))
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTsur(isalt),ng),r8)
          AVG(ng)%pioVar(idTsur(isalt))%dkind=PIO_FOUT
          AVG(ng)%pioVar(idTsur(isalt))%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idTsur(isalt))%vd,              &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4373, MyFile)) RETURN
        END IF
!
!  Define latent heat flux.
!
        IF (Aout(idLhea,ng)) THEN
          Vinfo( 1)=Vname(1,idLhea)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idLhea))
          Vinfo( 3)=Vname(3,idLhea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLhea)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idLhea)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idLhea,ng),r8)
          AVG(ng)%pioVar(idLhea)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idLhea)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idLhea)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4401, MyFile)) RETURN
        END IF
!
!  Define sensible heat flux.
!
        IF (Aout(idShea,ng)) THEN
          Vinfo( 1)=Vname(1,idShea)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idShea))
          Vinfo( 3)=Vname(3,idShea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idShea)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idShea)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idShea,ng),r8)
          AVG(ng)%pioVar(idShea)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idShea)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idShea)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4426, MyFile)) RETURN
        END IF
!
!  Define net longwave radiation flux.
!
        IF (Aout(idLrad,ng)) THEN
          Vinfo( 1)=Vname(1,idLrad)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idLrad))
          Vinfo( 3)=Vname(3,idLrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLrad)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idLrad)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idLrad,ng),r8)
          AVG(ng)%pioVar(idLrad)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idLrad)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idLrad)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4451, MyFile)) RETURN
        END IF
!
!  Define evaporation rate.
!
        IF (Aout(idevap,ng)) THEN
          Vinfo( 1)=Vname(1,idevap)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idevap))
          Vinfo( 3)=Vname(3,idevap)
          Vinfo(11)='downward flux, freshening (condensation)'
          Vinfo(12)='upward flux, salting (evaporation)'
          Vinfo(14)=Vname(4,idevap)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idevap)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idevap,ng),r8)
          AVG(ng)%pioVar(idevap)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idevap)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idevap)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4480, MyFile)) RETURN
        END IF
!
!  Define precipitation rate.
!
        IF (Aout(idrain,ng)) THEN
          Vinfo( 1)=Vname(1,idrain)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idrain))
          Vinfo( 3)=Vname(3,idrain)
          Vinfo(12)='downward flux, freshening (precipitation)'
          Vinfo(14)=Vname(4,idrain)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idrain)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idrain,ng),r8)
          AVG(ng)%pioVar(idrain)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idrain)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idrain)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4504, MyFile)) RETURN
        END IF
!
!  Define net shortwave radiation flux.
!
        IF (Aout(idSrad,ng)) THEN
          Vinfo( 1)=Vname(1,idSrad)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idSrad))
          Vinfo( 3)=Vname(3,idSrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idSrad)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idSrad)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSrad,ng),r8)
          AVG(ng)%pioVar(idSrad)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idSrad)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idSrad)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4532, MyFile)) RETURN
        END IF
!
!  Define surface u-momentum stress.
!
        IF (Aout(idUsms,ng)) THEN
          Vinfo( 1)=Vname(1,idUsms)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUsms))
          Vinfo( 3)=Vname(3,idUsms)
          Vinfo(14)=Vname(4,idUsms)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUsms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUsms,ng),r8)
          AVG(ng)%pioVar(idUsms)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idUsms)%gtype=u2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idUsms)%vd,                     &
     &                   PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4557, MyFile)) RETURN
        END IF
!
!  Define surface v-momentum stress.
!
        IF (Aout(idVsms,ng)) THEN
          Vinfo( 1)=Vname(1,idVsms)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVsms))
          Vinfo( 3)=Vname(3,idVsms)
          Vinfo(14)=Vname(4,idVsms)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVsms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVsms,ng),r8)
          AVG(ng)%pioVar(idVsms)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idVsms)%gtype=v2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idVsms)%vd,                     &
     &                   PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4580, MyFile)) RETURN
        END IF
!
!  Define bottom u-momentum stress.
!
        IF (Aout(idUbms,ng)) THEN
          Vinfo( 1)=Vname(1,idUbms)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUbms))
          Vinfo( 3)=Vname(3,idUbms)
          Vinfo(14)=Vname(4,idUbms)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUbms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbms,ng),r8)
          AVG(ng)%pioVar(idUbms)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idUbms)%gtype=u2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idUbms)%vd,                     &
     &                   PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4603, MyFile)) RETURN
        END IF
!
!  Define bottom v-momentum stress.
!
        IF (Aout(idVbms,ng)) THEN
          Vinfo( 1)=Vname(1,idVbms)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVbms))
          Vinfo( 3)=Vname(3,idVbms)
          Vinfo(14)=Vname(4,idVbms)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVbms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbms,ng),r8)
          AVG(ng)%pioVar(idVbms)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idVbms)%gtype=v2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idVbms)%vd,                     &
     &                   PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, 4626, MyFile)) RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL pio_netcdf_enddef (ng, model, ncname, AVG(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, 4680, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, model, AVG(ng)%pioFile, ncname)
        IF (FoundError(exit_flag, NoError, 4687, MyFile)) RETURN
      END IF DEFINE
!
!=======================================================================
!  Open an existing averages file, check its contents, and prepare
!  for appending data.
!=======================================================================
!
      QUERY : IF (.not.ldef) THEN
        ncname=AVG(ng)%name
!
!  Open averages file for read/write.
!
        CALL pio_netcdf_open (ng, model, ncname, 1, AVG(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, 4701, MyFile)) THEN
          WRITE (stdout,50) TRIM(ncname)
          RETURN
        END IF
!
!  Inquire about the dimensions and check for consistency.
!
        CALL pio_netcdf_check_dim (ng, model, ncname,                   &
     &                             pioFile = AVG(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, 4710, MyFile)) RETURN
!
!  Inquire about the variables.
!
        CALL pio_netcdf_inq_var (ng, model, ncname,                     &
     &                           pioFile = AVG(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, 4716, MyFile)) RETURN
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  average variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            AVG(ng)%pioVar(idtime)%vd=var_desc(i)
            AVG(ng)%pioVar(idtime)%dkind=PIO_TOUT
            AVG(ng)%pioVar(idtime)%gtype=0
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            AVG(ng)%pioVar(idFsur)%vd=var_desc(i)
            AVG(ng)%pioVar(idFsur)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idFsur)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            AVG(ng)%pioVar(idUbar)%vd=var_desc(i)
            AVG(ng)%pioVar(idUbar)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idUbar)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            AVG(ng)%pioVar(idVbar)%vd=var_desc(i)
            AVG(ng)%pioVar(idVbar)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idVbar)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu2dE))) THEN
            got_var(idu2dE)=.TRUE.
            AVG(ng)%pioVar(idu2dE)%vd=var_desc(i)
            AVG(ng)%pioVar(idu2dE)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idu2dE)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv2dN))) THEN
            got_var(idv2dN)=.TRUE.
            AVG(ng)%pioVar(idv2dN)%vd=var_desc(i)
            AVG(ng)%pioVar(idv2dN)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idv2dN)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            AVG(ng)%pioVar(idUvel)%vd=var_desc(i)
            AVG(ng)%pioVar(idUvel)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idUvel)%gtype=u3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            AVG(ng)%pioVar(idVvel)%vd=var_desc(i)
            AVG(ng)%pioVar(idVvel)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idVvel)%gtype=v3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu3dE))) THEN
            got_var(idu3dE)=.TRUE.
            AVG(ng)%pioVar(idu3dE)%vd=var_desc(i)
            AVG(ng)%pioVar(idu3dE)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idu3dE)%gtype=r3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv3dN))) THEN
            got_var(idv3dN)=.TRUE.
            AVG(ng)%pioVar(idv3dN)%vd=var_desc(i)
            AVG(ng)%pioVar(idv3dN)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idv3dN)%gtype=r3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOvel))) THEN
            got_var(idOvel)=.TRUE.
            AVG(ng)%pioVar(idOvel)%vd=var_desc(i)
            AVG(ng)%pioVar(idOvel)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idOvel)%gtype=w3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWvel))) THEN
            got_var(idWvel)=.TRUE.
            AVG(ng)%pioVar(idWvel)%vd=var_desc(i)
            AVG(ng)%pioVar(idWvel)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idWvel)%gtype=w3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            AVG(ng)%pioVar(idDano)%vd=var_desc(i)
            AVG(ng)%pioVar(idDano)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idDano)%gtype=r3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,id2dPV))) THEN
            got_var(id2dPV)=.TRUE.
            AVG(ng)%pioVar(id2dPV)%vd=var_desc(i)
            AVG(ng)%pioVar(id2dPV)%dkind=PIO_FOUT
            AVG(ng)%pioVar(id2dPV)%gtype=p2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,id2dRV))) THEN
            got_var(id2dRV)=.TRUE.
            AVG(ng)%pioVar(id2dRV)%vd=var_desc(i)
            AVG(ng)%pioVar(id2dRV)%dkind=PIO_FOUT
            AVG(ng)%pioVar(id2dRV)%gtype=p2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,id3dPV))) THEN
            got_var(id3dPV)=.TRUE.
            AVG(ng)%pioVar(id3dPV)%vd=var_desc(i)
            AVG(ng)%pioVar(id3dPV)%dkind=PIO_FOUT
            AVG(ng)%pioVar(id3dPV)%gtype=p3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,id3dRV))) THEN
            got_var(id3dRV)=.TRUE.
            AVG(ng)%pioVar(id3dRV)%vd=var_desc(i)
            AVG(ng)%pioVar(id3dRV)%dkind=PIO_FOUT
            AVG(ng)%pioVar(id3dRV)%gtype=p3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idZZav))) THEN
            got_var(idZZav)=.TRUE.
            AVG(ng)%pioVar(idZZav)%vd=var_desc(i)
            AVG(ng)%pioVar(idZZav)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idZZav)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU2av))) THEN
            got_var(idU2av)=.TRUE.
            AVG(ng)%pioVar(idU2av)%vd=var_desc(i)
            AVG(ng)%pioVar(idU2av)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idU2av)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV2av))) THEN
            got_var(idV2av)=.TRUE.
            AVG(ng)%pioVar(idV2av)%vd=var_desc(i)
            AVG(ng)%pioVar(idV2av)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idV2av)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHUav))) THEN
            got_var(idHUav)=.TRUE.
            AVG(ng)%pioVar(idHUav)%vd=var_desc(i)
            AVG(ng)%pioVar(idHUav)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idHUav)%gtype=u3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHVav))) THEN
            got_var(idHVav)=.TRUE.
            AVG(ng)%pioVar(idHVav)%vd=var_desc(i)
            AVG(ng)%pioVar(idHVav)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idHVav)%gtype=v3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUUav))) THEN
            got_var(idUUav)=.TRUE.
            AVG(ng)%pioVar(idUUav)%vd=var_desc(i)
            AVG(ng)%pioVar(idHVav)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idHVav)%gtype=v3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUVav))) THEN
            got_var(idUVav)=.TRUE.
            AVG(ng)%pioVar(idUVav)%vd=var_desc(i)
            AVG(ng)%pioVar(idUVav)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idUVav)%gtype=r3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVVav))) THEN
            got_var(idVVav)=.TRUE.
            AVG(ng)%pioVar(idVVav)%vd=var_desc(i)
            AVG(ng)%pioVar(idVVav)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idVVav)%gtype=v3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvis))) THEN
            got_var(idVvis)=.TRUE.
            AVG(ng)%pioVar(idVvis)%vd=var_desc(i)
            AVG(ng)%pioVar(idVvis)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idVvis)%gtype=w3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTdif))) THEN
            got_var(idTdif)=.TRUE.
            AVG(ng)%pioVar(idTdif)%vd=var_desc(i)
            AVG(ng)%pioVar(idTdif)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idTdif)%gtype=w3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSdif))) THEN
            got_var(idSdif)=.TRUE.
            AVG(ng)%pioVar(idSdif)%vd=var_desc(i)
            AVG(ng)%pioVar(idSdif)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idSdif)%gtype=w3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPair))) THEN
            got_var(idPair)=.TRUE.
            AVG(ng)%pioVar(idPair)%vd=var_desc(i)
            AVG(ng)%pioVar(idPair)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idPair)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUair))) THEN
            got_var(idUair)=.TRUE.
            AVG(ng)%pioVar(idUair)%vd=var_desc(i)
            AVG(ng)%pioVar(idUair)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idUair)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVair))) THEN
            got_var(idVair)=.TRUE.
            AVG(ng)%pioVar(idVair)%vd=var_desc(i)
            AVG(ng)%pioVar(idVair)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idVair)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUaiE))) THEN
            got_var(idUaiE)=.TRUE.
            AVG(ng)%pioVar(idUaiE)%vd=var_desc(i)
            AVG(ng)%pioVar(idUaiE)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idUaiE)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVaiN))) THEN
            got_var(idVaiN)=.TRUE.
            AVG(ng)%pioVar(idVaiN)%vd=var_desc(i)
            AVG(ng)%pioVar(idVaiN)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idVaiN)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idTsur(itemp)))) THEN
            got_var(idTsur(itemp))=.TRUE.
            AVG(ng)%pioVar(idTsur(itemp))%vd=var_desc(i)
            AVG(ng)%pioVar(idTsur(itemp))%dkind=PIO_FOUT
            AVG(ng)%pioVar(idTsur(itemp))%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idTsur(isalt)))) THEN
            got_var(idTsur(isalt))=.TRUE.
            AVG(ng)%pioVar(idTsur(isalt))%vd=var_desc(i)
            AVG(ng)%pioVar(idTsur(isalt))%dkind=PIO_FOUT
            AVG(ng)%pioVar(idTsur(isalt))%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLhea))) THEN
            got_var(idLhea)=.TRUE.
            AVG(ng)%pioVar(idLhea)%vd=var_desc(i)
            AVG(ng)%pioVar(idLhea)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idLhea)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idShea))) THEN
            got_var(idShea)=.TRUE.
            AVG(ng)%pioVar(idShea)%vd=var_desc(i)
            AVG(ng)%pioVar(idShea)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idShea)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLrad))) THEN
            got_var(idLrad)=.TRUE.
            AVG(ng)%pioVar(idLrad)%vd=var_desc(i)
            AVG(ng)%pioVar(idLrad)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idLrad)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTair))) THEN
            got_var(idTair)=.TRUE.
            AVG(ng)%pioVar(idTair)%vd=var_desc(i)
            AVG(ng)%pioVar(idTair)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idTair)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idevap))) THEN
            got_var(idevap)=.TRUE.
            AVG(ng)%pioVar(idevap)%vd=var_desc(i)
            AVG(ng)%pioVar(idevap)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idevap)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idrain))) THEN
            got_var(idrain)=.TRUE.
            AVG(ng)%pioVar(idrain)%vd=var_desc(i)
            AVG(ng)%pioVar(idrain)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idrain)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSrad))) THEN
            got_var(idSrad)=.TRUE.
            AVG(ng)%pioVar(idSrad)%vd=var_desc(i)
            AVG(ng)%pioVar(idSrad)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idSrad)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            AVG(ng)%pioVar(idUsms)%vd=var_desc(i)
            AVG(ng)%pioVar(idUsms)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idUsms)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            AVG(ng)%pioVar(idVsms)%vd=var_desc(i)
            AVG(ng)%pioVar(idVsms)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idVsms)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbms))) THEN
            got_var(idUbms)=.TRUE.
            AVG(ng)%pioVar(idUbms)%vd=var_desc(i)
            AVG(ng)%pioVar(idUbms)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idUbms)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbms))) THEN
            got_var(idVbms)=.TRUE.
            AVG(ng)%pioVar(idVbms)%vd=var_desc(i)
            AVG(ng)%pioVar(idVbms)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idVbms)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbrs))) THEN
            got_var(idUbrs)=.TRUE.
            AVG(ng)%pioVar(idUbrs)%vd=var_desc(i)
            AVG(ng)%pioVar(idUbrs)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idUbrs)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbrs))) THEN
            got_var(idVbrs)=.TRUE.
            AVG(ng)%pioVar(idVbrs)%vd=var_desc(i)
            AVG(ng)%pioVar(idVbrs)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idVbrs)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbws))) THEN
            got_var(idUbws)=.TRUE.
            AVG(ng)%pioVar(idUbws)%vd=var_desc(i)
            AVG(ng)%pioVar(idUbws)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idUbws)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbws))) THEN
            got_var(idVbws)=.TRUE.
            AVG(ng)%pioVar(idVbws)%vd=var_desc(i)
            AVG(ng)%pioVar(idVbws)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idVbws)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbcs))) THEN
            got_var(idUbcs)=.TRUE.
            AVG(ng)%pioVar(idUbcs)%vd=var_desc(i)
            AVG(ng)%pioVar(idUbcs)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idUbcs)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbcs))) THEN
            got_var(idVbcs)=.TRUE.
            AVG(ng)%pioVar(idVbcs)%vd=var_desc(i)
            AVG(ng)%pioVar(idVbcs)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idVbcs)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUVwc))) THEN
            got_var(idUVwc)=.TRUE.
            AVG(ng)%pioVar(idUVwc)%vd=var_desc(i)
            AVG(ng)%pioVar(idUVwc)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idUVwc)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbot))) THEN
            got_var(idUbot)=.TRUE.
            AVG(ng)%pioVar(idUbot)%vd=var_desc(i)
            AVG(ng)%pioVar(idUbot)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idUbot)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbot))) THEN
            got_var(idVbot)=.TRUE.
            AVG(ng)%pioVar(idVbot)%vd=var_desc(i)
            AVG(ng)%pioVar(idVbot)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idVbot)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbur))) THEN
            got_var(idUbur)=.TRUE.
            AVG(ng)%pioVar(idUbur)%vd=var_desc(i)
            AVG(ng)%pioVar(idUbur)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idUbur)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbvr))) THEN
            got_var(idVbvr)=.TRUE.
            AVG(ng)%pioVar(idVbvr)%vd=var_desc(i)
            AVG(ng)%pioVar(idVbvr)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idVbvr)%gtype=v2dvar
          END IF
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              AVG(ng)%pioTrc(itrc)%vd=var_desc(i)
              AVG(ng)%pioTrc(itrc)%dkind=PIO_FOUT
              AVG(ng)%pioTrc(itrc)%gtype=r3dvar
            END IF
          END DO
          DO itrc=1,NAT
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,iHUTav(itrc)))) THEN
              got_var(iHUTav(itrc))=.TRUE.
              AVG(ng)%pioVar(iHUTav(itrc))%vd=var_desc(i)
              AVG(ng)%pioVar(iHUTav(itrc))%dkind=PIO_FOUT
              AVG(ng)%pioVar(iHUTav(itrc))%gtype=u3dvar
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,iHVTav(itrc)))) THEN
              got_var(iHVTav(itrc))=.TRUE.
              AVG(ng)%pioVar(iHVTav(itrc))%vd=var_desc(i)
              AVG(ng)%pioVar(iHVTav(itrc))%dkind=PIO_FOUT
              AVG(ng)%pioVar(iHVTav(itrc))%gtype=v3dvar
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idUTav(itrc)))) THEN
              got_var(idUTav(itrc))=.TRUE.
              AVG(ng)%pioVar(idUTav(itrc))%vd=var_desc(i)
              AVG(ng)%pioVar(idUTav(itrc))%dkind=PIO_FOUT
              AVG(ng)%pioVar(idUTav(itrc))%gtype=u3dvar
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idVTav(itrc)))) THEN
              got_var(idVTav(itrc))=.TRUE.
              AVG(ng)%pioVar(idVTav(itrc))%vd=var_desc(i)
              AVG(ng)%pioVar(idVTav(itrc))%dkind=PIO_FOUT
              AVG(ng)%pioVar(idVTav(itrc))%gtype=v3dvar
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idTTav(itrc)))) THEN
              got_var(idTTav(itrc))=.TRUE.
              AVG(ng)%pioVar(idTTav(itrc))%vd=var_desc(i)
              AVG(ng)%pioVar(idTTav(itrc))%dkind=PIO_FOUT
              AVG(ng)%pioVar(idTTav(itrc))%gtype=r3dvar
            END IF
          END DO
        END DO
!
!  Check if averages variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsur).and.Aout(idFsur,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar).and.Aout(idUbar,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar).and.Aout(idVbar,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idu2dE).and.Aout(idu2dE,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idu2dE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv2dN).and.Aout(idv2dN,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idv2dN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUvel).and.Aout(idUvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel).and.Aout(idVvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idu3dE).and.Aout(idu3dE,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idu3dE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv3dN).and.Aout(idv3dN,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idv3dN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idOvel).and.Aout(idOvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idOvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWvel).and.Aout(idWvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idDano).and.Aout(idDano,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idDano)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(id2dPV).and.Aout(id2dPV,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,id2dPV)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(id2dRV).and.Aout(id2dRV,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,id2dRV)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(id3dPV).and.Aout(id3dPV,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,id3dPV)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(id3dRV).and.Aout(id3dRV,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,id3dRV)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idZZav).and.Aout(idZZav,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idZZav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idU2av).and.Aout(idU2av,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idU2av)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV2av).and.Aout(idV2av,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idV2av)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idHUav).and.Aout(idHUav,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHUav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idHVav).and.Aout(idHVav,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHVav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUUav).and.Aout(idUUav,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUUav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUVav).and.Aout(idUVav,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUVav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVVav).and.Aout(idVVav,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVVav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvis).and.Aout(idVvis,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVvis)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTdif).and.Aout(idTdif,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSdif).and.Aout(idSdif,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idPair).and.Aout(idPair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idPair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUair).and.Aout(idUair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVair).and.Aout(idVair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUaiE).and.Aout(idUaiE,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUaiE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVaiN).and.Aout(idVaiN,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVaiN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTsur(itemp)).and.Aout(idTsur(itemp),ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTsur(itemp))),   &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTsur(isalt)).and.Aout(idTsur(isalt),ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTsur(isalt))),   &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLhea).and.Aout(idLhea,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idLhea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idShea).and.Aout(idShea,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idShea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLrad).and.Aout(idLrad,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idLrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTair).and.Aout(idTair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idevap).and.Aout(idevap,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idevap)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idrain).and.Aout(idrain,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idrain)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSrad).and.Aout(idSrad,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUsms).and.Aout(idUsms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms).and.Aout(idVsms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbms).and.Aout(idUbms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbms).and.Aout(idVbms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc)).and.Aout(idTvar(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO itrc=1,NAT
          IF (.not.got_var(iHUTav(itrc)).and.Aout(iHUTav(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,iHUTav(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(iHVTav(itrc)).and.Aout(iHVTav(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,iHVTav(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idUTav(itrc)).and.Aout(idUTav(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUTav(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idVTav(itrc)).and.Aout(idVTav(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVTav(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idTTav(itrc)).and.Aout(idTTav(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTTav(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
!
!  Set unlimited time record dimension to the appropriate value.
!
        IF (nRST(ng).eq.nAVG(ng)) THEN
          IF (ndefAVG(ng).gt.0) THEN
            AVG(ng)%Rindex=((ntstart(ng)-1)-                            &
     &                      ndefAVG(ng)*((ntstart(ng)-1)/ndefAVG(ng)))/ &
     &                     nAVG(ng)
          ELSE
            AVG(ng)%Rindex=(ntstart(ng)-1)/nAVG(ng)
          END IF
        ELSE
          AVG(ng)%Rindex=rec_size
        END IF
      END IF QUERY
!
!  Set initial average time. Notice that the value is offset by half
!  nAVG*dt so there is not a special case when computing its value
!  in "set_avg".
!
      IF (ntsAVG(ng).eq.1) THEN
        AVGtime(ng)=time(ng)-0.5_r8*REAL(nAVG(ng),r8)*dt(ng)
      ELSE
        AVGtime(ng)=time(ng)+REAL(ntsAVG(ng),r8)*dt(ng)-                &
     &              0.5_r8*REAL(nAVG(ng),r8)*dt(ng)
      END IF
!
  10  FORMAT (2x,'DEF_AVG_PIO      - creating average file,',t56,       &
     &        'Grid ',i2.2,': ',a)
  20  FORMAT (2x,'DEF_AVG_PIO      - inquiring average file',t56,       &
     &        'Grid ',i2.2,': ',a)
  30  FORMAT (/,' DEF_AVG_PIO - unable to create averages NetCDF',      &
     &        ' file: ',a)
  40  FORMAT (1pe11.4,1x,'millimeter')
  50  FORMAT (/,' DEF_AVG_PIO - unable to open averages NetCDF'         &
     &        ' file: ',a)
  60  FORMAT (/,' DEF_AVG_PIO - unable to find variable: ',a,2x,        &
     &        ' in averages NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE def_avg_pio
      END MODULE def_avg_mod
