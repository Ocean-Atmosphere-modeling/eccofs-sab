      MODULE wrt_his_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module writes requested model fields into the HISTORY output   !
!  file using either the standard NetCDF library or the Parallel-IO    !
!  (PIO) library.                                                      !
!                                                                      !
!  Notice that only momentum is affected by the full time-averaged     !
!  masks.  If applicable, these mask contains information about        !
!  river runoff and time-dependent wetting and drying variations.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_coupling
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
      USE extract_slice_mod,     ONLY : extract_slice
      USE nf_fwrite2d_mod,       ONLY : nf_fwrite2d
      USE nf_fwrite3d_mod,       ONLY : nf_fwrite3d
      USE omega_mod,             ONLY : scale_omega
      USE strings_mod,           ONLY : FoundError
      USE uv_rotate_mod,         ONLY : uv_rotate2d
      USE uv_rotate_mod,         ONLY : uv_rotate3d
      USE vorticity_mod,         ONLY : pvorticity3d, rvorticity3d
!
      implicit none
!
      PUBLIC  :: wrt_his
      PRIVATE :: wrt_his_nf90
      PRIVATE :: wrt_his_pio
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wrt_his (ng, tile)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: IminS, ImaxS, JminS, JmaxS
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/wrt_his.F"
!
!-----------------------------------------------------------------------
!  Write out history fields according to IO type.
!-----------------------------------------------------------------------
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
      IminS=BOUNDS(ng)%Istr(tile)-3
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-3
      JmaxS=BOUNDS(ng)%Jend(tile)+3
!
      SELECT CASE (HIS(ng)%IOtype)
        CASE (io_nf90)
          CALL wrt_his_nf90 (ng, iNLM, tile,                            &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS)
        CASE (io_pio)
          CALL wrt_his_pio (ng, iNLM, tile,                             &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS)
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) HIS(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, 167, MyFile)) RETURN
!
  10  FORMAT (' WRT_HIS - Illegal output file type, io_type = ',i0,     &
     &        /,11x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE wrt_his
!
!***********************************************************************
      SUBROUTINE wrt_his_nf90 (ng, model, tile,                         &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      integer :: Fcount, gfactor, gtype, ifield, status
      integer :: i, itrc, j, k
!
      real(dp) :: scale
      real(r8), allocatable :: Ur2d(:,:)
      real(r8), allocatable :: Vr2d(:,:)
      real(r8), allocatable :: Fr3d(:,:,:)
      real(r8), allocatable :: Wr3d(:,:,:)
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/wrt_his.F"//", wrt_his_nf90"
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
!  Write out history fields.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, 222, MyFile)) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
      gfactor=1
!
!  Set time record index.
!
      HIS(ng)%Rindex=HIS(ng)%Rindex+1
      Fcount=HIS(ng)%load
      HIS(ng)%Nrec(Fcount)=HIS(ng)%Nrec(Fcount)+1
!
!  Report.
!
      IF (Master) WRITE (stdout,10) kstp(ng), nrhs(ng), HIS(ng)%Rindex
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, model, HIS(ng)%name,                    &
     &                      TRIM(Vname(1,idtime)), time(ng:),           &
     &                      (/HIS(ng)%Rindex/), (/1/),                  &
     &                      ncid = HIS(ng)%ncid,                        &
     &                      varid = HIS(ng)%Vid(idtime))
      IF (FoundError(exit_flag, NoError, 262, MyFile)) RETURN
!
!  Write time-varying depths of RHO-points.
!
      IF (Hout(idpthR,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idpthR,             &
     &                     HIS(ng)%Vid(idpthR),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % z_r,                              &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, 370, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthR)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write time-varying depths of U-points.
!
      IF (Hout(idpthU,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*u3dvar
        DO k=1,N(ng)
          DO j=Jstr-1,Jend+1
            DO i=IstrU-1,Iend+1
              GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i-1,j,k)+        &
     &                                    GRID(ng)%z_r(i  ,j,k))
            END DO
          END DO
        END DO
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idpthU,             &
     &                     HIS(ng)%Vid(idpthU),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % z_v,                              &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, 402, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthU)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write time-varying depths of V-points.
!
      IF (Hout(idpthV,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*v3dvar
        DO k=1,N(ng)
          DO j=JstrV-1,Jend+1
            DO i=Istr-1,Iend+1
              GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i,j-1,k)+        &
     &                                    GRID(ng)%z_r(i,j  ,k))
            END DO
          END DO
        END DO
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idpthV,             &
     &                     HIS(ng)%Vid(idpthV),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % vmask,                            &
     &                     GRID(ng) % z_v,                              &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, 434, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthV)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write time-varying depths of W-points.
!
      IF (Hout(idpthW,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idpthW,             &
     &                     HIS(ng)%Vid(idpthW),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % z_w,                              &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, 458, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthW)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out free-surface (m)
!
      IF (Hout(idFsur,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idFsur,             &
     &                     HIS(ng)%Vid(idFsur),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     OCEAN(ng) % zeta(:,:,kstp(ng)))
        IF (FoundError(status, nf90_noerr, 487, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idFsur)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D U-momentum component (m/s).
!
      IF (Hout(idUbar,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idUbar,             &
     &                     HIS(ng)%Vid(idUbar),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % umask_full,                       &
     &                     OCEAN(ng) % ubar(:,:,kstp(ng)))
        IF (FoundError(status, nf90_noerr, 552, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUbar)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D V-momentum component (m/s).
!
      IF (Hout(idVbar,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idVbar,             &
     &                     HIS(ng)%Vid(idVbar),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % vmask_full,                       &
     &                     OCEAN(ng) % vbar(:,:,kstp(ng)))
        IF (FoundError(status, nf90_noerr, 671, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVbar)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D Eastward and Northward momentum components (m/s) at
!  RHO-points.
!
      IF (Hout(idu2dE,ng).and.Hout(idv2dN,ng)) THEN
        IF (.not.allocated(Ur2d)) THEN
          allocate (Ur2d(LBi:UBi,LBj:UBj))
            Ur2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        IF (.not.allocated(Vr2d)) THEN
          allocate (Vr2d(LBi:UBi,LBj:UBj))
            Vr2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        CALL uv_rotate2d (ng, tile, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
     &                    GRID(ng) % rmask_full,                        &
     &                    OCEAN(ng) % ubar(:,:,kstp(ng)),               &
     &                    OCEAN(ng) % vbar(:,:,kstp(ng)),               &
     &                    Ur2d, Vr2d)
!
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idu2dE,             &
     &                     HIS(ng)%Vid(idu2dE),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_full,                       &
     &                     Ur2d)
        IF (FoundError(status, nf90_noerr, 810, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idu2dE)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idv2dN,             &
     &                     HIS(ng)%Vid(idv2dN),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_full,                       &
     &                     Vr2d)
        IF (FoundError(status, nf90_noerr, 827, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idv2dN)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF
!
!  Write out 3D U-momentum component (m/s).
!
      IF (Hout(idUvel,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idUvel,             &
     &                     HIS(ng)%Vid(idUvel),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % umask_full,                       &
     &                     OCEAN(ng) % u(:,:,:,nrhs(ng)))
        IF (FoundError(status, nf90_noerr, 854, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUvel)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D U-momentum component (m/s) at specified constant depth
!  slices.
!
      IF (Hout(idUzsl,ng).and.(Nslice.gt.0)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
          Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
        END IF
        scale=1.0_dp
        gtype=gfactor*u3dvar
!
        DO k=1,N(ng)
          DO j=Jstr-1,Jend+1
            DO i=IstrU-1,Iend+1
              GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i-1,j,k)+        &
     &                                    GRID(ng)%z_r(i  ,j,k))
            END DO
          END DO
        END DO
        CALL extract_slice (ng, model, tile, gtype,                     &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      OCEAN(ng) % u(:,:,:,nrhs(ng)),              &
     &                      GRID(ng) % z_v,                             &
     &                      GRID(ng) % umask_full,                      &
     &                      Zslice, Wr3d)
!
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idUzsl,             &
     &                     HIS(ng)%Vid(idUzsl),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, Nslice, scale,        &
     &                     GRID(ng) % umask_full,                       &
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, 918, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUzsl)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out 3D V-momentum component (m/s).
!
      IF (Hout(idVvel,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idVvel,             &
     &                     HIS(ng)%Vid(idVvel),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % vmask_full,                       &
     &                     OCEAN(ng) % v(:,:,:,nrhs(ng)))
        IF (FoundError(status, nf90_noerr, 967, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVvel)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D V-momentum component (m/s) at specified constant depth
!  slices.
!
      IF (Hout(idVzsl,ng).and.(Nslice.gt.0)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
          Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
        END IF
        scale=1.0_dp
        gtype=gfactor*v3dvar
!
        DO k=1,N(ng)
          DO j=JstrV-1,Jend+1
            DO i=Istr-1,Iend+1
              GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i,j-1,k)+        &
     &                                    GRID(ng)%z_r(i,j  ,k))
            END DO
          END DO
        END DO
        CALL extract_slice (ng, model, tile, gtype,                     &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      OCEAN(ng) % v(:,:,:,nrhs(ng)),              &
     &                      GRID(ng) % z_v,                             &
     &                      GRID(ng) % vmask_full,                      &
     &                      Zslice, Wr3d)
!
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idVzsl,             &
     &                     HIS(ng)%Vid(idVzsl),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, Nslice, scale,        &
     &                     GRID(ng) % vmask_full,                       &
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, 1031, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVvel)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out 3D Eastward momentum (m/s) at RHO-points, A-grid.
!
      IF (Hout(idu3dE,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idu3dE,             &
     &                     HIS(ng)%Vid(idu3dE),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % rmask_full,                       &
     &                     OCEAN(ng) % ua)
        IF (FoundError(status, nf90_noerr, 1080, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idu3dE)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Eastward momentum (m/s) at RHO-points, A-grid, at
!  specified constant depth slices.
!
      IF (Hout(idUzsE,ng).and.(Nslice.gt.0)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
          Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
        END IF
        scale=1.0_dp
        gtype=gfactor*r3dvar
        CALL extract_slice (ng, model, tile, gtype,                     &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      OCEAN(ng) % ua,                             &
     &                      GRID(ng) % z_r,                             &
     &                      GRID(ng) % rmask,                           &
     &                      Zslice, Wr3d)
!
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idUzsE,             &
     &                     HIS(ng)%Vid(idUzsE),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, Nslice, scale,        &
     &                     GRID(ng) % rmask_full,                       &
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, 1117, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idu3dE)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out 3D Northward momentum (m/s) at RHO-points, A-grid.
!
      IF (Hout(idv3dN,ng)) THEN
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idv3dN,             &
     &                     HIS(ng)%Vid(idv3dN),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % rmask_full,                       &
     &                     OCEAN(ng) % va)
        IF (FoundError(status, nf90_noerr, 1139, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idv3dN)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Northward momentum (m/s) at RHO-points, A-grid, at
!  specified constant depth slices.
!
      IF (Hout(idVzsN,ng).and.(Nslice.gt.0)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
          Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
        END IF
        scale=1.0_dp
        gtype=gfactor*r3dvar
        CALL extract_slice (ng, model, tile, gtype,                     &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      OCEAN(ng) % va,                             &
     &                      GRID(ng) % z_r,                             &
     &                      GRID(ng) % rmask,                           &
     &                      Zslice, Wr3d)
!
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idVzsN,             &
     &                     HIS(ng)%Vid(idVzsN),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, Nslice, scale,        &
     &                     GRID(ng) % rmask_full,                       &
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, 1176, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idv3dN)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out potential vorticity (m-1 s-1) at PSI-points at full column
!  or specified constant depth slices.
!
      IF (Hout(id3dPV,ng).or.(Hout(idPVzs,ng).and.(Nslice.gt.0))) THEN
        IF (.not.allocated(Fr3d)) THEN
          allocate ( Fr3d(LBi:UBi,LBj:UBj,N(ng)) )
          Fr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
        END IF
        scale=1.0_dp
        gtype=gfactor*p3dvar
!
        CALL pvorticity3d (ng, model, tile,                             &
                           LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS, nrhs(ng),        &
     &                     GRID(ng) % pmask,                            &
     &                     GRID(ng) % umask,   GRID(ng) % vmask,        &
                           GRID(ng) % f,                                &
                           GRID(ng) % om_u,    GRID(ng) % on_v,         &
                           GRID(ng) % pm,      GRID(ng) % pn,           &
                           GRID(ng) % z_r,     OCEAN(ng) % pden,        &
                           OCEAN(ng) % u,      OCEAN(ng) % v,           &
                           Fr3d)
!
        IF (Hout(id3dPV,ng)) THEN
          status=nf_fwrite3d(ng, model, HIS(ng)%ncid, id3dPV,           &
     &                       HIS(ng)%Vid(id3dPV),                       &
     &                       HIS(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
     &                       GRID(ng) % pmask,                          &
     &                       Fr3d)
          IF (FoundError(status, nf90_noerr, 1221, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,id3dPV)), HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
        IF (Hout(idPVzs,ng).and.(Nslice.gt.0)) THEN
          IF (.not.allocated(Wr3d)) THEN
            allocate ( Wr3d(LBi:UBi,LBj:UBj,N(ng)) )
            Wr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
          END IF
!
          DO k=1,N(ng)
            DO j=JstrP,JendT
              DO i=IstrP,IendT
                GRID(ng)%z_v(i,j,k)=0.25_r8*(GRID(ng)%z_r(i-1,j-1,k)+   &
     &                                       GRID(ng)%z_r(i-1,j  ,k)+   &
     &                                       GRID(ng)%z_r(i  ,j-1,k)+   &
     &                                       GRID(ng)%z_r(i ,j  ,k))
              END DO
            END DO
          END DO
!
          CALL extract_slice (ng, model, tile, p3dvar,                  &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        Fr3d,                                     &
     &                        GRID(ng) % z_v,                           &
     &                        GRID(ng) % pmask,                         &
     &                        Zslice, Wr3d)
!
          status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idPVzs,           &
     &                       HIS(ng)%Vid(idPVzs),                       &
     &                       HIS(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, Nslice, scale,      &
     &                       GRID(ng) % pmask,                          &
     &                       Wr3d)
          IF (FoundError(status, nf90_noerr, 1265, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idPVzs)), HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
          deallocate (Wr3d)
        END IF
        deallocate (Fr3d)
      END IF
!
!  Write out relative vorticity (s-1) at PSI-points at full column
!  or specified constant depth slices.
!
      IF (Hout(id3dRV,ng).or.(Hout(idRVzs,ng).and.(Nslice.gt.0))) THEN
        IF (.not.allocated(Fr3d)) THEN
          allocate ( Fr3d(LBi:UBi,LBj:UBj,N(ng)) )
          Fr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
        END IF
        scale=1.0_dp
        gtype=gfactor*p3dvar
!
        CALL rvorticity3d (ng, model, tile,                             &
                           LBi, UBi, LBj, UBj, nrhs(ng),                &
     &                     GRID(ng) % pmask,                            &
                           GRID(ng) % om_u, GRID(ng) % on_v,            &
                           GRID(ng) % pm,   GRID(ng) % pn,              &
                           OCEAN(ng) % u,   OCEAN(ng) % v,              &
                           Fr3d)
!
        IF (Hout(id3dRV,ng)) THEN
          status=nf_fwrite3d(ng, model, HIS(ng)%ncid, id3dRV,           &
     &                       HIS(ng)%Vid(id3dRV),                       &
     &                       HIS(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
     &                       GRID(ng) % pmask,                          &
     &                       Fr3d)
          IF (FoundError(status, nf90_noerr, 1308, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,id3dRV)), HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
        IF (Hout(idRVzs,ng).and.(Nslice.gt.0)) THEN
          IF (.not.allocated(Wr3d)) THEN
            allocate ( Wr3d(LBi:UBi,LBj:UBj,N(ng)) )
            Wr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
          END IF
!
          DO k=1,N(ng)
            DO j=JstrP,JendT
              DO i=IstrP,IendT
                GRID(ng)%z_v(i,j,k)=0.25_r8*(GRID(ng)%z_r(i-1,j-1,k)+   &
     &                                       GRID(ng)%z_r(i-1,j  ,k)+   &
     &                                       GRID(ng)%z_r(i  ,j-1,k)+   &
     &                                       GRID(ng)%z_r(i ,j  ,k))
              END DO
            END DO
          END DO
!
          CALL extract_slice (ng, model, tile, p3dvar,                  &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        Fr3d,                                     &
     &                        GRID(ng) % z_v,                           &
     &                        GRID(ng) % pmask,                         &
     &                        Zslice, Wr3d)
!
          status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idRVzs,           &
     &                       HIS(ng)%Vid(idRVzs),                       &
     &                       HIS(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, Nslice, scale,      &
     &                       GRID(ng) % pmask,                          &
     &                       Wr3d)
          IF (FoundError(status, nf90_noerr, 1352, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idRVzs)), HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
          deallocate (Wr3d)
        END IF
        deallocate (Fr3d)
      END IF
!
!  Write out S-coordinate omega vertical velocity (m/s).
!
      IF (Hout(idOvel,ng)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,0:N(ng)) )
          Wr3d(LBi:UBi,LBj:UBj,0:N(ng))=0.0_r8
        END IF
        scale=1.0_dp
        gtype=gfactor*w3dvar
        CALL scale_omega (ng, tile, LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    OCEAN(ng) % W,                                &
     &                    Wr3d)
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idOvel,             &
     &                     HIS(ng)%Vid(idOvel),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, 1387, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idOvel)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out S-coordinate implicit vertical "omega" momentum component.
!
      IF (Hout(idOvil,ng)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,0:N(ng)) )
          Wr3d(LBi:UBi,LBj:UBj,0:N(ng))=0.0_r8
        END IF
        scale=1.0_dp
        gtype=gfactor*w3dvar
        CALL scale_omega (ng, tile, LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    OCEAN(ng) % Wi,                               &
     &                    Wr3d(LBi:,LBj:,0:))
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idOvil,             &
     &                     HIS(ng)%Vid(idOvil),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     Wr3d(LBi:,LBj:,0:))
        IF (FoundError(status, nf90_noerr, 1422, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idOvil)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out vertical velocity (m/s).
!
      IF (Hout(idWvel,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idWvel,             &
     &                     HIS(ng)%Vid(idWvel),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     OCEAN(ng) % wvel)
        IF (FoundError(status, nf90_noerr, 1447, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idWvel)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        IF (Hout(idTvar(itrc),ng)) THEN
          scale=1.0_dp
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idTvar(itrc),     &
     &                       HIS(ng)%Tid(itrc),                         &
     &                       HIS(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
     &                       GRID(ng) % rmask,                          &
     &                       OCEAN(ng) % t(:,:,:,nrhs(ng),itrc))
          IF (FoundError(status, nf90_noerr, 1471, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idTvar(itrc))),            &
     &                          HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out tracer type variables at specified constant depth slices.
!
      DO itrc=1,NT(ng)
        IF (Hout(idzslT(itrc),ng).and.(Nslice.gt.0)) THEN
          IF (.not.allocated(Wr3d)) THEN
            allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
            Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
          END IF
          scale=1.0_dp
          gtype=gfactor*r3dvar
          CALL extract_slice (ng, model, tile, gtype,                   &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        OCEAN(ng) % t(:,:,:,nrhs(ng),itrc),       &
     &                        GRID(ng) % z_r,                           &
     &                        GRID(ng) % rmask,                         &
     &                        Zslice, Wr3d)
!
          status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idzslT(itrc),     &
     &                       HIS(ng)%Vid(idzslT(itrc)),                 &
     &                       HIS(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, Nslice, scale,      &
     &                       GRID(ng) % rmask,                          &
     &                       Wr3d)
          IF (FoundError(status, nf90_noerr, 1510, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idzslT(itrc))),            &
     &                          HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
          IF (itrc.eq.NT(ng)) deallocate (Wr3d)
        END IF
      END DO
!
!  Write out density anomaly.
!
      IF (Hout(idDano,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idDano,             &
     &                     HIS(ng)%Vid(idDano),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     OCEAN(ng) % rho)
        IF (FoundError(status, nf90_noerr, 1564, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idDano)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical viscosity coefficient.
!
      IF (Hout(idVvis,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idVvis,             &
     &                     HIS(ng)%Vid(idVvis),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     MIXING(ng) % Akv,                            &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, 1665, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVvis)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      IF (Hout(idTdif,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idTdif,             &
     &                     HIS(ng)%Vid(idTdif),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     MIXING(ng) % Akt(:,:,:,itemp),               &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, 1689, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTdif)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for salinity.
!
      IF (Hout(idSdif,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idSdif,             &
     &                     HIS(ng)%Vid(idSdif),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     MIXING(ng) % Akt(:,:,:,isalt),               &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, 1714, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSdif)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out turbulent kinetic energy.
!
      IF (Hout(idMtke,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idMtke,             &
     &                     HIS(ng)%Vid(idMtke),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     MIXING(ng) % tke(:,:,:,nrhs(ng)),            &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, 1740, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idMtke)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out turbulent length scale field.
!
      IF (Hout(idMtls,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, HIS(ng)%ncid, idMtls,             &
     &                     HIS(ng)%Vid(idMtls),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     MIXING(ng) % gls(:,:,:,nrhs(ng)),            &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, 1784, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idMtls)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface air pressure.
!
      IF (Hout(idPair,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idPair,             &
     &                     HIS(ng)%Vid(idPair),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     FORCES(ng) % Pair)
        IF (FoundError(status, nf90_noerr, 1849, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idPair)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface air temperature.
!
      IF (Hout(idTair,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idTair,             &
     &                     HIS(ng)%Vid(idTair),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     FORCES(ng) % Tair)
        IF (FoundError(status, nf90_noerr, 1874, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTair)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface winds.
!
      IF (Hout(idUair,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idUair,             &
     &                     HIS(ng)%Vid(idUair),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     FORCES(ng) % Uwind)
        IF (FoundError(status, nf90_noerr, 1899, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUair)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
      IF (Hout(idVair,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idVair,             &
     &                     HIS(ng)%Vid(idVair),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     FORCES(ng) % Vwind)
        IF (FoundError(status, nf90_noerr, 1920, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVair)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out Eastward/Northward surface wind (m/s) at RHO-points.
!
      IF (Hout(idUaiE,ng).and.Hout(idVaiN,ng)) THEN
        IF (.not.allocated(Ur2d)) THEN
          allocate (Ur2d(LBi:UBi,LBj:UBj))
          Ur2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        IF (.not.allocated(Vr2d)) THEN
          allocate (Vr2d(LBi:UBi,LBj:UBj))
          Vr2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        CALL uv_rotate2d (ng, tile, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
     &                    GRID(ng) % rmask_full,                        &
     &                    FORCES(ng) % Uwind,                           &
     &                    FORCES(ng) % Vwind,                           &
     &                    Ur2d, Vr2d)
!
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idUaiE,             &
     &                     HIS(ng)%Vid(idUaiE),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     Ur2d)
        IF (FoundError(status, nf90_noerr, 1962, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUaiE)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idVaiN,             &
     &                     HIS(ng)%Vid(idVaiN),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     Vr2d)
        IF (FoundError(status, nf90_noerr, 1981, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVaiN)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF
!
!  Write out surface active tracers fluxes.
!
      DO itrc=1,NAT
        IF (Hout(idTsur(itrc),ng)) THEN
          IF (itrc.eq.itemp) THEN
            scale=rho0*Cp                   ! Celsius m/s to W/m2
          ELSE IF (itrc.eq.isalt) THEN
            scale=1.0_dp
          END IF
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idTsur(itrc),     &
     &                       HIS(ng)%Vid(idTsur(itrc)),                 &
     &                       HIS(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       FORCES(ng) % stflx(:,:,itrc))
          IF (FoundError(status, nf90_noerr, 2016, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idTsur(itrc))),            &
     &                          HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out latent heat flux.
!
      IF (Hout(idLhea,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idLhea,             &
     &                     HIS(ng)%Vid(idLhea),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     FORCES(ng) % lhflx)
        IF (FoundError(status, nf90_noerr, 2043, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idLhea)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out sensible heat flux.
!
      IF (Hout(idShea,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idShea,             &
     &                     HIS(ng)%Vid(idShea),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     FORCES(ng) % shflx)
        IF (FoundError(status, nf90_noerr, 2066, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idShea)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out net longwave radiation flux.
!
      IF (Hout(idLrad,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idLrad,             &
     &                     HIS(ng)%Vid(idLrad),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     FORCES(ng) % lrflx)
        IF (FoundError(status, nf90_noerr, 2089, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idLrad)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out evaporation rate (kg/m2/s).
!
      IF (Hout(idevap,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idevap,             &
     &                     HIS(ng)%Vid(idevap),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     FORCES(ng) % evap)
        IF (FoundError(status, nf90_noerr, 2116, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idevap)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out precipitation rate (kg/m2/s).
!
      IF (Hout(idrain,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idrain,             &
     &                     HIS(ng)%Vid(idrain),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     FORCES(ng) % rain)
        IF (FoundError(status, nf90_noerr, 2139, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idrain)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out E-P (m/s).
!
      IF (Hout(idEmPf,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idEmPf,             &
     &                     HIS(ng)%Vid(idEmPf),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     FORCES(ng) % stflux(:,:,isalt))
        IF (FoundError(status, nf90_noerr, 2164, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idEmPf)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out net shortwave radiation flux.
!
      IF (Hout(idSrad,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idSrad,             &
     &                     HIS(ng)%Vid(idSrad),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     FORCES(ng) % srflx)
        IF (FoundError(status, nf90_noerr, 2188, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSrad)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface U-momentum stress.
!
      IF (Hout(idUsms,ng)) THEN
        scale=rho0                          ! m2/s2 to Pa
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idUsms,             &
     &                     HIS(ng)%Vid(idUsms),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % umask,                            &
     &                     FORCES(ng) % sustr)
        IF (FoundError(status, nf90_noerr, 2217, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUsms)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface V-momentum stress.
!
      IF (Hout(idVsms,ng)) THEN
        scale=rho0
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idVsms,             &
     &                     HIS(ng)%Vid(idVsms),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % vmask,                            &
     &                     FORCES(ng) % svstr)
        IF (FoundError(status, nf90_noerr, 2244, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVsms)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom U-momentum stress.
!
      IF (Hout(idUbms,ng)) THEN
        scale=-rho0
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idUbms,             &
     &                     HIS(ng)%Vid(idUbms),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % umask,                            &
     &                     FORCES(ng) % bustr)
        IF (FoundError(status, nf90_noerr, 2267, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUbms)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom V-momentum stress.
!
      IF (Hout(idVbms,ng)) THEN
        scale=-rho0
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, model, HIS(ng)%ncid, idVbms,             &
     &                     HIS(ng)%Vid(idVbms),                         &
     &                     HIS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % vmask,                            &
     &                     FORCES(ng) % bvstr)
        IF (FoundError(status, nf90_noerr, 2290, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVbms)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Synchronize history NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, model, HIS(ng)%name, HIS(ng)%ncid)
      IF (FoundError(exit_flag, NoError, 2368, MyFile)) RETURN
!
  10  FORMAT (2x,'WRT_HIS_NF90     - writing history', t42,             &
     &        'fields (Index=',i1,',',i1,') in record = ',i0)
  20  FORMAT (/,' WRT_HIS_NF90 - error while writing variable: ',a,     &
     &        /,16x,'into history NetCDF file for time record: ',i0)
!
      RETURN
      END SUBROUTINE wrt_his_nf90
!
!***********************************************************************
      SUBROUTINE wrt_his_pio (ng, model, tile,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      integer :: Fcount, ifield, status
      integer :: i, itrc, j, k
!
      real(dp) :: scale
      real(r8), allocatable :: Ur2d(:,:)
      real(r8), allocatable :: Vr2d(:,:)
      real(r8), allocatable :: Fr3d(:,:,:)
      real(r8), allocatable :: Wr3d(:,:,:)
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/wrt_his.F"//", wrt_his_pio"
!
      TYPE (IO_desc_t), pointer :: ioDesc
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
!  Write out history fields.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, 2442, MyFile)) RETURN
!
!  Set time record index.
!
      HIS(ng)%Rindex=HIS(ng)%Rindex+1
      Fcount=HIS(ng)%load
      HIS(ng)%Nrec(Fcount)=HIS(ng)%Nrec(Fcount)+1
!
!  Report.
!
      IF (Master) WRITE (stdout,10) kstp(ng), nrhs(ng), HIS(ng)%Rindex
!
!  Write out model time (s).
!
      CALL pio_netcdf_put_fvar (ng, model, HIS(ng)%name,                &
     &                          TRIM(Vname(1,idtime)), time(ng:),       &
     &                          (/HIS(ng)%Rindex/), (/1/),              &
     &                           pioFile = HIS(ng)%pioFile,             &
     &                           pioVar = HIS(ng)%pioVar(idtime)%vd)
      IF (FoundError(exit_flag, NoError, 2473, MyFile)) RETURN
!
!  Write time-varying depths of RHO-points.
!
      IF (Hout(idpthR,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idpthR)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idpthR,          &
     &                     HIS(ng)%pioVar(idpthR),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % z_r,                              &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, 2606, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthR)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write time-varying depths of U-points.
!
      IF (Hout(idpthU,ng)) THEN
        scale=1.0_dp
        DO k=1,N(ng)
          DO j=Jstr-1,Jend+1
            DO i=IstrU-1,Iend+1
              GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i-1,j,k)+        &
     &                                    GRID(ng)%z_r(i  ,j,k))
            END DO
          END DO
        END DO
        IF (HIS(ng)%pioVar(idpthU)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idpthU,          &
     &                     HIS(ng)%pioVar(idpthU),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % z_v,                              &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, 2643, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthU)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write time-varying depths of V-points.
!
      IF (Hout(idpthV,ng)) THEN
        scale=1.0_dp
        DO k=1,N(ng)
          DO j=JstrV-1,Jend+1
            DO i=Istr-1,Iend+1
              GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i,j-1,k)+        &
     &                                    GRID(ng)%z_r(i,j  ,k))
            END DO
          END DO
        END DO
        IF (HIS(ng)%pioVar(idpthV)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idpthV,          &
     &                     HIS(ng)%pioVar(idpthV),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % vmask,                            &
     &                     GRID(ng) % z_v,                              &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, 2680, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthV)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write time-varying depths of W-points.
!
      IF (Hout(idpthW,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idpthW)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idpthW,          &
     &                     HIS(ng)%pioVar(idpthW),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % z_w,                              &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, 2709, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthW)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out free-surface (m)
!
      IF (Hout(idFsur,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idFsur)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idFsur,          &
     &                     HIS(ng)%pioVar(idFsur),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     OCEAN(ng) % zeta(:,:,kstp(ng)))
        IF (FoundError(status, PIO_noerr, 2743, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idFsur)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D U-momentum component (m/s).
!
      IF (Hout(idUbar,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idUbar)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idUbar,          &
     &                     HIS(ng)%pioVar(idUbar),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % umask_full,                       &
     &                     OCEAN(ng) % ubar(:,:,kstp(ng)))
        IF (FoundError(status, PIO_noerr, 2823, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUbar)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D V-momentum component (m/s).
!
      IF (Hout(idVbar,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idVbar)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idVbar,          &
     &                     HIS(ng)%pioVar(idVbar),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % vmask_full,                       &
     &                     OCEAN(ng) % vbar(:,:,kstp(ng)))
        IF (FoundError(status, PIO_noerr, 2978, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVbar)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D Eastward and Northward momentum components (m/s) at
!  RHO-points.
!
      IF (Hout(idu2dE,ng).and.Hout(idv2dN,ng)) THEN
        IF (.not.allocated(Ur2d)) THEN
          allocate (Ur2d(LBi:UBi,LBj:UBj))
            Ur2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        IF (.not.allocated(Vr2d)) THEN
          allocate (Vr2d(LBi:UBi,LBj:UBj))
            Vr2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        CALL uv_rotate2d (ng, tile, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
     &                    GRID(ng) % rmask_full,                        &
     &                    OCEAN(ng) % ubar(:,:,kstp(ng)),               &
     &                    OCEAN(ng) % vbar(:,:,kstp(ng)),               &
     &                    Ur2d, Vr2d)
!
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idu2dE)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idu2dE,          &
     &                     HIS(ng)%pioVar(idu2dE),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_full,                       &
     &                     Ur2d)
        IF (FoundError(status, PIO_noerr, 3153, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idu2dE)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        IF (HIS(ng)%pioVar(idv2dN)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idv2dN,          &
     &                     HIS(ng)%pioVar(idv2dN),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_full,                       &
     &                     Vr2d)
        IF (FoundError(status, PIO_noerr, 3176, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idv2dN)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF
!
!  Write out 3D U-momentum component (m/s).
!
      IF (Hout(idUvel,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idUvel)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idUvel,          &
     &                     HIS(ng)%pioVar(idUvel),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % umask_full,                       &
     &                     OCEAN(ng) % u(:,:,:,nrhs(ng)))
        IF (FoundError(status, PIO_noerr, 3208, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUvel)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D U-momentum component (m/s) at specified constant depth
!  slices.
!
      IF (Hout(idUzsl,ng).and.(Nslice.gt.0)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
          Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
        END IF
!
        DO k=1,N(ng)
          DO j=Jstr-1,Jend+1
            DO i=IstrU-1,Iend+1
              GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i-1,j,k)+        &
     &                                    GRID(ng)%z_r(i  ,j,k))
            END DO
          END DO
        END DO
        CALL extract_slice (ng, model, tile, u3dvar,                    &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      OCEAN(ng) % u(:,:,:,nrhs(ng)),              &
     &                      GRID(ng) % z_v,                             &
     &                      GRID(ng) % umask_full,                      &
     &                      Zslice, Wr3d)
!
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idUzsl)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesZ_dp_u3dvar(ng)
        ELSE
          ioDesc => ioDesZ_sp_u3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idUvel,          &
     &                     HIS(ng)%pioVar(idUzsl),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, Nslice, scale,        &
     &                     GRID(ng) % umask_full,                       &
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, 3285, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUzsl)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out 3D V-momentum component (m/s).
!
      IF (Hout(idVvel,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idVvel)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idVvel,          &
     &                     HIS(ng)%pioVar(idVvel),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % vmask_full,                       &
     &                     OCEAN(ng) % v(:,:,:,nrhs(ng)))
        IF (FoundError(status, PIO_noerr, 3340, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVvel)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D V-momentum component (m/s) at specified constant depth
!  slices.
!
      IF (Hout(idVzsl,ng).and.(Nslice.gt.0)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
          Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
        END IF
!
        DO k=1,N(ng)
          DO j=JstrV-1,Jend+1
            DO i=Istr-1,Iend+1
              GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i,j-1,k)+        &
     &                                    GRID(ng)%z_r(i,j  ,k))
            END DO
          END DO
        END DO
        CALL extract_slice (ng, model, tile, v3dvar,                    &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      OCEAN(ng) % v(:,:,:,nrhs(ng)),              &
     &                      GRID(ng) % z_v,                             &
     &                      GRID(ng) % vmask_full,                      &
     &                      Zslice, Wr3d)
!
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idVzsl)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesZ_dp_v3dvar(ng)
        ELSE
          ioDesc => ioDesZ_sp_v3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idVzsl,          &
     &                     HIS(ng)%pioVar(idVzsl),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, Nslice, scale,        &
     &                     GRID(ng) % vmask_full,                       &
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, 3417, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVzsl)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out 3D Eastward momentum (m/s) at RHO-points, A-grid
!
      IF (Hout(idu3dE,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idu3dE)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idu3dE,          &
     &                     HIS(ng)%pioVar(idu3dE),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % rmask_full,                       &
     &                     OCEAN(ng) % ua)
        IF (FoundError(status, PIO_noerr, 3472, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idu3dE)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Eastward momentum (m/s) at RHO-points, A-grid, at
!  specified constant depth slices.
!
      IF (Hout(idUzsE,ng).and.(Nslice.gt.0)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
          Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
        END IF
        CALL extract_slice (ng, model, tile, r3dvar,                    &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      OCEAN(ng) % ua,                             &
     &                      GRID(ng) % z_r,                             &
     &                      GRID(ng) % rmask_full,                      &
     &                      Zslice, Wr3d)
!
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idUzsE)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesZ_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesZ_sp_r3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idUzsE,          &
     &                     HIS(ng)%pioVar(idUzsE),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, Nslice, scale,        &
     &                     GRID(ng) % rmask_full,                       &
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, 3514, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUzsE)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out 3D Nothward momentum (m/s) at RHO-points, A-grid
!
      IF (Hout(idv3dN,ng)) THEN
        IF (HIS(ng)%pioVar(idV3dN)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idv3dN,          &
     &                     HIS(ng)%pioVar(idv3dN),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % rmask_full,                       &
     &                     OCEAN(ng) % va)
        IF (FoundError(status, PIO_noerr, 3542, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idv3dN)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Northward momentum (m/s) at RHO-points, A-grid, at
!  specified constant depth slices.
!
      IF (Hout(idVzsN,ng).and.(Nslice.gt.0)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
          Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
        END IF
        CALL extract_slice (ng, model, tile, r3dvar,                    &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      OCEAN(ng) % va,                             &
     &                      GRID(ng) % z_r,                             &
     &                      GRID(ng) % rmask_full,                      &
     &                      Zslice, Wr3d)
!
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idv3dN)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesZ_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesZ_sp_r3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idVzsN,          &
     &                     HIS(ng)%pioVar(idVzsN),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, Nslice, scale,        &
     &                     GRID(ng) % rmask_full,                       &
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, 3584, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVzsN)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out potential vorticity (m-1 s-1) at PSI-points at full column
!  or specified constant depth slices.
!
      IF (Hout(id3dPV,ng).or.(Hout(idPVzs,ng).and.(Nslice.gt.0))) THEN
        IF (.not.allocated(Fr3d)) THEN
          allocate ( Fr3d(LBi:UBi,LBj:UBj,N(ng)) )
          Fr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
        END IF
        scale=1.0_dp
!
        CALL pvorticity3d (ng, model, tile,                             &
                           LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS, nrhs(ng),        &
     &                     GRID(ng) % pmask,                            &
     &                     GRID(ng) % umask,   GRID(ng) % vmask,        &
                           GRID(ng) % f,                                &
                           GRID(ng) % om_u,    GRID(ng) % on_v,         &
                           GRID(ng) % pm,      GRID(ng) % pn,           &
                           GRID(ng) % z_r,     OCEAN(ng) % pden,        &
                           OCEAN(ng) % u,      OCEAN(ng) % v,           &
                           Fr3d)
!
        IF (Hout(id3dPV,ng)) THEN
          IF (HIS(ng)%pioVar(id3dPV)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_p3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_p3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, id3dPV,        &
     &                       HIS(ng)%pioVar(id3dPV),                    &
     &                       HIS(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
     &                       GRID(ng) % pmask,                          &
     &                       Fr3d)
          IF (FoundError(status, PIO_noerr, 3634, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,id3dPV)), HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
        IF (Hout(idPVzs,ng).and.(Nslice.gt.0)) THEN
          IF (.not.allocated(Wr3d)) THEN
            allocate ( Wr3d(LBi:UBi,LBj:UBj,N(ng)) )
            Wr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
          END IF
!
          DO k=1,N(ng)
            DO j=JstrP,JendT
              DO i=IstrP,IendT
                GRID(ng)%z_v(i,j,k)=0.25_r8*(GRID(ng)%z_r(i-1,j-1,k)+   &
     &                                       GRID(ng)%z_r(i-1,j  ,k)+   &
     &                                       GRID(ng)%z_r(i  ,j-1,k)+   &
     &                                       GRID(ng)%z_r(i ,j  ,k))
              END DO
            END DO
          END DO
!
          CALL extract_slice (ng, model, tile, p3dvar,                  &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        Fr3d,                                     &
     &                        GRID(ng) % z_v,                           &
     &                        GRID(ng) % pmask,                         &
     &                        Zslice, Wr3d)
!
          IF (HIS(ng)%pioVar(idPVzs)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesZ_dp_p3dvar(ng)
          ELSE
            ioDesc => ioDesZ_sp_p3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idPVzs,        &
     &                       HIS(ng)%pioVar(idPVzs),                    &
     &                       HIS(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, Nslice, scale,      &
     &                       GRID(ng) % pmask,                          &
     &                       Wr3d)
          IF (FoundError(status, PIO_noerr, 3684, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idPVzs)), HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
          deallocate (Wr3d)
        END IF
        deallocate (Fr3d)
      END IF
!
!  Write out relative vorticity (s-1) at PSI-points at full column
!  or specified constant depth slices.
!
      IF (Hout(id3dRV,ng).or.(Hout(idRVzs,ng).and.(Nslice.gt.0))) THEN
        IF (.not.allocated(Fr3d)) THEN
          allocate ( Fr3d(LBi:UBi,LBj:UBj,N(ng)) )
          Fr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
        END IF
        scale=1.0_dp
!
        CALL rvorticity3d (ng, model, tile,                             &
                           LBi, UBi, LBj, UBj, nrhs(ng),                &
     &                     GRID(ng) % pmask,                            &
                           GRID(ng) % om_u, GRID(ng) % on_v,            &
                           GRID(ng) % pm,   GRID(ng) % pn,              &
                           OCEAN(ng) % u,   OCEAN(ng) % v,              &
                           Fr3d)
!
        IF (Hout(id3dRV,ng)) THEN
          IF (HIS(ng)%pioVar(id3dRV)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_p3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_p3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, id3dRV,        &
     &                       HIS(ng)%pioVar(id3dRV),                    &
     &                       HIS(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
     &                       GRID(ng) % pmask,                          &
     &                       Fr3d)
          IF (FoundError(status, PIO_noerr, 3732, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,id3dRV)), HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
        IF (Hout(idRVzs,ng).and.(Nslice.gt.0)) THEN
          IF (.not.allocated(Wr3d)) THEN
            allocate ( Wr3d(LBi:UBi,LBj:UBj,N(ng)) )
            Wr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
          END IF
!
          DO k=1,N(ng)
            DO j=JstrP,JendT
              DO i=IstrP,IendT
                GRID(ng)%z_v(i,j,k)=0.25_r8*(GRID(ng)%z_r(i-1,j-1,k)+   &
     &                                       GRID(ng)%z_r(i-1,j  ,k)+   &
     &                                       GRID(ng)%z_r(i  ,j-1,k)+   &
     &                                       GRID(ng)%z_r(i ,j  ,k))
              END DO
            END DO
          END DO
!
          CALL extract_slice (ng, model, tile, p3dvar,                  &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        Fr3d,                                     &
     &                        GRID(ng) % z_v,                           &
     &                        GRID(ng) % pmask,                         &
     &                        Zslice, Wr3d)
!
          IF (HIS(ng)%pioVar(idRVzs)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesZ_dp_p3dvar(ng)
          ELSE
            ioDesc => ioDesZ_sp_p3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idRVzs,        &
     &                       HIS(ng)%pioVar(idRVzs),                    &
     &                       HIS(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, Nslice, scale,      &
     &                       GRID(ng) % pmask,                          &
     &                       Wr3d)
          IF (FoundError(status, PIO_noerr, 3782, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idRVzs)), HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
          deallocate (Wr3d)
        END IF
        deallocate (Fr3d)
      END IF
!
!  Write out S-coordinate omega vertical velocity (m/s).
!
      IF (Hout(idOvel,ng)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,0:N(ng)) )
          Wr3d(LBi:UBi,LBj:UBj,0:N(ng))=0.0_r8
        END IF
        scale=1.0_dp
        CALL scale_omega (ng, tile, LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    OCEAN(ng) % W,                                &
     &                    Wr3d(LBi:,LBj:,0:))
!
        IF (HIS(ng)%pioVar(idOvel)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idOvel,          &
     &                     HIS(ng)%pioVar(idOvel),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     Wr3d)
        IF (FoundError(status, PIO_noerr, 3823, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idOvel)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out S-coordinate implicit omega vertical velocity (m/s).
!
      IF (Hout(idOvil,ng)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,0:N(ng)) )
          Wr3d(LBi:UBi,LBj:UBj,0:N(ng))=0.0_r8
        END IF
        scale=1.0_dp
        CALL scale_omega (ng, tile, LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    OCEAN(ng) % Wi,                               &
     &                    Wr3d(LBi:,LBj:,0:))
!
        IF (HIS(ng)%pioVar(idOvil)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idOvil,          &
     &                     HIS(ng)%pioVar(idOvil),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     Wr3d)
        IF (FoundError(status, PIO_noerr, 3864, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idOvil)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out vertical velocity (m/s).
!
      IF (Hout(idWvel,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idWvel)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idWvel,          &
     &                     HIS(ng)%pioVar(idWvel),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     OCEAN(ng) % wvel)
        IF (FoundError(status, PIO_noerr, 3894, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idWvel)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        IF (Hout(idTvar(itrc),ng)) THEN
          scale=1.0_dp
          IF (HIS(ng)%pioTrc(itrc)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_r3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idTvar(itrc),  &
     &                       HIS(ng)%pioTrc(itrc),                      &
     &                       HIS(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
     &                       GRID(ng) % rmask,                          &
     &                       OCEAN(ng) % t(:,:,:,nrhs(ng),itrc))
          IF (FoundError(status, PIO_noerr, 3923, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idTvar(itrc))),            &
     &                          HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out tracer type variables at specified constant depth slices.
!
      DO itrc=1,NT(ng)
        IF (Hout(idzslT(itrc),ng).and.(Nslice.gt.0)) THEN
          IF (.not.allocated(Wr3d)) THEN
            allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
            Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
          END IF
          CALL extract_slice (ng, model, tile, r3dvar,                  &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        OCEAN(ng) % t(:,:,:,nrhs(ng),itrc),       &
     &                        GRID(ng) % z_r,                           &
     &                        GRID(ng) % rmask,                         &
     &                        Zslice, Wr3d)
!
          scale=1.0_dp
          IF (HIS(ng)%pioVar(idzslT(itrc))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesZ_dp_r3dvar(ng)
          ELSE
            ioDesc => ioDesZ_sp_r3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idzslT(itrc),  &
     &                       HIS(ng)%pioVar(idzslT(itrc)),              &
     &                       HIS(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, Nslice, scale,      &
     &                       GRID(ng) % rmask,                          &
     &                       Wr3d)
          IF (FoundError(status, nf90_noerr, 3967, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idzslT(itrc))),            &
     &                          HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
          IF (itrc.eq.NT(ng)) deallocate (Wr3d)
        END IF
      END DO
!
!  Write out density anomaly.
!
      IF (Hout(idDano,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idDano)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idDano,          &
     &                     HIS(ng)%pioVar(idDano),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     OCEAN(ng) % rho)
        IF (FoundError(status, PIO_noerr, 4027, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idDano)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical viscosity coefficient.
!
      IF (Hout(idVvis,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idVvis)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idVvis,          &
     &                     HIS(ng)%pioVar(idVvis),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     MIXING(ng) % Akv,                            &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, 4149, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVvis)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      IF (Hout(idTdif,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idTdif)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idTdif,          &
     &                     HIS(ng)%pioVar(idTdif),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     MIXING(ng) % Akt(:,:,:,itemp),               &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, 4178, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTdif)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for salinity.
!
      IF (Hout(idSdif,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idSdif)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idSdif,          &
     &                     HIS(ng)%pioVar(idSdif),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     MIXING(ng) % Akt(:,:,:,isalt),               &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, 4209, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSdif)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out turbulent kinetic energy.
!
      IF (Hout(idMtke,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idMtke)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idMtke,          &
     &                     HIS(ng)%pioVar(idMtke),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     MIXING(ng) % tke(:,:,:,nrhs(ng)),            &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, 4240, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idMtke)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out turbulent length scale field.
!
      IF (Hout(idMtls,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idMtls)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, HIS(ng)%pioFile, idMtls,          &
     &                     HIS(ng)%pioVar(idMtls),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     MIXING(ng) % gls(:,:,:,nrhs(ng)),            &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, 4296, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idMtls)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface air pressure.
!
      IF (Hout(idPair,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idPair)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idPair,          &
     &                     HIS(ng)%pioVar(idPair),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     FORCES(ng) % Pair)
        IF (FoundError(status, PIO_noerr, 4379, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idPair)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface air temperature.
!
      IF (Hout(idTair,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idTair)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idTair,          &
     &                     HIS(ng)%pioVar(idTair),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     FORCES(ng) % Tair)
        IF (FoundError(status, PIO_noerr, 4409, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTair)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface winds.
!
      IF (Hout(idUair,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idUair)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idUair,          &
     &                     HIS(ng)%pioVar(idUair),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     FORCES(ng) % Uwind)
        IF (FoundError(status, PIO_noerr, 4439, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUair)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
      IF (Hout(idVair,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idVair)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idVair,          &
     &                     HIS(ng)%pioVar(idVair),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     FORCES(ng) % Vwind)
        IF (FoundError(status, PIO_noerr, 4465, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVair)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out Eastward/Northward surface wind (m/s) at RHO-points.
!
      IF (Hout(idUaiE,ng).and.Hout(idVaiN,ng)) THEN
        IF (.not.allocated(Ur2d)) THEN
          allocate (Ur2d(LBi:UBi,LBj:UBj))
          Ur2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        IF (.not.allocated(Vr2d)) THEN
          allocate (Vr2d(LBi:UBi,LBj:UBj))
          Vr2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        CALL uv_rotate2d (ng, tile, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
     &                    GRID(ng) % rmask_full,                        &
     &                    FORCES(ng) % Uwind,                           &
     &                    FORCES(ng) % Vwind,                           &
     &                    Ur2d, Vr2d)
!
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idUaiE)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idUaiE,          &
     &                     HIS(ng)%pioVar(idUaiE),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     Ur2d)
        IF (FoundError(status, PIO_noerr, 4512, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUaiE)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idVaiN)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idVaiN,          &
     &                     HIS(ng)%pioVar(idVaiN),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     Vr2d)
        IF (FoundError(status, PIO_noerr, 4536, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVaiN)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF
!
!  Write out surface active tracers fluxes.
!
      DO itrc=1,NAT
        IF (Hout(idTsur(itrc),ng)) THEN
          IF (itrc.eq.itemp) THEN
            scale=rho0*Cp                   ! Celsius m/s to W/m2
          ELSE IF (itrc.eq.isalt) THEN
            scale=1.0_dp
          END IF
          IF (HIS(ng)%pioVar(idTsur(itrc))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r2dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_r2dvar(ng)
          END IF
          status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idTsur(itrc),  &
     &                       HIS(ng)%pioVar(idTsur(itrc)),              &
     &                       HIS(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
     &                       GRID(ng) % rmask,                          &
     &                       FORCES(ng) % stflx(:,:,itrc))
          IF (FoundError(status, PIO_noerr, 4576, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idTsur(itrc))),            &
     &                          HIS(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out latent heat flux.
!
      IF (Hout(idLhea,ng)) THEN
        scale=rho0*Cp
        IF (HIS(ng)%pioVar(idLhea)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idLhea,          &
     &                     HIS(ng)%pioVar(idLhea),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     FORCES(ng) % lhflx)
        IF (FoundError(status, PIO_noerr, 4608, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idLhea)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out sensible heat flux.
!
      IF (Hout(idShea,ng)) THEN
        scale=rho0*Cp
        IF (HIS(ng)%pioVar(idShea)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idShea,          &
     &                     HIS(ng)%pioVar(idShea),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     FORCES(ng) % shflx)
        IF (FoundError(status, PIO_noerr, 4636, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idShea)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out net longwave radiation flux.
!
      IF (Hout(idLrad,ng)) THEN
        scale=rho0*Cp
        IF (HIS(ng)%pioVar(idLrad)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idLrad,          &
     &                     HIS(ng)%pioVar(idLrad),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     FORCES(ng) % lrflx)
        IF (FoundError(status, PIO_noerr, 4664, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idLrad)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out evaporation rate (kg/m2/s).
!
      IF (Hout(idevap,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idevap)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idevap,          &
     &                     HIS(ng)%pioVar(idevap),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     FORCES(ng) % evap)
        IF (FoundError(status, PIO_noerr, 4696, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idevap)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out precipitation rate (kg/m2/s).
!
      IF (Hout(idrain,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idrain)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idrain,          &
     &                     HIS(ng)%pioVar(idrain),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     FORCES(ng) % rain)
        IF (FoundError(status, PIO_noerr, 4724, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idrain)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out E-P (m/s).
!
      IF (Hout(idEmPf,ng)) THEN
        scale=1.0_dp
        IF (HIS(ng)%pioVar(idEmPf)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idEmPf,          &
     &                     HIS(ng)%pioVar(idEmPf),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     FORCES(ng) % stflux(:,:,isalt))
        IF (FoundError(status, PIO_noerr, 4754, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idEmPf)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out net shortwave radiation flux.
!
      IF (Hout(idSrad,ng)) THEN
        scale=rho0*Cp
        IF (HIS(ng)%pioVar(idSrad)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idSrad,          &
     &                     HIS(ng)%pioVar(idSrad),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     FORCES(ng) % srflx)
        IF (FoundError(status, PIO_noerr, 4784, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSrad)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface U-momentum stress.
!
      IF (Hout(idUsms,ng)) THEN
        scale=rho0                          ! m2/s2 to Pa
        IF (HIS(ng)%pioVar(idUsms)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idUsms,          &
     &                     HIS(ng)%pioVar(idUsms),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % umask,                            &
     &                     FORCES(ng) % sustr)
        IF (FoundError(status, PIO_noerr, 4818, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUsms)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface V-momentum stress.
!
      IF (Hout(idVsms,ng)) THEN
        scale=rho0
        IF (HIS(ng)%pioVar(idVsms)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idVsms,          &
     &                     HIS(ng)%pioVar(idVsms),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % vmask,                            &
     &                     FORCES(ng) % svstr)
        IF (FoundError(status, PIO_noerr, 4850, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVsms)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom U-momentum stress.
!
      IF (Hout(idUbms,ng)) THEN
        scale=-rho0
        IF (HIS(ng)%pioVar(idUbms)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idUbms,          &
     &                     HIS(ng)%pioVar(idUbms),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % umask,                            &
     &                     FORCES(ng) % bustr)
        IF (FoundError(status, PIO_noerr, 4878, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUbms)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom V-momentum stress.
!
      IF (Hout(idVbms,ng)) THEN
        scale=-rho0
        IF (HIS(ng)%pioVar(idVbms)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, HIS(ng)%pioFile, idVbms,          &
     &                     HIS(ng)%pioVar(idVbms),                      &
     &                     HIS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % vmask,                            &
     &                     FORCES(ng) % bvstr)
        IF (FoundError(status, PIO_noerr, 4906, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVbms)), HIS(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Synchronize history NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL pio_netcdf_sync (ng, model, HIS(ng)%name, HIS(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, 4984, MyFile)) RETURN
!
  10  FORMAT (2x,'WRT_HIS_PIO      - writing history', t42,             &
     &        'fields (Index=',i1,',',i1,') in record = ',i0)
  20  FORMAT (/,' WRT_HIS_PIO - error while writing variable: ',a,      &
     &        /,15x,'into history NetCDF file for time record: ',i0)
!
      RETURN
      END SUBROUTINE wrt_his_pio
!
      END MODULE wrt_his_mod
