      MODULE wrt_avg_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module writes model time-averaged fields into output file      !
!  using the standard NetCDF library or the Parallel-IO (PIO) library. !
!                                                                      !
!  Notice that only momentum is affected by the full time-averaged     !
!  masks.  If applicable, these mask contains information about        !
!  river runoff and time-dependent wetting and drying variations.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_average
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE nf_fwrite2d_mod,     ONLY : nf_fwrite2d
      USE nf_fwrite3d_mod,     ONLY : nf_fwrite3d
      USE strings_mod,         ONLY : FoundError
!
      implicit none
!
      PUBLIC  :: wrt_avg
      PRIVATE :: wrt_avg_nf90
      PRIVATE :: wrt_avg_pio
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wrt_avg (ng, tile)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/wrt_avg.F"
!
!-----------------------------------------------------------------------
!  Write out time-averaged fields according to IO type.
!-----------------------------------------------------------------------
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
      SELECT CASE (AVG(ng)%IOtype)
        CASE (io_nf90)
          CALL wrt_avg_nf90 (ng, iNLM, tile,                            &
     &                       LBi, UBi, LBj, UBj)
        CASE (io_pio)
          CALL wrt_avg_pio (ng, iNLM, tile,                             &
     &                      LBi, UBi, LBj, UBj)
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) AVG(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, 119, MyFile)) RETURN
!
  10  FORMAT (' WRT_AVG - Illegal output file type, io_type = ',i0,     &
     &        /,15x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE wrt_avg
!
!***********************************************************************
      SUBROUTINE wrt_avg_nf90 (ng, model, tile,                         &
     &                         LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: Fcount, gfactor, gtype, i, itrc, status
!
      real(dp) :: scale
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/wrt_avg.F"//", wrt_avg_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out time-averaged fields when appropriate.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, 154, MyFile)) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
      gfactor=1
!
!  Set time record index.
!
      AVG(ng)%Rindex=AVG(ng)%Rindex+1
      Fcount=AVG(ng)%load
      AVG(ng)%Nrec(Fcount)=AVG(ng)%Nrec(Fcount)+1
!
!  Report.
!
      IF (Master) WRITE (stdout,10) AVG(ng)%Rindex
!
!  Write out averaged time.
!
      CALL netcdf_put_fvar (ng, model, AVG(ng)%name,                    &
     &                      TRIM(Vname(1,idtime)), AVGtime(ng:),        &
     &                      (/AVG(ng)%Rindex/), (/1/),                  &
     &                      ncid = AVG(ng)%ncid,                        &
     &                      varid = AVG(ng)%Vid(idtime))
      IF (FoundError(exit_flag, NoError, 186, MyFile)) RETURN
!
!  Write out free-surface (m).
!
      IF (Aout(idFsur,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid, idFsur,             &
     &                     AVG(ng)%Vid(idFsur),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgzeta)
        IF (FoundError(status, nf90_noerr, 201, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idFsur)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      IF (Aout(idUbar,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid, idUbar,             &
     &                     AVG(ng)%Vid(idUbar),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % umask_avg,                        &
     &                     AVERAGE(ng) % avgu2d)
        IF (FoundError(status, nf90_noerr, 250, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUbar)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      IF (Aout(idVbar,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid, idVbar,             &
     &                     AVG(ng)%Vid(idVbar),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % vmask_avg,                        &
     &                     AVERAGE(ng) % avgv2d)
        IF (FoundError(status, nf90_noerr, 299, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVbar)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D Eastward momentum component (m/s) at RHO-points.
!
      IF (Aout(idu2dE,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid, idu2dE,             &
     &                     AVG(ng)%Vid(idu2dE),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_avg,                        &
     &                     AVERAGE(ng) % avgu2dE)
        IF (FoundError(status, nf90_noerr, 348, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idu2dE)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D Northward momentum component (m/s) at RHO-points.
!
      IF (Aout(idv2dN,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid, idv2dN,             &
     &                     AVG(ng)%Vid(idv2dN),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_avg,                        &
     &                     AVERAGE(ng) % avgv2dN)
        IF (FoundError(status, nf90_noerr, 371, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idv2dN)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      IF (Aout(idUvel,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid, idUvel,             &
     &                     AVG(ng)%Vid(idUvel),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % umask_avg,                        &
     &                     AVERAGE(ng) % avgu3d)
        IF (FoundError(status, nf90_noerr, 485, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUvel)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D momentum component (m/s) in the ETA-direction.
!
      IF (Aout(idVvel,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid, idVvel,             &
     &                     AVG(ng)%Vid(idVvel),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % vmask_avg,                        &
     &                     AVERAGE(ng) % avgv3d)
        IF (FoundError(status, nf90_noerr, 534, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVvel)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Eastward momentum component (m/s) at RHO-points.
!
      IF (Aout(idu3dE,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid, idu3dE,             &
     &                     AVG(ng)%Vid(idu3dE),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % rmask_avg,                        &
     &                     AVERAGE(ng) % avgu3dE)
        IF (FoundError(status, nf90_noerr, 583, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idu3dE)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Northward momentum component (m/s) at RHO-points.
!
      IF (Aout(idv3dN,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid, idv3dN,             &
     &                     AVG(ng)%Vid(idv3dN),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % rmask_avg,                        &
     &                     AVERAGE(ng) % avgv3dN)
        IF (FoundError(status, nf90_noerr, 606, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idv3dN)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out S-coordinate omega vertical velocity (m/s).
!
      IF (Aout(idOvel,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid, idOvel,             &
     &                     AVG(ng)%Vid(idOvel),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgw3d)
        IF (FoundError(status, nf90_noerr, 629, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idOvel)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out "true" vertical velocity (m/s).
!
      IF (Aout(idWvel,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid, idWvel,             &
     &                     AVG(ng)%Vid(idWvel),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgwvel)
        IF (FoundError(status, nf90_noerr, 652, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idOvel)), AVG(ng)%Rindex
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
        IF (Aout(idTvar(itrc),ng)) THEN
          scale=1.0_dp
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, model, AVG(ng)%ncid, idTvar(itrc),     &
     &                       AVG(ng)%Tid(itrc),                         &
     &                       AVG(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
     &                       GRID(ng) % rmask,                          &
     &                       AVERAGE(ng) % avgt(:,:,:,itrc))
          IF (FoundError(status, nf90_noerr, 676, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idTvar(itrc))),            &
     &                          AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out density anomaly.
!
      IF (Aout(idDano,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid, idDano,             &
     &                     AVG(ng)%Vid(idDano),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgrho)
        IF (FoundError(status, nf90_noerr, 730, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idDano)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D potential vorticity.
!
      IF (Aout(id2dPV,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*p2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid, id2dPV,             &
     &                     AVG(ng)%Vid(id2dPV),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % pmask,                            &
     &                     AVERAGE(ng) % avgpvor2d)
        IF (FoundError(status, nf90_noerr, 805, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,id2dPV)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D relative vorticity.
!
      IF (Aout(id2dRV,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*p2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid, id2dRV,             &
     &                     AVG(ng)%Vid(id2dRV),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % pmask,                            &
     &                     AVERAGE(ng) % avgrvor2d)
        IF (FoundError(status, nf90_noerr, 828, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,id2dRV)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D potential vorticity.
!
      IF (Aout(id3dPV,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*p3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid, id3dPV,             &
     &                     AVG(ng)%Vid(id3dPV),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % pmask,                            &
     &                     AVERAGE(ng) % avgpvor3d)
        IF (FoundError(status, nf90_noerr, 853, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,id3dPV)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D relative vorticity.
!
      IF (Aout(id3dRV,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*p3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid, id3dRV,             &
     &                     AVG(ng)%Vid(id3dRV),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % pmask,                            &
     &                     AVERAGE(ng) % avgrvor3d)
        IF (FoundError(status, nf90_noerr, 876, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,id3dRV)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out quadratic <zeta*zeta> term.
!
      IF (Aout(idZZav,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid, idZZav,             &
     &                     AVG(ng)%Vid(idZZav),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgZZ)
        IF (FoundError(status, nf90_noerr, 900, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idZZav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out quadratic <ubar*ubar> term.
!
      IF (Aout(idU2av,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid, idU2av,             &
     &                     AVG(ng)%Vid(idU2av),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % umask,                            &
     &                     AVERAGE(ng) % avgU2)
        IF (FoundError(status, nf90_noerr, 923, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idU2av)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out quadratic <vbar*vbar> term.
!
      IF (Aout(idV2av,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid, idV2av,             &
     &                     AVG(ng)%Vid(idV2av),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % vmask,                            &
     &                     AVERAGE(ng) % avgV2)
        IF (FoundError(status, nf90_noerr, 946, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idV2av)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out u-volume flux.
!
      IF (Aout(idHUav,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid, idHUav,             &
     &                     AVG(ng)%Vid(idHUav),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % umask,                            &
     &                     AVERAGE(ng) % avgHuon)
        IF (FoundError(status, nf90_noerr, 971, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idHUav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out v-volume flux.
!
      IF (Aout(idHVav,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid, idHVav,             &
     &                     AVG(ng)%Vid(idHVav),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % vmask,                            &
     &                     AVERAGE(ng) % avgHvom)
        IF (FoundError(status, nf90_noerr, 994, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idHVav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out quadratic <u*u> term.
!
      IF (Aout(idUUav,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid, idUUav,             &
     &                     AVG(ng)%Vid(idUUav),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % umask,                            &
     &                     AVERAGE(ng) % avgUU)
        IF (FoundError(status, nf90_noerr, 1017, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUUav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out quadratic <u*v> term.
!
      IF (Aout(idUVav,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid, idUVav,             &
     &                     AVG(ng)%Vid(idUVav),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgUV)
        IF (FoundError(status, nf90_noerr, 1040, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUVav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out quadratic <v*v> term.
!
      IF (Aout(idVVav,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid, idVVav,             &
     &                     AVG(ng)%Vid(idVVav),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % vmask,                            &
     &                     AVERAGE(ng) % avgVV)
        IF (FoundError(status, nf90_noerr, 1063, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVVav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out quadratic <t*t> term.
!
      DO i=1,NT(ng)
        IF (Aout(idTTav(i),ng)) THEN
          scale=1.0_dp
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, model, AVG(ng)%ncid, idTTav(i),        &
     &                       AVG(ng)%Vid(idTTav(i)),                    &
     &                       AVG(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
     &                       GRID(ng) % rmask,                          &
     &                       AVERAGE(ng) % avgTT(:,:,:,i))
          IF (FoundError(status, nf90_noerr, 1087, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idTTav(i))), AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out active tracer volume fluxes.
!
      DO i=1,NT(ng)
        IF (Aout(iHUTav(i),ng)) THEN
          scale=1.0_dp
          gtype=gfactor*u3dvar
          status=nf_fwrite3d(ng, model, AVG(ng)%ncid, iHUTav(i),        &
     &                       AVG(ng)%Vid(iHUTav(i)),                    &
     &                       AVG(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
     &                       GRID(ng) % umask,                          &
     &                       AVERAGE(ng) % avgHuonT(:,:,:,i))
          IF (FoundError(status, nf90_noerr, 1112, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,iHUTav(i))), AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
        IF (Aout(iHVTav(i),ng)) THEN
          scale=1.0_dp
          gtype=gfactor*v3dvar
          status=nf_fwrite3d(ng, model, AVG(ng)%ncid, iHVTav(i),        &
     &                       AVG(ng)%Vid(iHVTav(i)),                    &
     &                       AVG(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
     &                       GRID(ng) % vmask,                          &
     &                       AVERAGE(ng) % avgHvomT(:,:,:,i))
          IF (FoundError(status, nf90_noerr, 1133, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,iHVTav(i))), AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out quadratic <u*t> and <v*t> terms.
!
      DO i=1,NT(ng)
        IF (Aout(idUTav(i),ng)) THEN
          scale=1.0_dp
          gtype=gfactor*u3dvar
          status=nf_fwrite3d(ng, model, AVG(ng)%ncid, idUTav(i),        &
     &                       AVG(ng)%Vid(idUTav(i)),                    &
     &                       AVG(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
     &                       GRID(ng) % umask,                          &
     &                       AVERAGE(ng) % avgUT(:,:,:,i))
          IF (FoundError(status, nf90_noerr, 1158, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idUTav(i))), AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
        IF (Aout(idVTav(i),ng)) THEN
          scale=1.0_dp
          gtype=gfactor*v3dvar
          status=nf_fwrite3d(ng, model, AVG(ng)%ncid, idVTav(i),        &
     &                       AVG(ng)%Vid(idVTav(i)),                    &
     &                       AVG(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
     &                       GRID(ng) % vmask,                          &
     &                       AVERAGE(ng) % avgVT(:,:,:,i))
          IF (FoundError(status, nf90_noerr, 1179, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idVTav(i))), AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out vertical viscosity coefficient.
!
      IF (Aout(idVvis,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid, idVvis,             &
     &                     AVG(ng)%Vid(idVvis),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgAKv,                        &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, 1206, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVvis)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      IF (Aout(idTdif,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid, idTdif,             &
     &                     AVG(ng)%Vid(idTdif),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgAKt,                        &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, 1230, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTdif)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for salinity.
!
      IF (Aout(idSdif,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, AVG(ng)%ncid, idSdif,             &
     &                     AVG(ng)%Vid(idSdif),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgAKs,                        &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, 1256, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSdif)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface air pressure.
!
      IF (Aout(idPair,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid, idPair,             &
     &                     AVG(ng)%Vid(idPair),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgPair)
        IF (FoundError(status, nf90_noerr, 1282, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idPair)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface air temperature.
!
      IF (Aout(idTair,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid, idTair,             &
     &                     AVG(ng)%Vid(idTair),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgTair)
        IF (FoundError(status, nf90_noerr, 1307, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTair)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface winds.
!
      IF (Aout(idUair,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid, idUair,             &
     &                     AVG(ng)%Vid(idUair),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgUwind)
        IF (FoundError(status, nf90_noerr, 1332, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUair)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
      IF (Aout(idVair,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid, idVair,             &
     &                     AVG(ng)%Vid(idVair),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgVwind)
        IF (FoundError(status, nf90_noerr, 1353, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVair)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out Eastward/Northward surface wind (m/s) at RHO-points.
!
      IF (Aout(idUaiE,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid, idUaiE,             &
     &                     AVG(ng)%Vid(idUaiE),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_avg,                        &
     &                     AVERAGE(ng) % avguwindE)
        IF (FoundError(status, nf90_noerr, 1376, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUaiE)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
      IF (Aout(idVaiN,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid, idVaiN,             &
     &                     AVG(ng)%Vid(idVaiN),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_avg,                        &
     &                     AVERAGE(ng) % avgvwindN)
        IF (FoundError(status, nf90_noerr, 1397, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVaiN)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface net heat flux.
!
      IF (Aout(idTsur(itemp),ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid, idTsur(itemp),      &
     &                     AVG(ng)%Vid(idTsur(itemp)),                  &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgstf)
        IF (FoundError(status, nf90_noerr, 1430, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTsur(itemp))),             &
     &                        AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface salt flux  (PSU m/s = kg salt/m2/s).
!
      IF (Aout(idTsur(isalt),ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid, idTsur(isalt),      &
     &                     AVG(ng)%Vid(idTsur(isalt)),                  &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgswf)
        IF (FoundError(status, nf90_noerr, 1456, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTsur(isalt))),             &
     &                        AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out latent heat flux.
!
      IF (Aout(idLhea,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid, idLhea,             &
     &                     AVG(ng)%Vid(idLhea),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avglhf)
        IF (FoundError(status, nf90_noerr, 1487, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idLhea)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out sensible heat flux.
!
      IF (Aout(idShea,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid, idShea,             &
     &                     AVG(ng)%Vid(idShea),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgshf)
        IF (FoundError(status, nf90_noerr, 1514, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idShea)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out net longwave radiation flux.
!
      IF (Aout(idLrad,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid, idLrad,             &
     &                     AVG(ng)%Vid(idLrad),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avglrf)
        IF (FoundError(status, nf90_noerr, 1541, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idLrad)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out evaportaion rate (kg/m2/s).
!
      IF (Aout(idevap,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid, idevap,             &
     &                     AVG(ng)%Vid(idevap),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgevap)
        IF (FoundError(status, nf90_noerr, 1568, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idevap)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out precipitation rate (kg/m2/s).
!
      IF (Aout(idrain,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid, idrain,             &
     &                     AVG(ng)%Vid(idrain),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgrain)
        IF (FoundError(status, nf90_noerr, 1591, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idrain)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out net shortwave radiation flux.
!
      IF (Aout(idSrad,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid, idSrad,             &
     &                     AVG(ng)%Vid(idSrad),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgsrf)
        IF (FoundError(status, nf90_noerr, 1621, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSrad)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface u-momentum stress.
!
      IF (Aout(idUsms,ng)) THEN
        scale=rho0
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid, idUsms,             &
     &                     AVG(ng)%Vid(idUsms),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % umask,                            &
     &                     AVERAGE(ng) % avgsus)
        IF (FoundError(status, nf90_noerr, 1655, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUsms)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface v-momentum stress.
!
      IF (Aout(idVsms,ng)) THEN
        scale=rho0
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid, idVsms,             &
     &                     AVG(ng)%Vid(idVsms),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % vmask,                            &
     &                     AVERAGE(ng) % avgsvs)
        IF (FoundError(status, nf90_noerr, 1687, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVsms)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom u-momentum stress.
!
      IF (Aout(idUbms,ng)) THEN
        scale=rho0
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid, idUbms,             &
     &                     AVG(ng)%Vid(idUbms),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % umask,                            &
     &                     AVERAGE(ng) % avgbus)
        IF (FoundError(status, nf90_noerr, 1719, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUbms)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom v-momentum stress.
!
      IF (Aout(idVbms,ng)) THEN
        scale=rho0
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, model, AVG(ng)%ncid, idVbms,             &
     &                     AVG(ng)%Vid(idVbms),                         &
     &                     AVG(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % vmask,                            &
     &                     AVERAGE(ng) % avgbvs)
        IF (FoundError(status, nf90_noerr, 1751, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVbms)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Synchronize time-average NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, model, AVG(ng)%name, AVG(ng)%ncid)
      IF (FoundError(exit_flag, NoError, 1817, MyFile)) RETURN
!
  10  FORMAT (2x,'WRT_AVG_NF90     - writing averaged',t42,'fields',t61,&
     &        'in record = ',i0)
  20  FORMAT (/,' WRT_AVG_NF90 - error while writing variable: ',a,     &
     &        /,16x,'into averages NetCDF file for time record: ',i0)
!
      RETURN
      END SUBROUTINE wrt_avg_nf90
!
!***********************************************************************
      SUBROUTINE wrt_avg_pio (ng, model, tile,                          &
     &                        LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: Fcount, i, itrc, status
!
      real(dp) :: scale
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/wrt_avg.F"//", wrt_avg_pio"
!
      TYPE (IO_desc_t), pointer :: ioDesc
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out time-averaged fields when appropriate.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, 1863, MyFile)) RETURN
!
!  Set time record index.
!
      AVG(ng)%Rindex=AVG(ng)%Rindex+1
      Fcount=AVG(ng)%load
      AVG(ng)%Nrec(Fcount)=AVG(ng)%Nrec(Fcount)+1
!
!  Report.
!
      IF (Master) WRITE (stdout,10) AVG(ng)%Rindex
!
!  Write out averaged time.
!
      CALL pio_netcdf_put_fvar (ng, model, AVG(ng)%name,                &
     &                          TRIM(Vname(1,idtime)), AVGtime(ng:),    &
     &                          (/AVG(ng)%Rindex/), (/1/),              &
     &                          pioFile = AVG(ng)%pioFile,              &
     &                          pioVar = AVG(ng)%pioVar(idtime)%vd)
      IF (FoundError(exit_flag, NoError, 1886, MyFile)) RETURN
!
!  Write out free-surface (m).
!
      IF (Aout(idFsur,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idFsur)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile, idFsur,          &
     &                     AVG(ng)%pioVar(idFsur),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgzeta)
        IF (FoundError(status, nf90_noerr, 1906, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idFsur)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      IF (Aout(idUbar,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idUbar)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile, idUbar,          &
     &                     AVG(ng)%pioVar(idUbar),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % umask_avg,                        &
     &                     AVERAGE(ng) % avgu2d)
        IF (FoundError(status, nf90_noerr, 1965, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUbar)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      IF (Aout(idVbar,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idVbar)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile, idVbar,          &
     &                     AVG(ng)%pioVar(idVbar),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % vmask_avg,                        &
     &                     AVERAGE(ng) % avgv2d)
        IF (FoundError(status, nf90_noerr, 2024, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVbar)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D Eastward momentum component (m/s) at RHO-points.
!
      IF (Aout(idu2dE,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idu2dE)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile, idu2dE,          &
     &                     AVG(ng)%pioVar(idu2dE),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_avg,                        &
     &                     AVERAGE(ng) % avgu2dE)
        IF (FoundError(status, nf90_noerr, 2083, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idu2dE)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D Northward momentum component (m/s) at RHO-points.
!
      IF (Aout(idv2dN,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idv2dN)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile, idv2dN,          &
     &                     AVG(ng)%pioVar(idv2dN),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_avg,                        &
     &                     AVERAGE(ng) % avgv2dN)
        IF (FoundError(status, nf90_noerr, 2111, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idv2dN)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      IF (Aout(idUvel,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idUvel)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile, idUvel,          &
     &                     AVG(ng)%pioVar(idUvel),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % umask_avg,                        &
     &                     AVERAGE(ng) % avgu3d)
        IF (FoundError(status, nf90_noerr, 2250, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUvel)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D momentum component (m/s) in the ETA-direction.
!
      IF (Aout(idVvel,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idVvel)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile, idVvel,          &
     &                     AVG(ng)%pioVar(idVvel),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % vmask_avg,                        &
     &                     AVERAGE(ng) % avgv3d)
        IF (FoundError(status, nf90_noerr, 2309, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVvel)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Eastward momentum component (m/s) at RHO-points.
!
      IF (Aout(idu3dE,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idu3dE)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile, idu3dE,          &
     &                     AVG(ng)%pioVar(idu3dE),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % rmask_avg,                        &
     &                     AVERAGE(ng) % avgu3dE)
        IF (FoundError(status, nf90_noerr, 2368, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idu3dE)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Northward momentum component (m/s) at RHO-points.
!
      IF (Aout(idv3dN,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idv3dN)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile, idv3dN,          &
     &                     AVG(ng)%pioVar(idv3dN),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % rmask_avg,                        &
     &                     AVERAGE(ng) % avgv3dN)
        IF (FoundError(status, nf90_noerr, 2396, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idv3dN)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out S-coordinate omega vertical velocity (m/s).
!
      IF (Aout(idOvel,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idOvel)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile, idOvel,          &
     &                     AVG(ng)%pioVar(idOvel),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgw3d)
        IF (FoundError(status, nf90_noerr, 2424, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idOvel)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out "true" vertical velocity (m/s).
!
      IF (Aout(idWvel,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idWvel)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile, idWvel,          &
     &                     AVG(ng)%pioVar(idWvel),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgwvel)
        IF (FoundError(status, nf90_noerr, 2452, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idOvel)), AVG(ng)%Rindex
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
        IF (Aout(idTvar(itrc),ng)) THEN
          scale=1.0_dp
          IF (AVG(ng)%pioTrc(itrc)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_r3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, AVG(ng)%pioFile, idTvar(itrc),  &
     &                       AVG(ng)%pioTrc(itrc),                      &
     &                       AVG(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
     &                       GRID(ng) % rmask,                          &
     &                       AVERAGE(ng) % avgt(:,:,:,itrc))
          IF (FoundError(status, nf90_noerr, 2481, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idTvar(itrc))),            &
     &                          AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out density anomaly.
!
      IF (Aout(idDano,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idDano)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile, idDano,          &
     &                     AVG(ng)%pioVar(idDano),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgrho)
        IF (FoundError(status, nf90_noerr, 2545, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idDano)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D potential vorticity.
!
      IF (Aout(id2dPV,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(id2dPV)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_p2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_p2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile, id2dPV,          &
     &                     AVG(ng)%pioVar(id2dPV),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % pmask,                            &
     &                     AVERAGE(ng) % avgpvor2d)
        IF (FoundError(status, nf90_noerr, 2635, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,id2dPV)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D relative vorticity.
!
      IF (Aout(id2dRV,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(id2dRV)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_p2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_p2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile, id2dRV,          &
     &                     AVG(ng)%pioVar(id2dRV),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % pmask,                            &
     &                     AVERAGE(ng) % avgrvor2d)
        IF (FoundError(status, nf90_noerr, 2663, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,id2dRV)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D potential vorticity.
!
      IF (Aout(id3dPV,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(id3dPV)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_p3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_p3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile, id3dPV,          &
     &                     AVG(ng)%pioVar(id3dPV),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % pmask,                            &
     &                     AVERAGE(ng) % avgpvor3d)
        IF (FoundError(status, nf90_noerr, 2693, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,id3dPV)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D relative vorticity.
!
      IF (Aout(id3dRV,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(id3dRV)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_p3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_p3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile, id3dRV,          &
     &                     AVG(ng)%pioVar(id3dRV),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % pmask,                            &
     &                     AVERAGE(ng) % avgrvor3d)
        IF (FoundError(status, nf90_noerr, 2721, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,id3dRV)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out quadratic <zeta*zeta> term.
!
      IF (Aout(idZZav,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idZZav)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile, idZZav,          &
     &                     AVG(ng)%pioVar(idZZav),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgZZ)
        IF (FoundError(status, nf90_noerr, 2750, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idZZav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out quadratic <ubar*ubar> term.
!
      IF (Aout(idU2av,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idU2av)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile, idU2av,          &
     &                     AVG(ng)%pioVar(idU2av),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % umask,                            &
     &                     AVERAGE(ng) % avgU2)
        IF (FoundError(status, nf90_noerr, 2778, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idU2av)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out quadratic <vbar*vbar> term.
!
      IF (Aout(idV2av,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idV2av)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile, idV2av,          &
     &                     AVG(ng)%pioVar(idV2av),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % vmask,                            &
     &                     AVERAGE(ng) % avgV2)
        IF (FoundError(status, nf90_noerr, 2806, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idV2av)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out u-volume flux.
!
      IF (Aout(idHUav,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idHUav)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile, idHUav,          &
     &                     AVG(ng)%pioVar(idHUav),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % umask,                            &
     &                     AVERAGE(ng) % avgHuon)
        IF (FoundError(status, nf90_noerr, 2836, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idHUav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out v-volume flux.
!
      IF (Aout(idHVav,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idHVav)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile, idHVav,          &
     &                     AVG(ng)%pioVar(idHVav),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % vmask,                            &
     &                     AVERAGE(ng) % avgHvom)
        IF (FoundError(status, nf90_noerr, 2864, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idHVav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out quadratic <u*u> term.
!
      IF (Aout(idUUav,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idUUav)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile, idUUav,          &
     &                     AVG(ng)%pioVar(idUUav),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % umask,                            &
     &                     AVERAGE(ng) % avgUU)
        IF (FoundError(status, nf90_noerr, 2892, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUUav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out quadratic <u*v> term.
!
      IF (Aout(idUVav,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idUVav)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile, idUVav,          &
     &                     AVG(ng)%pioVar(idUVav),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgUV)
        IF (FoundError(status, nf90_noerr, 2920, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUVav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out quadratic <v*v> term.
!
      IF (Aout(idVVav,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idVVav)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile, idVVav,          &
     &                     AVG(ng)%pioVar(idVVav),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % vmask,                            &
     &                     AVERAGE(ng) % avgVV)
        IF (FoundError(status, nf90_noerr, 2948, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVVav)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out quadratic <t*t> term.
!
      DO i=1,NT(ng)
        IF (Aout(idTTav(i),ng)) THEN
          scale=1.0_dp
          IF (AVG(ng)%pioVar(idTTav(i))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_r3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, AVG(ng)%pioFile, idTTav(i),     &
     &                       AVG(ng)%pioVar(idTTav(i)),                 &
     &                       AVG(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
     &                       GRID(ng) % rmask,                          &
     &                       AVERAGE(ng) % avgTT(:,:,:,i))
          IF (FoundError(status, nf90_noerr, 2977, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idTTav(i))), AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out active tracer volume fluxes.
!
      DO i=1,NT(ng)
        IF (Aout(iHUTav(i),ng)) THEN
          scale=1.0_dp
          IF (AVG(ng)%pioVar(iHUTav(i))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_u3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_u3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, AVG(ng)%pioFile, iHUTav(i),     &
     &                       AVG(ng)%pioVar(iHUTav(i)),                 &
     &                       AVG(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
     &                       GRID(ng) % umask,                          &
     &                       AVERAGE(ng) % avgHuonT(:,:,:,i))
          IF (FoundError(status, nf90_noerr, 3007, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,iHUTav(i))), AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
        IF (Aout(iHVTav(i),ng)) THEN
          scale=1.0_dp
          IF (AVG(ng)%pioVar(iHVTav(i))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_v3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_v3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, AVG(ng)%pioFile, iHVTav(i),     &
     &                       AVG(ng)%pioVar(iHVTav(i)),                 &
     &                       AVG(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
     &                       GRID(ng) % vmask,                          &
     &                       AVERAGE(ng) % avgHvomT(:,:,:,i))
          IF (FoundError(status, nf90_noerr, 3033, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,iHVTav(i))), AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out quadratic <u*t> and <v*t> terms.
!
      DO i=1,NT(ng)
        IF (Aout(idUTav(i),ng)) THEN
          scale=1.0_dp
          IF (AVG(ng)%pioVar(idUTav(i))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_u3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_u3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, AVG(ng)%pioFile, idUTav(i),     &
     &                       AVG(ng)%pioVar(idUTav(i)),                 &
     &                       AVG(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
     &                       GRID(ng) % umask,                          &
     &                       AVERAGE(ng) % avgUT(:,:,:,i))
          IF (FoundError(status, nf90_noerr, 3063, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idUTav(i))), AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
        IF (Aout(idVTav(i),ng)) THEN
          scale=1.0_dp
          IF (AVG(ng)%pioVar(idVTav(i))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_v3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_v3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, AVG(ng)%pioFile, idVTav(i),     &
     &                       AVG(ng)%pioVar(idVTav(i)),                 &
     &                       AVG(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
     &                       GRID(ng) % vmask,                          &
     &                       AVERAGE(ng) % avgVT(:,:,:,i))
          IF (FoundError(status, nf90_noerr, 3089, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idVTav(i))), AVG(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out vertical viscosity coefficient.
!
      IF (Aout(idVvis,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idVvis)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile, idVvis,          &
     &                     AVG(ng)%pioVar(idVvis),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgAKv,                        &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, 3121, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVvis)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      IF (Aout(idTdif,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idTdif)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile, idTdif,          &
     &                     AVG(ng)%pioVar(idTdif),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgAKt,                        &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, 3150, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTdif)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for salinity.
!
      IF (Aout(idSdif,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idSdif)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, AVG(ng)%pioFile, idSdif,          &
     &                     AVG(ng)%pioVar(idSdif),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgAKs,                        &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, 3181, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSdif)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface air pressure.
!
      IF (Aout(idPair,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idPair)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile, idPair,          &
     &                     AVG(ng)%pioVar(idPair),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgPair)
        IF (FoundError(status, nf90_noerr, 3212, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idPair)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface air temperature.
!
      IF (Aout(idTair,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idTair)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile, idTair,          &
     &                     AVG(ng)%pioVar(idTair),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgTair)
        IF (FoundError(status, nf90_noerr, 3242, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTair)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface winds.
!
      IF (Aout(idUair,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idUair)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile, idUair,          &
     &                     AVG(ng)%pioVar(idUair),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgUwind)
        IF (FoundError(status, nf90_noerr, 3272, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUair)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
      IF (Aout(idVair,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idVair)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile, idVair,          &
     &                     AVG(ng)%pioVar(idVair),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgVwind)
        IF (FoundError(status, nf90_noerr, 3298, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVair)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out Eastward/Northward surface wind (m/s) at RHO-points.
!
      IF (Aout(idUaiE,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idUaiE)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile, idUaiE,          &
     &                     AVG(ng)%pioVar(idUaiE),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_avg,                        &
     &                     AVERAGE(ng) % avguwindE)
        IF (FoundError(status, nf90_noerr, 3326, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUaiE)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
      IF (Aout(idVaiN,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idVaiN)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile, idVaiN,          &
     &                     AVG(ng)%pioVar(idVaiN),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask_avg,                        &
     &                     AVERAGE(ng) % avgvwindN)
        IF (FoundError(status, nf90_noerr, 3352, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVaiN)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface net heat flux.
!
      IF (Aout(idTsur(itemp),ng)) THEN
        scale=rho0*Cp
        IF (AVG(ng)%pioVar(idTsur(itemp))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile, idTsur(itemp),   &
     &                     AVG(ng)%pioVar(idTsur(itemp)),               &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgstf)
        IF (FoundError(status, nf90_noerr, 3390, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTsur(itemp))),             &
     &                        AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface salt flux  (PSU m/s = kg salt/m2/s).
!
      IF (Aout(idTsur(isalt),ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idTsur(isalt))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile, idTsur(isalt),   &
     &                     AVG(ng)%pioVar(idTsur(isalt)),               &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgswf)
        IF (FoundError(status, nf90_noerr, 3421, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTsur(isalt))),             &
     &                        AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out latent heat flux.
!
      IF (Aout(idLhea,ng)) THEN
        scale=rho0*Cp
        IF (AVG(ng)%pioVar(idLhea)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile, idLhea,          &
     &                     AVG(ng)%pioVar(idLhea),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avglhf)
        IF (FoundError(status, nf90_noerr, 3457, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idLhea)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out sensible heat flux.
!
      IF (Aout(idShea,ng)) THEN
        scale=rho0*Cp
        IF (AVG(ng)%pioVar(idShea)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile, idShea,          &
     &                     AVG(ng)%pioVar(idShea),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgshf)
        IF (FoundError(status, nf90_noerr, 3489, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idShea)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out net longwave radiation flux.
!
      IF (Aout(idLrad,ng)) THEN
        scale=rho0*Cp
        IF (AVG(ng)%pioVar(idLrad)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile, idLrad,          &
     &                     AVG(ng)%pioVar(idLrad),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avglrf)
        IF (FoundError(status, nf90_noerr, 3521, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idLrad)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out evaportaion rate (kg/m2/s).
!
      IF (Aout(idevap,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idevap)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile, idevap,          &
     &                     AVG(ng)%pioVar(idevap),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgevap)
        IF (FoundError(status, nf90_noerr, 3553, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idevap)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out precipitation rate (kg/m2/s).
!
      IF (Aout(idrain,ng)) THEN
        scale=1.0_dp
        IF (AVG(ng)%pioVar(idrain)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile, idrain,          &
     &                     AVG(ng)%pioVar(idrain),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgrain)
        IF (FoundError(status, nf90_noerr, 3581, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idrain)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out net shortwave radiation flux.
!
      IF (Aout(idSrad,ng)) THEN
        scale=rho0*Cp
        IF (AVG(ng)%pioVar(idSrad)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile, idSrad,          &
     &                     AVG(ng)%pioVar(idSrad),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     AVERAGE(ng) % avgsrf)
        IF (FoundError(status, nf90_noerr, 3616, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSrad)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface u-momentum stress.
!
      IF (Aout(idUsms,ng)) THEN
        scale=rho0
        IF (AVG(ng)%pioVar(idUsms)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile, idUsms,          &
     &                     AVG(ng)%pioVar(idUsms),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % umask,                            &
     &                     AVERAGE(ng) % avgsus)
        IF (FoundError(status, nf90_noerr, 3655, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUsms)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface v-momentum stress.
!
      IF (Aout(idVsms,ng)) THEN
        scale=rho0
        IF (AVG(ng)%pioVar(idVsms)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile, idVsms,          &
     &                     AVG(ng)%pioVar(idVsms),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % vmask,                            &
     &                     AVERAGE(ng) % avgsvs)
        IF (FoundError(status, nf90_noerr, 3692, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVsms)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom u-momentum stress.
!
      IF (Aout(idUbms,ng)) THEN
        scale=rho0
        IF (AVG(ng)%pioVar(idUbms)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile, idUbms,          &
     &                     AVG(ng)%pioVar(idUbms),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % umask,                            &
     &                     AVERAGE(ng) % avgbus)
        IF (FoundError(status, nf90_noerr, 3729, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUbms)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom v-momentum stress.
!
      IF (Aout(idVbms,ng)) THEN
        scale=rho0
        IF (AVG(ng)%pioVar(idVbms)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, AVG(ng)%pioFile, idVbms,          &
     &                     AVG(ng)%pioVar(idVbms),                      &
     &                     AVG(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % vmask,                            &
     &                     AVERAGE(ng) % avgbvs)
        IF (FoundError(status, nf90_noerr, 3766, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVbms)), AVG(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Synchronize time-average NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL pio_netcdf_sync (ng, model, AVG(ng)%name, AVG(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, 3832, MyFile)) RETURN
!
  10  FORMAT (2x,'WRT_AVG_PIO      - writing averaged',t42,'fields',t61,&
     &        'in record = ',i0)
  20  FORMAT (/,' WRT_AVG_PIO - error while writing variable: ',a,      &
     &        /,16x,'into averages NetCDF file for time record: ',i0)
!
      RETURN
      END SUBROUTINE wrt_avg_pio
      END MODULE wrt_avg_mod
