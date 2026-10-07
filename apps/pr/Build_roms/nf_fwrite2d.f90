      MODULE nf_fwrite2d_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module writes out a generic floating point 2D array into       !
!  an output file using either the standard NetCDF library or the      !
!  Parallel-IO (PIO) library.                                          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number                                  !
!     model        Calling model identifier                            !
!     ncid         NetCDF file ID                                      !
!     ifield       Field metadata index (integer)                      !
!or   pioFile      PIO file descriptor structure, TYPE(file_desc_t)    !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
!     ncvarid      NetCDF variable ID                                  !
!or   pioVar       PIO variable descriptor structure, TYPE(My_VarDesc) !
!                    pioVar%vd     variable descriptor TYPE(Var_Desc_t)!
!                    pioVar%dkind  variable data kind                  !
!                    pioVar%gtype  variable C-gridtype                 !
!     tindex       NetCDF time record index to write                   !
!     gtype        Grid type. If negative, only write water points     !
!or   pioDesc      IO data decomposition descriptor, TYPE(IO_desc_t)   !
!     LBi          I-dimension Lower bound                             !
!     UBi          I-dimension Upper bound                             !
!     LBj          J-dimension Lower bound                             !
!     UBj          J-dimension Upper bound                             !
!     Amask        land/Sea mask, if any (real)                        !
!     Ascl         Factor to scale field before writing (real)         !
!     Adat         Field to write out (real)                           !
!     SetFillVal   Logical switch to set fill value in land areas      !
!                    (OPTIONAL)                                        !
!     ExtractField Field extraction flag (integer, OPTIONAL)           !
!                    ExtractField = 0   no extraction                  !
!                    ExtractField = 1   extraction by interpolation    !
!                    ExtractField > 1   extraction by decimation       !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     status       Error flag (integer)                                !
!     MinValue     Minimum value (real, OPTIONAL)                      !
!     MaxValue     Maximum value (real, OPTIONAL)                      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_scalars
!
      USE pack_field_mod,  ONLY : pack_field2d
!
      implicit none
!
      INTERFACE nf_fwrite2d
        MODULE PROCEDURE nf90_fwrite2d
        MODULE PROCEDURE pio_fwrite2d
      END INTERFACE nf_fwrite2d
!
      CONTAINS
!
!***********************************************************************
      FUNCTION nf90_fwrite2d (ng, model, ncid, ifield,                  &
     &                        ncvarid, tindex, gtype,                   &
     &                        LBi, UBi, LBj, UBj, Ascl,                 &
     &                        Amask,                                    &
     &                        Adat,                                     &
     &                        SetFillVal,                               &
     &                        ExtractField,                             &
     &                        MinValue, MaxValue) RESULT (status)
!***********************************************************************
!
      USE mod_netcdf
!
      USE distribute_mod, ONLY : mp_bcasti, mp_gather2d
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: SetFillVal
!
      integer, intent(in) :: ng, model, ncid, ncvarid, tindex, gtype
      integer, intent(in) :: ifield
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      integer, intent(in), optional :: ExtractField
!
      real(dp), intent(in) :: Ascl
!
      real(r8), intent(in) :: Amask(LBi:,LBj:)
      real(r8), intent(in) :: Adat(LBi:,LBj:)
      real(r8), intent(out), optional :: MinValue
      real(r8), intent(out), optional :: MaxValue
!
!  Local variable declarations.
!
      logical :: LandFill
!
      integer :: Extract_Flag
      integer :: Npts, i, tile
      integer :: status
      integer, dimension(3) :: start, total
!
      real(r8), dimension((Lm(ng)+2)*(Mm(ng)+2)) :: Awrk
!
!-----------------------------------------------------------------------
!  Initialize local variables.
!-----------------------------------------------------------------------
!
      status=nf90_noerr
!
!  Set parallel tile.
!
      tile=MyRank
!
!  Set switch to replace land areas with fill value, spval.
!
      IF (PRESENT(SetFillVal)) THEN
        LandFill=SetFillVal
      ELSE
        LandFill=tindex.gt.0
      END IF
!
!  If appropriate, set the field extraction flag to the provided grid
!  geometry through interpolation or decimation.
!
      IF (PRESENT(ExtractField)) THEN
        Extract_Flag=ExtractField
      ELSE
        Extract_Flag=0
      END IF
!
!  Initialize local array to avoid denormalized numbers. This
!  facilitates processing and debugging.
!
      Awrk=0.0_r8
!
!-----------------------------------------------------------------------
!  Pack 2D field data into 1D array.
!-----------------------------------------------------------------------
!
      CALL pack_field2d (ng, model, tile,                               &
     &                   gtype, ifield, tindex,                         &
     &                   LandFill, Extract_Flag,                        &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Amask,                                         &
     &                   Ascl, Adat,                                    &
     &                   start, total, Npts, Awrk)
!
!-----------------------------------------------------------------------
!  If applicable, compute output field minimum and maximum values.
!-----------------------------------------------------------------------
!
      IF (PRESENT(MinValue)) THEN
        IF (OutThread) THEN
          MinValue=spval
          MaxValue=-spval
          DO i=1,Npts
            IF (ABS(Awrk(i)).lt.spval) THEN
              MinValue=MIN(MinValue,Awrk(i))
              MaxValue=MAX(MaxValue,Awrk(i))
            END IF
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Write output buffer into NetCDF file.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
        status=nf90_put_var(ncid, ncvarid, Awrk, start, total)
      END IF
!
!-----------------------------------------------------------------------
!  Broadcast IO error flag to all nodes.
!-----------------------------------------------------------------------
!
      CALL mp_bcasti (ng, model, status)
!
      RETURN
      END FUNCTION nf90_fwrite2d
!
!***********************************************************************
      FUNCTION pio_fwrite2d (ng, model, pioFile, ifield,                &
     &                       pioVar, tindex, pioDesc,                   &
     &                       LBi, UBi, LBj, UBj, Ascl,                  &
     &                       Amask,                                     &
     &                       Adat,                                      &
     &                       SetFillVal,                                &
     &                       ExtractField,                              &
     &                       MinValue, MaxValue) RESULT (status)
!***********************************************************************
!
      USE mod_pio_netcdf
!
      USE distribute_mod, ONLY : mp_reduce
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: SetFillVal
!
      integer, intent(in) :: ng, model, tindex
      integer, intent(in) :: ifield
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      integer, intent(in), optional :: ExtractField
!
      real(dp), intent(in) :: Ascl
!
      real(r8), intent(in) :: Amask(LBi:,LBj:)
      real(r8), intent(in) :: Adat(LBi:,LBj:)
      real(r8), intent(out), optional :: MinValue
      real(r8), intent(out), optional :: MaxValue
!
      TYPE (File_desc_t), intent(inout) :: pioFile
      TYPE (IO_Desc_t),   intent(inout) :: pioDesc
      TYPE (My_VarDesc),  intent(inout) :: pioVar
!
!  Local variable declarations.
!
      logical :: LandFill, Lminmax
      logical,  pointer :: Lwater(:,:)
!
      integer :: Extract_Flag
      integer :: i, j, tile
      integer :: Imin, Imax, Jmin, Jmax
      integer :: Cgrid, dkind, ghost, gtype
      integer :: status
      integer, dimension(3) :: start, total
!
      real(r8), dimension(2) :: rbuffer
      real(r4), pointer :: Awrk4(:,:)
      real(r8), pointer :: Awrk8(:,:)
!
      character (len= 3), dimension(2) :: op_handle
!
!-----------------------------------------------------------------------
!  Set starting and ending indices to process.
!-----------------------------------------------------------------------
!
      status=PIO_noerr
!
      Awrk4 => NULL()
      Awrk8 => NULL()
!
!  Set first and last tile computational grid point according to the
!  staggered C-grid location. Ghost points are not included.
!
      ghost=0
      dkind=pioVar%dkind
      gtype=pioVar%gtype
!
      SELECT CASE (gtype)
        CASE (p2dvar, p3dvar)
          Cgrid=1                                         ! PSI-points
        CASE (r2dvar, r3dvar)
          Cgrid=2                                         ! RHO-points
        CASE (u2dvar, u3dvar)
          Cgrid=3                                         ! U-points
        CASE (v2dvar, v3dvar)
          Cgrid=4                                         ! V-points
        CASE DEFAULT
          Cgrid=2                                         ! RHO-points
      END SELECT
!
      Imin=BOUNDS(ng)%Imin(Cgrid,ghost,MyRank)
      Imax=BOUNDS(ng)%Imax(Cgrid,ghost,MyRank)
      Jmin=BOUNDS(ng)%Jmin(Cgrid,ghost,MyRank)
      Jmax=BOUNDS(ng)%Jmax(Cgrid,ghost,MyRank)
!
!  Set switch to compute minimum and maximum values.
!
      IF (PRESENT(MinValue)) THEN
        Lminmax=.TRUE.
        IF (.not.associated(Lwater)) THEN
           allocate ( Lwater(LBi:UBi,LBj:UBj) )
           Lwater=.TRUE.
        END IF
      ELSE
        Lminmax=.FALSE.
      END IF
!
!  Set switch to replace land areas with fill value, spval.
!
      IF (PRESENT(SetFillVal)) THEN
        LandFill=SetFillVal
      ELSE
        LandFill=tindex.gt.0
      END IF
!
!  If appropriate, set the field extraction flag to the provided grid
!  geometry through interpolation or decimation.
!
      IF (PRESENT(ExtractField)) THEN
        Extract_Flag=ExtractField
      ELSE
        Extract_Flag=0
      END IF
!
!-----------------------------------------------------------------------
!  Write out data into NetCDF file.
!-----------------------------------------------------------------------
!
!  Allocate, initialize and load data into local array used for
!  writing. Overwrite masked points with special value.
!
      IF (dkind.eq.PIO_double) THEN                 ! double precision
        IF (.not.associated(Awrk8)) THEN
          allocate ( Awrk8(LBi:UBi,LBj:UBj) )
          Awrk8=0.0_r8
        END IF
!
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            Awrk8(i,j)=Adat(i,j)*Ascl
            IF((Amask(i,j).eq.0.0_r8).and.LandFill) THEN
              Awrk8(i,j)=spval
              IF (Lminmax) Lwater(i,j)=.FALSE.
            END iF
          END DO
        END DO
        IF (Lminmax) THEN
          rbuffer(1)=MINVAL(Awrk8, MASK=Lwater)
          rbuffer(2)=MAXVAL(Awrk8, MASK=Lwater)
        END IF
      ELSE                                          ! single precision
        IF (.not.associated(Awrk4)) THEN
          allocate ( Awrk4(LBi:UBi,LBj:UBj) )
          Awrk4=0.0_r4
        END IF
!
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            Awrk4(i,j)=REAL(Adat(i,j)*Ascl, r4)
            IF((Amask(i,j).eq.0.0_r8).and.LandFill) THEN
              Awrk4(i,j)=REAL(spval, r4)
              IF (Lminmax) Lwater(i,j)=.FALSE.
            END iF
          END DO
        END DO
        IF (Lminmax) THEN
          rbuffer(1)=REAL(MINVAL(Awrk4, MASK=Lwater),r8)
          rbuffer(2)=REAL(MAXVAL(Awrk4, MASK=Lwater),r8)
        END IF
      END IF
!
!  Set unlimited time dimension record to write, if any.
!
      IF (tindex.gt.0) THEN
        CALL PIO_setframe (pioFile,                                     &
     &                     pioVar%vd,                                   &
     &                     INT(tindex, kind=PIO_OFFSET_KIND))
      END IF
!
!  Write out data into NetCDF.
!
      IF (dkind.eq.PIO_double) THEN                 ! double precision
        CALL PIO_write_darray (pioFile,                                 &
     &                         pioVar%vd,                               &
     &                         pioDesc,                                 &
     &                         Awrk8(Imin:Imax,Jmin:Jmax),              &
     &                         status)
      ELSE                                          ! single precision
        CALL PIO_write_darray (pioFile,                                 &
     &                         pioVar%vd,                               &
     &                         piodesc,                                 &
     &                         Awrk4(Imin:Imax,Jmin:Jmax),              &
     &                         status)
      END IF
!
!-----------------------------------------------------------------------
!  If applicable, compute global minimum and maximum values.
!-----------------------------------------------------------------------
!
      IF (Lminmax) THEN
        op_handle(1)='MIN'
        op_handle(2)='MAX'
        CALL mp_reduce (ng, model, 2, rbuffer, op_handle)
        MinValue=rbuffer(1)
        MaxValue=rbuffer(2)
        IF (associated(Lwater)) deallocate (Lwater)
      END IF
!
!  Deallocate local array.
!
      IF (dkind.eq.PIO_double) THEN
        IF (associated(Awrk8)) deallocate (Awrk8)
      ELSE
        IF (associated(Awrk4)) deallocate (Awrk4)
      END IF
!
      RETURN
      END FUNCTION pio_fwrite2d
!
      END MODULE nf_fwrite2d_mod
