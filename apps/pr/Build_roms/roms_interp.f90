      MODULE roms_interp_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module set the basic ROMS interpolation object for regridding  !
!  from source to destination grid geometries.                         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE distribute_mod,       ONLY : mp_reduce
      USE roms_interpolate_mod, ONLY : hindices
      USE stats_mod,            ONLY : stats_2dfld,                     &
     &                                 stats_3dfld
      USE strings_mod,          ONLY : assign_string,                   &
     &                                 find_string,                     &
     &                                 FoundError
!
      implicit none
!
!-----------------------------------------------------------------------
!  ROMS interpolation CLASS object.
!-----------------------------------------------------------------------
!
      TYPE, PUBLIC :: roms_interp
        logical :: rectangular = .FALSE.           ! plaid grid switch
        logical :: Have_geom = .FALSE.             ! geometry switch
!
        integer :: method                          ! interpolation type
        integer :: Cgrid                           ! C-grid location
!
        real(r8):: IJspv = 0                       ! unbounded value
!
        character (len=1) :: Cname                 ! C-grid name p,r,u,v
!
!     Source data NetCDF parameters.
!
        logical :: IsRomsFile                      ! Is a ROMS file?
!
        integer :: ncid = -1                       ! standard library ID
!
        TYPE (File_desc_t) :: pioFile              ! PIO library ID
!
        character (len=:), allocatable :: ncname   ! input filename
        character (len=:), allocatable :: Avname   ! rot. angle VarName
        character (len=:), allocatable :: Mvname   ! land/mask  VarName
        character (len=:), allocatable :: Xvname   ! X-location VarName
        character (len=:), allocatable :: Yvname   ! Y-location VarName
!
!     Source global grid geometry parameters and arrays.
!
        integer :: LBI_s,  UBI_s,  LBJ_s,  UBJ_s   ! declaration bounds
        integer :: Imin_s, Imax_s, Jmin_s, Jmax_s  ! data bounds
!
        real(r8), allocatable :: Asrc(:,:)         ! rotation angle
        real(r8), allocatable :: Msrc(:,:)         ! land/sea mask
        real(r8), allocatable :: Xsrc(:,:)         ! X-locations
        real(r8), allocatable :: Ysrc(:,:)         ! Y-locations
!
!     Destination tiled grid geometry parameters and arrays.
        integer :: LBI_d,  UBI_d,  LBJ_d,  UBJ_d   ! declaration bounds
        integer :: Imin_d, Imax_d, Jmin_d, Jmax_d  ! data bounds
!
        real(r8), allocatable :: Adst(:,:)         ! rotation angle
        real(r8), pointer :: Mdst(:,:)             ! land/sea mask
        real(r8), pointer :: Xdst(:,:)             ! X-locations
        real(r8), pointer :: Ydst(:,:)             ! Y-locations
        real(r8), allocatable :: Idst(:,:)         ! X-fractional coords
        real(r8), allocatable :: Jdst(:,:)         ! Y-fractional coords
!
        CONTAINS
!
        PROCEDURE :: create         => roms_interp_create
        PROCEDURE :: delete         => roms_interp_delete
        PROCEDURE :: frac_coords    => roms_interp_frac_coords
        PROCEDURE :: read_coords    => roms_interp_read_coords
        PROCEDURE :: regrid2d       => roms_interp_regrid2d
        PROCEDURE :: regrid_bylevel => roms_interp_regrid_bylevel
!
        PROCEDURE :: linear2d       => roms_interp_linear2d
!
      END TYPE roms_interp
!
!-----------------------------------------------------------------------
!  ROMS-to-ROMS interpolation CLASS object: regridding.
!-----------------------------------------------------------------------
!
      TYPE, PUBLIC :: roms2roms
        TYPE (roms_interp), allocatable :: Cvars(:)  ! C-type variables
!
        CONTAINS
!
        PROCEDURE :: create       => roms2roms_create
        PROCEDURE :: delete       => roms2roms_delete
      END TYPE roms2roms
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      PUBLIC  :: roms_interp_regrid2d
      PUBLIC  :: roms_interp_regrid_bylevel
      PUBLIC  :: roms2roms_create
      PRIVATE
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      CONTAINS
!
!<><><><><><><><><><><><><><><><><><><><><><><> CLASS ROMS_INTERP <><><>
!
!  It allocates and initializes field interpolation object.
!
      SUBROUTINE roms_interp_create (self, ng, tile, kernel, Ctype,     &
     &                               ncname)
!
!  Imported variable declarations.
!
      CLASS (roms_interp), intent(inout) :: self     ! grid geometries
      integer,             intent(in   ) :: ng       ! nested grid
      integer,             intent(in   ) :: tile     ! grid partition
      integer,             intent(in   ) :: kernel   ! kernel ID
      character (len=*),   intent(in   ) :: Ctype    ! C-type ID
      character (len=*),   intent(in   ) :: ncname   ! source file
!
!  Local variable declarations.
!
      integer                            :: LBi, UBi, LBj, UBj
      integer                            :: i, j, ghost
!
      character (len=*), parameter :: MyFile =                          &
     &    "ROMS/Utility/roms_interp.F"//", roms_interp_create"
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
!  Create intepolation object for the requested C-type variable.
!-----------------------------------------------------------------------
!
      self%Mdst => NULL()
      self%Xdst => NULL()
      self%Ydst => NULL()
!
!  Allocate and initialize destination interpolation parameters and
!  arrays. This arrays are available from the calling ROMS kernel.
!
      ghost = 0
      SELECT CASE (TRIM(Ctype))
        CASE ('r', 'R')
          self%Cname = 'r'
          self%Cgrid = 2
          self%Imin_d = BOUNDS(ng)%Imin(2, ghost, tile)
          self%Imax_d = BOUNDS(ng)%Imax(2, ghost, tile)
          self%Jmin_d = BOUNDS(ng)%Jmin(2, ghost, tile)
          self%Jmax_d = BOUNDS(ng)%Jmax(2, ghost, tile)
        CASE ('u', 'U')
          self%Cname = 'u'
          self%Cgrid = 3
          self%Imin_d = BOUNDS(ng)%Imin(3, ghost, tile)
          self%Imax_d = BOUNDS(ng)%Imax(3, ghost, tile)
          self%Jmin_d = BOUNDS(ng)%Jmin(3, ghost, tile)
          self%Jmax_d = BOUNDS(ng)%Jmax(3, ghost, tile)
        CASE ('v', 'V')
          self%Cname = 'v'
          self%Cgrid = 4
          self%Imin_d = BOUNDS(ng)%Imin(4, ghost, tile)
          self%Imax_d = BOUNDS(ng)%Imax(4, ghost, tile)
          self%Jmin_d = BOUNDS(ng)%Jmin(4, ghost, tile)
          self%Jmax_d = BOUNDS(ng)%Jmax(4, ghost, tile)
      END SELECT
!
      LBi = BOUNDS(ng)%LBi(tile)
      UBi = BOUNDS(ng)%UBi(tile)
      LBj = BOUNDS(ng)%LBj(tile)
      UBj = BOUNDS(ng)%UBj(tile)
!
      IF (.not.allocated(self%Adst))                                    &
        allocate ( self%Adst(LBi:UBi,LBj:UBj) )
      self%Adst = 0.0_r8
      IF (.not.allocated(self%Idst))                                    &
        allocate ( self%Idst(LBi:UBi,LBj:UBj) )
      self%Idst = 0.0_r8
      IF (.not.allocated(self%Jdst))                                    &
        allocate ( self%Jdst(LBi:UBi,LBj:UBj) )
      self%Jdst = 0.0_r8
!
      self%LBI_d = LBi
      self%UBI_d = UBi
      self%LBJ_d = LBj
      self%UBJ_d = UBj
!
      SELECT CASE (TRIM(Ctype))
        CASE ('r', 'R')
          self%Adst =  GRID(ng)%angler
          self%Mdst => GRID(ng)%rmask
          IF (spherical) THEN
            self%Xdst => GRID(ng)%lonr
            self%Ydst => GRID(ng)%latr
          ELSE
            self%Xdst => GRID(ng)%xr
            self%Ydst => GRID(ng)%yr
          END IF
        CASE ('u', 'U')
          DO j=Jstr-1,Jend+1
            DO i=IstrU-1, Iend+1
              self%Adst(i,j) = 0.5_r8*(GRID(ng)%angler(i-1,j)+          &
     &                                 GRID(ng)%angler(i  ,j))
            END DO
          END DO
          self%Mdst => GRID(ng)%umask
          IF (spherical) THEN
            self%Xdst => GRID(ng)%lonu
            self%Ydst => GRID(ng)%latu
          ELSE
            self%Xdst => GRID(ng)%xu
            self%Ydst => GRID(ng)%yu
          END IF
        CASE ('v', 'V')
          DO j=JstrV-1,Jend+1
            DO i=Istr-1, Iend+1
              self%Adst(i,j) = 0.5_r8*(GRID(ng)%angler(i,j-1)+          &
     &                                 GRID(ng)%angler(i,j  ))
            END DO
          END DO
          self%Mdst => GRID(ng)%vmask
          IF (spherical) THEN
            self%Xdst => GRID(ng)%lonv
            self%Ydst => GRID(ng)%latv
          ELSE
            self%Xdst => GRID(ng)%xv
            self%Ydst => GRID(ng)%yv
          END IF
      END SELECT
!
!  Get source interpolation parameters and arrays from input NetCDF
!  file.
!
      CALL self%read_coords (ng, kernel, Ctype, ncname)
      IF (FoundError(exit_flag, NoError, 290, MyFile)) RETURN
!
!  Compute destination grid fractional coordinates to accelerate the
!  interpolation in curvilinear, non-plaid grids. It identifies the
!  grid destination cell containing the interpolant.
!
      CALL self%frac_coords (ng)
!
      RETURN
      END SUBROUTINE roms_interp_create
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  It destroys/deallocates ROMS interpolation object
!
      SUBROUTINE roms_interp_delete (self)
!
!  Imported variable declarations.
!
      CLASS (roms_interp), intent(inout) :: self
!
!-----------------------------------------------------------------------
!  Deallocate source and destination arrays in object.
!-----------------------------------------------------------------------
!
      IF (allocated(self%Asrc)) deallocate (self%Asrc)
      IF (allocated(self%Msrc)) deallocate (self%Msrc)
      IF (allocated(self%Xsrc)) deallocate (self%Xsrc)
      IF (allocated(self%Ysrc)) deallocate (self%Ysrc)
!
      IF (allocated (self%Adst)) deallocate (self%Adst)
      IF (associated(self%Mdst)) nullify    (self%Mdst)
      IF (associated(self%Xdst)) nullify    (self%Xdst)
      IF (associated(self%Ydst)) nullify    (self%Ydst)
      IF (allocated (self%Idst)) deallocate (self%Idst)
      IF (allocated (self%Jdst)) deallocate (self%Jdst)
!
!  Turn off switches.
!
      self%rectangular = .FALSE.
      self%Have_geom = .FALSE.
!
      RETURN
      END SUBROUTINE roms_interp_delete
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  It computes the factional coordinates that identify the grid cells
!  for the horizontal interpolation from source to destination ROMS
!  grids. This task is challenging in curvilinear, non-plaid grids
!  with tiled arrays.
!
      SUBROUTINE roms_interp_frac_coords (self, ng)
!
!  Imported variable declarations.
!
      CLASS (roms_interp), intent(inout) :: self     ! grid geometry
      integer,             intent(in   ) :: ng       ! nested grid
!
!  Local variable declarations.
!
      character (len=*),       parameter :: MyFile =                    &
     &  "ROMS/Utility/roms_interp.F"//", roms_interp_frac_coords"
!
!-----------------------------------------------------------------------
!  Compute fractional grid cell coordinates.
!-----------------------------------------------------------------------
!
!  Compute fractional coarse grid cell indices containing the fine grid
!  points using geogrpahical coordinates. It uses a half binary search
!  and the inside/outside polygon algorithm, which affected by roundoff.
!
      CALL hindices (ng,                                                &
     &               self%LBI_s,  self%UBI_s,  self%LBJ_s,  self%UBJ_s, &
     &               self%Imin_s, self%Imax_s, self%Jmin_s, self%Jmax_s,&
     &               self%Asrc,   self%Xsrc,   self%Ysrc,               &
     &               self%LBI_d,  self%UBI_d,  self%LBJ_d,  self%UBJ_d, &
     &               self%Imin_d, self%Imax_d, self%Jmin_d, self%Jmax_d,&
     &               self%Xdst,   self%Ydst,   self%Idst,   self%Jdst,  &
     &               self%IJspv,  self%rectangular)
!
      RETURN
      END SUBROUTINE roms_interp_frac_coords
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  It inquires about source input NetCDF, its dimensions, and horizontal
!  geometry. Then, it allocates and initializes interpolation object
!  variables. Since mesh-to-mesh regridding in tiled distributed-memory
!  applications is challenging, the source geometry arrays and data are
!  processed serially instead of in parallel. Thus, the source arrays
!  have a global scope in all parallel tiles. We could use the optimized
!  FieldRegridStore in the ESMF library in the future.
!
      SUBROUTINE roms_interp_read_coords (self, ng, kernel, Ctype,      &
     &                                    ncname)
!
!  Imported variable declarations.
!
      CLASS (roms_interp), intent(inout) :: self     ! grid geometries
      integer,             intent(in   ) :: ng       ! nested grid
      integer,             intent(in   ) :: kernel   ! kernel ID
      character (len=*),   intent(in   ) :: Ctype    ! C-type ID
      character (len=*),   intent(in   ) :: ncname   ! source file
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &   "ROMS/Utility/roms_interp.F"//", roms_interp_read_coords"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Build source grid geometry by inquiring its NetCDF file.
!-----------------------------------------------------------------------
!
!  Invoke the requested NetCDF library.
!
      SELECT CASE (inp_lib)
        CASE (io_nf90)
          CALL roms_interp_read_nf90 (self, ng, kernel, Ctype, ncname)
        CASE (io_pio)
          CALL roms_interp_read_pio (self, ng, kernel, Ctype, ncname)
      END SELECT
      IF (FoundError(exit_flag, NoError, 441, MyFile)) RETURN
!
      RETURN
      END SUBROUTINE roms_interp_read_coords
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Read in source geometry using the standard NetCDF library. Then, load
!  themn to interpolation object for the reqested C-grid location.
!
      SUBROUTINE roms_interp_read_nf90 (self, ng, kernel, Ctype, ncname)
!
!  Imported variable declarations.
!
      CLASS (roms_interp), intent(inout) :: self     ! grid geometries
      integer,             intent(in   ) :: ng       ! nested grid
      integer,             intent(in   ) :: kernel   ! kernel ID
      character (len=*),   intent(in   ) :: Ctype    ! C-type ID
      character (len=*),   intent(in   ) :: ncname   ! source file
!
!  Local variable declarations.
!
      integer                            :: LBi, UBi, LBj, UBj
      integer                            :: Lp_src, Mp_src
      integer                            :: i, idim, ivar, j
      integer,              dimension(2) :: start, total
!
      real (r8),             allocatable :: wrk(:,:)
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/roms_interp.F"//", roms_interp_read_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set source grid geometry arrays (standard NetCDF library option).
!-----------------------------------------------------------------------
!
!  If appropriate, open source NetCDF, inquire about dimesions and
!  variables.
!
      IF (self%ncid.eq.-1) THEN
        IF (FoundError(assign_string(self%ncname, ncname),              &
     &                 NoError, 484, MyFile)) RETURN
!
        CALL netcdf_open (ng, kernel, self%ncname, 0, self%ncid)
        IF (FoundError(exit_flag, NoError, 487, MyFile)) THEN
          WRITE (stdout,10) TRIM(self%ncname)
          RETURN
        END IF
!
        CALL netcdf_get_dim (ng, kernel, self%ncname,                   &
     &                       ncid = self%ncid)
        IF (FoundError(exit_flag, NoError, 494, MyFile)) RETURN
!
        IF (find_string(dim_name, n_dim, 'xi_rho', idim)) THEN
          self%IsRomsFile = .TRUE.
          Lp_src = dim_size(idim)
        ELSE
          self%IsRomsFile = .FALSE.
        END IF
!
        IF (self%IsRomsFile) THEN
          IF (find_string(dim_name, n_dim, 'eta_rho', idim)) THEN
            Mp_src = dim_size(idim)
            IF ((Lp_src.eq.IOBOUNDS(ng)%xi_rho).and.                    &
     &          (Mp_src.eq.IOBOUNDS(ng)%eta_rho)) THEN
              IF (Master) WRITE (stdout,20) Lp_src, Mp_src,             &
     &                                      IOBOUNDS(ng)%xi_rho,        &
     &                                      IOBOUNDS(ng)%eta_rho,       &
     &                                      TRIM(ncname)
              exit_flag = 2
              RETURN
            END IF
          END IF
        END IF 
!
        CALL netcdf_inq_var (ng, kernel, self%ncname,                   &
     &                       ncid = self%ncid)
        IF (FoundError(exit_flag, NoError, 520, MyFile)) RETURN
      END IF
!
!  Determine source geometry arrays size (global scope).
!
      QUERY : IF (self%IsRomsFile) THEN
        SELECT CASE (TRIM(Ctype))
          CASE ('r', 'R')
            IF (find_string(dim_name, n_dim, 'xi_rho', idim)) THEN
              self%LBI_s  = 0
              self%Imin_s = 0
              self%UBI_s  = dim_size(idim)-1
              self%Imax_s = dim_size(idim)-1
              start(1) = 1
              total(1) = dim_size(idim)
            ELSE
              IF (Master) WRITE (stdout,30) 'xi_rho', self%ncname
              exit_flag=2
              RETURN
            END IF
!
            IF (find_string(dim_name, n_dim, 'eta_rho', idim)) THEN
              self%LBJ_s  = 0
              self%Jmin_s = 0
              self%UBJ_s  = dim_size(idim)-1
              self%Jmax_s = dim_size(idim)-1
              start(2) = 1
              total(2) = dim_size(idim)
            ELSE
              IF (Master) WRITE (stdout,30) 'eta_rho', self%ncname
              exit_flag=2
              RETURN
            END IF
!
            IF (FoundError(assign_string(self%Avname, 'angle'),         &
     &                     NoError, 556, MyFile)) RETURN
            IF (FoundError(assign_string(self%Mvname, 'mask_rho'),      &
     &                     NoError, 559, MyFile)) RETURN
            IF (spherical) THEN
              IF (FoundError(assign_string(self%Xvname, 'lon_rho'),     &
     &                       NoError, 563, MyFile)) RETURN
              IF (FoundError(assign_string(self%Yvname, 'lat_rho'),     &
     &                       NoError, 565, MyFile)) RETURN
            ELSE
              IF (FoundError(assign_string(self%Xvname, 'xr'),          &
     &                       NoError, 568, MyFile)) RETURN
              IF (FoundError(assign_string(self%Yvname, 'yr'),          &
     &                       NoError, 570, MyFile)) RETURN
            END IF
          CASE ('u', 'U')
            IF (find_string(dim_name, n_dim, 'xi_u', idim)) THEN
              self%LBI_s  = 1
              self%Imin_s = 1
              self%UBI_s  = dim_size(idim)
              self%Imax_s = dim_size(idim)
              start(1) = 1
              total(1) = dim_size(idim)
            ELSE
              IF (Master) WRITE (stdout,30) 'xi_u', self%ncname
              exit_flag=2
              RETURN
            END IF
!
            IF (find_string(dim_name, n_dim, 'eta_u', idim)) THEN
              self%LBJ_s  = 0
              self%Jmin_s = 0
              self%UBJ_s  = dim_size(idim)-1
              self%Jmax_s = dim_size(idim)-1
              start(2) = 1
              total(2) = dim_size(idim)
            ELSE
              IF (Master) WRITE (stdout,30) 'eta_u', self%ncname
              exit_flag=2
              RETURN
            END IF
!
            IF (FoundError(assign_string(self%Avname, 'angle'),         &
     &                     NoError, 601, MyFile)) RETURN
            IF (FoundError(assign_string(self%Mvname, 'mask_u'),        &
     &                     NoError, 604, MyFile)) RETURN
            IF (spherical) THEN
              IF (FoundError(assign_string(self%Xvname, 'lon_u'),       &
     &                       NoError, 608, MyFile)) RETURN
              IF (FoundError(assign_string(self%Yvname, 'lat_u'),       &
     &                       NoError, 610, MyFile)) RETURN
            ELSE
              IF (FoundError(assign_string(self%Xvname, 'xu'),          &
     &                       NoError, 613, MyFile)) RETURN
              IF (FoundError(assign_string(self%Yvname, 'yu'),          &
     &                       NoError, 615, MyFile)) RETURN
            END IF
          CASE ('v', 'V')
            IF (find_string(dim_name, n_dim, 'xi_v', idim)) THEN
              self%LBI_s  = 0
              self%Imin_s = 0
              self%UBI_s  = dim_size(idim)-1
              self%Imax_s = dim_size(idim)-1
              start(1) = 1
              total(1) = dim_size(idim)
            ELSE
              IF (Master) WRITE (stdout,30) 'xi_v', self%ncname
              exit_flag=2
              RETURN
            END IF
!
            IF (find_string(dim_name, n_dim, 'eta_v', idim)) THEN
              self%LBJ_s  = 1
              self%Jmin_s = 1
              self%UBJ_s  = dim_size(idim)
              self%Jmax_s = dim_size(idim)
              start(2) = 1
              total(2) = dim_size(idim)
            ELSE
              IF (Master) WRITE (stdout,30) 'eta_v', self%ncname
              exit_flag=2
              RETURN
            END IF
!
            IF (FoundError(assign_string(self%Avname, 'angle'),         &
     &                     NoError, 646, MyFile)) RETURN
            IF (FoundError(assign_string(self%Mvname, 'mask_v'),        &
     &                     NoError, 649, MyFile)) RETURN
            IF (spherical) THEN
              IF (FoundError(assign_string(self%Xvname, 'lon_v'),       &
     &                       NoError, 653, MyFile)) RETURN
              IF (FoundError(assign_string(self%Yvname, 'lat_v'),       &
     &                       NoError, 655, MyFile)) RETURN
            ELSE
              IF (FoundError(assign_string(self%Xvname, 'xv'),          &
     &                       NoError, 658, MyFile)) RETURN
              IF (FoundError(assign_string(self%Yvname, 'yv'),          &
     &                       NoError, 660, MyFile)) RETURN
            END IF
        END SELECT
!
!  Allocate and initialize.
!
        LBi = self%LBI_s
        UBi = self%UBI_s
        LBj = self%LBJ_s
        UBj = self%UBJ_s
!
        IF (.not.allocated(self%Asrc))                                  &
          allocate ( self%Asrc(LBi:UBi, LBj:UBj) )
        self%Asrc = 0.0_r8
        IF (.not.allocated(self%Msrc))                                  &
          allocate ( self%Msrc(LBi:UBi, LBj:UBj) )
        self%Msrc = 1.0_r8                           ! no land mask
        IF (.not.allocated(self%Xsrc))                                  &
          allocate ( self%Xsrc(LBi:UBi, LBj:UBj) )
        self%Xsrc = 0.0_r8
        IF (.not.allocated(self%Ysrc))                                  &
          allocate ( self%Ysrc(LBi:UBi, LBj:UBj) )
        self%Ysrc = 0.0_r8
        IF (self%Cname.eq.'u') THEN
          allocate ( wrk(LBi:UBi+1, LBj:UBj) )
          wrk = 0.0_r8
        ELSE IF (self%Cname.eq.'v') THEN
          allocate ( wrk(LBi:UBi+1, LBj:UBj) )
          wrk = 0.0_r8
        END IF
!
!  Read in geometry arrays needed for interpolation.
!
        IF (find_string(var_name, n_var, self%Avname, ivar)) THEN
          IF ((self%Cname.eq.'u').or.(self%Cname.eq.'v')) THEN
            CALL netcdf_get_fvar (ng, kernel, self%ncname,              &
     &                            self%Avname, wrk(LBi:,LBj:),          &
     &                            ncid = self%ncid,                     &
     &                            start = start,                        &
     &                            total = total)
          ELSE
            CALL netcdf_get_fvar (ng, kernel, self%ncname,              &
     &                            self%Avname, self%Asrc,               &
     &                            ncid = self%ncid,                     &
     &                            start = start,                        &
     &                            total = total)
          END IF
          IF (FoundError(exit_flag, NoError, 711, MyFile)) RETURN
!
          IF (self%Cname.eq.'u') THEN
            DO j = self%Jmin_s,self%Jmax_s
              DO i = self%Imin_s,self%Imax_s-1
                self%Asrc(i,j) = 0.5_r8*(wrk(i,j)+wrk(i+1,j))
              END DO
            END DO
          ELSE IF (self%Cname.eq.'v') THEN
            DO j = self%Jmin_s,self%Jmax_s-1
              DO i = self%Imin_s,self%Imax_s
                self%Asrc(i,j) = 0.5_r8*(wrk(i,j)+wrk(i,j+1))
              END DO
            END DO
          END IF
          IF (allocated(wrk)) deallocate (wrk)
        END IF
!
        IF (find_string(var_name, n_var, self%Mvname, ivar)) THEN
          CALL netcdf_get_fvar (ng, kernel, self%ncname,                &
     &                          self%Mvname, self%Msrc,                 &
     &                          ncid = self%ncid,                       &
     &                          start = start,                          &
     &                          total = total)
          IF (FoundError(exit_flag, NoError, 737, MyFile)) RETURN
        END IF
!
        IF (find_string(var_name, n_var, self%Xvname, ivar)) THEN
          CALL netcdf_get_fvar (ng, kernel, self%ncname,                &
     &                          self%Xvname, self%Xsrc,                 &
     &                          ncid = self%ncid,                       &
     &                          start = start,                          &
     &                          total = total)
          IF (FoundError(exit_flag, NoError, 747, MyFile)) RETURN
        END IF
!
        IF (find_string(var_name, n_var, self%Yvname, ivar)) THEN
          CALL netcdf_get_fvar (ng, kernel, self%ncname,                &
     &                          self%Yvname, self%Ysrc,                 &
     &                          ncid = self%ncid,                       &
     &                          start = start,                          &
     &                          total = total)
          IF (FoundError(exit_flag, NoError, 756, MyFile)) RETURN
        END IF
      ELSE
        IF (Master) WRITE (stdout,40) self%ncname
        exit_flag=5
        RETURN
      END IF QUERY
!  
  10  FORMAT (/,' ROMS_INTERP::read_nf90 - unable to open file: ',      &
              /,26x,a)
  20  FORMAT (/,' ROMS_INTERP::read_nf90 - equal source and ',          &
     &        'destination grids:',2(2x,i0,'x',i0),/,26x,'file: ',a,    &
     &        /,26x,'Must have coarse and fine resolution grids!')
  30  FORMAT (/,' ROMS_INTERP::read_nf90 - find dimension: ',a,         &
     &        'in file: ',/,26x,a)
  40  FORMAT (/,' ROMS_INTERP::read_nf90 - not a ROMS file: ',a,        &
     &        /,26x,a)
!
      RETURN
      END SUBROUTINE roms_interp_read_nf90
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Read in source geometry using the PIO NetCDF library. Then, load
!  themn to interpolation object for the reqested C-grid location.
!
      SUBROUTINE roms_interp_read_pio (self, ng, kernel, Ctype, ncname)
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      CLASS (roms_interp), intent(inout) :: self     ! grid geometries
      integer,             intent(in   ) :: ng       ! nested grid
      integer,             intent(in   ) :: kernel   ! kernel ID
      character (len=*),   intent(in   ) :: Ctype    ! C-type ID
      character (len=*),   intent(in   ) :: ncname   ! source file
!
!  Local variable declarations.
!
      integer                            :: LBi, UBi, LBj, UBj
      integer                            :: Lp_src, Mp_src
      integer                            :: i, idim, ivar, j
      integer,              dimension(2) :: start, total
!
      real (r8),             allocatable :: wrk(:,:)
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/roms_interp.F"//", roms_interp_read_pio"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set source grid geometry arrays (PIO NetCDF library option).
!-----------------------------------------------------------------------
!
!  If appropriate, open source NetCDF, inquire about dimesions and
!  variables.
!
      IF (self%pioFile%fh.eq.-1) THEN
        IF (FoundError(assign_string(self%ncname, ncname),              &
     &                 NoError, 822, MyFile)) RETURN
!
        CALL pio_netcdf_open (ng, kernel, self%ncname, 0, self%pioFile)
        IF (FoundError(exit_flag, NoError, 825, MyFile)) THEN
          WRITE (stdout,10) TRIM(self%ncname)
          RETURN
        END IF
!
        CALL pio_netcdf_get_dim (ng, kernel, self%ncname,               &
     &                           pioFile = self%pioFile)
        IF (FoundError(exit_flag, NoError, 832, MyFile)) RETURN
!
        IF (find_string(dim_name, n_dim, 'xi_rho', idim)) THEN
          self%IsRomsFile = .TRUE.
          Lp_src = dim_size(idim)
        ELSE
          self%IsRomsFile = .FALSE.
        END IF
!
        IF (self%IsRomsFile) THEN
          IF (find_string(dim_name, n_dim, 'eta_rho', idim)) THEN
            Mp_src = dim_size(idim)
            IF ((Lp_src.eq.IOBOUNDS(ng)%xi_rho).and.                    &
     &          (Mp_src.eq.IOBOUNDS(ng)%eta_rho)) THEN
              IF (Master) WRITE (stdout,20) Lp_src, Mp_src,             &
     &                                      IOBOUNDS(ng)%xi_rho,        &
     &                                      IOBOUNDS(ng)%eta_rho,       &
     &                                      TRIM(ncname)
              exit_flag = 2
              RETURN
            END IF
          END IF
        END IF 
!
        CALL pio_netcdf_inq_var (ng, kernel, self%ncname,               &
     &                           pioFile = self%pioFile)
        IF (FoundError(exit_flag, NoError, 858, MyFile)) RETURN
      END IF
!
!  Determine source geometry arrays size (global scope).
!
      QUERY : IF (self%IsRomsFile) THEN
        SELECT CASE (TRIM(Ctype))
          CASE ('r', 'R')
            IF (find_string(dim_name, n_dim, 'xi_rho', idim)) THEN
              self%LBI_s  = 0
              self%Imin_s = 0
              self%UBI_s  = dim_size(idim)-1
              self%Imax_s = dim_size(idim)-1
              start(1) = 1
              total(1) = dim_size(idim)
            ELSE
              IF (Master) WRITE (stdout,30) 'xi_rho', self%ncname
              exit_flag=2
              RETURN
            END IF
!
            IF (find_string(dim_name, n_dim, 'eta_rho', idim)) THEN
              self%LBJ_s  = 0
              self%Jmin_s = 0
              self%UBJ_s  = dim_size(idim)-1
              self%Jmax_s = dim_size(idim)-1
              start(2) = 1
              total(2) = dim_size(idim)
            ELSE
              IF (Master) WRITE (stdout,30) 'eta_rho', self%ncname
              exit_flag=2
              RETURN
            END IF
!
            IF (FoundError(assign_string(self%Avname, 'angle'),         &
     &                     NoError, 893, MyFile)) RETURN
            IF (FoundError(assign_string(self%Mvname, 'mask_rho'),      &
     &                     NoError, 896, MyFile)) RETURN
            IF (spherical) THEN
              IF (FoundError(assign_string(self%Xvname, 'lon_rho'),     &
     &                       NoError, 900, MyFile)) RETURN
              IF (FoundError(assign_string(self%Yvname, 'lat_rho'),     &
     &                       NoError, 902, MyFile)) RETURN
            ELSE
              IF (FoundError(assign_string(self%Xvname, 'xr'),          &
     &                       NoError, 905, MyFile)) RETURN
              IF (FoundError(assign_string(self%Yvname, 'yr'),          &
     &                       NoError, 907, MyFile)) RETURN
            END IF
          CASE ('u', 'U')
            IF (find_string(dim_name, n_dim, 'xi_u', idim)) THEN
              self%LBI_s  = 1
              self%Imin_s = 1
              self%UBI_s  = dim_size(idim)
              self%Imax_s = dim_size(idim)
              start(1) = 1
              total(1) = dim_size(idim)
            ELSE
              IF (Master) WRITE (stdout,30) 'xi_u', self%ncname
              exit_flag=2
              RETURN
            END IF
!
            IF (find_string(dim_name, n_dim, 'eta_u', idim)) THEN
              self%LBJ_s  = 0
              self%Jmin_s = 0
              self%UBJ_s  = dim_size(idim)-1
              self%Jmax_s = dim_size(idim)-1
              start(2) = 1
              total(2) = dim_size(idim)
            ELSE
              IF (Master) WRITE (stdout,30) 'eta_u', self%ncname
              exit_flag=2
              RETURN
            END IF
!
            IF (FoundError(assign_string(self%Avname, 'angle'),         &
     &                     NoError, 938, MyFile)) RETURN
            IF (FoundError(assign_string(self%Mvname, 'mask_u'),        &
     &                     NoError, 941, MyFile)) RETURN
            IF (spherical) THEN
              IF (FoundError(assign_string(self%Xvname, 'lon_u'),       &
     &                       NoError, 945, MyFile)) RETURN
              IF (FoundError(assign_string(self%Yvname, 'lat_u'),       &
     &                       NoError, 947, MyFile)) RETURN
            ELSE
              IF (FoundError(assign_string(self%Xvname, 'xu'),          &
     &                       NoError, 950, MyFile)) RETURN
              IF (FoundError(assign_string(self%Yvname, 'yu'),          &
     &                       NoError, 952, MyFile)) RETURN
            END IF
          CASE ('v', 'V')
            IF (find_string(dim_name, n_dim, 'xi_v', idim)) THEN
              self%LBI_s  = 0
              self%Imin_s = 0
              self%UBI_s  = dim_size(idim)-1
              self%Imax_s = dim_size(idim)-1
              start(1) = 1
              total(1) = dim_size(idim)
            ELSE
              IF (Master) WRITE (stdout,30) 'xi_v', self%ncname
              exit_flag=2
              RETURN
            END IF
!
            IF (find_string(dim_name, n_dim, 'eta_v', idim)) THEN
              self%LBJ_s  = 1
              self%Jmin_s = 1
              self%UBJ_s  = dim_size(idim)
              self%Jmax_s = dim_size(idim)
              start(2) = 1
              total(2) = dim_size(idim)
            ELSE
              IF (Master) WRITE (stdout,30) 'eta_v', self%ncname
              exit_flag=2
              RETURN
            END IF
!
            IF (FoundError(assign_string(self%Avname, 'angle'),         &
     &                     NoError, 983, MyFile)) RETURN
            IF (FoundError(assign_string(self%Mvname, 'mask_v'),        &
     &                     NoError, 986, MyFile)) RETURN
            IF (spherical) THEN
              IF (FoundError(assign_string(self%Xvname, 'lon_v'),       &
     &                       NoError, 990, MyFile)) RETURN
              IF (FoundError(assign_string(self%Yvname, 'lat_v'),       &
     &                       NoError, 992, MyFile)) RETURN
            ELSE
              IF (FoundError(assign_string(self%Xvname, 'xv'),          &
     &                       NoError, 995, MyFile)) RETURN
              IF (FoundError(assign_string(self%Yvname, 'yv'),          &
     &                       NoError, 997, MyFile)) RETURN
            END IF
        END SELECT
!
!  Allocate and initialize.
!
        LBi = self%LBI_s
        UBi = self%UBI_s
        LBj = self%LBJ_s
        UBj = self%UBJ_s
!
        IF (.not.allocated(self%Asrc))                                  &
          allocate ( self%Asrc(LBi:UBi, LBj:UBj) )
        self%Asrc = 0.0_r8
        IF (.not.allocated(self%Msrc))                                  &
          allocate ( self%Msrc(LBi:UBi, LBj:UBj) )
        self%Msrc = 1.0_r8                           ! no land mask
        IF (.not.allocated(self%Xsrc))                                  &
          allocate ( self%Xsrc(LBi:UBi, LBj:UBj) )
        self%Xsrc = 0.0_r8
        IF (.not.allocated(self%Ysrc))                                  &
          allocate ( self%Ysrc(LBi:UBi, LBj:UBj) )
        self%Ysrc = 0.0_r8
        IF (self%Cname.eq.'u') THEN
          allocate ( wrk(LBi:UBi+1, LBj:UBj) )
          wrk = 0.0_r8
        ELSE IF (self%Cname.eq.'v') THEN
          allocate ( wrk(LBi:UBi, LBj:UBj+1) )
          wrk = 0.0_r8
        END IF
!
!  Read in geometry arrays needed for interpolation.
!
        IF (find_string(var_name, n_var, self%Avname, ivar)) THEN
          IF ((self%Cname.eq.'u').or.(self%Cname.eq.'v')) THEN
            CALL pio_netcdf_get_fvar (ng, kernel, self%ncname,          &
     &                                self%Avname, wrk(LBi:,LBj:),      &
     &                                pioFile = self%pioFile)
          ELSE
            CALL pio_netcdf_get_fvar (ng, kernel, self%ncname,          &
     &                                self%Avname, self%Asrc,           &
     &                                pioFile = self%pioFile,           &
     &                                start = start,                    &
     &                                total = total)
          END IF
          IF (FoundError(exit_flag, NoError, 1046, MyFile)) RETURN
!
          IF (self%Cname.eq.'u') THEN
            DO j = self%Jmin_s,self%Jmax_s
              DO i = self%Imin_s,self%Imax_s-1
                self%Asrc(i,j) = 0.5_r8*(wrk(i,j)+wrk(i+1,j))
              END DO
            END DO
          ELSE IF (self%Cname.eq.'v') THEN
            DO j = self%Jmin_s,self%Jmax_s-1
              DO i = self%Imin_s,self%Imax_s
                self%Asrc(i,j) = 0.5_r8*(wrk(i,j)+wrk(i,j+1))
              END DO
            END DO
          END IF
        END IF
!
        IF (find_string(var_name, n_var, self%Mvname, ivar)) THEN
          CALL pio_netcdf_get_fvar (ng, kernel, self%ncname,            &
     &                              self%Mvname, self%Msrc,             &
     &                              pioFile = self%pioFile,             &
     &                              start = start,                      &
     &                              total = total)
          IF (FoundError(exit_flag, NoError, 1071, MyFile)) RETURN
        END IF
!
        IF (find_string(var_name, n_var, self%Xvname, ivar)) THEN
          CALL pio_netcdf_get_fvar (ng, kernel, self%ncname,            &
     &                              self%Xvname, self%Xsrc,             &
     &                              pioFile = self%pioFile,             &
     &                              start = start,                      &
     &                              total = total)
          IF (FoundError(exit_flag, NoError, 1081, MyFile)) RETURN
        END IF
!
        IF (find_string(var_name, n_var, self%Yvname, ivar)) THEN
          CALL pio_netcdf_get_fvar (ng, kernel, self%ncname,            &
     &                              self%Yvname, self%Ysrc,             &
     &                              pioFile = self%pioFile,             &
     &                              start = start,                      &
     &                              total = total)
          IF (FoundError(exit_flag, NoError, 1090, MyFile)) RETURN
        END IF
        IF (allocated(wrk)) deallocate (wrk)
      ELSE
        IF (Master) WRITE (stdout,40) self%ncname
        exit_flag=5
        RETURN
      END IF QUERY
!
  10  FORMAT (/,' ROMS_INTERP::read_pio - unable to open file: ',       &
     &        /,25x,a)
  20  FORMAT (/,' ROMS_INTERP::read_pio - equal source and ',           &
     &        'destination grids:',2(2x,i0,'x',i0),/,25x,'file: ',a,    &
     &          /,25x,'Must have coarse and fine resolution grids!')
  30  FORMAT (/,' ROMS_INTERP::read_pio - unable to find dimension:',   &
     &        1x,a,' in file: ',/,25x,a)
  40  FORMAT (/,' ROMS_INTERP::read_pio - not a ROMS file: ',a,         &
     &        /,25x,a)
!
      RETURN
      END SUBROUTINE roms_interp_read_pio
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Read in requested source 2D variable and regrid it to the destination
!  grid.
!
      SUBROUTINE roms_interp_regrid2d (self,                            &
     &                                 ng, tile, kernel, method,        &
     &                                 ifield, Linp,                    &
     &                                 ncname, ncvname,                 &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 Imin, Imax, Jmin, Jmax,          &
     &                                 Fmask,                           &
     &                                 Fout)
!
      USE mod_pio_netcdf, ONLY : pio_netcdf_get_fvar
!
!  Imported variable declarations.
!
      CLASS (roms_interp), intent(inout) :: self     ! grid geometries
      integer,             intent(in   ) :: ng       ! nested grid
      integer,             intent(in   ) :: tile     ! domain partition
      integer,             intent(in   ) :: kernel   ! kernel ID
      integer,             intent(in   ) :: ifield   ! metadata index
      integer,             intent(in   ) :: method   ! regrid method
      integer,             intent(in   ) :: Linp     ! time record
      character (len=*),   intent(in   ) :: ncname   ! source file
      character (len=*),   intent(in   ) :: ncvname  ! variable name
      integer,             intent(in   ) :: LBi, UBi, LBj, UBj
      integer,             intent(in   ) :: Imin, Imax, Jmin, Jmax
      real (r8),           intent(in   ) :: Fmask(LBi:,LBj:)
      real (r8),           intent(inout) :: Fout(LBi:,LBj:)
!
!  Local variable declarations.
!
      TYPE (T_STATS)                     :: Stats
!
      integer                            :: LBx, UBx, LBy, UBy
      integer                            :: ivar
      integer,              dimension(3) :: start, total
!
      real (r8),             allocatable :: Finp(:,:)
!
      character (len=*),       parameter :: MyFile =                    &
     &  "ROMS/Utility/roms_interp.F"//", roms_interp_regrid2d"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Interpolate 2D field from source to destination grid.
!-----------------------------------------------------------------------
!
!  Read in requested 2D source variable.
!
      IF (find_string(var_name, n_var, TRIM(ncvname), ivar)) THEN
        LBx = self%Imin_s
        UBx = self%Imax_s
        LBy = self%Jmin_s
        UBy = self%Jmax_s
!
        start(1) = 1
        start(2) = 1
        start(3) = Linp
        total(1) = UBx-LBx+1
        total(2) = UBy-LBy+1
        total(3) = 1
!
        allocate ( Finp(LBx:UBx, LBy:UBy) )
!
        SELECT CASE (inp_lib)
          CASE (io_nf90)
            CALL netcdf_get_fvar (ng, kernel, ncname,                   &
     &                            TRIM(ncvname), Finp(LBx:,LBy:),       &
     &                            ncid = self%ncid,                     &
     &                            start = start,                        &
     &                            total = total)
            CASE (io_pio)
              CALL pio_netcdf_get_fvar (ng, kernel, ncname,             &
     &                                  TRIM(ncvname), Finp(LBx:,LBy:), &
     &                                  pioFile = self%pioFile,         &
     &                                  start = start,                  &
     &                                  total = total)
        END SELECT
        IF (FoundError(exit_flag, NoError, 1209, MyFile)) RETURN
      ELSE
        IF (Master) WRITE (stdout,10) TRIM(ncvname), TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Regrid requested 2D variable.
!
      Fout = 0.0_r8
!
      IF (method.eq.linear) THEN
        CALL self%linear2d (LBx, UBx, LBy, UBy, Finp(LBx:,LBy:),        &
     &                      LBi, UBi, LBj, UBj, Fout(LBi:,LBj:))
      ELSE
      END IF
      IF (allocated(Finp)) deallocate (Finp)
!
!-----------------------------------------------------------------------
!  Compute interpolated field statistics.
!-----------------------------------------------------------------------
!
!  Initialize statistics object.
!
      Stats % checksum=0_i8b
      Stats % count=0
      Stats % min=spval
      Stats % max=-spval
      Stats % avg=0.0_r8
      Stats % rms=0.0_r8
!
      CALL stats_2dfld (ng, tile, kernel, self%Cgrid, Stats, 0,         &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  Fout,                                           &
     &                  Fmask = Fmask,                                  &
     &                  debug = .FALSE.)
      IF (Master) THEN
        WRITE (stdout,20) TRIM(Vname(2,ifield)), TRIM(ncvname),         &
                          Stats%min, Stats%max, Stats%checksum
      END IF
!
  10  FORMAT (/,' ROMS_INTERP::regrid2d - unable to find variable: ',   &
     &        a,2x,'in file: ',/,24x,a)
  20  FORMAT (19x,'- ',a,': ',a,/,22x,'(Min = ',1p,e15.8,',  Max = ',   &
     &        1p,e15.8,')  CheckSum = ',i0)
!
      RETURN
      END SUBROUTINE roms_interp_regrid2d
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Read in requested source 3D variable and regrid it to the destination
!  grid. The regridding in done level-by-level.
!
      SUBROUTINE roms_interp_regrid_bylevel (self,                      &
     &                                       ng, tile, kernel, method,  &
     &                                       ifield, Linp,              &
     &                                       ncname, ncvname,           &
     &                                       LBi, UBi, LBj, UBj,        &
     &                                       LBk, UBk,                  &
     &                                       Imin, Imax, Jmin, Jmax,    &
     &                                       Fmask,                     &
     &                                       Fout)
!
      USE mod_pio_netcdf, ONLY : pio_netcdf_get_fvar
!
!  Imported variable declarations.
!
      CLASS (roms_interp), intent(inout) :: self     ! grid geometries
      integer,             intent(in   ) :: ng       ! nested grid
      integer,             intent(in   ) :: tile     ! domain partition
      integer,             intent(in   ) :: kernel   ! kernel ID
      integer,             intent(in   ) :: method   ! regrid method
      integer,             intent(in   ) :: ifield   ! metadata index
      integer,             intent(in   ) :: Linp     ! time record
      character (len=*),   intent(in   ) :: ncname   ! source file
      character (len=*),   intent(in   ) :: ncvname  ! variable name
      integer,             intent(in   ) :: LBi, UBi, LBj, UBj
      integer,             intent(in   ) :: LBk, UBk
      integer,             intent(in   ) :: Imin, Imax, Jmin, Jmax
      real (r8),           intent(in   ) :: Fmask(LBi:,LBj:)
      real (r8),           intent(inout) :: Fout(LBi:,LBj:,LBk:)
!
!  Local variable declarations.
!
      TYPE (T_STATS)                     :: Stats
!
      integer                            :: LBx, UBx, LBy, UBy
      integer                            :: ivar, k
      integer,              dimension(4) :: start, total
      real (r8),             allocatable :: Finp(:,:,:)
!
      character (len=*),       parameter :: MyFile =                    &
     &  "ROMS/Utility/roms_interp.F"//", roms_interp_regrid_bylevel"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Interpolate 3D field from source to destination grid.
!-----------------------------------------------------------------------
!
!  Read in requested 2D source variable.
!
      IF (find_string(var_name, n_var, TRIM(ncvname), ivar)) THEN
        LBx = self%Imin_s
        UBx = self%Imax_s
        LBy = self%Jmin_s
        UBy = self%Jmax_s
!
        start(1) = 1
        start(2) = 1
        start(3) = 1
        start(4) = Linp
        total(1) = UBx-LBx+1
        total(2) = UBy-LBy+1
        total(3) = UBk-LBk+1
        total(4) = 1
!
        allocate ( Finp(LBx:UBx, LBy:UBy, LBk:UBk) )
!
        SELECT CASE (inp_lib)
          CASE (io_nf90)
            CALL netcdf_get_fvar (ng, kernel, ncname,                   &
     &                            TRIM(ncvname), Finp(LBx:,LBy:,:),     &
     &                            ncid = self%ncid,                     &
     &                            start = start,                        &
     &                            total = total)
          CASE (io_pio)
            CALL pio_netcdf_get_fvar (ng, kernel, ncname,               &
     &                                TRIM(ncvname), Finp(LBx:,LBy:,:), &
     &                                pioFile = self%pioFile,           &
     &                                start = start,                    &
     &                                total = total)
        END SELECT
        IF (FoundError(exit_flag, NoError, 1355, MyFile)) RETURN
      ELSE
        IF (Master) WRITE (stdout,10) TRIM(ncvname), TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Regrid requested 3D variable, level-by-level.
!
      Fout = 0.0_r8
!
      IF (method.eq.linear) THEN
        DO k = LBk,UBk
          CALL self%linear2d (LBx, UBx, LBy, UBy, Finp(LBx:,LBy:,k),    &
     &                        LBi, UBi, LBj, UBj, Fout(LBi:,LBj:,k))
        END DO
      ELSE
      END IF
      IF (allocated(Finp)) deallocate (Finp)
!
!-----------------------------------------------------------------------
!  Compute interpolated field statistics.
!-----------------------------------------------------------------------
!
!  Initialize statistics object.
!
      Stats % checksum=0_i8b
      Stats % count=0
      Stats % min=spval
      Stats % max=-spval
      Stats % avg=0.0_r8
      Stats % rms=0.0_r8
!
      CALL stats_3dfld (ng, tile, kernel, self%Cgrid, Stats, 0,         &
     &                  LBi, UBi, LBj, UBj, LBk, UBk,                   &
     &                  Fout,                                           &
     &                  Fmask = Fmask,                                  &
     &                  debug = .FALSE.)
      IF (Master) THEN
        WRITE (stdout,20) TRIM(Vname(2,ifield)), TRIM(ncvname),         &
                          Stats%min, Stats%max, Stats%checksum
      END IF
!
  10  FORMAT (/,' ROMS_INTERP::regrid_bylevel - unable to find ',       &
              'variable: ',a,2x,'in file: ',/,30x,a)
  20  FORMAT (19x,'- ',a,': ',a,/,22x,'(Min = ',1p,e15.8,',  Max = ',   &
     &        1p,e15.8,')  CheckSum = ',i0)
!
      RETURN
      END SUBROUTINE roms_interp_regrid_bylevel
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Linear interpolation for any 2D field that includes land/sea mask
!  in source (SRC) and destination (DST) grids. The interplation is done
!  in fractional mesh cell indices in SRC and DST grids.
!
!  The input source field (Finp) has a globlal scope where as the output
!  field (Fout) is a tiled array of the calling ROMS kernel. Thus, the
!  interpolation is in parallel.
!
      SUBROUTINE roms_interp_linear2d (self,                            &
     &                                 LBx, UBx, LBy, UBy, Finp,        &
     &                                 LBi, UBi, LBj, UBj, Fout)
!
!  Imported variable declarations.
!
      CLASS (roms_interp), intent(in ) :: self            ! geometries
      integer,             intent(in ) :: LBx, UBx        ! SRC x-bounds
      integer,             intent(in ) :: LBy, UBy        ! SRC y-bounds
      integer,             intent(in ) :: LBi, UBi        ! DST i-bounds
      integer,             intent(in ) :: LBj, UBj        ! DST j-bounds
      real (r8),           intent(in ) :: Finp(LBx:,LBy:) ! SRC 2D field
      real (r8),           intent(out) :: Fout(LBi:,LBj:) ! DST 2D field
!
!  Local variable declarations.
!
      integer                          :: i, ic, i1, i2, j, j1, j2
      real (r8)                        :: Hmat(4)
      real (r8)                        :: p1, p2, q1, q2, wsum
!
!-----------------------------------------------------------------------
!  2D Linear interpolation..
!-----------------------------------------------------------------------
!
!  Notice that is possible for the source and destination have
!  coinciding domain boundaries because of decimation. Thus, we
!  to add limiters at the Eastern (i2=i1+1) and Northern (j2=j1+1)
!  boundaries out-range array indices.
!
      DO j = self%Jmin_d, self%Jmax_d
        DO i = self%Imin_d, self%Imax_d
          i1 = INT(self%Idst(i,j))
          i2 = MIN(i1 + 1, self%Imax_s)            ! Eastern limiter
          j1 = INT(self%Jdst(i,j))
          j2 = MIN(j1 + 1, self%Jmax_s)            ! Northern limiter
          IF (self%Mdst(i,j).gt.0.5_r8) THEN       ! DST seawater cell
            IF (((LBx.le.i1).and.(i1.le.UBx)).and.                      &
     &          ((LBy.le.j1).and.(j1.le.UBy))) THEN
              p2 = REAL(i2-i1, r8) * (self%Idst(i,j) - REAL(i1, r8))
              q2 = REAL(j2-j1, r8) * (self%Jdst(i,j) - REAL(j1, r8))
              p1 = 1.0_r8 - p2
              q1 = 1.0_r8 - q2
!
              Hmat(1) = p1 * q1 * self%Msrc(i1,j1)
              Hmat(2) = p2 * q1 * self%Msrc(i2,j1)
              Hmat(3) = p2 * q2 * self%Msrc(i2,j2)
              Hmat(3) = p1 * q2 * self%Msrc(i1,j2)
!
              wsum = 0.0_r8                        ! Interpolation
              DO ic=1,4                            ! weights
                wsum = wsum + Hmat(ic)             ! check
              END DO
              IF (wsum.gt.0.0_r8) THEN             ! Scale interpolation
                wsum = 1.0_r8 / wsum               ! weights to account
                DO ic=1,4                          ! land masked cells
                  Hmat(ic) = Hmat(ic) * wsum       ! in SRC field
                END DO
              END IF
!
              Fout(i,j) = Hmat(1) * Finp(i1,j1)+                        &
     &                    Hmat(2) * Finp(i2,j1)+                        &
     &                    Hmat(3) * Finp(i2,j2)+                        &
     &                    Hmat(4) * Finp(i1,j2)
            END IF
          ELSE
            Fout(i,j) = 0.0_r8                     ! DST land cell
          END IF
        END DO
      END DO
!
      RETURN
      END SUBROUTINE roms_interp_linear2d
!
!<><><><><><><><><><><><><><><><><><><><><><><><> CLASS ROMS2ROMS <><><>
!
!  It creates ROMS-to-ROMS interpolation object.
!
      SUBROUTINE roms2roms_create (self, ng, tile, kernel, Ctype,       &
     &                             ncname)
!
!  Imported variable declarations.
!
      CLASS (roms2roms), intent(inout) :: self      ! grid geometries
      integer,           intent(in   ) :: ng        ! nested grid
      integer,           intent(in   ) :: tile      ! grid partition
      integer,           intent(in   ) :: kernel    ! kernel ID
      character (len=*), intent(in   ) :: Ctype(:)  ! C-type ID
      character (len=*), intent(in   ) :: ncname    ! source file
!
!  Local variable declarations.
!
      integer                          :: i, nvars
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/roms_interp.F"//", roms2roms_create"
!
!-----------------------------------------------------------------------
!  Create ROMS-to-ROMS interpolation objects according to C-type grid.
!-----------------------------------------------------------------------
!
!  Allocate ROMS-to-ROMS object.
!
      nvars = SIZE(Ctype)
      allocate (self%Cvars(nvars))
!
!  Allocate source and destination interpolation arrays for each C-type
!  variable requested.
!
      DO i=1,nvars
        IF (i.gt.1) THEN           ! avoid inquiring about source file
          SELECT CASE (inp_lib)    
            CASE (io_nf90)
              self%Cvars(i)%ncid = self%Cvars(1)%ncid
            CASE (io_pio)
              self%Cvars(i)%pioFile = self%Cvars(1)%pioFile
          END SELECT
          self%Cvars(i)%IsRomsFile = self%Cvars(1)%IsRomsFile
          IF (FoundError(assign_string(self%Cvars(i)%ncname, ncname),   &
     &                 NoError, 1538, MyFile)) RETURN
        ELSE
          self%Cvars(1)%ncid = -1        ! force open file on first pass
          self%Cvars(1)%pioFile%fh = -1
        END IF
!
        CALL self%Cvars(i)%create (ng, tile, kernel, Ctype(i), ncname)
      END DO
!
      RETURN
      END SUBROUTINE roms2roms_create
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  It deallocates ROMS-to-ROMS interpolation object.
!
      SUBROUTINE roms2roms_delete (self)
!
!  Imported variable declarations.
!
      CLASS (roms2roms), intent(inout) :: self      ! grid geometries
!
!  Local variable declarations.
!
      integer                          :: i
!
!-----------------------------------------------------------------------
!  Deallocate source and destination arrays in ROMS-to-ROMS object.
!-----------------------------------------------------------------------
!
      DO i=1,SIZE(self%Cvars)
         CALL self%Cvars(i)%delete ()
      END DO
!
!  Destroy ROMS-to-ROMS object.
!
      deallocate (self%Cvars)
!
      RETURN
      END SUBROUTINE roms2roms_delete
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      END MODULE roms_interp_mod
