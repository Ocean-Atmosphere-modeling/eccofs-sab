      MODULE mod_arrays
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module is used to allocate, initialize, and deallocate ROMS    !
!  state arrays for each nested and/or multiple connected grids.       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars
!
      USE mod_average,         ONLY : allocate_average,                 &
     &                                deallocate_average,               &
     &                                initialize_average
      USE mod_boundary,        ONLY : allocate_boundary,                &
     &                                deallocate_boundary,              &
     &                                initialize_boundary
      USE mod_clima,           ONLY : allocate_clima,                   &
     &                                deallocate_clima,                 &
     &                                initialize_clima
      USE mod_coupling,        ONLY : allocate_coupling,                &
     &                                deallocate_coupling,              &
     &                                initialize_coupling
      USE mod_forces,          ONLY : allocate_forces,                  &
     &                                deallocate_forces,                &
     &                                initialize_forces
      USE mod_grid,            ONLY : allocate_grid,                    &
     &                                deallocate_grid,                  &
     &                                initialize_grid
      USE mod_iounits,         ONLY : deallocate_iounits
      USE mod_mixing,          ONLY : allocate_mixing,                  &
     &                                deallocate_mixing,                &
     &                                initialize_mixing
      USE mod_ocean,           ONLY : allocate_ocean,                   &
     &                                deallocate_ocean,                 &
     &                                initialize_ocean
      USE mod_sources,         ONLY : allocate_sources,                 &
     &                                deallocate_sources
      USE set_pio_mod,         ONLY : set_iodecomp
!
      implicit none
!
      PUBLIC :: ROMS_allocate_arrays
      PUBLIC :: ROMS_deallocate_arrays
      PUBLIC :: ROMS_initialize_arrays
!
      logical :: LallocateClima = .FALSE.
!
      CONTAINS
!
      SUBROUTINE ROMS_allocate_arrays (allocate_vars)
!
!=======================================================================
!                                                                      !
!  This routine allocates ROMS state variables.                        !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations
!
      logical, intent(in) :: allocate_vars
!
!  Local variable declarations.
!
      integer :: ng, thread, tile
      integer :: IminS, ImaxS, JminS, JmaxS
      integer :: LBi, UBi, LBj, UBj, LBij, UBij
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Modules/mod_arrays.F"//", ROMS_allocate_arrays"
!
!-----------------------------------------------------------------------
!  Turn on allocation time wall clock.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO thread=MyRank,MyRank
          CALL wclock_on (ng, iNLM, 1, 156, MyFile)
        END DO
!$OMP BARRIER
      END DO
!
!-----------------------------------------------------------------------
!  Allocate model type-derived structures and its variables.
!-----------------------------------------------------------------------
!
      IF (allocate_vars) then
        tile=MyRank
        DO ng=1,Ngrids
!$OMP MASTER
          LBi=BOUNDS(ng)%LBi(tile)
          UBi=BOUNDS(ng)%UBi(tile)
          LBj=BOUNDS(ng)%LBj(tile)
          UBj=BOUNDS(ng)%UBj(tile)
          LBij=BOUNDS(ng)%LBij
          UBij=BOUNDS(ng)%UBij
          CALL allocate_average (ng, LBi, UBi, LBj, UBj)
          CALL allocate_boundary (ng)
          IF (LallocateClima.or.Lclimatology(ng)) THEN
            CALL allocate_clima (ng, LBi, UBi, LBj, UBj)
          END IF
          CALL allocate_coupling (ng, LBi, UBi, LBj, UBj)
          CALL allocate_forces (ng, LBi, UBi, LBj, UBj)
          CALL allocate_grid (ng, ExtractFlag(ng),                      &
     &                        LBi, UBi, LBj, UBj, LBij, UBij)
          CALL allocate_mixing (ng, LBi, UBi, LBj, UBj)
          CALL allocate_ocean (ng, LBi, UBi, LBj, UBj)
          IF (LuvSrc(ng).or.LwSrc(ng).or.ANY(LtracerSrc(:,ng))) THEN
            CALL allocate_sources (ng)
          END IF
!$OMP END MASTER
!$OMP BARRIER
        END DO
!
!  Allocate and initialize PIO single and double precision data
!  decomposition for mapping between computational and I/O processes.
!  It needs to done after all parameters have been allocated and
!  initialize. In particular, we need tidal constituents and state
!  propagator  parameters.
!
        CALL set_iodecomp
        LallocatedMemory=.TRUE.
      END IF
!
!-----------------------------------------------------------------------
!  Turn off allocation time wall clock.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO thread=MyRank,MyRank
          CALL wclock_off (ng, iNLM, 1, 268, MyFile)
        END DO
!$OMP BARRIER
      END DO
!
      RETURN
      END SUBROUTINE ROMS_allocate_arrays
!
      SUBROUTINE ROMS_deallocate_arrays
!
!=======================================================================
!                                                                      !
!  This routine deallocates ROMS state objects, arrays and vectors.    !
!                                                                      !
!=======================================================================
!
      USE mod_param, ONLY : Ngrids
!
!  Local variable declarations.
!
      integer :: ng
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Modules/mod_arrays.F"//", ROMS_deallocate_arrays"
!
!-----------------------------------------------------------------------
!  Deallocate all structures.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
!$OMP MASTER
        CALL deallocate_average (ng)
        CALL deallocate_boundary (ng)
        IF (LallocateClima.or.Lclimatology(ng)) THEN
          CALL deallocate_clima (ng)
        END IF
        CALL deallocate_coupling (ng)
        CALL deallocate_forces (ng)
        CALL deallocate_grid (ng)
        CALL deallocate_mixing (ng)
        CALL deallocate_ocean (ng)
        IF (LuvSrc(ng).or.LwSrc(ng).or.ANY(LtracerSrc(:,ng))) THEN
          CALL deallocate_sources (ng)
        END IF
!$OMP END MASTER
!$OMP BARRIER
      END DO
!
!  Deallocate I/O derived-type structures.
!
      CALL deallocate_iounits
!
!  Deallocate main configuration dimensions and associated parameters.
!
      CALL deallocate_param
!
      RETURN
      END SUBROUTINE ROMS_deallocate_arrays
!
      SUBROUTINE ROMS_initialize_arrays
!
!=======================================================================
!                                                                      !
!  This routine initialize ROMS state variables. In shared-memory it   !
!  important for the first-touch policy in memory.                     !
!                                                                      !
!=======================================================================
!
!  Local variable declarations.
!
      integer :: ng, thread, tile
!
      integer, parameter :: model = 0
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Modules/mod_arrays.F"//", ROMS_initialize_arrays"
!
!-----------------------------------------------------------------------
!  Turn on allocation time wall clock.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO thread=MyRank,MyRank
          CALL wclock_on (ng, iNLM, 1, 405, MyFile)
        END DO
!$OMP BARRIER
      END DO
!
!-----------------------------------------------------------------------
!  Intialize variables within structures for each grid.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL initialize_average (ng, tile)
          CALL initialize_boundary (ng, tile, model)
          IF (LallocateClima.or.Lclimatology(ng)) THEN
            CALL initialize_clima (ng, tile)
          END IF
          CALL initialize_coupling (ng, tile, model)
          CALL initialize_forces (ng, tile, model)
          CALL initialize_grid (ng, tile, model)
          CALL initialize_mixing (ng, tile, model)
          CALL initialize_ocean (ng, tile, model)
        END DO
!$OMP BARRIER
      END DO
!
!-----------------------------------------------------------------------
!  Turn off allocation time wall clock.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO thread=MyRank,MyRank
          CALL wclock_off (ng, iNLM, 1, 488, MyFile)
        END DO
!$OMP BARRIER
      END DO
!
      RETURN
      END SUBROUTINE ROMS_initialize_arrays
!
      END MODULE mod_arrays
