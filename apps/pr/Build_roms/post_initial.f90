      MODULE post_initial_mod
!
!git $Id$
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  On the first timestep, it computes the initial depths and level     !
!  thicknesses from the initial free-surface field. Additionally, it   !
!  initializes the nonlinear state variables for all time levels and   !
!  applies lateral boundary conditions.                                !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
!
      USE ini_fields_mod, ONLY : ini_fields, ini_zeta
      USE set_depth_mod,  ONLY : set_depth
!
      implicit none
!
      PUBLIC  :: post_initial
      PRIVATE
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE post_initial (ng, model)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      integer :: tile
!
!-----------------------------------------------------------------------
!  Initialize free-surface and compute initial level thicknesses and
!  depths.
!-----------------------------------------------------------------------
!
      DO tile=first_tile(ng),last_tile(ng),+1
        CALL ini_zeta (ng, tile, model)
        CALL set_depth (ng, tile, model)
      END DO
!$OMP BARRIER
!
!-----------------------------------------------------------------------
!  Initialize other state variables.
!-----------------------------------------------------------------------
!
      DO tile=last_tile(ng),first_tile(ng),-1
        CALL ini_fields (ng, tile, model)
      END DO
!$OMP BARRIER
!
      RETURN
      END SUBROUTINE post_initial
      END MODULE post_initial_mod
