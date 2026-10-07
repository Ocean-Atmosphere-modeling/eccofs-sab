      MODULE omega_mod
!
!git $Id$
!=======================================================================
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                            Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine computes S-coordinate vertical velocity (m^3/s),       !
!                                                                      !
!                  W=[Hz/(m*n)]*omega,                                 !
!                                                                      !
!  diagnostically at horizontal RHO-points and vertical W-points.      !
!                                                                      !
!  Added implicit vertical advection following Shchepetkin (2015).     !
!  The vertical velocity is split into explicit (W) and implicit (Wi)  !
!  parts that adjust automatically to the local flow conditions based  !
!  on the Courant number for stability, allowing larger time steps.    !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Shchepetkin, A.F., 2015: An adaptive, Courant-number-dependent    !
!      implicit scheme for vertical advection in oceanic modeling,     !
!      Ocean Modelling, 91, 38-69, doi: 10.1016/j.ocemod.2015.03.006   !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: omega, scale_omega
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE omega (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Nonlinear/omega.F"
!
      integer :: IminS, ImaxS, JminS, JmaxS
      integer :: LBi, UBi, LBj, UBj, LBij, UBij
!
!  Set horizontal starting and ending indices for automatic private
!  storage arrays.
!
      IminS=BOUNDS(ng)%Istr(tile)-3
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-3
      JmaxS=BOUNDS(ng)%Jend(tile)+3
!
!  Determine array lower and upper bounds in the I- and J-directions.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
!  Set array lower and upper bounds for MIN(I,J) directions and
!  MAX(I,J) directions.
!
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
!
      CALL wclock_on (ng, model, 13, 64, MyFile)
      CALL omega_tile (ng, tile, model,                                 &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 GRID(ng) % Huon,                                 &
     &                 GRID(ng) % Hvom,                                 &
     &                 GRID(ng) % pm,                                   &
     &                 GRID(ng) % pn,                                   &
     &                 GRID(ng) % z_w,                                  &
     &                 OCEAN(ng) % Wi,                                  &
     &                 OCEAN(ng) % W)
      CALL wclock_off (ng, model, 13, 89, MyFile)
!
      RETURN
      END SUBROUTINE omega
!
!***********************************************************************
      SUBROUTINE omega_tile (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       Huon, Hvom,                                &
     &                       pm, pn,                                    &
     &                       z_w,                                       &
     &                       Wi,                                        &
     &                       W)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_sources
!
      USE bc_3d_mod, ONLY : bc_w3d_tile
      USE mp_exchange_mod, ONLY : mp_exchange3d
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
      real(r8), intent(in) :: Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:,:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(out) :: Wi(LBi:,LBj:,0:)
      real(r8), intent(out) :: W(LBi:,LBj:,0:)
!
!  Local variable declarations.
!
      integer :: i, ii, is, j, jj, k
      real(r8), dimension(IminS:ImaxS) :: wrk
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: Cu_adv
      real(r8) :: cff
      real(r8) :: cw, c2d, dh, cw_max, cw_max2, cw_min
!
      real(r8), parameter :: amax = 0.75_r8    ! Maximum Courant number
      real(r8), parameter :: amin = 0.60_r8    ! Minimum Courant number
      real(r8), parameter :: cmnx_ratio = amin/amax
      real(r8), parameter :: cutoff = 2.0_r8-amin/amax
      real(r8), parameter :: r4cmx = 1.0_r8/(4.0_r8-4.0_r8*amin/amax)
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
!------------------------------------------------------------------------
!  Vertically integrate horizontal mass flux divergence:
!
!  Starting with zero vertical velocity at the bottom, integrate
!  from the bottom (k=0) to the free-surface (k=N).  The w(:,:,N(ng))
!  contains the vertical velocity at the free-surface, d(zeta)/d(t).
!  Notice that barotropic mass flux divergence is not used directly.
!------------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          W(i,j,0)=0.0_r8
        END DO
        DO k=1,N(ng)
          DO i=Istr,Iend
            W(i,j,k)=W(i,j,k-1)-                                        &
     &               (Huon(i+1,j,k)-Huon(i,j,k)+                        &
     &                Hvom(i,j+1,k)-Hvom(i,j,k))
!
!  Compute the horizontal Courant number.
!
            Cu_adv(i,k)=MAX(Huon(i+1,j  ,k),0.0_r8)-                    &
     &                  MIN(Huon(i  ,j  ,k),0.0_r8)+                    &
     &                  MAX(Hvom(i  ,j+1,k),0.0_r8)-                    &
     &                  MIN(Hvom(i  ,j  ,k),0.0_r8)
          END DO
        END DO
!
!  Apply mass point sources (volume vertical influx), if any.
!
!  Overwrite W(Isrc,Jsrc,k) with the same divergence of Huon,Hvom as
!  above but add in point source Qsrc(k) and reaccumulate the vertical
!  sum to obtain the correct net Qbar given in user input - J. Levin
!  (Jupiter Intelligence Inc.) and J. Wilkin
!
!    Dsrc(is) = 2,  flow across grid cell w-face (positive or negative)
!
        IF (LwSrc(ng)) THEN
          DO is=1,Nsrc(ng)
            IF (INT(SOURCES(ng)%Dsrc(is)).eq.2) THEN
              ii=SOURCES(ng)%Isrc(is)
              jj=SOURCES(ng)%Jsrc(is)
              IF (((IstrR.le.ii).and.(ii.le.IendR)).and.                &
     &            ((JstrR.le.jj).and.(jj.le.JendR)).and.                &
     &            (j.eq.jj)) THEN
                DO k=1,N(ng)
                  W(ii,jj,k)=W(ii,jj,k-1)-                              &
     &                       (Huon(ii+1,jj,k)-Huon(ii,jj,k)+            &
     &                        Hvom(ii,jj+1,k)-Hvom(ii,jj,k))+           &
     &                       SOURCES(ng)%Qsrc(is,k)
                END DO
              END IF
            END IF
          END DO
        END IF
!
        DO i=Istr,Iend
          wrk(i)=W(i,j,N(ng))/(z_w(i,j,N(ng))-z_w(i,j,0))
          Cu_adv(i,0)=dt(ng)*pm(i,j)*pn(i,j)
        END DO
!
!  In order to insure zero vertical velocity at the free-surface,
!  subtract the vertical velocities of the moving S-coordinates
!  isosurfaces. These isosurfaces are proportional to d(zeta)/d(t).
!  The proportionally coefficients are a linear function of the
!  S-coordinate with zero value at the bottom (k=0) and unity at
!  the free-surface (k=N).
!
        DO k=N(ng)-1,1,-1
          DO i=Istr,Iend
            W(i,j,k)=W(i,j,k)-                                          &
     &               wrk(i)*(z_w(i,j,k)-z_w(i,j,0))
!
!  Determine implicit part "Wi" of vertical advection: "W"  becomes the
!  explicit part "We".
!
            Wi(i,j,k)=W(i,j,k)
            IF (Wi(i,j,k).ge.0.0_r8) THEN     ! Three different variants
              c2d=Cu_adv(i,k)                 ! for computing 2D Courant
              dh=z_w(i,j,k)-z_w(i,j,k-1)      ! number at the interface:
            ELSE                              ! (1) use value from the
              c2d=Cu_adv(i,k+1)               !     grid box upstream in
              dh=z_w(i,j,k+1)-z_w(i,j,k)      !     vertical direction;
            END IF
!
            cw_max=amax*dh-c2d*Cu_adv(i,0)    ! compare the vertical
            IF (cw_max.ge.0.0_r8) THEN        ! displacement to dz*amax.
              cw_max2=cw_max*cw_max           ! Partition W into Wi and
              cw_min=cw_max*cmnx_ratio        ! We.
              cw=ABS(Wi(i,j,k))*Cu_adv(i,0)
              IF (cw.le.cw_min) THEN
                cff=cw_max2
              ELSE IF (cw.le.cutoff*cw_max) THEN
                cff=cw_max2+r4cmx*(cw-cw_min)**2
              ELSE
                cff=cw_max*cw
              END IF
!
              W(i,j,k)=cw_max2*Wi(i,j,k)/cff
              Wi(i,j,k)=Wi(i,j,k)-W(i,j,k)
            ELSE                              ! All the displacement is
              W(i,j,k)=0.0_r8                 ! greater than amax*dz, so
            END IF                            ! keep it all into Wi.
          END DO
        END DO
        DO i=Istr,Iend
          W(i,j,N(ng))=0.0_r8
        END DO
      END DO
!
!  Set lateral boundary conditions.
!
      CALL bc_w3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 0, N(ng),                   &
     &                  W)
      CALL bc_w3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 0, N(ng),                   &
     &                  Wi)
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    W)
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Wi)
!
      RETURN
      END SUBROUTINE omega_tile
!
!***********************************************************************
      SUBROUTINE scale_omega (ng, tile, LBi, UBi, LBj, UBj, LBk, UBk,   &
     &                        pm, pn, W, Wscl)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_3d_mod, ONLY : exchange_w3d_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
!
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: W(LBi:,LBj:,LBk:)
      real(r8), intent(out) :: Wscl(LBi:,LBj:,LBk:)
!
!  Local variable declarations.
!
      integer :: i, j, k
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
!-----------------------------------------------------------------------
!  Scale omega vertical velocity to m/s.
!-----------------------------------------------------------------------
!
      DO k=LBk,UBk
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            Wscl(i,j,k)=W(i,j,k)*pm(i,j)*pn(i,j)
          END DO
        END DO
      END DO
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, LBk, UBk,           &
     &                          Wscl)
      END IF
!
      RETURN
      END SUBROUTINE scale_omega
      END MODULE omega_mod
