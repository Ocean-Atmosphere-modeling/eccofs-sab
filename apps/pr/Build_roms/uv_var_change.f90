      MODULE uv_var_change_mod
!
!git $Id$
!=======================================================================
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                            Hernan G. Arango   !
!=======================================================================
!                                                                      !
!  These routines are used for ocean current variable changes from     !
!  C-grid to A-grid and vice versa. It is done for output purposes     !
!  and data assimilation where the state vector is located at the      !
!  cell-center (Arakawa A-grid).                                       !
!                                                                      !
!  * If transforming vector components from C-grid to A-grid, rotate   !
!    to geographical Eastward and Northward directions.                !
!                                                                      !
!  Ur(i,j,k) = 0.5 * [u(i,j,k,ninp) + u(i+1,j,k,ninp)]   i=Istr:Iend   !
!  Vr(i,j,k) = 0.5 * [v(i,j,k,ninp) + v(i,j+1,k,ninp)]   j=Jstr:Jend   !
!                                                                      !
!  Apply lateral boundary conditions (gradient) via 'bc_r3d_tile'      !
!                                                                      !
!  ua(i,j,k) = Ur(i,j,k) * CosAngler(i,j) - Vr(i,j,k) * SinAngler(i,j) !
!  va(i,j,k) = Vr(i,j,k) * CosAngler(i,j) + Ur(i,j,k) * SinAngler(i,j) !
!                                                                      !
!  PUBLIC:  uv_C2A_grid, ad_uv_C2A_grid, tl_uv_C2A_grid                !
!  PRIVATE: uv_C2A_grid_tile, ad_uv_C2A_grid_tile, tl_uv_C2A_grid_tile !
!                                                                      !
!  * If transforming vector components from A-grid to C-grid, rotate   !
!    to computational XI and ETA directions.                           !
!                                                                      !
!  Ur(i,j,k) = ua(i,j,k) * CosAngler(i,j) + va(i,j,k) * SinAngler(i,j) !
!  Vr(i,j,k) = va(i,j,k) * CosAngler(i,j) - ua(i,j,k) * SinAngler(i,j) !
!                                                                      !
!  u(i,j,k,nout) = 0.5 * [Ur(i-1,j,k) + Ur(i,j,k)]      i=Istr:IendR   !
!  v(i,j,k<nout) = 0.5 * [Vr(i,j-1,k) + Vr(i,j,k)]      j=Jstr:JendR   !
!                                                                      !
!  PUBLIC:  uv_A2C_grid, ad_uv_A2C_grid, tl_uv_A2C_grid                !
!  PRIVATE: uv_A2C_grid_tile, ad_uv_A2C_grid_tile, tl_uv_A2C_grid_tile !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_scalars
!
      USE mp_exchange_mod, ONLY : mp_exchange3d
!
      implicit none
!
      PRIVATE
      PUBLIC  :: uv_A2C_grid
      PUBLIC  :: uv_C2A_grid
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE uv_C2A_grid (ng, tile, model, ninp)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, ninp
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/uv_var_change.F"
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
      CALL wclock_on (ng, model, 34, 89, MyFile)
      CALL uv_C2A_grid_tile (ng, tile, model, ninp,                     &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       GRID(ng) % rmask_full,                     &
     &                       GRID(ng) % CosAngler,                      &
     &                       GRID(ng) % SinAngler,                      &
     &                       OCEAN(ng) % u,                             &
     &                       OCEAN(ng) % v,                             &
     &                       OCEAN(ng) % ua,                            &
     &                       OCEAN(ng) % va)
      CALL wclock_off (ng, model, 34, 104, MyFile)
!
      RETURN
      END SUBROUTINE uv_C2A_grid
!
!***********************************************************************
      SUBROUTINE uv_C2A_grid_tile (ng, tile, model, ninp,               &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             rmask,                               &
     &                             CosAngler, SinAngler,                &
     &                             u, v, ua, va)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, ninp
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
      real(r8), intent(in) :: CosAngler(LBi:,LBj:)
      real(r8), intent(in) :: SinAngler(LBi:,LBj:)
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ua(LBi:,LBj:,:)
      real(r8), intent(inout) :: va(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      integer :: i, j, k
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Urho, Vrho
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
!  Tranform vector components from C-grid to A-grid.
!-----------------------------------------------------------------------
!
      K_LOOP : DO k=1,N(ng)
!
!  Compute A-grid (cell center) components. Apply gradient condition.
!
        DO j=JstrR,JendR
          DO i=Istr,Iend
            Urho(i,j)=0.5_r8*(u(i,j,k,ninp)+u(i+1,j,k,ninp))
            IF (.not.EWperiodic(ng)) THEN
              IF (DOMAIN(ng)%Western_Edge(tile)) THEN
                Urho(Istr-1,j)=Urho(Istr,j)
              END IF
              IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
                Urho(Iend+1,j)=Urho(Iend,j)
              END IF
            END IF
          END DO
        END DO
!
        DO j=Jstr,Jend
          DO i=IstrR,IendR
            Vrho(i,j)=0.5_r8*(v(i,j,k,ninp)+v(i,j+1,k,ninp))
            IF (.not.NSperiodic(ng)) THEN
              IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
                Vrho(i,Jstr-1) = Vrho(i,Jstr)
              END IF
              IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
                Vrho(i,Jend+1) = Vrho(i,Jend)
              END IF
            END IF
          END DO
        END DO
!
!  Rotate from computational to gegraphical Eastward and Northward
!  directions.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            ua(i,j,k)=Urho(i,j)*CosAngler(i,j)-                         &
     &                Vrho(i,j)*SinAngler(i,j)
            va(i,j,k)=Vrho(i,j)*CosAngler(i,j)+                         &
     &                Urho(i,j)*SinAngler(i,j)
            ua(i,j,k)=ua(i,j,k)*rmask(i,j)
            va(i,j,k)=va(i,j,k)*rmask(i,j)
          END DO
        END DO
      END DO K_LOOP
!
      CALL mp_exchange3d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    ua, va)
!
      RETURN
      END SUBROUTINE uv_C2A_grid_tile    
!
!***********************************************************************
      SUBROUTINE uv_A2C_grid (ng, tile, model, nout)
!***********************************************************************
!
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, nout
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/uv_var_change.F"
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
      CALL wclock_on (ng, model, 34, 242, MyFile)
      CALL uv_A2C_grid_tile (ng, tile, model, nout,                     &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       GRID(ng) % umask_full,                     &
     &                       GRID(ng) % vmask_full,                     &
     &                       GRID(ng) % CosAngler,                      &
     &                       GRID(ng) % SinAngler,                      &
     &                       OCEAN(ng) % ua,                            &
     &                       OCEAN(ng) % va,                            &
     &                       OCEAN(ng) % u,                             &
     &                       OCEAN(ng) % v)
      CALL wclock_off (ng, model, 34, 258, MyFile)
!
      RETURN
      END SUBROUTINE uv_A2C_grid
!
!***********************************************************************
      SUBROUTINE uv_A2C_grid_tile (ng, tile, model, nout,               &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             umask, vmask,                        &
     &                             CosAngler, SinAngler,                &
     &                             ua, va, u, v)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, nout
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
      real(r8), intent(in) :: CosAngler(LBi:,LBj:)
      real(r8), intent(in) :: SinAngler(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
      real(r8), intent(in) :: ua(LBi:,LBj:,:)
      real(r8), intent(in) :: va(LBi:,LBj:,:)
      real(r8), intent(inout) :: u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: v(LBi:,LBj:,:,:)
!
!  Local variable declarations.
!
      integer :: i, j, k
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Urho, Vrho
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
!  Transform vector components from A-grid to C-grid.
!-----------------------------------------------------------------------
!
      K_LOOP : DO k=1,N(ng)
!
!  Rotate vector components to computational (XI,ETA) directions.
!
        DO j=Jstr-1,JendR
          DO i=Istr-1,IendR
            Urho(i,j)=ua(i,j,k)*CosAngler(i,j)+                         &
     &                va(i,j,k)*SinAngler(i,j)
            Vrho(i,j)=va(i,j,k)*CosAngler(i,j)-                         &
     &                ua(i,j,k)*SinAngler(i,j)
          END DO
        END DO
!
!  Compute staggered C-grid components.
!
        DO j=JstrR,JendR
          DO i=Istr,IendR
            u(i,j,k,nout)=0.5_r8*(Urho(i-1,j)+Urho(i,j))
            u(i,j,k,nout)=u(i,j,k,nout)*umask(i,j)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            v(i,j,k,nout)=0.5_r8*(Vrho(i,j-1)+Vrho(i,j))
            v(i,j,k,nout)=v(i,j,k,nout)*vmask(i,j)
          END DO
        END DO
      END DO K_LOOP
!
!  Exchange boundary data.
!
      CALL mp_exchange3d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    u(:,:,:,nout), v(:,:,:,nout))
!
      RETURN
      END SUBROUTINE uv_A2C_grid_tile    
      END MODULE uv_var_change_mod
