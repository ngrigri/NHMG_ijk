module mg_relax

  use mg_cst
  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids
  use mg_mpi_exchange

  implicit none

contains

  !----------------------------------------
  subroutine relax(lev,nsweeps)
    integer(kind=ip), intent(in):: lev
    integer(kind=ip), intent(in):: nsweeps

    real(kind=rp),dimension(:,:,:), pointer:: p
    real(kind=rp),dimension(:,:,:), pointer:: b
    real(kind=rp),dimension(:,:,:,:), pointer:: cA

    integer(kind=ip) :: nx, ny, nz, nd

    p  => grid(lev)%p
    b  => grid(lev)%b
    cA => grid(lev)%cA

    nx = grid(lev)%nx
    ny = grid(lev)%ny
    nz = grid(lev)%nz
    nd = size(cA(:,:,:,:),dim=1)

    if (grid(lev)%nz == 1) then

       call relax_2D_5(lev,p,b,cA,nsweeps,nx,ny)

    else

       select case(trim(relax_method))

       case('Gauss-Seidel','GS')
          call relax_3D_8_GS(lev,p,b,cA,nsweeps,nx,ny,nz)

       case('Red-Black','RB')
          call relax_3D_8_RB(lev,p,b,cA,nsweeps,nx,ny,nz)

       case('Four-Color','FC')
          call relax_3D_8_FC(lev,p,b,cA,nsweeps,nx,ny,nz)

       end select

    end if

  end subroutine relax

  !----------------------------------------
  subroutine relax_2D_5(lev,p,b,cA,nsweeps,nx,ny)

    integer(kind=ip)                         , intent(in)   :: lev
    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout):: p
    real(kind=rp),dimension(:,:,:)  , pointer, intent(in)   :: b
    real(kind=rp),dimension(:,:,:,:), pointer, intent(in)   :: cA
    integer(kind=ip)                        , intent(in)   :: nsweeps
    integer(kind=ip)                        , intent(in)   :: nx, ny

    integer(kind=ip)           :: i,j,k, it,rb
    integer(kind=ip)            :: ib,ie,jb,je,rbb,rbe,rbi
    real(kind=rp) :: z,gamma,g1,g2

    gamma = one
    g1 = gamma
    g2 = one - gamma

    k=1

    do it = 1,nsweeps
       if (mod(it,1) == 0) then
          ib = 1 
          ie = nx
          jb = 1
          je = ny
       else
          ib = 0
          ie = nx+1
          jb = 0
          je = ny+1
       endif

       if (red_black) then
          rbb = 1
          rbe = 2
          rbi = 2
       else
          rbb = 0
          rbe = 0
          rbi = 1
       endif

       do rb = rbb,rbe
          do i = ib,ie
             do j = jb+mod(i+rb,rbi),je,rbi

                z =    b(i,j,k)                                                    &
                     - cA(2,i,j,k) * p(i  ,j-1,k) - cA(2,i  ,j+1,k) * p(i  ,j+1,k) &
                     - cA(3,i,j,k) * p(i-1,j  ,k) - cA(3,i+1,j  ,k) * p(i+1,j  ,k) &
                     - cA(4,i,j,k) * p(i-1,j-1,k) - cA(4,i+1,j+1,k) * p(i+1,j+1,k) &
                     - cA(5,i,j,k) * p(i-1,j+1,k) - cA(5,i+1,j-1,k) * p(i+1,j-1,k)

                p(i,j,k) = z / cA(1,i,j,k)

             enddo
          enddo
       enddo

       call fill_halo(lev,p)

    enddo

  end subroutine relax_2D_5

 !----------------------------------------
  subroutine relax_3D_8_GS(lev,p,b,cA,nsweeps,nx,ny,nz)

    integer(kind=ip)                        , intent(in)   :: lev
    integer(kind=ip)                        , intent(in)   :: nsweeps
    integer(kind=ip)                        , intent(in)   :: nx, ny, nz

    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout):: p
    real(kind=rp),dimension(:,:,:)  , pointer, intent(in)   :: b
    real(kind=rp),dimension(:,:,:,:), pointer, intent(in)   :: cA

    integer(kind=ip)            :: i,j,it

    call tic(lev,'relax_3D_8_GS')

    ! add a loop on smoothing
    do it = 1,nsweeps

       do j = 1, ny
          do i = 1, nx

             call relax_3D_8_heart(p,b,cA,i,j,nz)

          enddo ! i
       enddo    ! j

       call fill_halo(lev,p)

    enddo  !it

    call toc(lev,'relax_3D_8_GS')

  end subroutine relax_3D_8_GS

 !----------------------------------------
  subroutine relax_3D_8_RB(lev,p,b,cA,nsweeps,nx,ny,nz)

    integer(kind=ip)                        , intent(in)   :: lev
    integer(kind=ip)                        , intent(in)   :: nsweeps
    integer(kind=ip)                        , intent(in)   :: nx, ny, nz

    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout):: p
    real(kind=rp),dimension(:,:,:)  , pointer, intent(in)   :: b
    real(kind=rp),dimension(:,:,:,:), pointer, intent(in)   :: cA

    !- Local -!
    integer(kind=ip) :: i,j,it
    integer(kind=ip) :: rb

    real(kind=rp),dimension(nx,ny,nz) :: gam
    real(kind=rp),dimension(nx,ny,nz)    :: bet

    integer(kind=ip) :: k
    real(kind=rp)    :: rhs
    real(kind=rp), dimension(nx,ny,nz) :: p2

    call tic(lev,'relax_3D_8_RB')

    ! Coefficients are stored in order of diagonals
    ! cA(1,:,:,:)      -> p(i  ,j  ,k  )
    ! cA(2,:,:,:)      -> p(i  ,j  ,k-1)
    ! cA(3,:,:,:)      -> p(i  ,j-1,k+1)
    ! cA(4,:,:,:)      -> p(i  ,j-1,k  )
    ! cA(5,:,:,:)      -> p(i  ,j-1,k-1)
    ! cA(6,:,:,:)      -> p(i-1,j  ,k+1)
    ! cA(7,:,:,:)      -> p(i-1,j  ,k  )
    ! cA(8,:,:,:)      -> p(i-1,j  ,k-1)

    ! add a loop on smoothing
    do it = 1,nsweeps

       do rb = 1, 2 ! Red black loop

          k=1 !lower level
          do j = 1, ny
             do i = 1+mod(j+rb,2),nx,2

                rhs  = b(i,j,k)                                                          & ! ijk
                     - cA(3,i  ,j  ,k  ) * p(i,j-1,k+1)                                  & ! ijk
                     - cA(4,i  ,j  ,k  ) * p(i,j-1,k  ) - cA(4,i  ,j+1,k) * p(i  ,j+1,k) & ! ijk
                     - cA(5,i  ,j+1,k+1) * p(i,j+1,k+1)                                  & ! ijk
                     - cA(6,i  ,j  ,k  ) * p(i-1,j,k+1)                                  & ! ijk
                     - cA(7,i  ,j  ,k  ) * p(i-1  ,j,k) - cA(7,i+1,j  ,k) * p(i+1,j  ,k) & ! ijk
                     - cA(8,i+1,j  ,k+1) * p(i+1,j,k+1)                                    ! ijk

                if (cmatrix == 'real') then
                   !- Exception for the redefinition of the coef for the bottom level
                   rhs   = rhs   &
                        - cA(5,i,j,k) * p(i-1,j+1,k) - cA(5,i+1,j-1,k) * p(i+1,j-1,k) & ! ijk
                        - cA(8,i,j,k) * p(i-1,j-1,k) - cA(8,i+1,j+1,k) * p(i+1,j+1,k)    ! ijk
                endif

                bet(i,j,k) = one / cA(1,i,j,k)! k=1
                p2(i,j,k)  = rhs * bet(i,j,k) ! k=1

             enddo ! i
          enddo  ! j

          do k=2,nz-1
             do j = 1, ny
                do i = 1+mod(j+rb,2),nx,2

                   rhs  = b(i,j,k)                                                          & ! ijk
                        - cA(3,i,j,k) * p(i  ,j-1,k+1) - cA(3,i  ,j+1,k-1) * p(i  ,j+1,k-1) & ! ijk
                        - cA(4,i,j,k) * p(i  ,j-1,k  ) - cA(4,i  ,j+1,k  ) * p(i  ,j+1,k  ) & ! ijk
                        - cA(5,i,j,k) * p(i  ,j-1,k-1) - cA(5,i  ,j+1,k+1) * p(i  ,j+1,k+1) & ! ijk
                        - cA(6,i,j,k) * p(i-1,j  ,k+1) - cA(6,i+1,j  ,k-1) * p(i+1,j  ,k-1) & ! ijk
                        - cA(7,i,j,k) * p(i-1,j  ,k  ) - cA(7,i+1,j  ,k  ) * p(i+1,j  ,k  ) & ! ijk
                        - cA(8,i,j,k) * p(i-1,j,  k-1) - cA(8,i+1,j  ,k+1) * p(i+1,j  ,k+1)   ! ijk

                   gam(i,j,k) = cA(2,i,j,k) * bet(i,j,k-1)
                   bet(i,j,k) = one / ( cA(1,i,j,k) - cA(2,i,j,k) * gam(i,j,k) )
                   p2(i,j,k)  = ( rhs - cA(2,i,j,k) * p2(i,j,k-1) ) * bet(i,j,k)

                enddo ! i
             enddo ! j
          enddo ! k=2,nz

          k=nz !upper level
          do j = 1, ny
             do i = 1+mod(j+rb,2),nx,2

                rhs  = b(i,j,k)                                                            & ! ijk
                     - cA(3,i  ,j+1,k-1) * p(i  ,j+1,k-1)                                  & ! ijk
                     - cA(4,i  ,j  ,k  ) * p(i  ,j-1,k  ) - cA(4,i  ,j+1,k)  *p(i  ,j+1,k) & ! ijk
                     - cA(5,i  ,j  ,k  ) * p(i  ,j-1,k-1)                                  & ! ijk
                     - cA(6,i+1,j  ,k-1) * p(i+1,j  ,k-1)                                  & ! ijk
                     - cA(7,i  ,j  ,k  ) * p(i-1,j  ,k  ) - cA(7,i+1,j  ,k) * p(i+1,j  ,k) & ! ijk
                     - cA(8,i  ,j  ,k  ) * p(i-1,j  ,k-1)                                    ! ijk

!                d(i,j,k)   = cA(1,i,j,k) ! ijk

                gam(i,j,k) = cA(2,i,j,k) * bet(i,j,k-1)
                bet(i,j,k) = one / ( cA(1,i,j,k) - cA(2,i,j,k) * gam(i,j,k) )
                p2(i,j,k)  = ( rhs  - cA(2,i,j,k) * p2(i,j,k-1) ) * bet(i,j,k)

             enddo ! i
          enddo ! j

          do k=nz-1,1,-1
             do j = 1, ny
                do i = 1+mod(j+rb,2),nx,2
                   p2(i,j,k) = p2(i,j,k) - gam(i,j,k+1) * p2(i,j,k+1)
                enddo ! i
             enddo ! j
          enddo ! k = nz-1,1,-1

          do k=1,nz
             do j = 1, ny
                do i = 1+mod(j+rb,2),nx,2
                   p(i,j,k) = p2(i,j,k)
                enddo ! i
             enddo ! j
          enddo ! k=1,nz

          call fill_halo(lev,p)

       enddo       ! red-black

    enddo  !it

    call toc(lev,'relax_3D_8_RB')

  end subroutine relax_3D_8_RB

 !----------------------------------------
  subroutine relax_3D_8_FC(lev,p,b,cA,nsweeps,nx,ny,nz)

    integer(kind=ip)                        , intent(in)   :: lev
    integer(kind=ip)                        , intent(in)   :: nsweeps
    integer(kind=ip)                        , intent(in)   :: nx, ny, nz

    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout):: p
    real(kind=rp),dimension(:,:,:)  , pointer, intent(in)   :: b
    real(kind=rp),dimension(:,:,:,:), pointer, intent(in)   :: cA

    integer(kind=ip)            :: i,j,it
    integer(kind=ip)            :: fc1,fc2

    call tic(lev,'relax_3D_8_FC')

    ! add a loop on smoothing
    do it = 1,nsweeps

       do fc1 = 1, 2 ! 
          do fc2 = 1, 2 ! 
             do j = 1 + mod(fc1-1,2), ny, 2
                do i = 1 + mod(fc2-1,2), nx, 2

                   call relax_3D_8_heart(p,b,cA,i,j,nz)

                enddo ! j
             enddo    ! i

             call fill_halo(lev,p)

          enddo  ! fc2
       enddo  ! fc1

    enddo  !it

    call toc(lev,'relax_3D_8_FC')

  end subroutine relax_3D_8_FC

  !---------------------------------------------------------------------
  subroutine relax_3D_8_heart(p,b,cA,i,j,nz)

    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout):: p
    real(kind=rp),dimension(:,:,:)  , pointer, intent(in)   :: b
    real(kind=rp),dimension(:,:,:,:), pointer, intent(in)   :: cA

    integer(kind=ip)                          , intent(in)  :: i
    integer(kind=ip)                          , intent(in)  :: j
    integer(kind=ip)                          , intent(in)  :: nz

    !- Local -!
    integer(kind=ip) :: k
    real(kind=rp), dimension(nz) :: rhs, d, ud

    real(kind=rp),dimension(nz):: gam
    real(kind=rp)              :: bet


    ! Coefficients are stored in order of diagonals
    ! cA(1,:,:,:)      -> p(i  ,j  ,k  )
    ! cA(2,:,:,:)      -> p(i  ,j  ,k-1)
    ! cA(3,:,:,:)      -> p(i  ,j-1,k+1)
    ! cA(4,:,:,:)      -> p(i  ,j-1,k  )
    ! cA(5,:,:,:)      -> p(i  ,j-1,k-1)
    ! cA(6,:,:,:)      -> p(i-1,j  ,k+1)
    ! cA(7,:,:,:)      -> p(i-1,j  ,k  )
    ! cA(8,:,:,:)      -> p(i-1,j  ,k-1)

    k=1 !lower level
    rhs(k) = b(i,j,k)                                                        & ! ijk
         - cA(3,i  ,j  ,k  ) * p(i,j-1,k+1)                                  & ! ijk
         - cA(4,i  ,j  ,k  ) * p(i,j-1,k  ) - cA(4,i  ,j+1,k) * p(i  ,j+1,k) & ! ijk
         - cA(5,i  ,j+1,k+1) * p(i,j+1,k+1)                                  & ! ijk
         - cA(6,i  ,j  ,k  ) * p(i-1,j,k+1)                                  & ! ijk
         - cA(7,i  ,j  ,k  ) * p(i-1  ,j,k) - cA(7,i+1,j  ,k) * p(i+1,j  ,k) & ! ijk
         - cA(8,i+1,j  ,k+1) * p(i+1,j,k+1)                                    ! ijk

    if (cmatrix == 'real') then
       !- Exception for the redefinition of the coef for the bottom level
       rhs(k) = rhs(k) &
            - cA(5,i,j,k) * p(i-1,j+1,k) - cA(5,i+1,j-1,k) * p(i+1,j-1,k) & ! ijk
            - cA(8,i,j,k) * p(i-1,j-1,k) - cA(8,i+1,j+1,k) * p(i+1,j+1,k)    ! ijk
    endif

    d(k)   = cA(1,i,j,k)   ! ijk
    ud(k)  = cA(2,i,j,k+1) ! ijk 

    do k = 2,nz-1 !interior levels
       rhs(k) = b(i,j,k)                                                        & ! ijk
            - cA(3,i,j,k) * p(i  ,j-1,k+1) - cA(3,i  ,j+1,k-1) * p(i  ,j+1,k-1) & ! ijk
            - cA(4,i,j,k) * p(i  ,j-1,k  ) - cA(4,i  ,j+1,k  ) * p(i  ,j+1,k  ) & ! ijk
            - cA(5,i,j,k) * p(i  ,j-1,k-1) - cA(5,i  ,j+1,k+1) * p(i  ,j+1,k+1) & ! ijk
            - cA(6,i,j,k) * p(i-1,j  ,k+1) - cA(6,i+1,j  ,k-1) * p(i+1,j  ,k-1) & ! ijk
            - cA(7,i,j,k) * p(i-1,j  ,k  ) - cA(7,i+1,j  ,k  ) * p(i+1,j  ,k  ) & ! ijk
            - cA(8,i,j,k) * p(i-1,j,  k-1) - cA(8,i+1,j  ,k+1) * p(i+1,j  ,k+1)   ! ijk
       d(k)   = cA(1,i,j,k  ) ! ijk
       ud(k)  = cA(2,i,j,k+1) ! ijk
    enddo

    k=nz !upper level
    rhs(k) = b(i,j,k)                                                          & ! ijk
         - cA(3,i  ,j+1,k-1) * p(i  ,j+1,k-1)                                  & ! ijk
         - cA(4,i  ,j  ,k  ) * p(i  ,j-1,k  ) - cA(4,i  ,j+1,k)  *p(i  ,j+1,k) & ! ijk
         - cA(5,i  ,j  ,k  ) * p(i  ,j-1,k-1)                                  & ! ijk
         - cA(6,i+1,j  ,k-1) * p(i+1,j  ,k-1)                                  & ! ijk
         - cA(7,i  ,j  ,k  ) * p(i-1,j  ,k  ) - cA(7,i+1,j  ,k) * p(i+1,j  ,k) & ! ijk
         - cA(8,i  ,j  ,k  ) * p(i-1,j  ,k-1)                                    ! ijk
    d(k)   = cA(1,i,j,k) ! ijk

    bet      = one / d(1)
    p(i,j,1) = rhs(1) * bet

    do k=2,nz
       gam(k)   = ud(k-1) * bet
       bet      = one / ( d(k) - ud(k-1) * gam(k) )
       p(i,j,k) = ( rhs(k) - ud(k-1) * p(i,j,k-1) ) * bet
    enddo

    do k=nz-1,1,-1
       p(i,j,k) = p(i,j,k) - gam(k+1) * p(i,j,k+1)
    enddo

  end subroutine relax_3D_8_heart

  !----------------------------------------
  subroutine compute_residual(lev,res)
    integer(kind=ip), intent(in) :: lev
    real(kind=rp)   , intent(out):: res

    real(kind=rp),dimension(:,:,:)  , pointer:: p
    real(kind=rp),dimension(:,:,:)  , pointer:: b
    real(kind=rp),dimension(:,:,:)  , pointer:: r
    real(kind=rp),dimension(:,:,:,:), pointer:: cA

    integer(kind=ip) :: nx, ny, nz, nd
    real(kind=rp) ::resloc

    p  => grid(lev)%p
    b  => grid(lev)%b
    r  => grid(lev)%r
    cA => grid(lev)%cA

    nx = grid(lev)%nx
    ny = grid(lev)%ny
    nz = grid(lev)%nz
    nd = size(cA(:,:,:,:),dim=1)

    if (grid(lev)%nz == 1) then

       call compute_residual_2D_5(res,p,b,r,cA,nx,ny)

    else

       call compute_residual_3D_8(res,p,b,r,cA,nx,ny,nz)

    end if

    call fill_halo(lev,r)

    if (lev >-1) then
       resloc=res
       call global_sum(lev,resloc,res)
       res = sqrt(res)
    else
       res = -999._rp
    endif

  end subroutine compute_residual

  !----------------------------------------
  subroutine compute_residual_2D_5(res,p,b,r,cA,nx,ny)

    real(kind=rp)                            , intent(out)  :: res
    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout):: p
    real(kind=rp),dimension(:,:,:)  , pointer, intent(in)   :: b
    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout)   :: r
    real(kind=rp),dimension(:,:,:,:), pointer, intent(in)   :: cA
    integer(kind=ip)                        , intent(in)   :: nx, ny

    integer(kind=ip) :: i,j,k
    real(kind=rp)  :: z

    res = zero

    k=1

    do j = 1,ny
       do i = 1,nx
          z = b(i,j,k) - cA(1,i,j,k)*p(i,j,k)                               &
               - cA(2,i,j,k) * p(i  ,j-1,k) - cA(2,i  ,j+1,k) * p(i  ,j+1,k)&
               - cA(3,i,j,k) * p(i-1,j  ,k) - cA(3,i+1,j  ,k) * p(i+1,j  ,k)&
               - cA(4,i,j,k) * p(i-1,j-1,k) - cA(4,i+1,j+1,k) * p(i+1,j+1,k)&
               - cA(5,i,j,k) * p(i-1,j+1,k) - cA(5,i+1,j-1,k) * p(i+1,j-1,k)

          r(i,j,k) = z
          !          res = max(res,abs(r(i,j,k)))
          res = res+z*z

       enddo
    enddo

  end subroutine compute_residual_2D_5

  !----------------------------------------
  subroutine compute_residual_3D_8(res,p,b,r,cA,nx,ny,nz)

    real(kind=rp)                            , intent(out)   :: res
    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout) :: p
    real(kind=rp),dimension(:,:,:)  , pointer, intent(in)    :: b
    real(kind=rp),dimension(:,:,:)  , pointer, intent(inout) :: r
    real(kind=rp),dimension(:,:,:,:), pointer, intent(in)    :: cA
    integer(kind=ip)                        , intent(in)     :: nx, ny, nz

    ! Coefficients are stored in order of diagonals
    ! cA(1,:,:,:)      -> p(i  ,j  ,k  )
    ! cA(2,:,:,:)      -> p(i  ,j  ,k-1)
    ! cA(3,:,:,:)      -> p(i  ,j-1,k+1)
    ! cA(4,:,:,:)      -> p(i  ,j-1,k  )
    ! cA(5,:,:,:)      -> p(i  ,j-1,k-1)
    ! cA(6,:,:,:)      -> p(i-1,j  ,k+1)
    ! cA(7,:,:,:)      -> p(k-1,j  ,i  )
    ! cA(8,:,:,:)      -> p(i-1,j  ,k-1)

    integer(kind=ip)           :: i,j,k

    res = zero

    k=1 !lower level
    do j = 1,ny    ! ijk
       do i = 1,nx ! ijk

          r(i,j,k) = b(i,j,k)                                                        & ! ijk
               - cA(1,i  ,j  ,k  ) * p(i  ,j  ,k  )                                  & ! ijk
               - cA(2,i  ,j  ,k+1) * p(i  ,j  ,k+1)                                  & ! ijk
               - cA(3,i  ,j  ,k  ) * p(i  ,j-1,k+1)                                  & ! ijk
               - cA(4,i  ,j  ,k  ) * p(i  ,j-1,k  ) - cA(4,i  ,j+1,k) * p(i  ,j+1,k) & ! ijk
               - cA(5,i  ,j+1,k+1) * p(i  ,j+1,k+1)                                  & ! ijk
               - cA(6,i  ,j  ,k  ) * p(i-1,j  ,k+1)                                  & ! ijk
               - cA(7,i  ,j  ,k  ) * p(i-1,j  ,k  ) - cA(7,i+1,j  ,k) * p(i+1,j  ,k) & ! ijk
               - cA(8,i+1,j  ,k+1) * p(i+1,j  ,k+1)                                    ! ijk

          if (cmatrix == 'real') then
             !- Exception for the redefinition of the coef for the bottom level
             r(i,j,k) = r(i,j,k) &                                                ! ijk
                  - cA(5,i,j,k) * p(i-1,j+1,k) - cA(5,i+1,j-1,k) * p(i+1,j-1,k) & ! ijk
                  - cA(8,i,j,k) * p(i-1,j-1,k) - cA(8,i+1,j+1,k) * p(i+1,j+1,k)   ! ijk
          endif

          res = res+r(i,j,k)*r(i,j,k)  ! ijk
       enddo
    enddo

    do k = 2,nz-1 !interior levels  ! ijk
       do j = 1,ny    ! ijk
          do i = 1,nx ! ijk

             r(i,j,k) = b(i,j,k)                                                    & ! ijk
                  - cA(1,i,j,k) * p(i,j  ,k  )                                      & ! ijk
                  - cA(2,i,j,k) * p(i,j  ,k-1) - cA(2,i  ,j  ,k+1) * p(i  ,j  ,k+1) & ! ijk
                  - cA(3,i,j,k) * p(i,j-1,k+1) - cA(3,i  ,j+1,k-1) * p(i  ,j+1,k-1) & ! ijk
                  - cA(4,i,j,k) * p(i,j-1,k  ) - cA(4,i  ,j+1,k  ) * p(i  ,j+1,k  ) & ! ijk
                  - cA(5,i,j,k) * p(i,j-1,k-1) - cA(5,i  ,j+1,k+1) * p(i  ,j+1,k+1) & ! ijk
                  - cA(6,i,j,k) * p(i-1,j,k+1) - cA(6,i+1,j  ,k-1) * p(i+1,j  ,k-1) & ! ijk
                  - cA(7,i,j,k) * p(i-1,j,k  ) - cA(7,i+1,j  ,k  ) * p(i+1,j  ,k  ) & ! ijk
                  - cA(8,i,j,k) * p(i-1,j,k-1) - cA(8,i+1,j  ,k+1) * p(i+1,j  ,k+1)   ! ijk

             res = res+r(i,j,k)*r(i,j,k) ! ijk
          enddo

       enddo
    enddo

    k=nz !upper level
    do j = 1,ny    ! ijk
       do i = 1,nx ! ijk
          r(i,j,k) = b(i,j,k)                                                        & ! ijk
               - cA(1,i  ,j  ,k  ) * p(i  ,j  ,k  )                                  & ! ijk
               - cA(2,i  ,j  ,k  ) * p(i  ,j  ,k-1)                                  & ! ijk
               - cA(3,i  ,j+1,k-1) * p(i  ,j+1,k-1)                                  & ! ijk
               - cA(4,i  ,j  ,k  ) * p(i  ,j-1,k  ) - cA(4,i  ,j+1,k) * p(i  ,j+1,k) & ! ijk
               - cA(5,i  ,j  ,k  ) * p(i  ,j-1,k-1)                                  & ! ijk
               - cA(6,i+1,j  ,k-1) * p(i+1,j  ,k-1)                                  & ! ijk
               - cA(7,i  ,j  ,k  ) * p(i-1,j  ,k  ) - cA(7,i+1,j  ,k) * p(i+1,j  ,k) & ! ijk
               - cA(8,i  ,j  ,k  ) * p(i-1,j  ,k-1)                                    ! ijk

          res = res+r(i,j,k)*r(i,j,k)  ! ijk
       enddo
    enddo

  end subroutine compute_residual_3D_8

end module mg_relax
