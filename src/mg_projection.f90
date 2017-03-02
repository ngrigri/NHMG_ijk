module mg_projection

  use mg_cst
  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids
  use mg_mpi_exchange
  use mg_gather
  use mg_netcdf_out

  implicit none

contains
  !-------------------------------------------------------------------------     
  subroutine set_rhs

    integer(kind=ip):: k,j,i
    integer(kind=ip):: nx,ny,nz

    real(kind=rp), dimension(:,:,:), pointer :: u,v,w,rhs

    !    if (myrank==0) write(*,*)'   - set rhs:'

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    u => grid(1)%u
    v => grid(1)%v
    w => grid(1)%w

    rhs => grid(1)%b

    !! What comes into nhmg_solve are area integrated u,v,w.
    do k = 1,nz       ! ijk
       do j = 1,ny    ! ijk
          do i = 1,nx ! ijk
             rhs(i,j,k) = &                     ! ijk
                  + u(i+1,j  ,k  ) - u(i,j,k) & ! ijk
                  + v(i  ,j+1,k  ) - v(i,j,k) & ! ijk
                  + w(i  ,j  ,k+1) - w(i,j,k)   ! ijk
          enddo
       enddo
    enddo

  end subroutine set_rhs

  !-----------------------------------------------------------------------------------
  subroutine set_matrices()

    ! Define matrix coefficients cA
    ! Coefficients are stored in order of diagonals
    ! cA(1,:,:,:)      -> p(k,j,i)
    ! cA(2,:,:,:)      -> p(k-1,j,i)
    ! cA(3,:,:,:)      -> p(k+1,j-1,i)
    ! cA(4,:,:,:)      -> p(k,j-1,i)
    ! cA(5,:,:,:)      -> p(k-1,j-1,i)
    ! cA(6,:,:,:)      -> p(k+1,j,i-1)
    ! cA(7,:,:,:)      -> p(k,j,i-1)
    ! cA(8,:,:,:)      -> p(k-1,j,i-1)

    integer(kind=ip) :: lev
    integer(kind=ip) :: k,j,i
    integer(kind=ip) :: nx,ny,nz

    real(kind=rp), dimension(:,:),     pointer :: dx,dy
    real(kind=rp), dimension(:,:),     pointer :: dxu,dyv
    real(kind=rp), dimension(:,:,:),   pointer :: dzw
    real(kind=rp), dimension(:,:,:),   pointer :: Arx,Ary
    real(kind=rp), dimension(:,:)  ,   pointer :: Arz
    real(kind=rp), dimension(:,:,:),   pointer :: zxdy,zydx
    real(kind=rp), dimension(:,:,:),   pointer :: alpha
    real(kind=rp), dimension(:,:),     pointer :: beta
    real(kind=rp), dimension(:,:),     pointer :: gamu
    real(kind=rp), dimension(:,:),     pointer :: gamv
    real(kind=rp), dimension(:,:,:,:), pointer :: cA

    !    if (myrank==0) write(*,*)'   - set matrices:'

    do lev = 1, nlevs

       nx = grid(lev)%nx
       ny = grid(lev)%ny
       nz = grid(lev)%nz

       dx    => grid(lev)%dx    !
       dy    => grid(lev)%dy    !
       dxu   => grid(lev)%dxu   !
       dyv   => grid(lev)%dyv   !
       dzw   => grid(lev)%dzw   !
       Arx   => grid(lev)%Arx   !
       Ary   => grid(lev)%Ary   !
       Arz   => grid(lev)%Arz   !
       alpha => grid(lev)%alpha !
       beta  => grid(lev)%beta  !
       gamu  => grid(lev)%gamu  !
       gamv  => grid(lev)%gamv  !
       zxdy  => grid(lev)%zxdy  !
       zydx  => grid(lev)%zydx  !
       cA    => grid(lev)%cA    !

       !! interaction coeff with neighbours !!

       !---------------!
       !- lower level -!
       !---------------!
       k = 1
       do j = 1,ny+1  ! ijk
          do i = 1,nx ! ijk
             ! couples with k+1,j-1
             cA(3,i,j,k) = qrt * ( zydx(i,j,k+1) + zydx(i,j-1,k) ) ! ijk
             ! couples with j-1
             cA(4,i,j,k) = hlf * (gamv(i,j) + gamv(i,j-1)) * Ary(i,j,k) / dyv(i,j) & ! ijk
                  - qrt * ( zydx(i,j-1,k) - zydx(i,j,k) ) ! ijk
          enddo
       enddo

       do j = 1,ny      ! ijk
          do i = 1,nx+1 ! ijk
             ! couples with k+1,i-1
             cA(6,i,j,k) = qrt * ( zxdy(i,j,k+1) + zxdy(i-1,j,k) ) ! ijk
             ! couples with i-1
             cA(7,i,j,k) = hlf * (gamu(i,j) + gamu(i-1,j))  * Arx(i,j,k) / dxu(i,j) & ! ijk
                  - qrt * ( zxdy(i-1,j,k) - zxdy(i,j,k) ) ! ijk
          enddo
       enddo

       do j = 0,ny      ! ijk
          do i = 1,nx+1 ! ijk
             ! only for k==1, couples with j+1,i-1
             cA(5,i,j,k) = beta(i-1,j) + beta(i,j+1) ! ijk
          enddo
       enddo

       do j = 1,ny+1
          do i = 1,nx+1
             ! only for k==1, couples with j-1,i-1
             cA(8,i,j,k) = - beta(i-1,j) - beta(i,j-1) ! ijk
          enddo
       enddo

       !-------------------!
       !- interior levels -!
       !-------------------!
       do k = 2,nz-1     ! ijk
          do j = 1,ny    ! ijk
             do i = 1,nx ! ijk
                ! couples with k-1
                cA(2,i,j,k) =                           & ! ijk
                     ( Arz(i,j) / dzw(i,j,k) ) *        & ! ijk
                     hlf * (alpha(i,j,k)+alpha(i,j,k-1))  ! ijk
             enddo
          enddo
       enddo

       do k = 2,nz-1      ! ijk
          do j = 1,ny+1   ! ijk
             do i = 1,nx  ! ijk
                ! couples with k+1,j-1
                cA(3,i,j,k) = qrt * ( zydx(i,j,k+1) + zydx(i,j-1,k) ) ! ijk
                ! couples with j-1
                cA(4,i,j,k) =  Ary(i,j,k) / dyv(i,j)  ! ijk
                ! couples with k-1,j-1
                cA(5,i,j,k) = - qrt * ( zydx(i,j,k-1) + zydx(i,j-1,k) ) ! ijk
             enddo
          enddo
       enddo

       do k = 2,nz-1       ! ijk
          do j = 1,ny      ! ijk
             do i = 1,nx+1 ! ijk
                ! couples with k+1,i-1
                cA(6,i,j,k) = qrt * ( zxdy(i,j,k+1) + zxdy(i-1,j,k) ) ! ijk
                ! couples with i-1
                cA(7,i,j,k) = Arx(i,j,k) / dxu(i,j) ! ijk
                ! couples with k-1,i-1
                cA(8,i,j,k) = - qrt * ( zxdy(i,j,k-1) + zxdy(i-1,j,k) ) ! ijk
             enddo
          enddo
       enddo

       !---------------!
       !- upper level -!
       !---------------!
       k = nz         ! ijk
       do j = 1,ny    ! ijk
          do i = 1,nx ! ijk    
             ! couples with k-1
             cA(2,i,j,k) = (Arz(i,j)/dzw(i,j,k)) * hlf * ( alpha(i,j,k) + alpha(i,j,k-1) ) ! ijk
          enddo
       enddo

       do j = 1,ny+1  ! ijk
          do i = 1,nx ! ijk
             ! couples with j-1
             cA(4,i,j,k) = Ary(i,j,k) / dyv(i,j) &  ! ijk
                  + qrt * ( - zydx(i,j-1,k) + zydx(i,j,k) ) ! ijk
             ! couples with k-1,j-1
             cA(5,i,j,k) = - qrt * ( zydx(i,j,k-1) + zydx(i,j-1,k) ) ! ijk
          enddo
       enddo

       do j = 1,ny      ! ijk
          do i = 1,nx+1 ! ijk
             ! couples with i-1
             cA(7,i,j,k) = Arx(i,j,k) / dxu(i,j) &          ! ijk
                  + qrt * ( -zxdy(i-1,j,k) + zxdy(i,j,k) )  ! ijk
             ! couples with k-1,i-1
             cA(8,i,j,k) = - qrt * ( zxdy(i,j,k-1) + zxdy(i-1,j,k) ) ! ijk
          enddo
       enddo

       call fill_halo(lev,cA)

       !! self-interaction coeff !!

       !- lower level
       k = 1
       do j = 1,ny    ! ijk
          do i = 1,nx ! ijk
             cA(1,i,j,k) = &                                                          ! ijk
                  - Arz(i  ,j)    /dzw(i,j,k+1) * hlf * (alpha(i,j,k+1) + alpha(i,j,k)) & ! ijk
                  - Arx(i  ,j  ,k)/dxu(i  ,j  ) * hlf * (gamu(i,j) + gamu(i-1,j  )) & ! ijk
                  - Arx(i+1,j  ,k)/dxu(i+1,j  ) * hlf * (gamu(i,j) + gamu(i+1,j  )) & ! ijk
                  - Ary(i  ,j  ,k)/dyv(i  ,j  ) * hlf * (gamv(i,j) + gamv(i  ,j-1)) & ! ijk
                  - Ary(i  ,j+1,k)/dyv(i  ,j+1) * hlf * (gamv(i,j) + gamv(i  ,j+1))   ! ijk
          enddo
       enddo

       !- interior levels
       do k = 2,nz-1 
          do j = 1,ny    ! ijk
             do i = 1,nx ! ijk
                cA(1,i,j,k) = &                                                              ! ijk
                     - Arz(i  ,j    )/dzw(i,j,k+1) * hlf * (alpha(i,j,k+1) + alpha(i,j,k)) & ! ijk
                     - Arz(i  ,j    )/dzw(i,j,k  ) * hlf * (alpha(i,j,k-1) + alpha(i,j,k)) & ! ijk
                     - Arx(i  ,j  ,k)/dxu(i  ,j  )  & ! ijk
                     - Arx(i+1,j  ,k)/dxu(i+1,j  )  & ! ijk
                     - Ary(i  ,j  ,k)/dyv(i  ,j  )  & ! ijk
                     - Ary(i  ,j+1,k)/dyv(i  ,j+1)    ! ijk
             enddo
          enddo
       enddo

       !- upper level
       k=nz 
       do j = 1,ny    ! ijk
          do i = 1,nx ! ijk
             cA(1,i,j,k) = &                                                              ! ijk
                  - Arz(i  ,j    )/dzw(i,j,k+1) * alpha(i,j,k) &                          ! ijk
                  - Arz(i  ,j    )/dzw(i,j,k  ) * hlf * (alpha(i,j,k-1) + alpha(i,j,k)) & ! ijk
                  - Arx(i  ,j  ,k)/dxu(i  ,j  )  & ! ijk
                  - Arx(i+1,j  ,k)/dxu(i+1,j  )  & ! ijk
                  - Ary(i  ,j  ,k)/dyv(i  ,j  )  & ! ijk
                  - Ary(i  ,j+1,k)/dyv(i  ,j+1)    ! ijk

          enddo
       enddo

       if (netcdf_output) then
          if (myrank==0) write(*,*)'       write cA in a netcdf file'
          call write_netcdf(grid(lev)%cA,vname='ca',netcdf_file_name='cA.nc',rank=myrank,iter=lev)
       endif

    enddo

  end subroutine set_matrices

  !-------------------------------------------------------------------------     
  subroutine correct_uvw()

    !! u,v,w are fluxes, the correction is T*grad(p)

    integer(kind=ip):: i,j,k
    integer(kind=ip):: nx, ny, nz

    real(kind=rp) :: gamma
    real(kind=rp), dimension(:,:)  , pointer :: dx,dy
    real(kind=rp), dimension(:,:)  , pointer :: dxu,dyv
    real(kind=rp), dimension(:,:)  , pointer :: Arz
    real(kind=rp), dimension(:,:,:), pointer :: dzw
    real(kind=rp), dimension(:,:,:), pointer :: Arx,Ary
    real(kind=rp), dimension(:,:,:), pointer :: zxdy,zydx
    real(kind=rp), dimension(:,:,:), pointer :: alpha
    real(kind=rp), dimension(:,:)  , pointer :: beta
    real(kind=rp), dimension(:,:,:), pointer :: du,dv,dw
    real(kind=rp), dimension(:,:,:), pointer :: p

    real(kind=rp), dimension(:,:,:), pointer :: px,py,pz

    !    if (myrank==0) write(*,*)'   - compute pressure gradient and translate to fluxes'

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    dx    => grid(1)%dx    !
    dy    => grid(1)%dy    !
    dxu   => grid(1)%dxu   !
    dyv   => grid(1)%dyv   !
    dzw   => grid(1)%dzw   !
    Arx   => grid(1)%Arx   !
    Ary   => grid(1)%Ary   !
    Arz   => grid(1)%Arz   !
    alpha => grid(1)%alpha !
    beta  => grid(1)%beta  !
    zxdy  => grid(1)%zxdy  !
    zydx  => grid(1)%zydx  !
    p     => grid(1)%p

    px => grid(1)%u
    py => grid(1)%v
    pz => grid(1)%w

    !! Pressure gradient -
    do k = 1,nz         ! ijk
       do j = 0,ny+1    ! ijk
          do i = 1,nx+1 ! ijk
             px(i,j,k) = -one / dxu(i,j) * (p(i,j,k)-p(i-1,j,k)) ! ijk
          enddo
       enddo
    enddo

    do k = 1,nz         ! ijk
       do j = 1,ny+1    ! ijk
          do i = 0,nx+1 ! ijk
             py(i,j,k) = -one / dyv(i,j) * (p(i,j,k)-p(i,j-1,k)) ! ijk
          enddo
       enddo
    enddo

    k = 1 !bottom pressure gradient is undefined
    do j = 0,ny+1    ! ijk
       do i = 0,nx+1 ! ijk
          pz(i,j,k) = 9999999999. ! ijk
       enddo
    enddo

    do k = 2,nz !interior levels
       do j = 0,ny+1    ! ijk
          do i = 0,nx+1 ! ijk
             pz(i,j,k) = -one / dzw(i,j,k) * (p(i,j,k)-p(i,j,k-1)) ! ijk
          enddo
       enddo
    enddo

    k = nz+1 !surface
    do j = 0,ny+1    ! ijk
       do i = 0,nx+1 ! ijk
          pz(i,j,k) =  -one / dzw(i,j,k) * (-p(i,j,k-1)) ! ijk
       enddo
    enddo

    !! Correct U -

    du => grid(1)%du

    k = 1
    do j = 1,ny      ! ijk
       do i = 1,nx+1 ! ijk  
          gamma = one - qrt * ( &
               (zxdy(i  ,j,k)/dy(i  ,j))**2/alpha(i  ,j,k) + & ! ijk
               (zxdy(i-1,j,k)/dy(i-1,j))**2/alpha(i-1,j,k) )   ! ijk
          du(i,j,k) = gamma * Arx(i,j,k) * px(i,j,k) & ! ijk
               - qrt * ( &
               + zxdy(i  ,j,k) * dzw(i  ,j  ,k+1) * pz(i  ,j  ,k+1) & ! ijk
               + zxdy(i-1,j,k) * dzw(i-1,j  ,k+1) * pz(i-1,j  ,k+1) )  & ! ijk
               - beta(i-1,j)   * dyv(i-1,j  )     * py(i-1,j  ,k  ) & ! ijk
               - beta(i-1,j)   * dyv(i-1,j+1)     * py(i-1,j+1,k  ) & ! ijk
               - beta(i  ,j)   * dyv(i  ,j  )     * py(i  ,j  ,k  ) & ! ijk
               - beta(i  ,j)   * dyv(i  ,j+1)     * py(i  ,j+1,k  )   ! ijk
       enddo
    enddo

    do k = 2,nz-1 
       do j = 1,ny      ! ijk
          do i = 1,nx+1 ! ijk  
             du(i,j,k) = Arx(i,j,k) * px(i,j,k) & ! ijk
                  - qrt * ( &
                  + zxdy(i  ,j,k) * dzw(i  ,j,k  ) * pz(i  ,j,k  ) & ! ijk
                  + zxdy(i  ,j,k) * dzw(i  ,j,k+1) * pz(i  ,j,k+1) & ! ijk
                  + zxdy(i-1,j,k) * dzw(i-1,j,k  ) * pz(i-1,j,k  ) & ! ijk
                  + zxdy(i-1,j,k) * dzw(i-1,j,k+1) * pz(i-1,j,k+1) ) ! ijk
          enddo
       enddo
    enddo

    k = nz
    do j = 1,ny      ! ijk
       do i = 1,nx+1 ! ijk  
          du(i,j,k) = Arx(i,j,k) * px(i,j,k) & ! ijk
               - qrt * ( &
               + zxdy(i  ,j,k) *       dzw(i  ,j,k  ) * pz(i  ,j,k  ) & ! ijk
               + zxdy(i  ,j,k) * two * dzw(i  ,j,k+1) * pz(i  ,j,k+1) & ! ijk
               + zxdy(i-1,j,k) *       dzw(i-1,j,k  ) * pz(i-1,j,k  ) & ! ijk
               + zxdy(i-1,j,k) * two * dzw(i-1,j,k+1) * pz(i-1,j,k+1) ) ! ijk
       enddo
    enddo

    !! Correct V - 

    dv => grid(1)%dv

    k = 1
    do j = 1,ny+1  ! ijk
       do i = 1,nx ! ijk
          gamma = one - qrt * (  &
               (zydx(i,j  ,k)/dx(i,j  ))**2/alpha(i,j  ,k) + & ! ijk
               (zydx(i,j-1,k)/dx(i,j-1))**2/alpha(i,j-1,k) ) ! ijk
          dv(i,j,k) = gamma * Ary(i,j,k) * py(i,j,k) & ! ijk
               - qrt * ( &
               + zydx(i,j  ,k) * dzw(i  ,j  ,k+1) * pz(i  ,j  ,k+1) & ! ijk
               + zydx(i,j-1,k) * dzw(i  ,j-1,k+1) * pz(i  ,j-1,k+1) ) & ! ijk
               - beta(i,j-1)   * dxu(i  ,j-1)     * px(i  ,j-1,k  ) & ! ijk
               - beta(i,j-1)   * dxu(i+1,j-1)     * px(i+1,j-1,k  ) & ! ijk
               - beta(i,j  )   * dxu(i  ,j  )     * px(i  ,j  ,k  ) & ! ijk
               - beta(i,j  )   * dxu(i+1,j  )     * px(i+1,j  ,k  )   ! ijk
       enddo
    enddo

    do k = 2,nz-1
       do j = 1,ny+1  ! ijk
          do i = 1,nx ! ijk
             dv(i,j,k) =  Ary(i,j,k) * py(i,j,k) & ! ijk
                  - qrt * ( &
                  + zydx(i,j  ,k) * dzw(i,j  ,k  ) * pz(i,j  ,k  ) & ! ijk
                  + zydx(i,j  ,k) * dzw(i,j  ,k+1) * pz(i,j  ,k+1) & ! ijk
                  + zydx(i,j-1,k) * dzw(i,j-1,k  ) * pz(i,j-1,k  ) & ! ijk
                  + zydx(i,j-1,k) * dzw(i,j-1,k+1) * pz(i,j-1,k+1) ) ! ijk
          enddo
       enddo
    enddo

    k = nz
    do j = 1,ny+1  ! ijk
       do i = 1,nx ! ijk
          dv(i,j,k) = Ary(i,j,k) * py(i,j,k) & ! ijk
               - qrt * ( &
               + zydx(i,j  ,k)       * dzw(i,j  ,k  ) * pz(i,j  ,k  ) & ! ijk
               + zydx(i,j  ,k) * two * dzw(i,j  ,k+1) * pz(i,j  ,k+1) & ! ijk
               + zydx(i,j-1,k)       * dzw(i,j-1,k  ) * pz(i,j-1,k  ) & ! ijk
               + zydx(i,j-1,k) * two * dzw(i,j-1,k+1) * pz(i,j-1,k+1) ) ! ijk

       enddo
    enddo

    !! Correct W -

    dw => grid(1)%dw

    do k = 2,nz
       do j = 1,ny    ! ijk
          do i = 1,nx ! ijk
             dw(i,j,k) =  hlf * (alpha(i,j,k-1) + alpha(i,j,k)) * Arz(i,j) * pz(i,j,k) & ! ijk
                  - qrt * ( &
                  + zxdy(i,j,k  ) * dxu(i  ,j) * px(i  ,j,k  ) & ! ijk
                  + zxdy(i,j,k  ) * dxu(i+1,j) * px(i+1,j,k  ) & ! ijk
                  + zxdy(i,j,k-1) * dxu(i  ,j) * px(i  ,j,k-1) & ! ijk
                  + zxdy(i,j,k-1) * dxu(i+1,j) * px(i+1,j,k-1) ) & ! ijk
                  - qrt * ( &
                  + zydx(i,j,k  ) * dyv(i,j  ) * py(i,j  ,k  ) & ! ijk
                  + zydx(i,j,k  ) * dyv(i,j+1) * py(i,j+1,k  ) & ! ijk
                  + zydx(i,j,k-1) * dyv(i,j  ) * py(i,j  ,k-1) & ! ijk
                  + zydx(i,j,k-1) * dyv(i,j+1) * py(i,j+1,k-1) ) ! ijk
          enddo
       enddo
    enddo

    k = nz+1 
    do j = 1,ny    ! ijk
       do i = 1,nx ! ijk
          dw(i,j,k) = alpha(i,j,k-1) * Arz(i,j) * pz(i,j,k) & ! ijk
               - hlf * ( &
               + zxdy(i,j,k-1) * dxu(i  ,j) * px(i  ,j,k-1) & ! ijk
               + zxdy(i,j,k-1) * dxu(i+1,j) * px(i+1,j,k-1) ) & ! ijk
               - hlf * ( &
               + zydx(i,j,k-1) * dyv(i,j  ) * py(i,j  ,k-1) & ! ijk
               + zydx(i,j,k-1) * dyv(i,j+1) * py(i,j+1,k-1) ) ! ijk
       enddo
    enddo

  end subroutine correct_uvw

end module mg_projection

