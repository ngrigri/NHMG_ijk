module mg_vert_grids

  use mg_cst
  use mg_mpi
  use mg_tictoc
  use mg_grids
  use mg_mpi_exchange
  use mg_gather
  use mg_namelist
  use mg_netcdf_out

  implicit none

contains

  !----------------------------------------
  subroutine set_vert_grids(z_r,Hz)

    real(kind=rp), dimension(:,:,:), pointer, intent(in)  :: z_r,Hz

    integer(kind=ip) :: lev
    integer(kind=ip) :: nx,ny,nz
    integer(kind=ip) :: nxf,nxc
    integer(kind=ip) :: nyf,nyc
    integer(kind=ip) :: nzf,nzc

    real(kind=rp), dimension(:,:,:), pointer :: zr,dz

    real(kind=rp), dimension(:,:,:), pointer :: zrf,zrc 
    real(kind=rp), dimension(:,:,:), pointer :: dzf,dzc

    real(kind=rp), dimension(:,:)  , pointer :: dx,dy
    real(kind=rp), dimension(:,:,:), pointer :: dzw
    real(kind=rp), dimension(:,:,:), pointer :: Arx,Ary
    real(kind=rp), dimension(:,:,:), pointer :: zydx,zxdy
    real(kind=rp), dimension(:,:,:), pointer :: alpha
    real(kind=rp), dimension(:,:)  , pointer :: beta
    real(kind=rp), dimension(:,:)  , pointer :: gamu
    real(kind=rp), dimension(:,:)  , pointer :: gamv
    real(kind=rp), dimension(:,:)  , pointer :: Arz

    integer(kind=ip) :: i,j,k

    !    if (myrank==0) write(*,*)'   - set vertical grids:'

    do lev = 1, nlevs

       !! fill and coarsen zr and dz

       nx=grid(lev)%nx
       ny=grid(lev)%ny
       nz=grid(lev)%nz

       if (lev == 1) then ! zr,dz from croco

          grid(lev)%zr(0:nx+1,0:ny+1,1:nz) = z_r !only 1 extra point ?????????? ! ijk
          grid(lev)%dz(0:nx+1,0:ny+1,1:nz) = Hz  !only 1 extra point ?????????? ! ijk

       else               ! coarsen zr,dz (needed when directly discretizing on coarser grids)

          nxf = grid(lev-1)%nx
          nyf = grid(lev-1)%ny
          nzf = grid(lev-1)%nz

          zrf => grid(lev-1)%zr
          dzf => grid(lev-1)%dz

          if (grid(lev)%gather == 1) then
             nxc = nx/grid(lev)%ngx
             nyc = ny/grid(lev)%ngy
             nzc = nz
             allocate(zrc(0:nxc+1,0:nyc+1,1:nzc)) ! ijk
             allocate(dzc(0:nxc+1,0:nyc+1,1:nzc)) ! ijk
          else
             nxc = nx
             nyc = ny
             nzc = nz
             zrc => grid(lev)%zr
             dzc => grid(lev)%dz
          endif

          ! Call fine2coarse
          zrc(1:nxc,1:nyc,1:nzc) = eighth * (       & !only interior points ??????????  ! ijk
               zrf(1:nxf  :2,1:nyf  :2,1:nzf  :2) + & ! ijk
               zrf(1:nxf  :2,2:nyf+1:2,1:nzf  :2) + & ! ijk
               zrf(2:nxf+1:2,1:nyf  :2,1:nzf  :2) + & ! ijk
               zrf(2:nxf+1:2,2:nyf+1:2,1:nzf  :2) + & ! ijk
               zrf(1:nxf  :2,1:nyf  :2,2:nzf+1:2) + & ! ijk
               zrf(1:nxf  :2,2:nyf+1:2,2:nzf+1:2) + & ! ijk
               zrf(2:nxf+1:2,1:nyf  :2,2:nzf+1:2) + & ! ijk
               zrf(2:nxf+1:2,2:nyf+1:2,2:nzf+1:2) )   ! ijk

          ! Call fine2coarse
          dzc(1:nxc,1:nyc,1:nzc) = 2._rp * eighth * (        & !only interior points ?????????? ! ijk
               dzf(1:nxf  :2,1:nyf  :2,1:nzf  :2) + & ! ijk
               dzf(1:nxf  :2,2:nyf+1:2,1:nzf  :2) + & ! ijk
               dzf(2:nxf+1:2,1:nyf  :2,1:nzf  :2) + & ! ijk
               dzf(2:nxf+1:2,2:nyf+1:2,1:nzf  :2) + & ! ijk
               dzf(1:nxf  :2,1:nyf  :2,2:nzf+1:2) + & ! ijk
               dzf(1:nxf  :2,2:nyf+1:2,2:nzf+1:2) + & ! ijk
               dzf(2:nxf+1:2,1:nyf  :2,2:nzf+1:2) + & ! ijk
               dzf(2:nxf+1:2,2:nyf+1:2,2:nzf+1:2) )   ! ijk

          if (grid(lev)%gather == 1) then
             call gather(lev,zrc,grid(lev)%zr)
             call gather(lev,dzc,grid(lev)%dz)
             deallocate(zrc)
             deallocate(dzc)
          endif

       end if

       call fill_halo(lev,grid(lev)%zr) ! special fill_halo of zr (nh=2)
       call fill_halo(lev,grid(lev)%dz) ! special fill_halo of dz (nh=2)

       if (netcdf_output) then
          call write_netcdf(grid(lev)%zr,vname='zr',netcdf_file_name='zr.nc',rank=myrank,iter=lev)
          call write_netcdf(grid(lev)%dz,vname='dz',netcdf_file_name='dz.nc',rank=myrank,iter=lev)
       endif

       !! compute derived qties

       dx    => grid(lev)%dx
       dy    => grid(lev)%dy
       zr    => grid(lev)%zr
       dz    => grid(lev)%dz

       dzw   => grid(lev)%dzw
       Arx   => grid(lev)%Arx
       Ary   => grid(lev)%Ary
       Arz   => grid(lev)%Arz
       zxdy  => grid(lev)%zxdy
       zydx  => grid(lev)%zydx
       alpha => grid(lev)%alpha
       beta  => grid(lev)%beta
       gamu  => grid(lev)%gamu
       gamv  => grid(lev)%gamv

       !! Cell height 
       do j = 0,ny+1                      ! ijk
          do i = 0,nx+1                   ! ijk
             dzw(i,j,1) = hlf * dz(i,j,1) ! ijk
          enddo                           ! ijk
       enddo                              ! ijk
       do k = 2,nz
          do j = 0,ny+1                   ! ijk
             do i = 0,nx+1                ! ijk
                dzw(i,j,k) = hlf * (dz(i,j,k-1) + dz(i,j,k)) ! ijk
             enddo                        ! ijk
          enddo                           ! ijk
       enddo
       do j = 0,ny+1                      ! ijk
          do i = 0,nx+1                   ! ijk
             dzw(i,j,nz+1) = hlf * dz(i,j,nz) ! ijk
          enddo
       enddo

       !!  Cell faces area

       do k = 1,nz         ! ijk
          do j = 0,ny+1    ! ijk
             do i = 1,nx+1 ! ijk
                Arx(i,j,k) = hlf * ( dy(i,j) * dz(i,j,k) + dy(i-1,j) * dz(i-1,j,k) ) ! ijk
             enddo
          enddo
       enddo

       do k = 1,nz         ! ijk
          do j = 1,ny+1    ! ijk
             do i = 0,nx+1 ! ijk
                Ary(i,j,k) = hlf * ( dx(i,j) * dz(i,j,k) + dx(i,j-1) * dz(i,j-1,k) ) ! ijk
             enddo
          enddo
       enddo

       do j = 1,ny+1    ! ijk
          do i = 0,nx+1 ! ijk
             Arz(i,j) = dx(i,j) * dy(i,j) ! ijk
          enddo
       enddo

       !! Slopes! We need zr with 2 halo points !
       do k = 1, nz        ! ijk
          do j = 0,ny+1    ! ijk 
             do i = 0,nx+1 ! ijk      
                zxdy(i,j,k) = hlf * (( zr(i+1,j  ,k) - zr(i-1,j  ,k) ) / dx(i,j) ) * dy(i,j) ! ijk
                zydx(i,j,k) = hlf * (( zr(i  ,j+1,k) - zr(i  ,j-1,k) ) / dy(i,j) ) * dx(i,j) ! ijk
             enddo
          enddo
       enddo

       call set_phybound2zero(lev,zxdy,gt='u')
       call set_phybound2zero(lev,zydx,gt='v')

       !!- Used in set_matrices and fluxes
       do k = 1, nz        ! ijk
          do j = 0,ny+1    ! ijk 
             do i = 0,nx+1 ! ijk    
                alpha(i,j,k) = one + (zxdy(i,j,k)/dy(i,j))**2 + (zydx(i,j,k)/dx(i,j))**2 ! ijk
             enddo
          enddo
       enddo

       do j = 0,ny+1    ! ijk
          do i = 0,nx+1 ! ijk
             gamu(i,j) = one - hlf * ( zxdy(i,j,1) / dy(i,j) )**2 / alpha(i,j,1) ! ijk
          enddo
       enddo

       do j = 0,ny+1    ! ijk
          do i = 0,nx+1 ! ijk
             gamv(i,j) = one - hlf * ( zydx(i,j,1) / dx(i,j) )**2 / alpha(i,j,1) ! ijk
          enddo
       enddo

       do j = 0,ny+1    ! ijk
          do i = 0,nx+1 ! ijk
             beta(i,j) = eighth * zxdy(i,j,1)/dy(i,j) * zydx(i,j,1)/dx(i,j) * dz(i,j,1) / alpha(i,j,1) ! ijk
          enddo
       end do

       if (netcdf_output) then
          call write_netcdf(grid(lev)%dzw,vname='dzw',netcdf_file_name='dzw.nc',rank=myrank,iter=lev)
          call write_netcdf(grid(lev)%zxdy,vname='zxdy',netcdf_file_name='zxdy.nc',rank=myrank,iter=lev)
          call write_netcdf(grid(lev)%zydx,vname='zydx',netcdf_file_name='zydx.nc',rank=myrank,iter=lev)
          call write_netcdf(grid(lev)%alpha,vname='alpha',netcdf_file_name='alpha.nc',rank=myrank,iter=lev)
       endif

    enddo

  end subroutine set_vert_grids

end module mg_vert_grids
