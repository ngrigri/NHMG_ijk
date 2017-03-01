module mg_horiz_grids

  use mg_cst
  use mg_mpi
  use mg_tictoc
  use mg_grids
  use mg_mpi_exchange
  use mg_gather
  use mg_namelist

  implicit none

contains

  !----------------------------------------
  subroutine set_horiz_grids(dx,dy)

    real(kind=rp), dimension(:,:), pointer, intent(in) :: dx,dy

    integer(kind=ip) :: lev
    integer(kind=ip) :: nx,ny
    integer(kind=ip) :: nxf,nxc
    integer(kind=ip) :: nyf,nyc

    real(kind=rp), dimension(:,:), pointer :: dxf,dxc
    real(kind=rp), dimension(:,:), pointer :: dyf,dyc
    real(kind=rp), dimension(:,:), pointer :: dxuf,dxuc
    real(kind=rp), dimension(:,:), pointer :: dyvf,dyvc

    if (myrank==0) write(*,*)'   - set horizontal grids:' 

    do lev = 1, nlevs

       nx=grid(lev)%nx
       ny=grid(lev)%ny

       if (lev == 1) then ! dx,dy from croco

          grid(lev)%dx(0:nx+1,0:ny+1) = dx ! ijk
          grid(lev)%dy(0:nx+1,0:ny+1) = dy ! ijk

          grid(lev)%dxu(1:nx+1,0:ny+1) = hlf*(dx(0:nx  ,0:ny+1)+dx(1:nx+1,0:ny+1)) ! ijk
          grid(lev)%dyv(0:nx+1,1:ny+1) = hlf*(dy(0:nx+1,0:ny  )+dy(0:nx+1,1:ny+1)) ! ijk

       else               ! coarsen dx,dy 
          ! (needed when directly discretizing on coarser grids)

          nxf =grid(lev-1)%nx
          nyf =grid(lev-1)%ny

          dxf => grid(lev-1)%dx
          dyf => grid(lev-1)%dy
          dxuf => grid(lev-1)%dxu
          dyvf => grid(lev-1)%dyv

          if (grid(lev)%gather == 1) then
             nxc= nx/grid(lev)%ngx
             nyc= ny/grid(lev)%ngy
             allocate( dxc(0:nxc+1,0:nyc+1)) ! ijk
             allocate( dyc(0:nxc+1,0:nyc+1)) ! ijk
             allocate(dxuc(0:nxc+1,0:nyc+1)) ! ijk
             allocate(dyvc(0:nxc+1,0:nyc+1)) ! ijk
          else
             nxc = nx
             nyc = ny
             dxc => grid(lev)%dx
             dyc => grid(lev)%dy
             dxuc => grid(lev)%dxu
             dyvc => grid(lev)%dyv
          endif

          dxc(1:nxc,1:nyc) = hlf      * ( & ! only interior points ! ijk
               dxf(1:nxf  :2,1:nyf  :2) + & ! ijk
               dxf(1:nxf  :2,2:nyf+1:2) + & ! ijk
               dxf(2:nxf+1:2,1:nyf  :2) + & ! ijk
               dxf(2:nxf+1:2,2:nyf+1:2) )   ! ijk

          dyc(1:nxc,1:nyc) = hlf      * ( & ! ijk
               dyf(1:nxf  :2,1:nyf  :2) + & ! ijk
               dyf(1:nxf  :2,2:nyf+1:2) + & ! ijk
               dyf(2:nxf+1:2,1:nyf  :2) + & ! ijk
               dyf(2:nxf+1:2,2:nyf+1:2) )   ! ijk

          dxuc(1:nxc,1:nyc) = hlf      * ( & ! only interior points ! ijk
               dxuf(1:nxf  :2,1:nyf  :2) + & ! ijk
               dxuf(1:nxf  :2,2:nyf+1:2) + & ! ijk
               dxuf(2:nxf+1:2,1:nyf  :2) + & ! ijk
               dxuf(2:nxf+1:2,2:nyf+1:2) )   ! ijk

          dyvc(1:nxc,1:nyc) = hlf      * ( & ! ijk
               dyvf(1:nxf  :2,1:nyf  :2) + & ! ijk
               dyvf(1:nxf  :2,2:nyf+1:2) + & ! ijk
               dyvf(2:nxf+1:2,1:nyf  :2) + & ! ijk
               dyvf(2:nxf+1:2,2:nyf+1:2) )   ! ijk

          if (grid(lev)%gather == 1) then
             call gather(lev,dxc,grid(lev)%dx)
             call gather(lev,dyc,grid(lev)%dy)
             call gather(lev,dxuc,grid(lev)%dxu)
             call gather(lev,dyvc,grid(lev)%dyv)
             deallocate(dxc)
             deallocate(dyc)
             deallocate(dxuc)
             deallocate(dyvc)
          endif

       endif

       call fill_halo(lev,grid(lev)%dx)
       call fill_halo(lev,grid(lev)%dy)
       call fill_halo(lev,grid(lev)%dxu)
       call fill_halo(lev,grid(lev)%dyv)

    enddo

  end subroutine set_horiz_grids

end module mg_horiz_grids
