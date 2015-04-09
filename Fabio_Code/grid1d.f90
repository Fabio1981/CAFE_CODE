
subroutine grid1d

  use arrays 
  use global_numbers

  implicit none

  integer i

  dx = (xmax - xmin)/dble(Nx)

  do i=0,Nx
     x(i)  = xmin + dble(i) * dx !- half * dx  
  end do

  dt = courant * dx

end subroutine grid1d
