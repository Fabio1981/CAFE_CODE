
subroutine boundaries

  use arrays
  use global_numbers

  implicit none

  u(1,0)  = zero
  u(1,Nx) = zero

  u(2,0)  = (one/10.0d0)*(15.0d0*u(2,1) - 6.0d0*u(2,2) + u(2,3))
  u(2,Nx) = zero

end subroutine boundaries
