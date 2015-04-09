
subroutine allocate_memory

  use arrays
  use global_numbers

  implicit none

  !! Grid

  allocate(x(0:Nx))

  !! Scalar Field Variables

  allocate(u(1:nvars,0:Nx))
  allocate(u_p(1:nvars,0:Nx))
  allocate(rhs_u(1:nvars,0:Nx))
  allocate(du(1:nvars,0:Nx))

!!$ contraint

  allocate(Mom(0:Nx))

end subroutine allocate_memory
