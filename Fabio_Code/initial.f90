
subroutine initial

  use arrays
  use global_numbers
  use ode

  implicit none 

  integer i

!!$ ***************************************************************
!!$ THE RK stuff   

  integer, parameter            :: neq = 2 !!$ # of eqs.
  real(kind=8), dimension(neq)  :: w       !!$ solution vector
  real(kind=8), dimension(neq)  :: rhs     !!$ rhs vector
  integer                       :: step
  logical                       :: done = .false.
  character(len=42)             :: method
  real(kind=8) :: r, dr, r_final

!!$ ***************************************************************

!!$ u(1,:) = PHI
!!$ u(2,:) = PI
!!$ u(3,:) = A
!!$ u(4,:) = Delta

!!$ ***************************************************************

  u(1,:) = zero

  u(2,:) = (two*amplitude/pii)*exp(-4.0d0*(tan(x-x0))**2/(pii**2*sigma**2)) &
       + 1.0D-20

!!$ ***************************************************************
!!$ Integration of A and Delta

  w(1) = 1.0D0  !!$ A
  w(2) = 0.0D0  !!$ Delta

  step = 0
  r    = x(0)
  r_final = x(Nx)
  method = 'rk4'

  u(3,0) = w(1)
  u(4,0) = w(2)

  do i=0,Nx

     u(3,i) = w(1)
     u(4,i) = w(2)

     call odestep(w, r, dx, method)

  end do

!!$ ***************************************************************

  u(3,Nx) = 1.0d0

!!$ ***************************************************************

end subroutine initial

!!$ ***************************************************************
!!$ Subroutine that ODE subroutine calls. Sources
!!$ ***************************************************************

subroutine calcrhs(w, r, rhs)

  use global_numbers

  implicit none

  real(kind=8), dimension(:),        intent(in)  :: w
  real(kind=8),                      intent(in)  :: r
  real(kind=8), dimension(size(w)),  intent(out) :: rhs

  Real(kind=8) :: PHI, PI

  PHI = zero

  PI = (two*amplitude/pii)*exp(-4.0d0*(tan(r-x0))**2/(pii**2*sigma**2)) &
       + 1.0D-20

  rhs(1) = (one + two*(sin(r))**2)*(1 - w(1))/(sin(r)*cos(r)) &
         -sin(r)*cos(r)*w(1)*(PHI**2 + PI**2)

  rhs(2) = -sin(r)*cos(r)*(PHI**2 + PI**2)

end subroutine calcrhs
