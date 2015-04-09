
subroutine metric_evolve

use arrays
use global_numbers
use ode_fly

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
!!$ Integration of Delta

  w(1) = 1.0D0  !!$ A
  w(2) = 0.0D0  !!$ Delta

  step = 0
  r    = x(0)
  r_final = x(Nx)
  method = 'rk4'

  do i=0,Nx

     u(3,i) = w(1)     
     u(4,i) = w(2)

     call odestep_fly(w, r, dx, method)

  end do

  u(3,Nx) =  1.0d0

!!$ ***************************************************************

end subroutine metric_evolve

subroutine calcrhs_fly(w, r, rhs)

  use arrays
  use global_numbers

  implicit none

  real(kind=8), dimension(:),        intent(in)  :: w
  real(kind=8),                      intent(in)  :: r
  real(kind=8), dimension(size(w)),  intent(out) :: rhs

  Real(kind=8) :: PHI, PI, temp
  integer  aux, i

  aux = r/dx

  PHI = u(1,aux)

  PI = u(2,aux)

  rhs(1) = (one + two*(sin(r))**2)*(1 - w(1))/(sin(r)*cos(r)) &
         -sin(r)*cos(r)*w(1)*(PHI**2 + PI**2)

  rhs(2) = -sin(r)*cos(r)*(PHI**2 + PI**2)

end subroutine calcrhs_fly
