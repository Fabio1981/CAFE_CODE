
  module ode

  implicit none

  real(kind=8), dimension(:,:), private, allocatable, save :: tmp

  contains

  subroutine odestep(statevect, t, dt, method)

  implicit none

  interface
  subroutine calcrhs(statevect, t, rhs)
  real(kind=8), dimension(:),               intent(in)  :: statevect
  real(kind=8),                             intent(in)  :: t
  real(kind=8), dimension(size(statevect)), intent(out) :: rhs
  end subroutine calcrhs
  end interface

  real(kind=8), dimension(:),        intent(inout) :: statevect
  real(kind=8),                      intent(inout) :: t
  real(kind=8),                      intent(in)    :: dt
  character(len=*),                  intent(in)    :: method

  real(kind=8), dimension(size(statevect)) :: rhs
  integer                                  :: tmpsteps

  real(kind=8), parameter :: sixth = 1.0d0 / 6.0d0

  select case(method)

  case ('euler')
  tmpsteps = 0
  case ('rk2')
  tmpsteps = 1
  case ('rk3')
  tmpsteps = 2
  case ('rk4')
  tmpsteps = 3
  case DEFAULT
  tmpsteps = 2
  end select

  if (.NOT.(allocated(tmp)) .AND.  (tmpsteps > 0))  then
  allocate(tmp(size(statevect), tmpsteps))
  endif

  select case(method)

  case ('euler')
  call calcrhs(statevect, t, rhs)
  statevect = statevect + dt * rhs
  t = t + dt

  case ('rk2')
  call calcrhs(statevect, t, rhs)
  tmp(:,1) = statevect + dt * rhs
  t = t + dt

  call calcrhs(tmp(:,1), t, rhs)
  statevect = 0.5d0 * (statevect + tmp(:,1) + dt * rhs)

  case ('rk4')  
  call calcrhs(statevect, t, rhs)
  tmp(:, 1) = statevect + rhs * dt * 0.5d0
  t = t + dt / 2.0d0
  tmp(:, 3) = rhs
  tmp(:, 3) = rhs

  call calcrhs(tmp(:, 1), t, rhs)
  tmp(:, 2) = statevect + rhs * dt * 0.5d0
  tmp(:, 3) =  tmp(:, 3) + 2.0d0 * rhs
  call calcrhs(tmp(:, 2), t, rhs)
  tmp(:, 1) = statevect  + rhs * dt
  t = t + dt / 2.0d0
  tmp(:, 3) =  tmp(:, 3) + 2.0d0 * rhs
  call calcrhs(tmp(:, 1), t, rhs)
  tmp(:, 3) =  tmp(:, 3) + rhs
  statevect = statevect + dt / (6.0d0) *  tmp(:, 3)

  case DEFAULT
  stop 'method not supported'

  end select
  end subroutine odestep
  end module ode


