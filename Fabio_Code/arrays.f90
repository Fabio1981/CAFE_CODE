
module arrays

  implicit none

!! Grid

  real(kind=8), allocatable, dimension (:) :: x

!! Scalar Field Variables

  real(kind=8), allocatable, dimension (:,:) :: u
  real(kind=8), allocatable, dimension (:,:) :: u_p
  real(kind=8), allocatable, dimension (:,:) :: rhs_u
  real(kind=8), allocatable, dimension (:,:) :: du

!! Constraint

  real(kind=8), allocatable, dimension (:) :: Mom

end module arrays
