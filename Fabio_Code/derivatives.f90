
module derivatives

  implicit none

contains

  subroutine derivs(u,du,dx,Nx)
    real(kind=8), dimension(0:Nx), intent(in) :: u
    real(kind=8), dimension(0:Nx), intent(inout) :: du
    real(kind=8) :: dx,idx
    integer :: i,Nx

    idx = 1.0d0/dx

    do i = 1, Nx-1
       du(i) = 0.5d0*(u(i+1) - u(i-1))*idx
    end do

    du(0) = -0.5d0*idx*(3.0D0*u(0) - 4.0D0*u(1) + u(2))
    du(Nx) = 0.5d0*idx*(3.0D0*u(Nx) - 4.0D0*u(Nx-1) + u(Nx-2))

  end subroutine derivs

  ! -------------------------------------------------------
  ! -----------    DISSIPATION  STARTS   ------------------
  ! -------------------------------------------------------
  
  !dissipation operator Q that satisfies (u,Qu) <= 0
  !we take Q = -h^3 (D_+ D_-)^2 in the interior
  !with appropriate changes near the boundary points
  
  subroutine diss(u,du,dx, Nx)
    real(kind=8), dimension(0:Nx), intent(in) :: u
    real(kind=8), dimension(0:Nx), intent(inout) :: du
    real(kind=8) :: dx
    integer :: i, Nx


    do i = 2, Nx-2
       du(i) = -( u(i+2) + u(i-2) - 4.0*(u(i+1) + u(i-1)) + 6.0*u(i) )/dx
    end do

    !something special for the first 2, last 2 points
    du(0) = 0.0
    du(1) = -( u(3) - 4.0*u(2) + 5.0*u(1) - 2.0*u(0) )/dx

    du(Nx-1) = -( u(Nx-3) - 4.0*u(Nx-2) + 5.0*u(Nx-1) - 2.0*u(Nx) )/dx
    du(Nx) = -2.0*( u(Nx-2) - 2.0*u(Nx-1) + u(Nx) )/dx

  end subroutine diss

end module derivatives
