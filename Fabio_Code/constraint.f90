
subroutine constraint

  use arrays
  use global_numbers
  use derivatives

  implicit none 

  integer i 

  call derivs(u(3,:),du(3,:),dx,Nx)

  do i=1,Nx-1

     Mom(i) = du(3,i) &
            - (one + two*(sin(x(i)))**2)*(1 - u(3,i))/(sin(x(i))*cos(x(i))) &
            + sin(x(i))*cos(x(i))*u(3,i)*(u(1,i)**2 + u(2,i)**2)

  end do

  Mom(0)  = 3.0D0*Mom(1) - 3.0D0*Mom(2) + Mom(3)

  Mom(Nx) = 3.0D0*Mom(Nx-1) - 3.0D0*Mom(Nx-2) + Mom(Nx-3)


end subroutine constraint
