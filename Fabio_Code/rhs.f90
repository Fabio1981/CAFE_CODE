
subroutine rhs

  use arrays
  use global_numbers
  use derivatives

  implicit none

  integer i,k

  do k=1,nvars
     call derivs(u(k,:),du(k,:),dx,Nx)
  end do

  do i=1,Nx-1

     rhs_u(1,i) = u(3,i)*exp(-u(4,i))*du(2,i) &
                - u(3,i)*exp(-u(4,i))*du(4,i)*u(2,i) &
                + du(3,i)*exp(-u(4,i))*u(2,i) 

     rhs_u(2,i) = u(3,i)*exp(-u(4,i))*du(1,i) &
                - u(3,i)*exp(-u(4,i))*du(4,i)*u(1,i) &
                + du(3,i)*exp(-u(4,i))*u(1,i) &
                + two*u(3,i)*exp(-u(4,i))*u(1,i)/(sin(x(i))*cos(x(i))) 

  end do

!     rhs_u(3,:) = -two*sin(x)*cos(x)*u(3,:)**2*exp(-u(4,:))*u(1,:)*u(2,:)    

  do k=1,nvars
      call diss(u(k,:),du(k,:),dx,Nx)
   end do

  do k=1,2
   rhs_u(k,Nx-1) = rhs_u(k,Nx-1) + diss_coef*du(k,Nx-1)
  end do

end subroutine rhs
