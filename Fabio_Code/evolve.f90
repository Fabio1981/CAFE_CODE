
  subroutine evolve

  use arrays 
  use global_numbers

  implicit none 

  integer l,k,i
  
  real(kind=8) dt_temp

!!$ *****************

!!$ u(1,:) = PHI
!!$ u(2,:) = PI
!!$ u(3,:) = A
!!$ u(4,:) = Delta
!!$ Hola Mundo 

!!$ *****************

  call allocate_memory
  call initialize_variables
  call grid1d  

  t = zero

  print *,'----------------------------'
  print *,'|  Time step  |    Time    |'
  print *,'----------------------------'

  write(*,"(A5,I6,A6,ES9.2,A3)") ' |   ',0,'    | ',t,'  |'

  call initial
  call constraint

  call save1Ddata(u(1,:),'u1',0)
  call save1Ddata(u(2,:),'u2',0)
  call save1Ddata(u(3,:),'u3',0)
  call save1Ddata(u(4,:),'u4',0)
  call save1Ddata(Mom,'Mom',0)

!!$ Main loop 

  do l=1,Nt

     t = t + dt

     if (mod(l,every_1D).eq.0) then
        write(*,"(A5,I6,A6,ES9.2,A3)") ' |   ',l,'    | ',t,'  |'
     end if

     u_p = u

     do k=1,3

        call rhs

        if (k.eq.1) then
           dt_temp = dt
           do i=1,Nx-1
              u(1,i) = u_p(1,i) + dt_temp*rhs_u(1,i)
              u(2,i) = u_p(2,i) + dt_temp*rhs_u(2,i)
           end do
!              u(3,:) = u_p(3,:) + dt_temp*rhs_u(3,:)
        else if (k.eq.2) then
           dt_temp = 0.25*dt
           do i=1,Nx-1
              u(1,i) = 0.75*u_p(1,i) + 0.25*u(1,i) + dt_temp*rhs_u(1,i)
              u(2,i) = 0.75*u_p(2,i) + 0.25*u(2,i) + dt_temp*rhs_u(2,i)
           end do
!              u(3,:) = 0.75*u_p(3,:) + 0.25*u(3,:) + dt_temp*rhs_u(3,:)
        else
           dt_temp = 2.0D0*dt/3.0D0
           do i=1,Nx-1
              u(1,i) = u_p(1,i)/3.0D0 + 2.0D0*u(1,i)/3.0D0 + dt_temp*rhs_u(1,i)
              u(2,i) = u_p(2,i)/3.0D0 + 2.0D0*u(2,i)/3.0D0 + dt_temp*rhs_u(2,i)
           end do
!              u(3,:) = u_p(3,:)/3.0D0 + 2.0D0*u(3,:)/3.0D0 + dt_temp*rhs_u(3,:)
        end if

        call boundaries
        call metric_evolve

     end do

     call constraint


     if (mod(l,every_0D).eq.0) then

     end if

     if (mod(l,every_1D).eq.0) then

        call save1Ddata(u(1,:),'u1',1)
        call save1Ddata(u(2,:),'u2',1)
        call save1Ddata(u(3,:),'u3',1)
        call save1Ddata(u(4,:),'u4',1)
        call save1Ddata(Mom,'Mom',1)

     end if

  end do

  end subroutine evolve
