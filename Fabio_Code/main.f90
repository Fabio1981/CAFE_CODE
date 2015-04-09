
program main

  use global_numbers

  implicit none

  integer Nxx,Ntt
  integer every_0Dt, every_1Dt

  Namelist /ads_Input/ xmin, xmax, &
                       res_num, nvars, &
                       Nxx, courant, Ntt, &
                       amplitude, sigma, x0, &
                       every_0Dt, every_1Dt, diss_coef

  open (3, file='iii.par', status = 'old' )
  read (3, nml = ads_Input)
  close(3)

  zero  = 0.0D0
  half  = 0.5D0
  one   = 1.0D0
  two   = 2.0D0
  pii = 4.0d0*atan(1.0d0)

  Nx = 2**(res_num-1)*Nxx
  Nt = 2**(res_num-1)*Ntt
  every_0D = 2**(res_num-1)*every_0Dt
  every_1D = 2**(res_num-1)*every_1Dt

  call evolve

  print *
  print *, 'PROGRAM CAFE-ads HAS FINISHED'
  print *

end program main

