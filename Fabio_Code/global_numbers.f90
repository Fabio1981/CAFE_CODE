
module global_numbers

  ! Some numbers

  real(kind=8) zero,half,one,two,pii

  ! Driver stuff

  real(kind=8) xmin, xmax, dx, courant, t, dt
  integer  res_num, nvars
  integer  Nx, Nt
  integer  every_0D, every_1D

  ! Initial data parameters

  real(kind=8) amplitude, sigma, x0, diss_coef

end module global_numbers
