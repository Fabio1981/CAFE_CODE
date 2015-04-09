#FC = mpif90
FC=gfortran
#FC=/opt/intel/fc/10.1.012/bin/ifort

FFLAGS = -O3
#FFLAGS = -g

#LNK = mpif90
LNK=gfortran
#LNK=/opt/intel/fc/10.1.012/bin/ifort

OBJS = main.o evolve.o allocate_memory.o initialize_variables.o grid1d.o initial.o rhs.o boundaries.o constraint.o metric_evolve.o save1Ddata.o save0Ddata.o 

MODS = arrays.o global_numbers.o derivatives.o ode.o ode_fly.o

$(OBJS):	$(MODS)

ads_wave:	$(OBJS) $(MODS)
		$(LNK) $(FFLAGS) -o xADS $(OBJS) $(MODS) 
	@ mkdir -p xxx
	@ mv xADS xxx

.PHONY:	clean

clean:
	-rm -f *.o *.mod xxx/xADS xxx/*xl xxx/*tl

%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@
