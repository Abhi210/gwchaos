# Makefile for SUN -- Ga Tech architecture
# Makefile for NONLINEAR-SMALL
FFLAGS= -O
LDFLAGS= -llapack
# change previous line to location of blas and lapack
# and comment out lines below as appropriate
                                                                                              
all : driv.e 

driv.e : driv.o lesnls.o utils.o 
	gfortran $(FFLAGS) -o driv.e driv.o lesnls.o utils.o $(LDFLAGS)
# uncomment previous two lines if have linked to blas/lapack
# comment next two lines if have not linked to blas/lapack
#driv.e : driv.o lesnls.o utils.o blalap.obj 
#       g77 $(FFLAGS) -o driv.e driv.o lesnls.o utils.o blalap.o $(LDFLAGS)

# DRIVERS and lesnls  
driv.o : driv.f 
	gfortran -c $(FFLAGS) driv.f
                                                                                              
lesnls.o : lesnls.f
	gfortran -c $(FFLAGS) lesnls.f

# COMMON CODES
# common routines we wrote
utils.o : utils.f 
	gfortran -c $(FFLAGS) utils.f
# routines from BLAS and LAPACK
# (un)comment next two lines if have (not) linked to blas and lapack
#blalap.o : blalap.f 
	gfortran -c $(FFLAGS) blalap.f
# CLEAN
clean : 
	rm -f a.out core *.o *.e *~
