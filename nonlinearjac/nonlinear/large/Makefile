# Makefile for SUN -- Ga Tech architecture
# Makefile for NONLINEAR-LARGE
FFLAGS= -O
LDFLAGS= -llapack
# change previous line to location of blas and lapack
# and comment out lines below as appropriate
                                                                                              
all : driv.e 

driv.e : driv.o lesnll.o utils.o 
	g77 $(FFLAGS) -o driv.e driv.o lesnll.o utils.o $(LDFLAGS)
# uncomment previous two lines if have linked to blas/lapack
# comment next two lines if have not linked to blas/lapack
#driv.e : driv.o lesnll.o utils.o blalap.obj 
#       g77 $(FFLAGS) -o driv.e driv.o lesnll.o utils.o blalap.o $(LDFLAGS)

# DRIVERS and lesnll  
driv.o : driv.f 
	g77 -c $(FFLAGS) driv.f
                                                                                              
lesnll.o : lesnll.f
	g77 -c $(FFLAGS) lesnll.f

# COMMON CODES
# common routines we wrote
utils.o : utils.f 
	g77 -c $(FFLAGS) utils.f
# routines from BLAS and LAPACK
# (un)comment next two lines if have (not) linked to blas and lapack
#blalap.o : blalap.f 
	g77 -c $(FFLAGS) blalap.f
# CLEAN
clean : 
	rm -f a.out core *.o *.e *~
