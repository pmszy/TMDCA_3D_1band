# Makefile for AVEDISDCA - library maintenance
# Author: C.E. Ekuma
# Please note that the compiler flags for optimization will
# have to be changed for machines other than linus alphas (fort) 
# and intels (pgi).
RM = rm -f *.mod
TMDCA = mod_tprof.o module_Global.o dateregupg.o geod.o ginit.o main.o meas.o\
	put.o  readin.o sumup.o tables.o broyden.o 
SPEC = mod_tprof.o module_Global.o spectra_standalone.o tables.o  analyze.o
EXETMDCA=AVEDISDCA_intel
EXESPEC=spectra
#
#################################
#       Compiler Flags          #
#################################
#INTEL_F90 = ifort
# Flags for debugging
#FLAGS = -132 -check all -c
# Flags for pentium, athlon, centrino
FLAGS = -O3 -132 -r8 -funroll-loops -inline all -g -traceback 
#FLAGS = -O5 -132 -r8 -traceback
#FLAGS = -O5 -xHost -ipo -no-prec-div -132 -r8 -funroll-loops -inline all -g -traceback
# Flags for recent Intel P4
#FLAGS = -132 -ipo -O3 -xP -c
F90=ifort 
#################################
#       Linker Flags            #
#################################
LINK = ifort
MKLPATH= /opt/intel/mkl/lib/intel64/ 
MKL_LIBS=-L$(MKLPATH) -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread_nonshared $(FFLAGS) $(LFLAGS)
FFLAGS = -openmp -fpp -i4
#################################

all: $(TMDCA) $(SPEC)
	$(LINK) $(TMDCA) $(FLAGS) $(MKL_LIBS) $(FFLAGS) $(LFLAGS) -parallel -o AVEDISDCA_intel;
	$(LINK)  $(SPEC) $(FLAGS) $(MKL_LIBS) $(FFLAGS) $(LFLAGS) -o spectra 

AVEDCA: $(TMDCA)
	$(LINK) $(TMDCA) $(FLAGS) $(MKL_LIBS) $(FFLAGS) $(LFLAGS) -parallel -o AVEDISDCA_intel

SPECT: $(SPEC)
	$(LINK) $(SPEC) $(FLAGS) $(MKL_LIBS) $(FFLAGS) $(LFLAGS) -o spectra 
  
mod_tprof.o: mod_tprof.f
	$(F90) $(FLAGS) $(FFLAGS) $(link) -c mod_tprof.f

module_Global.o: module_Global.for
	$(F90) $(FLAGS) $(FFLAGS) $(link) -c module_Global.for

dateregupg.o: module_Global.for dateregupg.f
	$(F90) $(FLAGS) $(FFLAGS) $(link) -c dateregupg.f

geod.o: module_Global.for geod.f
	$(F90) $(FLAGS) $(FFLAGS) $(link) -c geod.f

ginit.o: module_Global.for ginit.f
	$(F90) $(FLAGS) $(FFLAGS) $(link) -c ginit.f

main.o: module_Global.for main.f
	$(F90) $(FLAGS) $(FFLAGS) $(link) -c main.f

meas.o: module_Global.for meas.f
	$(F90) $(FLAGS) $(FFLAGS) $(link) -c meas.f

broyden.o : module_Global.for broyden.f
	$(F90) $(FLAGS) $(FFLAGS) $(link) -c broyden.f

put.o: module_Global.for put.f
	$(F90) $(FLAGS) $(FFLAGS) $(link) -c put.f


readin.o: module_Global.for readin.f
	$(F90) $(FLAGS) $(FFLAGS) $(link) -c readin.f

sumup.o: module_Global.for sumup.f
	$(F90) $(FLAGS) $(FFLAGS) $(link) -c sumup.f

tables.o: module_Global.for tables.f
	$(F90) $(FLAGS) $(FFLAGS) $(link) -c tables.f

spectra_standalone.o: module_Global.for spectra_standalone.f90
	$(F90) $(FLAGS) $(FFLAGS) $(link) -c spectra_standalone.f90

analyze.o: module_Global.for analyze.f90
	$(F90) $(FLAGS) $(FFLAGS) $(link) -c analyze.f90

clean:
	$(RM) $(EXETMDCA)  $(EXESPEC) $(TMDCA) $(SPEC)
realclean:
	 $(RM) *.o core dump.F global.mod module_Global.d work.pc *.F *~ *.f 

