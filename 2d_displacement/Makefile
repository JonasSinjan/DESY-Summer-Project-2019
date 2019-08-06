DEBUG = Y
DOUBLE = Y

FC     = gfortran
FFLAGS = -cpp -ffree-line-length-none -fbackslash -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1
LIBS   = 

FFTW3_PATH = /lustre/fs23/group/that/jonas/lib/fftwlib
#FFTW3_PATH = /usr
#FFTW3_PATH = /lustre/fs17/group/that/rs/mylibs/fftw-3.3.4

FFLAGS += -I${FFTW3_PATH}/lib -I${FFTW3_PATH}/include
ifeq ($(DOUBLE), Y)
  FFLAGS += -DDP
  LIBS   += -L${FFTW3_PATH}/lib -lfftw3
else
  LIBS   += -L${FFTW3_PATH}/lib -lfftw3f
endif

ifeq ($(DEBUG), Y)
  #FFLAGS += -O0 -g -fcheck=bounds -fbacktrace -ffpe-trap=invalid,zero,overflow
  FFLAGS += -O0 -g -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow
else
  FFLAGS += -O2 -march=native -mtune=native -funroll-loops
endif

SOURCES = $(wildcard *.f08)
OBJECTS = $(subst .f08,.o,$(SOURCES))

.PHONY: clean test

run.x: $(OBJECTS)
	$(FC) $(OBJECTS) $(LIBS) -o -fopenmp run.x

.SUFFIXES:
.SUFFIXES: .f08 .o
.f08.o:
	$(FC) -c $(FFLAGS) $<

clean:
	rm -f *.x *.o *.mod *.log *.DAT *~

include Makefile.depend
