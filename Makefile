# Fortran compiler
FC = gfortran

# Optimization flags
#OPT = -O3 -fdefault-real-8 -fdefault-double-8 -fbounds-check -fbacktrace -Wall -floop-nest-optimize
OPT = -O3 -floop-nest-optimize

ifeq ($(BUILD), parallel)
	PAR = -fopenmp
endif

# Glueing stuff together
FFLAGS = $(OPT) $(PAR)

# Source directory
SOURCE_DIR = sources

# Main source file
SRC_MAIN = $(SOURCE_DIR)/particle_code.f90

# Main executable
MAIN_EXE = particle_code

# Object files
OBJ_FILES = nrecip.o interpolation.o constants.o variables.o aerodynamics.o io.o memory.o
OBJ_MAIN  = particle_code.o

all: $(MAIN_EXE) clean

$(MAIN_EXE): $(OBJ_FILES) $(OBJ_MAIN)
	$(FC) $(FFLAGS) -o $(@)  $(OBJ_FILES) $(OBJ_MAIN)

$(OBJ_FILES): $(SOURCE_DIR)/*.f90 Makefile
	$(FC) $(FFLAGS) -c $(SOURCE_DIR)/$(@:%.o=%.f90)

$(OBJ_MAIN): $(SOURCE_DIR)/*.f90 Makefile
	$(FC) $(FFLAGS) -c  $(SRC_MAIN)

clean:
	@rm -rf *.o *~ core *.mod

clobber:
	@rm -rf *.o *~ core *.mod $(MAIN_EXE)
