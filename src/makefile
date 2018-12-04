PROGRAM = run
OBJECTS = variables_v3s.f90 \
       inputdata_v3s.f90 \
       outputdata_v3s.f90 \
       excesscalc_v3s.f90 \
       qchannel_mc_v3s.f90 \
       riverrouting_v3s.f90 \
       main_v3s.f90

FC = gfortran
FFLAGS = 

RM = rm -f


$(PROGRAM): $(OBJECTS)
	$(FC) -o $(PROGRAM) $(OBJECTS)

.SUFFIXES: .f90 .o
.f90.o:
	$(FC) -c $(FFLAGS) $*.f90 -o $*.o

all: $(PROGRAM)
.PHONY: clean
clean:
	$(RM) $(PROGRAM) *.o *.mod
