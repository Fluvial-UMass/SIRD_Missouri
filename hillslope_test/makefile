PROGRAM = run
OBJECTS = variables_v3s.f90 \
       inputdata_v3s.f90 \
       outputdata_v3s.f90 \
       qsurface_kw_v1.f90 \
       hillsloperouting_v3s.f90 \
       excesscalc_v3s.f90 \
       qchannel_mc_v3s.f90 \
       riverrouting_v3s.f90 \
       main_v3s.f90

FC = ifort
FFLAGS = -O3

RM = rm -f


$(PROGRAM): $(OBJECTS)
	$(FC) -o $(PROGRAM) $(OBJECTS) $(FFLAGS)

.SUFFIXES: .f90 .o
.f90.o:
	$(FC) -c $(FFLAGS) $*.f90 -o $*.o

all: $(PROGRAM) \
     $(SUBPROGRAM)
.PHONY: clean
clean:
	$(RM) $(PROGRAM) *.o *.mod
