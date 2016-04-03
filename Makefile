FC=gfortran
FFLAGS=-finit-local-zero -Og -ffixed-line-length-132 -Wall -g -fbounds-check -fbacktrace
F90FLAGS = -finit-local-zero -Og -ffree-line-length-none -x f95-cpp-input -Wall -g -fbounds-check -fbacktrace
#FFLAGS=-finit-local-zero -O3 -ffixed-line-length-132 -Wall
#F90FLAGS = -finit-local-zero -O0 -ffree-line-length-none -x f95-cpp-input -Wall
QUENCH  = -L$(HOME)/SimulatedAnnealing/quench_anneal/lib -lquench -lquench_seq
LINPACK = -L$(HOME)/lib2/linpack -llinpack
BLAS = -L$(HOME)/lib2/blas -lblas
#LIBS = $(LINPACK) $(QUENCH)
#LIBS = $(QUENCH) $(LINPACK) $(BLAS)
LIBS = $(BLAS)
LIBS = -llapack-3
#LIBS = -L/usr/lib/libblas -lblas -L/usr/lib -llapack_atlas
LIBS = -L/home/cyrus/lib2/lapack -llapack -L/usr/lib/libblas -lblas -L/home/cyrus/lib -lcyrus

.SUFFIXES:
.SUFFIXES: .f90 .f95 .o .f .c

.f90.o:
	$(FC) $(F90FLAGS) -o $@ -c $<

OBJS = types.o rannyu.o tools.o more_tools.o fixed_node.o
OBJS2 = types.o rannyu.o tools.o more_tools.o fixed_node2.o
OBJSIMPSAMP = types.o tools.o more_tools.o imp_samp.o
OBJSIMPSAMP2 = types.o rannyu.o tools.o more_tools.o imp_samp2.o
OBJSTEST = types.o test.o

all: fixed_node imp_samp
#all: fixed_node fixed_node2

fixed_node: $(OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) $(OBJS) $(LIBS) -o $@

fixed_node2: $(OBJS2)
	$(FC) $(FFLAGS) $(LDFLAGS) $(OBJS2) $(LIBS) -o $@

imp_samp: $(OBJSIMPSAMP)
	$(FC) $(FFLAGS) $(LDFLAGS) $(OBJSIMPSAMP) $(LIBS) -o $@

imp_samp2: $(OBJSIMPSAMP2)
	$(FC) $(FFLAGS) $(LDFLAGS) $(OBJSIMPSAMP2) $(LIBS) -o $@

test: $(OBJSTEST)
	$(FC) $(FFLAGS) $(LDFLAGS) $(OBJSTEST) $(LIBS) -o $@
