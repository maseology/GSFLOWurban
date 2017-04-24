# mf2k5 libraray makefile
# $Id: Makefile 8266 2016-05-09 16:13:58Z rsregan $

include ../makelist

MF2K5LIB = $(LIBDIR)/libmf2k5.a

LIBOBJS = \
        utl7.o \
        gwf2bas7_NWT.o \
        pcg7_NWT.o \
        sip7_NWT.o \
        mf2kgmg.o \
	mhc7.o \
        gmg7.o \
        de47_NWT.o \
        gwf2bcf7.o \
        gwf2lpf7.o \
        gwf2huf7.o \
        hufutl7.o \
        mach_mod.o \
        NWT1_module.o \
        gwf2upw1.o \
        gwfuzfmodule_NWT.o \
        gwfsfrmodule_NWT.o \
        gwflakmodule_NWT.o \
        gwf2lak7_NWT.o \
        gwf2sfr7_NWT.o \
        gwf2uzf1_NWT.o \
        gwf2wel7_NWT.o \
        gwf2ghb7_NWT.o \
        gwf2fhb7.o \
        gwf2mnw17_NWT.o \
        gwf2mnw27_NWT.o \
        gwf2mnw2i7.o \
        NWT1_ilupc_mod.o \
        NWT1_xmdlib.o \
        NWT1_xmd.o \
        NWT1_gmres.o \
        NWT1_solver.o \
        gwf2chd7.o \
        gwf2hfb7_NWT.o \
        gwf2gag7.o \
        Irestart.o \
        parutl7.o \
        r_vector.o \
        solvers.o \
        ccfd.o \
	obs2bas7.o \
	obs2chd7.o \
	obs2ghb7.o

all: lib
	echo made $(MF2K5LIB)

install: lib

lib: $(MF2K5LIB)

$(MF2K5LIB): $(LIBOBJS)
	$(RM) $(MF2K5LIB)
	$(AR) $(MF2K5LIB) $(LIBOBJS)
	$(RANLIB) $(MF2K5LIB)

clean:
	$(RM) *.o *.mod *~ $(MF2K5LIB)

# on Windows machines comment out next 3 lines
openspec.inc:
	$(RM) openspec.inc
	$(CP) openspec_linux.inc openspec.inc

modules.o: modules.f90 mach_mod.o
	$(FC) -c $(FFLAGS) modules.f90

mach_mod.o: mach_mod.f90
	$(FC) -c $(FFLAGS) mach_mod.f90

NWT1_gmres.o: NWT1_gmres.f90
	$(FC) -c $(FFLAGS) NWT1_gmres.f90

NWT1_ilupc_mod.o: NWT1_ilupc_mod.f90
	$(FC) -c $(FFLAGS) NWT1_ilupc_mod.f90

ccfd.o: ccfd.c ccfd.h
	$(FC) -c $(CFLAGS) ccfd.c

mhc7.o: mhc7.f
	$(FC) -c $(FFLAGS) mhc7.f

mf2kgmg.o: mf2kgmg.c mf2kgmg.h
	$(CC) -c $(CFLAGS) mf2kgmg.c

mf2kgmg_wrap.o: mf2kgmg_wrap.c mf2kgmg.h
	$(CC) -c $(CFLAGS) mf2kgmg_wrap.c

r_vector.o: r_vector.c r_vector.h
	$(CC) -c $(CFLAGS) r_vector.c

NWT1_module.o: NWT1_module.f
	$(FC) -c $(FFLAGS) NWT1_module.f

solvers.o: solvers.c solvers.h
	$(CC) -c $(CFLAGS) solvers.c

NWT1_xmd.o: NWT1_xmd.f NWT1_xmdlib.o gwf2bas7_NWT.o NWT1_module.o
	$(FC) -c $(FFLAGS) NWT1_xmd.f

NWT1_xmdlib.o: NWT1_xmdlib.f
	$(FC) -c $(FFLAGS) NWT1_xmdlib.f

NWT1_solver.o: NWT1_solver.f gwf2bas7_NWT.o NWT1_module.o gwf2upw1.o NWT1_xmd.o NWT1_gmres.o NWT1_ilupc_mod.o mach_mod.o
	$(FC) -c $(FFLAGS) NWT1_solver.f

utl7.o: utl7.f openspec.inc
	$(FC) $(FFLAGS) -c utl7.f

gwf2bas7_NWT.o: gwf2bas7_NWT.f openspec.inc
	$(FC) $(FFLAGS) -c gwf2bas7_NWT.f

gwf2bcf7.o: gwf2bcf7.f gwf2bas7_NWT.o
	$(FC) $(FFLAGS) -c gwf2bcf7.f

hufutl7.o: hufutl7.f gwf2bas7_NWT.o gwf2huf7.o
	$(FC) $(FFLAGS) -c hufutl7.f

gwf2upw7.o: gwf2upw1.f gwf2bas7_NWT.o NWT1_module.o
	$(FC) $(FFLAGS) -c gwf2upw1.f

gwf2wel7.o: gwf2wel7.f gwf2bas7_NWT.o NWT1_module.o gwf2upw1.o
	$(FC) $(FFLAGS) -c gwf2wel7.f

gwf2sfr7_NWT.o: gwf2sfr7_NWT.f gwf2bas7_NWT.o gwfsfrmodule_NWT.o gwf2lpf7.o gwf2bcf7.o gwf2huf7.o gwf2upw1.o gwfuzfmodule_NWT.o NWT1_module.o gwflakmodule_NWT.o
	$(FC) $(FFLAGS) -c gwf2sfr7_NWT.f

gwf2uzf1_NWT.o: gwf2uzf1_NWT.f gwf2bas7_NWT.o gwfsfrmodule_NWT.o gwf2lpf7.o gwf2bcf7.o gwf2huf7.o gwf2upw1.o gwfuzfmodule_NWT.o NWT1_module.o gwflakmodule_NWT.o
	$(FC) $(FFLAGS) -c gwf2uzf1_NWT.f

gwf2mnw17_NWT.o: gwf2mnw17_NWT.f gwf2bas7_NWT.o gmg7.o sip7_NWT.o de47_NWT.o pcg7_NWT.o NWT1_module.o gwf2bcf7.o gwf2lpf7.o gwf2huf7.o gwf2upw1.o
	$(FC) $(FFLAGS) -c gwf2mnw17_NWT.f

gwf2mnw27_NWT.o: gwf2mnw27_NWT.f gwf2bas7_NWT.o gmg7.o sip7_NWT.o de47_NWT.o pcg7_NWT.o NWT1_module.o gwf2bcf7.o gwf2lpf7.o gwf2huf7.o gwf2upw1.o
	$(FC) $(FFLAGS) -c gwf2mnw27_NWT.f

gwf2mnw2i7.o: gwf2mnw2i7.f gwf2bas7_NWT.o gwf2mnw27_NWT.o
	$(FC) $(FFLAGS) -c gwf2mnw2i7.f

gwf2chd7.o: gwf2chd7.f gwf2bas7_NWT.o
	$(FC) $(FFLAGS) -c gwf2chd7.f

gwf2fbh7.o: gwf2fhb7.f gwf2bas7_NWT.o
	$(FC) $(FFLAGS) -c gwf2fhb7.f

gwf2gag7.o: gwf2gag7.f gwf2bas7_NWT.o gwflakmodule_NWT.o gwf2sfr7_NWT.o
	$(FC) $(FFLAGS) -c gwf2gag7.f

gwf2ghb7_NWT.o: gwf2ghb7_NWT.f gwf2bas7_NWT.o
	$(FC) $(FFLAGS) -c gwf2ghb7_NWT.f

gwf2lpf7.o: gwf2lpf7.f gwf2bas7_NWT.o
	$(FC) $(FFLAGS) -c gwf2lpf7.f

gwf2hfb7_NWT.o: gwf2hfb7_NWT.f openspec.inc gwf2bas7_NWT.o
	$(FC) $(FFLAGS) -c gwf2hfb7_NWT.f

pcg7_NWT.o: pcg7_NWT.f gwf2bas7_NWT.o
	$(FC) $(FFLAGS) -c pcg7_NWT.f

sip7_NWT.o: sip7_NWT.f gwf2bas7_NWT.o
	$(FC) $(FFLAGS) -c sip7_NWT.f

gmg7.o: gmg7.f gwf2bas7_NWT.o mhc7.o
	$(FC) $(FFLAGS) -c gmg7.f

Irestart.o: Irestart.f gwfuzfmodule_NWT.o gwfsfrmodule_NWT.o gsfmodflow.mod
	$(FC) $(FFLAGS) -c Irestart.f

gsfmodflow.mod:
	$(CD) ../gsflow;make gsfmodflow.mod
	$(CP) ../gsflow/gsfmodflow.mod .

de47_NWT.o: de47_NWT.f gwf2bas7_NWT.o
	$(FC) $(FFLAGS) -c de47_NWT.f

parutl7.o: parutl7.f gwf2bas7_NWT.o
	$(FC) $(FFLAGS) -c parutl7.f

obs2bas7.o: obs2bas7.f gwf2bas7_NWT.o
	$(FC) $(FFLAGS) -c obs2bas7.f

obs2chd7.o: obs2chd7.f gwf2bas7_NWT.o
	$(FC) $(FFLAGS) -c obs2chd7.f

gwf2lak7_NWT.o: gwf2lak7_NWT.f gwflakmodule_NWT.o gwf2bas7_NWT.o gwfsfrmodule_NWT.o gwfuzfmodule_NWT.o
	$(FC) $(FFLAGS) -c gwf2lak7_NWT.f

obs2ghb7.o: obs2ghb7.f gwf2bas7_NWT.o gwf2ghb7_NWT.o
	$(FC) $(FFLAGS) -c obs2ghb7.f

gwfsfrmodule_NWT.o: gwfsfrmodule_NWT.f
	$(FC) $(FFLAGS) -c gwfsfrmodule_NWT.f

gwfuzfmodule_NWT.o: gwfuzfmodule_NWT.f
	$(FC) $(FFLAGS) -c gwfuzfmodule_NWT.f

gwflakmodule_NWT.o: gwflakmodule_NWT.f
	$(FC) $(FFLAGS) -c gwflakmodule_NWT.f

global.mod: gwf2bas7_NWT.o
	$(FC) $(FFLAGS) -c gwf2bas7_NWT.f
#  end
