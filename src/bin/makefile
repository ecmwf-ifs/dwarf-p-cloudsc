COMP?=UNKNOWN
GPROF?=0

MOD = mod
FC = error
FCOPTS =
FCFLAGS = 

CC = error
COPTS = 
CFLAGS = 

LD = $(FC)
LDFLAGS = 
LIBS = 

ifeq ($(COMP),UNKNOWN)
error:
	@echo "Invalid compiler COMP=$(COMP). Currently supported ones are COMP=GNU or COMP=INTEL or COMP=PGI"
	@exit 1 
endif

ifeq ($(COMP),GNU)
FC = gfortran -g -fPIC -fopenmp -fconvert=big-endian -fbacktrace
#FCOPTS = -O0 -finit-real=snan -ffpe-trap=invalid,zero,overflow
FCOPTS = -O2 -fstack-arrays -ffpe-trap=invalid,zero,overflow
FCFLAGS  = -J .module -I .module -D$(COMP) -DECMWF -DLINUX -I include/41r2/ifs -I include/41r2/ifsaux -I ifs/function
FCFLAGS += -fno-second-underscore -ffree-form 

CC = gcc -g -fPIC -fopenmp
COPTS = -O2
CFLAGS = -D$(COMP) -DECMWF -DLINUX -I include/41r2/ifs -I include/41r2/ifsaux -I ifs/function

LD = $(FC)
LDFLAGS = -Wl,--as-needed -ffast-math

FCVERS = gfortran -v
CCVERS = gcc -v

ifeq ($(GPROF),1)
FCFLAFS += -pg -DDEBUG_OUT
CFLAGS += -pg -DDEBUG_OUT
LDFLAGS += -pg
endif
endif

ifeq ($(COMP),INTEL)
FC = ifort -g -fpic -qopenmp -qopenmp-threadprivate compat -convert big_endian -traceback
#FCOPTS = -O0 -check uninit -fpe-all=0 -ftz -qno-openmp
#FCOPTS = -O0 -check stack,uninit,pointers,bounds -fpe-all=0 -ftz
FCOPTS = -O2 # -fpe-all=0 -ftz
FCFLAGS  = -module .module -I .module -D$(COMP) -DECMWF -DLINUX -I include/41r2/ifs -I include/41r2/ifsaux -I ifs/function
FCFLAGS += -assume byterecl -align array64byte
FCFLAGS += -warn nounused,nouncalled
FCFLAGS += -xHost
FCFLAGS += -finline-functions -finline-limit=500 -Winline
#FCFLAGS += -heap-arrays 64

CC = icc -g -qopenmp -qopenmp-threadprivate compat -fpic
COPTS = -O2
CFLAGS  = -D$(COMP) -DECMWF -DLINUX -I include/41r2/ifs -I include/41r2/ifsaux -I ifs/function
CFLAGS += -xHost

LD = $(FC)
LDFLAGS = -Wl,--as-needed

FCVERS = ifort -V
CCVERS = icc -V

ifeq ($(GPROF),1)
FCFLAFS += -p -DDEBUG_OUT
CFLAGS += -p -DDEBUG_OUT
LDFLAGS += -p
endif
endif

ifeq ($(COMP),PGI)
FC = pgfortran -g -fpic -mp=bind,allcores,numa -traceback -Mbyteswapio
FCOPTS = -O2 -Mstack_arrays -Mrecursive -Ktrap=fp -Kieee -Mdaz
#FCOPTS = -O0 -Mstack_arrays -Mrecursive -Mbounds -Mchkptr -Ktrap=fp -Kieee
FCFLAGS  = -module .module -I .module -D$(COMP) -DECMWF -DLINUX -I include/41r2/ifs -I include/41r2/ifsaux -I ifs/function

CC = pgcc -g -fpic -mp -traceback
COPTS = -O2
CFLAGS = -D$(COMP) -DECMWF -DLINUX -I include/41r2/ifs -I include/41r2/ifsaux -I ifs/function

LD = $(FC)
LDFLAGS = -Wl,--as-needed

FCVERS = pgfortran -V
CCVERS = pgcc -V

ifeq ($(GPROF),1)
FCFLAFS += -pg -DDEBUG_OUT
CFLAGS += -pg -DDEBUG_OUT
LDFLAGS += -pg
endif
endif

ifeq ($(COMP),CCE)
FC = ftn 
FCOPTS = -hflex_mp=conservative -hfp1 -hadd_paren -hcontiguous
FCFLAGS  = -J .module -I .module -D$(COMP) -DECMWF -DLINUX -I include/41r2/ifs -I include/41r2/ifsaux -I ifs/function

CC = cc
COPTS = -O0
CFLAGS = -D$(COMP) -DECMWF -DLINUX -I include/41r2/ifs -I include/41r2/ifsaux -I ifs/function

LD = $(FC)
LDFLAGS = -hbyteswapio -Ktrap=fp -Wl,--as-needed

FCVERS = ftn -V
CCVERS = cc -V

ifeq ($(GPROF),1)
# gprof not available ?
FCFLAFS += -DDEBUG_OUT
CFLAGS += -DDEBUG_OUT
LDFLAGS += 
endif
endif

LIB = libcloudsc.a

EXE = main.x.$(COMP)

OBJS = \
	dummy/abor1.o \
	dummy/routines.o \
	ifs/module/diag_mod.o \
	ifs/module/diff_mod.o \
	ifs/module/expand_mod.o \
	ifs/module/surface_fields_mix.o \
	ifs/module/timer_mod.o \
	ifs/module/yoecldp.o \
	ifs/module/yoephli.o \
	ifs/module/yoethf.o \
	ifs/module/yomcst.o \
	ifs/module/yomct0.o \
	ifs/module/yomct3.o \
	ifs/module/yomdim.o \
	ifs/module/yomdyn.o \
	ifs/module/yomjfh.o \
	ifs/module/yomlun.o \
	ifs/module/yommp0.o \
	ifs/module/yomphy2.o \
	ifs/module/yomphyder.o \
	ifs/phys_ec/cloud_layer.o \
	ifs/phys_ec/cloudsc.o \
	ifs/phys_ec/cloudsc_driver.o \
	ifs/phys_ec/cloudsc_in.o \
	ifs/phys_ec/cloudsc_out.o \
	ifs/phys_ec/cuadjtq.o \
	ifs/programs/cloudsc_dwarf.o \
	ifs/programs/main.o \
	ifsaux/module/parkind1.o \
	ifsaux/module/yomlun_ifsaux.o \
	ifsaux/utilities/mycpu.o

$(EXE): .module $(LIB)
	$(LD) $(LDFLAGS) ifs/programs/main.o -o $@ -L. -lcloudsc $(LIBS)
	-ldd $@

.module:
	mkdir -p .module || :

version:
	@echo " "
	-$(FCVERS)
	@echo " "
	-$(CCVERS)
	@echo " "

options:
	@echo " "
	@echo "=== Fortran compiler options:"
	@echo "$(FC) -c $(FCFLAGS) $(FCOPTS)"
	@echo " "
	@echo "=== C compiler options:"
	@echo "$(CC) -c $(CFLAGS) $(COPTS)"
	@echo " "
	@echo "=== Linking options:"
	@echo "$(LD) $(LDFLAGS) -L. -lcloudsc $(LIBS)"
	@echo " "

$(LIB): $(OBJS)
	/bin/rm -f $@
	ar -crv $@ $(OBJS)
	ar -dv $@ main.o
	ranlib $@

.SUFFIXES:
.SUFFIXES: .F90 .c

dummy/abor1.o: dummy/abor1.F90

dummy/routines.o: dummy/routines.F90 \
	include/41r2/ifsaux/abor1.intfb.h

ifs/module/diag_mod.o: ifs/module/diag_mod.F90 \
	ifsaux/module/parkind1.o

ifs/module/diff_mod.o: ifs/module/diff_mod.F90 \
	ifsaux/module/parkind1.o \
	ifs/module/diag_mod.o \
	ifs/module/expand_mod.o

ifs/module/expand_mod.o: ifs/module/expand_mod.F90 \
	ifs/module/diag_mod.o \
	ifs/module/yomphyder.o \
	ifsaux/module/parkind1.o

ifs/module/surface_fields_mix.o: ifs/module/surface_fields_mix.F90 \
	ifs/module/yomct0.o \
	ifs/module/yomdim.o \
	ifs/module/yomdyn.o \
	ifs/module/yomlun.o \
	ifsaux/module/parkind1.o \
	ifsaux/module/yomlun_ifsaux.o \
	include/41r2/ifs/dotprod2.intfb.h \
	include/41r2/ifs/dotprod3.intfb.h \
	include/41r2/ifsaux/abor1.intfb.h

ifs/module/timer_mod.o: ifs/module/timer_mod.F90 \
	ifsaux/module/parkind1.o

ifs/module/yoecldp.o: ifs/module/yoecldp.F90 \
	ifsaux/module/parkind1.o

ifs/module/yoephli.o: ifs/module/yoephli.F90 \
	ifsaux/module/parkind1.o

ifs/module/yoethf.o: ifs/module/yoethf.F90 \
	ifsaux/module/parkind1.o

ifs/module/yomcst.o: ifs/module/yomcst.F90 \
	ifsaux/module/parkind1.o

ifs/module/yomct0.o: ifs/module/yomct0.F90 \
	ifsaux/module/parkind1.o

ifs/module/yomct3.o: ifs/module/yomct3.F90 \
	ifsaux/module/parkind1.o

ifs/module/yomdim.o: ifs/module/yomdim.F90 \
	ifsaux/module/parkind1.o

ifs/module/yomdyn.o: ifs/module/yomdyn.F90 \
	ifsaux/module/parkind1.o

ifs/module/yomjfh.o: ifs/module/yomjfh.F90 \
	ifsaux/module/parkind1.o

ifs/module/yomlun.o: ifs/module/yomlun.F90 \
	ifsaux/module/parkind1.o \
	ifsaux/module/yomlun_ifsaux.o

ifs/module/yommp0.o: ifs/module/yommp0.F90 \
	ifsaux/module/parkind1.o

ifs/module/yomphy2.o: ifs/module/yomphy2.F90 \
	ifsaux/module/parkind1.o

ifs/module/yomphyder.o: ifs/module/yomphyder.F90 \
	ifsaux/module/parkind1.o

ifs/phys_ec/cloud_layer.o: ifs/phys_ec/cloud_layer.F90 \
	ifs/module/surface_fields_mix.o \
	ifs/module/yoecldp.o \
	ifs/module/yomcst.o \
	ifs/module/yomct0.o \
	ifs/module/yomct3.o \
	ifs/module/yomdim.o \
	ifs/module/yomdyn.o \
	ifs/module/yomlun.o \
	ifs/module/yommp0.o \
	ifs/module/yomphy2.o \
	ifs/module/yomphyder.o \
	ifsaux/module/parkind1.o \
	ifsaux/module/yomlun_ifsaux.o \
	include/41r2/ifs/cloudsc.intfb.h

ifs/phys_ec/cloudsc.o: ifs/phys_ec/cloudsc.F90 \
	ifs/function/fccld.func.h \
	ifs/function/fcttre.func.h \
	ifs/module/yoecldp.o \
	ifs/module/yoethf.o \
	ifs/module/yomcst.o \
	ifs/module/yomjfh.o \
	ifs/module/yommp0.o \
	ifs/module/yomphyder.o \
	ifsaux/module/parkind1.o \
	include/41r2/ifs/cuadjtq.intfb.h \
	include/41r2/ifsaux/abor1.intfb.h

ifs/phys_ec/cloudsc_driver.o: ifs/phys_ec/cloudsc_driver.F90 \
	ifs/module/diag_mod.o \
	ifs/module/diff_mod.o \
	ifs/module/expand_mod.o \
	ifs/module/timer_mod.o \
	ifs/module/yoecldp.o \
	ifs/module/yomphyder.o \
	ifsaux/module/parkind1.o \
	include/41r2/ifs/cloudsc.intfb.h \
	include/41r2/ifs/cloudsc_in.intfb.h \
	include/41r2/ifsaux/mycpu.intfb.h

ifs/phys_ec/cloudsc_in.o: ifs/phys_ec/cloudsc_in.F90 \
	ifs/module/yoecldp.o \
	ifs/module/yoephli.o \
	ifs/module/yoethf.o \
	ifs/module/yomcst.o \
	ifs/module/yomjfh.o \
	ifs/module/yommp0.o \
	ifs/module/yomphyder.o \
	ifsaux/module/parkind1.o \
	include/41r2/ifsaux/abor1.intfb.h

ifs/phys_ec/cloudsc_out.o: ifs/phys_ec/cloudsc_out.F90 \
	ifs/module/yoecldp.o \
	ifs/module/yoephli.o \
	ifs/module/yoethf.o \
	ifs/module/yomcst.o \
	ifs/module/yomjfh.o \
	ifs/module/yommp0.o \
	ifs/module/yomphyder.o \
	ifsaux/module/parkind1.o \
	include/41r2/ifsaux/abor1.intfb.h

ifs/phys_ec/cuadjtq.o: ifs/phys_ec/cuadjtq.F90 \
	ifs/function/fcttre.func.h \
	ifs/module/yoephli.o \
	ifs/module/yoethf.o \
	ifs/module/yomcst.o \
	ifs/module/yomjfh.o \
	ifsaux/module/parkind1.o

ifs/programs/cloudsc_dwarf.o: ifs/programs/cloudsc_dwarf.F90 \
	ifs/module/diag_mod.o \
	ifsaux/module/parkind1.o \
	include/41r2/ifs/cloudsc_driver.intfb.h

ifs/programs/main.o: ifs/programs/main.F90

ifsaux/module/parkind1.o: ifsaux/module/parkind1.F90

ifsaux/module/yomlun_ifsaux.o: ifsaux/module/yomlun_ifsaux.F90 \
	ifsaux/module/parkind1.o

ifsaux/utilities/mycpu.o: ifsaux/utilities/mycpu.c

%.o:	%.F90
	$(FC) -c $(FCFLAGS) $(FCOPTS) $< -o $@

%.o:	%.c
	$(CC) -c $(CFLAGS) $(COPTS) $< -o $@

clean:
	/bin/rm -f $(LIB)
	/bin/rm -f $(OBJS)
	/bin/rm -f .module/*.$(MOD)
	/bin/rm -rf .module

veryclean: clean
#	/bin/rm -f $(EXE)
	/bin/rm -f *.x*
	/bin/rm -f gmon.*
