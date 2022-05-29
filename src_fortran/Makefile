#!/usr/bin/env gmake

EXEC     = featfem

FC       = gfortran
ifeq ($(FC),ifort)
 OPTFLAGS = -O3 -xHost -ipo
 FCFLAGS  =
else
 OPTFLAGS = -Ofast -fno-strict-overflow
 FCFLAGS  =
endif
LDFLAGS  =
INCDIR   = -I.
FCOMP    = $(FC) $(FCFLAGS) $(OPTFLAGS) $(INCDIR)
OBJDIR   = obj
OBJ      = $(filter %.o,$(SRC:%.f=$(OBJDIR)/%.o))
F90COMP  = $(FCOMP)
OBJDIRF90= $(OBJDIR)
OBJF90   = $(filter %.o,$(SRC:%.f90=$(OBJDIRF90)/%.o))

SRC= featfem.f usqgrd.f gmvout.f assemble.f bndry.f util.f mgrout.f timer.f90 \
     \
     lcl1.f lax1.f lcp1.f lcp2.f lcp3.f ll21.f llc1.f lcl2.f lcl3.f lsc1.f lsp1.f lvm.f \
     lli1.f ia11.f ia21.f ia01.f ib01.f ic21.f ic01.f id11.f id21.f ie01.f \
     ea00.f e011.f ap7.f ab07.f vb0.f ndfl.f ndfg.f ndfgl.f cb2q.f \
     sb0.f s2c.f s2v.f s2m.f s2db.f xsb0x.f sbd00.f sveb.f s2di.f s2a.f \
     zinit.f znew.f zdisp.f zclear.f zctype.f zfree.f zcpy.f ztype.f zlen.f of0.f omsg.f oerr.f otrc.f werr.f \
     dcopy.f scopy.f dscal.f daxpy.f ddot.f idamax.f

vpath %.f src src_featlib
vpath %.f90 src


.PHONY: all
all: greet $(EXEC)
	@echo "Done," $(EXEC) "is ready."
	@echo

.PHONY: greet
greet:
	@echo "Compiling module" $(EXEC) "in"
	@pwd


$(EXEC): $(OBJDIR) $(OBJ) $(OBJDIRF90) $(OBJF90)
	$(FCOMP) $(OBJ) $(OBJF90) $(LDFLAGS) -o $@

$(OBJDIR)/%.o : %.f
	$(FCOMP) -c -o $@ $<

$(OBJDIRF90)/%.o : %.f90
	$(F90COMP) -c -o $@ $<

.PHONY: clean
clean:
	-rm -f $(OBJDIR)/*.o

.PHONY: clean_exec
clean_exec:
	-rm -f $(EXEC)

.PHONY: purge
purge: clean clean_exec

.PHONY: purgeall
purge_all: purge
	@(for i in *-*-*-* ; do if [ -x $$i ] ; then echo rm $$i ; rm -f $$i ; fi ; done )
	-rm -f -r obj/*
	-rm -f -r *~ \src/*~ \src_featlib/*~
	-rm -f -r \data/*~ \data/*.out \data/*.err \data/*.prt \data/*.sys \data/*.trc

# Force remake
re: purge $(EXEC)
	$(EXEC)
