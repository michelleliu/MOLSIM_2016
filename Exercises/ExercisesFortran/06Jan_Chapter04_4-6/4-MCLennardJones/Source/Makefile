FCOMPFLAGS    =	-O -ffixed-line-length-132
CCOMPFLAGS    =	-O

FFLAGS        =	$(FCOMPFLAGS) -fsecond-underscore
CFLAGS        = $(CCOMPFLAGS)
LDFLAGS       =	$(FCOMPFLAGS)

LD            =	gfortran
FC            =	gfortran
CC            =	cc

MAKEFILE      =	Makefile
PROGRAM       =	mc_nvt

INCS	      =	system.inc conf.inc parameter.inc potential.inc

OBJS	      =	adjust.o          \
		energyparticle.o  \
		lattice.o         \
		mc_nvt.o          \
		mcmove.o          \
		readdat.o         \
		sample.o          \
		sstmm.o           \
		store.o           \
		energysystem.o    \
		writepdb.o        \
		ran_uniform.o

all:		$(PROGRAM)

$(PROGRAM)::	$(INCS)
		@/bin/rm -f $(OBJS) core

$(PROGRAM)::	$(OBJS) $(MAKEFILE)
		@$(LD) $(LDFLAGS) $(OBJS) -o $(PROGRAM)

clean:;		@rm -f $(OBJS) core
