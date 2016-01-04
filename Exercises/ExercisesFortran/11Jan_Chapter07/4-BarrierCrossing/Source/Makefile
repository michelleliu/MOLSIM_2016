FCOMPFLAGS    =	-O -v
CCOMPFLAGS    =	-O -v

FFLAGS        =	$(FCOMPFLAGS) -fsecond-underscore
CFLAGS        = $(CCOMPFLAGS)
LDFLAGS       =	$(FCOMPFLAGS)

LD            =	gfortran
FC            =	gfortran
CC            =	cc

MAKEFILE      =	Makefile
PROGRAM       =	barrier2

INCS	      =	system.inc

OBJS	      =	barrier2.o          \
		force.o             \
		integrate.o         \
		mdloop.o            \
		ran_uniform.o       \
		ran_vel.o           \
		ran_gauss.o         \
		readdat.o           \
		sample.o            \
		sstmm.o

all:		$(PROGRAM)

$(PROGRAM)::	$(INCS)
		@/bin/rm -f $(OBJS) core

$(PROGRAM)::	$(OBJS) $(MAKEFILE)
		@$(LD) $(LDFLAGS) $(OBJS) -o $(PROGRAM)

clean:;		@rm -f $(OBJS) core
