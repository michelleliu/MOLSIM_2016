FCOMPFLAGS    =	-Wall -O
CCOMPFLAGS    = -Wall -O

FFLAGS        =	$(FCOMPFLAGS) -fsecond-underscore
CFLAGS        = $(CCOMPFLAGS)
LDFLAGS       =	$(FCOMPFLAGS)

LD            =	gfortran
FC            =	gfortran
CC            =	cc

MAKEFILE      =	Makefile
PROG1         =	hard2

INCS	      =	system.inc

OBJ1	      =	hard2.o writepdb.o ran_uniform.o sample.o sstmm.o

all:		$(PROG1)

$(PROG1)::	$(INCS)
		@/bin/rm -f $(OBJS) core

$(PROG1)::	$(OBJ1) $(OBJS) $(MAKEFILE)
		@$(LD) $(LDFLAGS) $(OBJ1) -o $(PROG1)

clean:;		@rm -f $(OBJ1) core
