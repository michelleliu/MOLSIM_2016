CC          = gcc
CLINKER     = gcc

CFLAGS      =   -Wall -O4 -g 
LIBS        = -lm
DEPEND= makedepend

SRC     = system.c ewald.c fourierspace.c realspace.c lattice.c ran_uniform.c
OBJS    = system.o ewald.o fourierspace.o realspace.o lattice.o ran_uniform.o

default: ewald

all: $(EXECS)

ewald:$(OBJS)
	$(CLINKER) $(OPTFLAGS) -o ewald $(OBJS) $(LIBS)

clean:
	/bin/rm -f *.o *~ $(EXECS)

.c.o:
	$(CC) $(CFLAGS) -c $*.c

ewald.o: system.h
fourierspace.o: system.h  
lattice.o: system.h  
realspace.o: system.h  
system.o: system.h
