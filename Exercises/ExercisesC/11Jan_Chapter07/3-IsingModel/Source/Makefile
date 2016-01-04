CC          = gcc
CLINKER     = gcc

CFLAGS      =   -Wall -O4 -g 
LIBS        = -lm
DEPEND= makedepend
EXECS = ising

SRC        = ran_uniform.c ising.c
OBJS        = ran_uniform.o ising.o

default: ising

all: $(EXECS)

ising:$(OBJS)
	$(CLINKER) $(OPTFLAGS) -o ising $(OBJS) $(LIBS)

clean:
	/bin/rm -f *.o *~ $(EXECS)

.c.o:
	$(CC) $(CFLAGS) -c $*.c
