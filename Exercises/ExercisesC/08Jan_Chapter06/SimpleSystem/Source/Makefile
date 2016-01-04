CC          = gcc
CLINKER     = gcc

CFLAGS      =   -Wall -O4 -g 
LIBS        = -lm
DEPEND= makedepend

SRC     = cbmc.c ran_uniform.c
OBJS    = cbmc.o ran_uniform.o
EXECS   = cbmc

default: cbmc

all: $(EXECS)

cbmc:$(OBJS)
	$(CLINKER) $(OPTFLAGS) -o cbmc $(OBJS) $(LIBS)

clean:
	/bin/rm -f *.o *~ $(EXECS)

.c.o:
	$(CC) $(CFLAGS) -c $*.c
