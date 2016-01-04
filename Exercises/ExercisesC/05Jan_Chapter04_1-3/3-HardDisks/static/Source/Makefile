CC          = gcc
CLINKER     = gcc

CFLAGS      =   -Wall -O4 -g 
LIBS        = -lm
DEPEND= makedepend

SRC     = hard1.c sample.c system.c writepdb.c ran_uniform.c
OBJS    = hard1.o sample.o system.o writepdb.o ran_uniform.o

default: hard1

all: $(EXECS)

hard1:$(OBJS)
	$(CLINKER) $(OPTFLAGS) -o hard1 $(OBJS) $(LIBS)

clean:
	/bin/rm -f *.o *~ $(EXECS)

.c.o:
	$(CC) $(CFLAGS) -c $*.c

hard1.o: system.h
sample.o: system.h
system.o: system.h
writepdb.o: system.h

