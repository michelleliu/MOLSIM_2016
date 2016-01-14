CC          = gcc
CLINKER     = gcc

CFLAGS      =   -Wall -O4 -g -DDANGLEMAX=${angle}
LIBS        = -lm
DEPEND= makedepend

SRC     = tps.c utils.c read.c write.c integrate.c force.c init.c
OBJS    = tps.o utils.o read.o write.o integrate.o force.o init.o
EXECS   = tps.x

default: tps.x

all: $(EXECS)

tps.x:$(OBJS)
	$(CLINKER) $(OPTFLAGS) -o tps.x $(OBJS) $(LIBS)

clean:
	/bin/rm -f *.o *~ $(EXECS)

.c.o:
	$(CC) $(CFLAGS) -c $*.c

tps.o: tps.h
utils.o: tps.h
read.o: tps.h
write.o: tps.h
integrate.o: tps.h
force.o: tps.h
init.o: tps.h
