# lib name
PROGRAM=lattice2D

# sources for lib
SOURCES= Chain2D.cpp Lattice2D.cpp Stats.cpp main.cpp Simulation2D.cpp\
			Native2D.cpp AA.cpp Pos2D.cpp Design2D.cpp 
OBJECTS=$(SOURCES:.cpp=.o)

# non-exported headers
HEADERS=Lattice2D.hpp Chain2D.hpp Residue2D.hpp Stats.hpp \
	Native2D.hpp AA.hpp Pos2D.hpp StatsMethods2D.hpp Design2D.hpp

CFLAGS=-Wall -O3

CC=g++

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS) 
	$(CC) -o $@ $(OBJECTS)

%.o:%.cpp
	$(CC) $(CFLAGS) -c $<


clean:
	rm $(OBJECTS) $(PROGRAM)


tidy:
	rm *~

run:
	./$(PROGRAM)
