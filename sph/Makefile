OBJECTS := calculators.o kernel.o main.o setup.o smoothing_length.o sph_simulation.o ghost_particles.o

CXX := g++
CXXFLAGS := -std=c++17 -Wall -Wno-unknown-pragmas -Ofast
LDFLAGS := -lgsl -lgslcblas -lm -lstdc++fs

all: $(OBJECTS)
	${CXX} ${LDFLAGS} -o sph ${OBJECTS}

$(OBJECTS): %.o: %.cpp

clean:
	rm -f sph
	rm -f *.o
