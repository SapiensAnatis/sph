SOURCES = $(wildcard ./src/*.cpp)
OBJECTS = $(SOURCES:.c=.o)
BINARYNAME = sph

HEADERS = $(wildcard ./src/*.hpp)

TEST_SOURCES = $(wildcard ./unittest/*.cpp)
TEST_OBJECTS = $(TEST_SOURCES:.c=.o)
TEST_BINARYNAME = sph-test

default: main

%.o: %.c $(HEADERS)
	g++ -c -g -o $@ $<

main: $(OBJECTS)
	g++ $(OBJECTS) -o $(BINARYNAME) -Wall -Wextra -Og

test: $(TEST_OBJECTS)
	g++ $(TEST_OBJECTS) -o $(TEST_BINARYNAME) -Wall -Wextra -Og

clean:
	rm -f *.o
	rm -f $(BINARYNAME)
	rm -f $(TEST_BINARYNAME)