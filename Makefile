SOURCES = $(wildcard *.cpp)
HEADERS = $(wildcard *.hpp)
OBJECTS = $(SOURCES:.c=.o)
BINARYNAME = sph


default: main

%.o: %.c $(HEADERS)
	g++ -c -g -o $@ $<

main: $(OBJECTS)
	echo $(SOURCES)
	g++ $(OBJECTS) -o $(BINARYNAME) -Wall -Wpedantic -lm -Og

clean:
	rm -f *.o
	rm -f $(BINARYNAME)