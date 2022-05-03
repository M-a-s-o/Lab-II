# -*- MakeFile -*-
CPPC = g++
ROOTFLAGS = `root-config --cflags --glibs`
CPPDEPS = functions.hpp
CPPOBJ = main.o
OBJFLAGS = -c -o

all: main $(CPPOBJ)

main: $(CPPOBJ) $(CPPDEPS)
	$(CPPC) $^ $(ROOTFLAGS) -o $@

main.o: main.cpp
	$(CPPC) $(OBJFLAGS) $@ $< $(ROOTFLAGS) -O3

.PHONY: all clean

clean:
	rm *.o