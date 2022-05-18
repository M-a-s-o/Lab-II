# -*- MakeFile -*-
CPPC = g++
ROOTFLAGS = `root-config --cflags --glibs`
CPPDEPS = main.cpp functions.cpp
CPPOBJ = main.o functions.o
OBJFLAGS = -c -o

all: main

main: $(CPPOBJ)
	$(CPPC) $^ $(ROOTFLAGS) -o $@

%.o: %.cpp
	$(CPPC) $(OBJFLAGS) $@ $< $(ROOTFLAGS) -O3

.PHONY: all clean

clean:
	rm *.o