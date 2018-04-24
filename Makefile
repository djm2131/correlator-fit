TARGET := ./build/correlator-fit.x
LIBS := -lm -lgomp -lxml2 -lgsl -lgslcblas
CXX := g++
CXXFLAGS := -g -Wall -O3 -I./include -I/usr/include/libxml2 -fopenmp -std=c++11

.PHONY: default all clean

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.cc, %.o, $(wildcard ./src/*.cc))
HEADERS = $(wildcard ./include/*.h)

%.o: %.c $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -Wall -o $@ $(LIBS)

clean:
	-rm -f ./src/*.o
	-rm -f $(TARGET)
