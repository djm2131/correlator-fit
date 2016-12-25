#!/bin/bash

g++ -c ../src/correlator.cc -g -Wall -O3 -lgsl -lm -lgslcblas -std=c++11 -o ../src/correlator.o -I../include -I/usr/include -I/usr/include/libxml2
g++ -c ../src/fitter.cc -g -Wall -O3 -lgsl -lm -lgslcblas -std=c++11 -o ../src/fitter.o -I../include -I/usr/include -I/usr/include/libxml2
g++ -c ../src/main.cc -g -Wall -O3 -lgsl -lm -lgslcblas -std=c++11 -o ../src/main.o -I../include -I/usr/include -I/usr/include/libxml2
g++ ../src/correlator.o ../src/fitter.o ../src/main.o -o ../build/correlator_fit -L/usr/local/lib -lxml2 -lz -llzma -lm -ldl -lgsl -lgslcblas
