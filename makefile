CCC=g++
CC=gcc
CFLAGS=-O3 #-g -Wall
LDFLAGS=-lglpk -lgmp
all: scabble

scabcc: scabble.cc word.cc rational.cc lp.cc
	$(CCC) $(CFLAGS) -c scabble.cc word.cc rational.cc lp.cc

scabc: matrix.c
	$(CC) $(CFLAGS) -c matrix.c

scabble: scabcc scabc
	cd exlp-package; make
	$(CCC) $(CFLAGS) -o scabble scabble.o word.o rational.o matrix.o lp.o exlp-package/*.o $(LDFLAGS)


clean: 
	rm scabble
	rm *.o
	cd exlp-package; rm *.o
