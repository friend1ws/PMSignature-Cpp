CC=g++
CFLAGS=-g -Wall

PMSignature: pmsig.o pmsig_independent.o pmsig_full.o main.o 
	$(CC) $(CFLAGS) -o PMSignature pmsig.o pmsig_independent.o pmsig_full.o main.o

pmsig.o: pmsig.cpp
	$(CC) $(CFLAGS) -c pmsig.cpp

pmsig_independent.o: pmsig_independent.cpp
	$(CC) $(CFLAGS) -c pmsig_independent.cpp

pmsig_full.o: pmsig_full.cpp
	$(CC) $(CFLAGS) -c pmsig_full.cpp

main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp

