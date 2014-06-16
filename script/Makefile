CC=g++
CFLAGS=-g -Wall -O2

PMSignature: pmsig.o main.o 
	$(CC) $(CFLAGS) -o PMSignature pmsig.o main.o

pmsig.o: pmsig.cpp
	$(CC) $(CFLAGS) -c pmsig.cpp

main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp

