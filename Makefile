CC = g++
CFLAGS = -c -std=c++11 -Ofast -march=native

all: ran

ran: main.o particle.o particlesystem.o vector3d.o
	$(CC) main.o particle.o particlesystem.o vector3d.o -fopenmp -o ran

main.o : main.cpp
	$(CC) $(CFLAGS) main.cpp

particle.o: particle.cpp
	$(CC) $(CFLAGS) particle.cpp

particlesystem.o: particlesystem.cpp
	$(CC) $(CFLAGS) particlesystem.cpp

vector3d.o: vector3d.cpp
	$(CC) $(CFLAGS) vector3d.cpp

clean: 
	rm -rf *.o ran