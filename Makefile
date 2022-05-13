CXX=g++
NVCC=nvcc

all: main

main: main.o gpu.o
	${CXX} -o main main.o gpu.o -lSDL2 -lcudart

main.o: main.cpp
	${CXX} -c -o main.o main.cpp

gpu.o: gpu.cu gpu.h
	${NVCC} -c -o gpu.o gpu.cu

clean: 
	rm -f *.o main *.ppm
