CXX=g++
NVCC=nvcc
MYPROGRAM=cundelbrot

all: ${MYPROGRAM}

cundelbrot: main.o gpu.o
	${CXX} -o ${MYPROGRAM} main.o gpu.o -lSDL2 -lglfw -lcudart

main.o: main.cpp
	${CXX} -c -o main.o main.cpp

gpu.o: gpu.cu gpu.h
	${NVCC} -c -o gpu.o gpu.cu

clean: 
	rm -f *.o ${MYPROGRAM} *.ppm
