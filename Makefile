#------------------------------------------------------------------------------
SOURCE=main.cpp
MYPROGRAM=mb
LIBS=-lsfml-graphics -lsfml-window -lsfml-system -fopenmp
CC=g++
FLAGS=-O3 -mavx2 -march=native -fopt-info-vec
DEBUG_FLAGS= -g -O0 -Wall

#------------------------------------------------------------------------------

all: $(MYPROGRAM)

$(MYPROGRAM): main.cpp mandelbrot.h color.h cIterations.h
	$(CC) $(FLAGS) main.cpp mandelbrot.h color.h cIterations.h -o $(MYPROGRAM) $(LIBS)

clean:
	rm -f $(MYPROGRAM) *.ppm
