#------------------------------------------------------------------------------
SOURCE=./src/main.cpp
MYPROGRAM=mandelbrot
LIBS=
CC=g++

#------------------------------------------------------------------------------

all: $(MYPROGRAM)

$(MYPROGRAM): $(SOURCE)
	$(CC) $(SOURCE) -o $(MYPROGRAM) $(LIBS)

clean:
	rm -f $(MYPROGRAM) *.ppm
