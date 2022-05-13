#------------------------------------------------------------------------------
SOURCE=./src/main.cu
MYPROGRAM=mandelbrot
LIBS= 
CC=nvcc


#------------------------------------------------------------------------------

all: $(MYPROGRAM)

$(MYPROGRAM): $(SOURCE)
	$(CC) -O3 $(SOURCE) -o $(MYPROGRAM) $(LIBS)

clean:
	rm -f $(MYPROGRAM) *.ppm
