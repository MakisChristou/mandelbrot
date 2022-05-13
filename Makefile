#------------------------------------------------------------------------------
SOURCE=./src/main.cu
MYPROGRAM=cudelbrot
LIBS= 
CC=nvcc


#------------------------------------------------------------------------------

all: $(MYPROGRAM)

$(MYPROGRAM): $(SOURCE)
	$(CC) -O3 $(SOURCE) -o $(MYPROGRAM) $(LIBS)

clean:
	rm -f $(MYPROGRAM) *.ppm
