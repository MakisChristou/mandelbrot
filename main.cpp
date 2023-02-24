#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <SFML/Graphics.hpp>
#include <omp.h>

#include "mandelbrot.h"

int main(int argc, char* argv[])
{
    int image_width = 1000;
    int image_height = 1000; 

    long double output_start = -2.0f;
    long double output_end = 2.0f;
    long double factor = 1.0f;

    int n_max = 64; // 4096
    int s_max = 4; // prefer to be a power of 2

    Mandelbrot mandelbrot = Mandelbrot{image_width, image_height, output_start, output_end, factor, n_max, s_max};
    
    // for(int i = 0; i < 10; i++)
    // {
    //     mandelbrot.zoom();
    // }

    mandelbrot.render();
    // mandelbrot.writePPM();    
    
    


    return 0;
}