#include <iostream>
#include <vector>
#include <cstdlib>
#include <bits/stdc++.h>
#include <math.h>

// Global Declarations
int image_width = 1000;
int image_height = 1000; 

long double output_start = -2.0f;
long double output_end = 2.0f;

// long double output_start = 0.2f;
// long double output_end = 0.5f;

// long double output_start = 0.2f;
// long double output_end = 0.5f;


int n_max = 4096; // 4096

int s_max = 1; // anti-aliasing



typedef struct Color
{
    int R;
    int G;
    int B;
};

typedef struct Complex
{
    long double Re;
    long double Im;
};


typedef struct cIterations
{
    Complex c;
    int n;
};


// Convert pixel coords to complex coords
long double map(long double input, long double output_start, long double output_end, long double input_start, long double input_end)
{
    return (output_start + ((output_end - output_start) / (input_end - input_start)) * (input - input_start));
}

// Returns smooth colour based on iteration and C value when escape
long double smoothColor(int n, Complex c)
{
    long double Zr = c.Re;
    long double Zi = c.Im; 
    long double N = (long double) n;

    // std::cerr << "Zr: " << Zr << " Zi: " << Zi << " N: " << N << std::endl;

    return 1.0 + N - log(log(sqrt(Zr*Zr + Zi*Zi)));
}

// Prints PPM in std
void writePPM(std::vector<std::vector<int>> IterationCounts, std::vector<std::vector<Complex>> IterationValues, std::vector<Color> colorPalete)
{
    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";


    for(int j = image_height-1; j >=0; --j)
    {    
        for(int i = 0; i < image_width; ++i)
        {

            int n = IterationCounts[i][j];
            Complex c = IterationValues[i][j];

            // Map iterations to brightness
            // int color = map(n, 0, n_max, 0, 255);

            // long double color = smoothColor(n,c);

            Color color = colorPalete[n];
    
            // if(n == n_max || n < 20)
            if(n == n_max)
            {
               std::cout << 0 << " " << 0 << " " << 0 << std::endl;
            }
            else
            {
                std::cout << color.R << " " << color.G << " " << color.B << "\n";
            }
        }

    }
}


cIterations iterateMandelbrot(int i, int j)
{
    long double complex_i = map(i, output_start, output_end, 0, image_width);
    long double complex_j = map(j, output_start, output_end, 0, image_height);

    long double x0 = 0;
    long double y0 = 0;
    
    long double x = 0;
    long double y = 0;

    int n = 0; // iterations

    while(x*x + y*y <= 4 && n < n_max)
    {
        x = x0*x0 - y0*y0 + complex_i;
        y = 2*x0*y0 + complex_j;

        x0 = x;
        y0 = y;
        n++;
    }

    cIterations citerations;
    citerations.n = n;
    citerations.c.Re = x;
    citerations.c.Re = y;

    return citerations;
}



int main(int argc, char* argv[])
{
    std::vector<std::vector<int>> IterationCounts; // Number of iterations per pixel
    std::vector<std::vector<Complex>> IterationValues; // Value of function after n iterations
    std::vector<int> NumPixelsPerItteration(n_max, 0); // Position i stores number of pixels that have i iterations

    std::vector<Color> colorPalete;

    // Color 1 #0048BA // R = 0x00, G = 0x48, B = 0xBA
    // Color 2 #BA7300 // R = 0xBA, G = 0x73, B = 0x00

    Color color1;
    color1.R = 0xFF;
    color1.G = 0xFF;
    color1.B = 0xFF;

    Color color2;
    color2.R = 0x00;
    color2.G = 0x00;
    color2.B = 0x00;

    Color black;
    black.R = 0x00;
    black.G = 0x00;
    black.B = 0x00;

    // double rStep = -1.0 * (double) (color1.R - color2.R) / n_max;
    // double gStep = -1.0 * (double) (color1.G - color2.G) / n_max;
    // double bStep = -1.0 * (double) (color1.B - color2.B) / n_max;

    long double sumR = 0;
    long double sumG = 0;
    long double sumB = 0;



    for(int i = 0; i < n_max; i++)
    {
        if(i == 0)
        {
            colorPalete.push_back(color1);
        }
        else if(i == n_max-1)
        {
            colorPalete.push_back(color2);
        }
        else if(i == n_max)
        {
            colorPalete.push_back(black);
        }
        else
        {
            Color tempColor;

            // Calculate log step
            double long rStep = log2(1 + (color1.R * log2(i)) / 10) / 8;
            double long gStep = log2(1 + (color1.G * log2(i)) / 10) / 8;
            double long bStep = log2(1 + (color1.B * log2(i)) / 10) / 8;

            sumR += rStep;
            sumG += gStep;
            sumB += bStep;

            tempColor.R = sumR;
            tempColor.G = sumG;
            tempColor.B = sumB;


            colorPalete.push_back(tempColor);

            // std::cerr << tempColor.R << " " << tempColor.G << " " << tempColor.B << "\n";
        }
    }



    IterationCounts.reserve(image_width);
    IterationValues.reserve(image_width);

    srand(time(0));

    for(int i = 0; i < image_width; i++)
    {
        std::cerr << "\rScanlines remaining: " << image_width - i << " " <<  std::flush;

        IterationCounts.push_back(std::vector<int>());
        IterationValues.push_back(std::vector<Complex>());

        for(int j = 0; j < image_height; j++)
        {

            int sum_n = 0;
            int ii = i;
            int jj = j;

            cIterations citerations;
            Complex c;
            int n = 0;

            for(int k = 0; k < s_max; k++)
            {
                long double r = (double)(rand()%100)/1000;

                ii+=r;
                jj+=r;

                citerations = iterateMandelbrot(ii,jj);

                n = citerations.n;
                c = citerations.c;

                sum_n = sum_n + citerations.n;
            }

            sum_n = sum_n / s_max;


            // int n = iterateMandelbrot(i,j);

            IterationCounts[i].push_back(sum_n);
            IterationValues[i].push_back(c);

            
        }
    }

    writePPM(IterationCounts, IterationValues, colorPalete);

    return 0;
}