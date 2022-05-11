#include <iostream>
#include <vector>
#include <cstdlib>
#include <bits/stdc++.h>
#include <math.h>

// Global Declarations
int image_width = 5000;
int image_height = 5000; 

// long double output_start = -2.0f;
// long double output_end = 2.0f;

// long double output_start = 0.2f;
// long double output_end = 0.5f;

long double output_start = 0.35f;
long double output_end = 0.36f;



int n_max = 512; // 4096

int s_max = 8; // prefer to be a power of 2



typedef struct Color
{
    int R;
    int G;
    int B;

    Color()
    {
        R = 0;
        G = 0;
        B = 0;
    }

    Color(int r, int g, int b)
    {
        R = r;
        G = g;
        B = b;
    }
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

// Stolen code from https://github.com/sevity/mandelbrot
inline Color linearInterpolation(const Color& v, const Color& u, double a)
{
	auto const b = 1 - a;
	return Color(b * v.R + a * u.R, b*v.G + a * u.G, b*v.B + a * u.B);
}

// Prints PPM in std
void writePPM(std::vector<std::vector<int>> IterationCounts, std::vector<std::vector<Complex>> IterationValues, std::vector<Color> colors)
{
    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    for(int j = image_height-1; j >=0; --j)
    {
        std::cerr << "\rScanlines remaining: " << j << " " <<  std::flush;

        for(int i = 0; i < image_width; ++i)
        {
            int iter = IterationCounts[i][j];
            
            // Stolen Code from https://github.com/sevity/mandelbrot
            static const auto max_color = colors.size() - 1;
			// if (iter == n_max)iter = 0;
			double mu = 1.0*iter / n_max;
			//scale mu to be in the range of colors
			mu *= max_color;
			auto i_mu = static_cast<size_t>(mu);
			auto color1 = colors[i_mu];
			auto color2 = colors[std::min(i_mu + 1, max_color)];
			Color c = linearInterpolation(color1, color2, mu - i_mu);

            if(iter == n_max)
            {
                std::cout << 0 << " " << 0 << " " << 0 << "\n";
            }
            else
            {
                std::cout << c.R << " " << c.G << " " << c.B << "\n";
            }
        }

    }


    // return;

    // for(int j = image_height-1; j >=0; --j)
    // {    
    //     for(int i = 0; i < image_width; ++i)
    //     {

    //         int n = IterationCounts[i][j];
    //         Complex c = IterationValues[i][j];

    //         // Map iterations to brightness
    //         // int color = map(n, 0, n_max, 0, 255);

    //         // long double color = smoothColor(n,c);

    //         Color color = colorPalete[n];
    
    //         // if(n == n_max || n < 20)
    //         if(n == n_max)
    //         {
    //            std::cout << 0 << " " << 0 << " " << 0 << std::endl;
    //         }
    //         else
    //         {
    //             std::cout << color.R << " " << color.G << " " << color.B << "\n";
    //         }
    //     }

    // }
}


cIterations iterateMandelbrot(long double i, long double j)
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


std::vector<Color> generateColorPalete()
{
    std::vector<Color> colorPalete;

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

            // Calculate log step for current iteration
            double long rStep = log2(1 + (color1.R * log2(i)) / 10) / 8;
            double long gStep = log2(1 + (color1.G * log2(i)) / 10) / 8;
            double long bStep = log2(1 + (color1.B * log2(i)) / 10) / 8;

            // 
            sumR += rStep;
            sumG += gStep;
            sumB += bStep;

            tempColor.R = sumR;
            tempColor.G = sumG;
            tempColor.B = sumB;

            colorPalete.push_back(tempColor);
        }
    }

    return colorPalete;
}


int main(int argc, char* argv[])
{
    std::vector<std::vector<int>> IterationCounts; // Number of iterations per pixel
    std::vector<std::vector<Complex>> IterationValues; // Value of function after n iterations
    std::vector<int> NumPixelsPerItteration(n_max, 0); // Position i stores number of pixels that have i iterations


    // std::vector<Color> colorPalete;
    // colorPalete = generateColorPalete();

    // Stolen color palette
    std::vector<Color> colors{
        {0,7,100},
        {32,107,203},
        {237,255,255},
        {255,170,0},
        {0,2,0},
    };


    // Preallocation
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

            cIterations citerations;
            Complex c;
            int n = 0;

            int sum = 0;


            for(double k = 0.0; k < 1.0; k+=1.0/s_max)
            {
                // std::cerr << k << std::endl;

                double ii  = i+k;
                double jj = j+k;

                // std::cerr << ii << " " << jj << std::endl;

                citerations = iterateMandelbrot(ii,jj);

                n = citerations.n;
                c = citerations.c;
                
                sum+=n;
            }

            sum = sum / s_max;

            n = sum;

            // return 0;



            // citerations = iterateMandelbrot(i,j);

            // n = citerations.n;
            // c = citerations.c;

            IterationCounts[i].push_back(n);
            IterationValues[i].push_back(c);
        }
    }

    writePPM(IterationCounts, IterationValues, colors);

    return 0;
}