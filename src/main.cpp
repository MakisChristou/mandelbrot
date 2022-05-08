#include <iostream>
#include <vector>

// Global Declarations
int image_width = 1000;
int image_height = 1000;

long double output_start = -2.0f;
long double output_end = 2.0f;

int n_max = 1024;

// Convert pixel coords to complex coords
long double map(long double input, long double output_start, long double output_end, long double input_start, long double input_end)
{
    return (output_start + ((output_end - output_start) / (input_end - input_start)) * (input - input_start));
}

// Prints PPM in std
void writePPM(std::vector<std::vector<int>> IterationCounts)
{
    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    for(int j = image_height-1; j >=0; --j)
    {    
        for(int i = 0; i < image_width; ++i)
        {

            int n = IterationCounts[i][j];

            // Map iterations to brightness
            int bright = map(n, 0, n_max, 0, 255);
    
            if(n == n_max || n < 20)
            {
               std::cout << 0 << " " << 0 << " " << 0 << std::endl;
            }
            else
            {
                std::cout << bright << " " << bright << " " << bright << "\n";
            }
        }

    }
}


int main(int argc, char* argv[])
{
    std::vector<std::vector<int>> IterationCounts;

    IterationCounts.reserve(image_width);


    for(int i = 0; i < image_width; i++)
    {

        IterationCounts.push_back(std::vector<int>());

        for(int j = 0; j < image_height; j++)
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

            IterationCounts[i].push_back(n);
        }
    }

    writePPM(IterationCounts);

    return 0;
}