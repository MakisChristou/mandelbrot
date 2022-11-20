#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <SFML/Graphics.hpp>

class Color
{   
    public:
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

class Complex
{
    public:
        long double Re;
        long double Im;

        Complex()
        {

        }

        Complex(long double Re, long double Im)
        {
            this->Re = Re;
            this->Im = Im;
        }
};

class cIterations
{
    public:
        Complex c;
        int n;

        cIterations()
        {

        }

        cIterations(int n, Complex c)
        {
            this->n = n;
            this->c = c;
        }

};

class Mandelbrot
{
    private:
        int image_width = 500;
        int image_height = 500; 

        long double output_start = -2.0f;
        long double output_end = 2.0f;

        long double factor = 1.0f;

        int n_max = 64; // 4096
        int s_max = 4; // prefer to be a power of 2

        std::vector<std::vector<int>> IterationCounts; // Number of iterations per pixel
        std::vector<std::vector<Complex>> IterationValues; // Value of function after n iterations
        std::vector<std::vector<Color>> pixelColours;
        
        // Stolen color Pallete
        std::vector<Color> colorPallete{
            {0,7,100},
            {32,107,203},
            {237,255,255},
            {255,170,0},
            {0,2,0},
        };

        // Convert pixel coords to complex coords
        long double map(long double input, long double output_start, long double output_end, long double input_start, long double input_end)
        {
            return (output_start + ((output_end - output_start) / (input_end - input_start)) * (input - input_start));
        }

        // Returns smooth colour based on iteration and C value when escape
        long double smoothColor(int n, const Complex& c)
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

        Color getColor(int iter, const std::vector<Color>& colorPallete)
        {
            // Stolen Code from https://github.com/sevity/mandelbrot
            static const auto max_color = colorPallete.size() - 1;
            // if (iter == n_max)iter = 0;
            double mu = 1.0*iter / n_max;
            //scale mu to be in the range of colors
            mu *= max_color;
            auto i_mu = static_cast<size_t>(mu);
            auto color1 = colorPallete[i_mu];
            auto color2 = colorPallete[std::min(i_mu + 1, max_color)];
            Color c = linearInterpolation(color1, color2, mu - i_mu);

            if(iter == n_max)
            {
                c.R = c.G = c.B = 0;
            }

            return c;
        }

        cIterations iterateMandelbrot(const long double& i,const long double& j)
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

            return cIterations{n, Complex{x,y}};
        }

        std::vector<Color> generateColorPalete()
        {
            std::vector<Color> colorPalete;

            Color color1{0xFF, 0xFF, 0xFF};

            Color color2{0x00, 0x00, 0x00};

            Color black{0x00, 0x00, 0x00};

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
                    // Calculate log step for current iteration
                    double long rStep = log2(1 + (color1.R * log2(i)) / 10) / 8;
                    double long gStep = log2(1 + (color1.G * log2(i)) / 10) / 8;
                    double long bStep = log2(1 + (color1.B * log2(i)) / 10) / 8;

                    sumR += rStep;
                    sumG += gStep;
                    sumB += bStep;

                    Color tempColor{sumR, sumG, sumB};
                    colorPalete.emplace_back(tempColor);
                }
            }

            return colorPalete;
        }

    public:

        Mandelbrot()
        {
            IterationCounts.reserve(image_width);
            IterationValues.reserve(image_width);
        }

        Mandelbrot(int image_width, int image_height, long double output_start, long double output_end, long double factor, int n_max, int s_max)
        {
            IterationCounts.reserve(image_width);
            IterationValues.reserve(image_width);


            this->image_width = image_width;
            this->image_height = image_height;

            this->output_start = output_start;
            this->output_end = output_end;

            this->factor = factor;

            this->n_max = n_max;
            this->s_max = s_max;
        }

        void writePPM(const std::vector<std::vector<Color>>& pixelColors)
        {
            std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

            for(int j = image_height-1; j >=0; --j)
            {
                for(int i = 0; i < image_width; ++i)
                {
                    Color c = pixelColors[i][j];
                    std::cout << c.R << " " << c.G << " " << c.B << "\n";
                }
            }
        }

        void writePPM()
        {
            std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

            for(int j = image_height-1; j >=0; --j)
            {
                for(int i = 0; i < image_width; ++i)
                {
                    Color c = pixelColours[i][j];
                    std::cout << c.R << " " << c.G << " " << c.B << "\n";
                }
            }
        }

        void render()
        {
            srand(time(0));

            // Main Loop
            for(int i = 0; i < image_width; i++)
            {
                // std::cerr << "\rScanlines remaining: " << image_width - i << " " <<  std::flush;

                IterationCounts.emplace_back(std::vector<int>());
                IterationValues.emplace_back(std::vector<Complex>());
                pixelColours.emplace_back(std::vector<Color>());

                for(int j = 0; j < image_height; j++)
                {
                    cIterations citerations;
                    Complex c;
                    int n = 0;
                    int sum = 0;


                    for(double k = 0.0; k < 1.0; k+=1.0/s_max)
                    {
                        double ii  = i+k;
                        double jj = j+k;

                        citerations = iterateMandelbrot(ii,jj);

                        n = citerations.n;
                        c = citerations.c;
                        
                        sum+=n;
                    }

                    sum = sum / s_max;
                    n = sum;

                    Color color = getColor(n, colorPallete);

                    IterationCounts[i].emplace_back(n);
                    IterationValues[i].emplace_back(c);
                    pixelColours[i].emplace_back(color);
                }
            }
        }
};

int main(int argc, char* argv[])
{
    int image_width = 500;
    int image_height = 500; 

    long double output_start = -2.0f;
    long double output_end = 2.0f;

    long double factor = 1.0f;

    int n_max = 64; // 4096
    int s_max = 4; // prefer to be a power of 2

    Mandelbrot mandelbrot = Mandelbrot{image_width, image_height, output_start, output_end, factor, n_max, s_max};
    mandelbrot.render();

    return 0;
}