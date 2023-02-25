#ifndef MANDELBROT_H
#define MANDELBROT_H

#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <immintrin.h>

#include "complex.h"
#include "cIterations.h"
#include "color.h"


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
        inline long double map(long double input, long double output_start, long double output_end, long double input_start, long double input_end)
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

        inline Color getColor(int iter, const std::vector<Color>& colorPallete)
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

        inline cIterations iterateMandelbrot(const long double& i,const long double& j)
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
        }

        Mandelbrot(int image_width, int image_height, long double output_start, long double output_end, long double factor, int n_max, int s_max)
        {
            this->image_width = image_width;
            this->image_height = image_height;

            this->output_start = output_start;
            this->output_end = output_end;

            this->factor = factor;

            this->n_max = n_max;
            this->s_max = s_max;

            for(int i = 0; i < image_width; i++)
            {
                pixelColours.emplace_back(std::vector<Color>());

                for(int j = 0; j < image_height; j++)
                {                    
                    pixelColours[i].emplace_back(Color{0,0,0});
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

        void vectorized_loop(double* x0_array, double* y0_array, long double* complex_i_array, long double* complex_j_array, int* n_max_array, int array_size, int* sum_array)
        {
            int i = 0;
            while (i < array_size)
            {
                int n = 0;
                __m256d x_vec = _mm256_set_pd(x0_array[i], x0_array[i+1], x0_array[i+2], x0_array[i+3]);
                __m256d y_vec = _mm256_set_pd(y0_array[i], y0_array[i+1], y0_array[i+2], y0_array[i+3]);
                __m256d complex_i_vec = _mm256_set1_pd(*complex_i_array);
                __m256d complex_j_vec = _mm256_set1_pd(*complex_j_array);
                __m256d threshold_vec = _mm256_set1_pd(4);
                __m256i n_max_vec = _mm256_set1_epi32(*n_max_array);

                while (true)
                {
                    __m256d x_squared = _mm256_mul_pd(x_vec, x_vec);
                    __m256d y_squared = _mm256_mul_pd(y_vec, y_vec);
                    __m256d sum_squared = _mm256_add_pd(x_squared, y_squared);
                    __m256d mask = _mm256_cmp_pd(sum_squared, threshold_vec, _CMP_LE_OQ);
                    __m256i n_vec = _mm256_set1_epi32(n);
                    __m256i n_max_mask = _mm256_cmpgt_epi32(n_max_vec, n_vec);
                    mask = _mm256_and_pd(mask, _mm256_castsi256_pd(n_max_mask));

                    int active_mask = _mm256_movemask_pd(mask);
                    if (active_mask == 0)
                    {
                        break;
                    }

                    n += 4 - active_mask;

                    __m256d x0_vec = x_vec;
                    __m256d y0_vec = y_vec;

                    x_vec = _mm256_sub_pd(_mm256_add_pd(_mm256_mul_pd(x_vec, x_vec), _mm256_mul_pd(y_vec, y_vec)), _mm256_sub_pd(complex_i_vec, complex_i_vec));
                    y_vec = _mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(_mm256_set1_pd(2.0), x0_vec), y0_vec), complex_j_vec);

                    x0_vec = _mm256_blendv_pd(x0_vec, x_vec, mask);
                    y0_vec = _mm256_blendv_pd(y0_vec, y_vec, mask);

                    _mm256_store_pd(&x0_array[i], x0_vec);
                    _mm256_store_pd(&y0_array[i], y0_vec);
                }

                sum_array[i] += n;
                i += 4;
            }
        }

        void render()
        {   
            for(int i = 0; i < image_width; i++)
            {
                for(int j = 0; j < image_height; j++)
                {
                    uint64_t sum = 0;

                    for(double k = 0.0; k < 1.0; k+=1.0/s_max)
                    {
                        double ii  = i+k;
                        double jj = j+k;

                        long double complex_i = map(ii, output_start, output_end, 0, image_width);
                        long double complex_j = map(jj, output_start, output_end, 0, image_height);
                        
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
                        sum+=n;                    
                    }

                    sum/=s_max;
                    pixelColours[i][j] = getColor(sum, colorPallete);
                }
            }
        }

        void zoom()
        {
            // Automatic Zoom 
            output_start+=0.15*factor;
            output_end-=0.1*factor;
            factor *= 0.9550;
        }
};

#endif