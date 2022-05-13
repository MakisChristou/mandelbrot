#include "gpu.h"

#include <stdio.h>


// Convert pixel coords to complex coords
__device__ inline long double map(long double input, long double output_start, long double output_end, long double input_start, long double input_end)
{
    return (output_start + ((output_end - output_start) / (input_end - input_start)) * (input - input_start));
}

// Stolen code from https://github.com/sevity/mandelbrot
__device__ inline Color linearInterpolation(const Color& v, const Color& u, double a)
{
	auto const b = 1 - a;

    Color color;

    color.R = b * v.R + a * u.R;
    color.G = b*v.G + a * u.G;
    color.B = b*v.B + a * u.B;

	return color;
}

__device__ inline cIterations iterateMandelbrot(long double i, long double j, long double output_start, long double output_end, int image_width, int image_height, int n_max)
{

    // printf("i: %lf\n", i);
    // printf("j: %lf\n", j);
    

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

    // printf("complex_i: %lf\n", complex_i);
    // printf("complex_j: %lf\n", complex_j);
    

    cIterations citerations;
    citerations.n = n;
    citerations.c.Re = x;
    citerations.c.Re = y;

    return citerations;
}

__device__ inline Color getColor(int iter, Color* colorPallete, int palleteSize, int n_max)
{
    // Stolen Code from https://github.com/sevity/mandelbrot
    

    size_t max_color = palleteSize - 1;

    // if (iter == n_max)iter = 0;
    double mu = 1.0*iter / n_max;


    //scale mu to be in the range of colors
    mu *= max_color;
    size_t i_mu = static_cast<size_t>(mu);

    Color color1 = colorPallete[i_mu];
    Color color2;

    

    if(i_mu + 1 < max_color)
    {
        color2 = colorPallete[i_mu + 1];
    }
    else
    {
        color2 = colorPallete[max_color];
    }

    Color c = linearInterpolation(color1, color2, mu - i_mu);

    if(iter == n_max)
    {
        c.R = c.G = c.B = 0;
    }

    return c;
}

// Kernel
__global__ void mandelbortKernel(Color *pixelColours, Color* colorPallete, int palleteSize, int image_height, int image_width, double* output_start, double* output_end, int n_max, int s_max)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x; 
    int j = blockIdx.y * blockDim.y + threadIdx.y;


    if (i < image_width && j < image_height)
    {
        cIterations citerations;
        Complex c;
        int n = 0;
        int sum = 0;


        for(double k = 0.0; k < 1.0; k+=1.0/s_max)
        {
            double ii  = i+k;
            double jj = j+k;

            // Pass value of bounds
            citerations = iterateMandelbrot(ii,jj, *output_start, *output_end, image_width, image_height, n_max);

            n = citerations.n;
            c = citerations.c;
            
            sum+=n;
        }

        sum = sum / s_max;
        n = sum;
        
        // bug here, cannot access colorPallete (and probably not other arrays either)
        Color color = getColor(n, colorPallete, palleteSize, n_max);

        // printf("%d%d%d\n",color.R,color.G,color.B);

        // printf("*output_start: %lf\n", *output_start);
        // printf("*output_end: %lf\n", *output_end);
        // printf("n_max: %d\n", n_max);
        // printf("s_max: %d\n", s_max);
        // printf("image_width: %d\n", image_width);
        // printf("image_height: %d\n", image_height);

        // printf("Hello1 thread row=%d, col=%d, n = %d \n", i, j, n);

        // iterationCounts[i * image_height + j] = n;
        pixelColours[i * image_height + j] = color;
    } 
}


Color* gpuAllocColor(int N, int M, int bytes)
{
    Color *d_B;
    cudaMalloc(&d_B, bytes);
    return d_B;
}

double* gpuAllocDouble(int N, int M, int bytes)
{
    double* d_output_start;
    cudaMalloc(&d_output_start, bytes);
    return d_output_start;
}

void gpuFree(Color* d_B, Color* d_P, double* d_output_start, double* d_output_end)
{
        cudaFree(d_B);
        cudaFree(d_P);
        cudaFree(d_output_start);
        cudaFree(d_output_end);
}

void gpuCopyToDevice(int N, int M, int palleteSize, Color* d_B, Color* d_P, Color* B, Color* P, double* d_output_start, double* d_output_end, double* output_start_host, double* output_end_host)
{
    size_t pixelBytes = N*M*sizeof(Color);
    size_t palleteBytes = palleteSize*sizeof(Color);

    cudaMemcpy(d_B, B, pixelBytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_P, P, palleteBytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_output_start, output_start_host, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_output_end, output_end_host, sizeof(double), cudaMemcpyHostToDevice);

}

void gpuCopyFromDevice(int N, int M, int palleteSize, Color* d_B, Color* d_P, Color* B, Color* P)
{   
    size_t pixelBytes = N*M*sizeof(Color);
    size_t palleteBytes = palleteSize*sizeof(Color);


    cudaMemcpy(B, d_B, pixelBytes, cudaMemcpyDeviceToHost);
}

void gpuRender(Color* d_B, Color* d_P, int palleteSize, int N, int M, double* d_output_start, double* d_output_end, int n_max, int s_max)
{
    int thr_per_blk = 16;
    int blk_in_grid = (N + thr_per_blk -1 )/ thr_per_blk ;

    dim3 threads(thr_per_blk, thr_per_blk);
    dim3 blocks(blk_in_grid, blk_in_grid);

    mandelbortKernel<<< blocks, threads >>>(d_B, d_P, palleteSize, N, M, d_output_start, d_output_end, n_max, s_max);

    cudaThreadSynchronize();
}