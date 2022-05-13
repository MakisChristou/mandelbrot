#include <iostream>
#include <vector>
#include <cstdlib>
#include <math.h>
#include <stdio.h>
#include <cuda.h>
#include <chrono>


using namespace std::chrono;

// #include <SDL2/SDL.h>

// Global Declarations
int image_width = 5000;
int image_height = 5000; 



// long double output_start = 0.2f;
// long double output_end = 0.5f;

// long double output_start = 0.35f;
// long double output_end = 0.36f;

int n_max = 64; // 4096
int s_max = 8; // prefer to be a power of 2

typedef struct gpuColor
{
    int R;
    int G;
    int B;

    __device__ __host__ gpuColor()
    {

    }

    __device__ __host__ gpuColor(int r, int g, int b)
    {
        R = r;
        G = g;
        B = b;
    }
    __device__ __host__ gpuColor(double r, double g, double b)
    {
        R = (int)r;
        G = (int)g;
        B = (int)b;
    }
};

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
__device__ inline long double map(long double input, long double output_start, long double output_end, long double input_start, long double input_end)
{
    return (output_start + ((output_end - output_start) / (input_end - input_start)) * (input - input_start));
}

// Stolen code from https://github.com/sevity/mandelbrot
__device__ inline gpuColor linearInterpolation(const gpuColor& v, const gpuColor& u, double a)
{
	auto const b = 1 - a;
	return gpuColor(b * v.R + a * u.R, b*v.G + a * u.G, b*v.B + a * u.B);
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

__device__ inline gpuColor getColor(int iter, gpuColor* colorPallete, int palleteSize, int n_max)
{
    // Stolen Code from https://github.com/sevity/mandelbrot
    

    size_t max_color = palleteSize - 1;

    // if (iter == n_max)iter = 0;
    double mu = 1.0*iter / n_max;


    //scale mu to be in the range of colors
    mu *= max_color;
    size_t i_mu = static_cast<size_t>(mu);

    gpuColor color1 = colorPallete[i_mu];
    gpuColor color2;

    

    if(i_mu + 1 < max_color)
    {
        color2 = colorPallete[i_mu + 1];
    }
    else
    {
        color2 = colorPallete[max_color];
    }

    gpuColor c = linearInterpolation(color1, color2, mu - i_mu);

    if(iter == n_max)
    {
        c.R = c.G = c.B = 0;
    }

    return c;
}

// Kernel
__global__ void mandelbortKernel(gpuColor *pixelColours, gpuColor* colorPallete, int palleteSize, int image_height, int image_width, double* output_start, double* output_end, int n_max, int s_max)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x; 
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    // For some reason
    // double output_start = *output_startt;
    // double output_end = *output_endd;
    // double factor = *factorr;


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

            citerations = iterateMandelbrot(ii,jj, *output_start, *output_end, image_width, image_height, n_max);

            n = citerations.n;
            c = citerations.c;
            
            sum+=n;
        }

        sum = sum / s_max;
        n = sum;
        
        // bug here, cannot access colorPallete (and probably not other arrays either)
        gpuColor color = getColor(n, colorPallete, palleteSize, n_max);

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


// Returns smooth colour based on iteration and C value when escape
inline long double smoothColor(int n, Complex c)
{
    long double Zr = c.Re;
    long double Zi = c.Im; 
    long double N = (long double) n;

    // std::cerr << "Zr: " << Zr << " Zi: " << Zi << " N: " << N << std::endl;

    return 1.0 + N - log(log(sqrt(Zr*Zr + Zi*Zi)));
}

// Prints PPM in std
void writePPM(gpuColor* pixelColors)
{
    printf("P3\n");
    printf("%d %d",image_width,image_height);
    printf("\n255\n");

    for(int j = image_height-1; j >=0; --j)
    {
        for(int i = 0; i < image_width; ++i)
        {
            gpuColor c = pixelColors[i*image_height + j];
            printf("%d %d %d\n", c.R, c.G, c.B);
        }
    }
}


inline std::vector<Color> generateColorPalete()
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
    // SDL_Init(SDL_INIT_EVERYTHING);
    // SDL_Window* window;
    // SDL_Renderer* renderer;
    // SDL_Event event;
    // SDL_CreateWindowAndRenderer(image_width, image_height, 0, &window, &renderer);
    // SDL_RenderSetLogicalSize(renderer, image_width, image_height);


    int N = image_width;
    int M = image_height; 

    double output_start = 0.35f;
    double output_end = 0.36f;
    double factor = 1.0f;

     int palleteSize = 5;

    int *iterationCounts = (int*) malloc(N*M*sizeof(int));
    Color *pixelColours = (Color*) malloc(N*M*sizeof(Color));

    while(1)
    {
        // SDL_RenderPresent(renderer);
        // SDL_Event event;
        // while (SDL_PollEvent(&event))
        // {
        //     if (event.type == SDL_QUIT) {
        //         SDL_Quit();
        //         return 0;
        //     }
        // }
    



        // Number of bytes to allocate for N doubles
        // size_t iterationBytes = N*M*sizeof(int);
        size_t pixelBytes = N*M*sizeof(gpuColor);
        size_t palleteBytes = palleteSize*sizeof(gpuColor);

        // Allocate memory for arrays A, B, and C on host
        // int *A = (int*)malloc(iterationBytes); // iterationCounts
        gpuColor *B = (gpuColor*)malloc(pixelBytes); // pixelColours
        gpuColor *P = (gpuColor*)malloc(palleteBytes); // colorPalete
        double* output_start_host = (double*)malloc(sizeof(double));
        double* output_end_host = (double*)malloc(sizeof(double));
        double* factor_host = (double*)malloc(sizeof(double));


        // Allocate memory for arrays d_A, d_B, and d_C on device
        int *d_A;
        gpuColor *d_B;
        gpuColor *d_P;
        double* d_output_start;
        double* d_output_end;

        // cudaMalloc(&d_A, iterationBytes);
        cudaMalloc(&d_B, pixelBytes);
        cudaMalloc(&d_P, palleteBytes);
        cudaMalloc(&d_output_start, sizeof(double));
        cudaMalloc(&d_output_end, sizeof(double));

        // Fill host arrays data structures
        //  for(int i = 0; i < image_width; i++)
        // {
        //     for(int j = 0; j < image_height; j++)
        //     {
        //         B[i*image_height+j].R = 0;
        //         B[i*image_height+j].G = 0;
        //         B[i*image_height+j].B = 0;
        //         // A[i*image_height+j] = 0;
        //     }
        // }

        P[0] = gpuColor(0,7,100);
        P[1] = gpuColor(32,107,203);
        P[2] = gpuColor(237,255,255);
        P[3] = gpuColor(255,170,0);
        P[4] = gpuColor(0,2,0);

        // Fill host arrays data structures
        output_start_host[0] = output_start;
        output_end_host[0] = output_end;


        // Copy data from host arrays A and B to device arrays d_A and d_B
        // cudaMemcpy(d_A, A, iterationBytes, cudaMemcpyHostToDevice);
        cudaMemcpy(d_B, B, pixelBytes, cudaMemcpyHostToDevice);
        cudaMemcpy(d_P, P, palleteBytes, cudaMemcpyHostToDevice);
        cudaMemcpy(d_output_start, output_start_host, sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(d_output_end, output_end_host, sizeof(double), cudaMemcpyHostToDevice);
        

        // Set execution configuration parameters
        //      thr_per_blk: number of CUDA threads per grid block
        //      blk_in_grid: number of blocks in grid
        int thr_per_blk = 16;
        int blk_in_grid = (N + thr_per_blk -1 )/ thr_per_blk ;

        dim3 threads(thr_per_blk, thr_per_blk);
        dim3 blocks(blk_in_grid, blk_in_grid);

        auto start = high_resolution_clock::now();

        mandelbortKernel<<< blocks, threads >>>(d_B, d_P, palleteSize, N, N, d_output_start, d_output_end, n_max, s_max);

        cudaThreadSynchronize();

        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);


        std::cerr << "Mandel frame computed in: " << duration.count() << "us" << "\n";

        // Copy data from device array d_C to host array C
        cudaMemcpy(B, d_B, pixelBytes, cudaMemcpyDeviceToHost);


        writePPM(B);

        // Free CPU memory
        // free(A);
        free(B);
        free(P);
        free(output_start_host);
        free(output_end_host);

        // Free GPU memory
        // cudaFree(d_A);
        cudaFree(d_B);
        cudaFree(d_P);
        cudaFree(d_output_start);
        cudaFree(d_output_end);

        // printf("\n---------------------------\n");
        // printf("__SUCCESS__\n");
        // printf("---------------------------\n");
        // printf("N                 = %d\n", K);
        // printf("Threads Per Block = %d\n", thr_per_blk);
        // printf("Blocks In Grid    = %d\n", blk_in_grid);
        // printf("---------------------------\n\n");


        // printf("Rendered! \n");

        return 0;



        // // Main Loop
        // for(int i = 0; i < image_width; i++)
        // {
        //     // std::cerr << "\rScanlines remaining: " << image_width - i << " " <<  std::flush;

        //     for(int j = 0; j < image_height; j++)
        //     {

        //         cIterations citerations;
        //         Complex c;
        //         int n = 0;
        //         int sum = 0;
  


        //         for(double k = 0.0; k < 1.0; k+=1.0/s_max)
        //         {
        //             double ii  = i+k;
        //             double jj = j+k;

        //             citerations = iterateMandelbrot(ii,jj);

        //             n = citerations.n;
        //             c = citerations.c;
                    
        //             sum+=n;
        //         }

        //         sum = sum / s_max;
        //         n = sum;

        //         Color color = getColor(n, colorPallete, palleteSize);

        //         iterationCounts[i * image_height + j] = n;
        //         pixelColours[i * image_height + j] = color;
                

        //         // SDL Draw
        //         // SDL_SetRenderDrawColor(renderer, color.R, color.G, color.B, 255);
        //         // SDL_RenderDrawPoint(renderer, i, j);
        //     }
        // }

        // Zoom in code by https://www.youtube.com/watch?v=KnCNfBb2ODQ
        // output_start+=0.15*factor;
        // output_end-=0.1*factor;
        // factor *= 0.9349;
        // n_max+=5;

        

        return 0;

    }

    free(iterationCounts);
    free(pixelColours);

    return 0;
}