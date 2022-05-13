#include <iostream>
#include <vector>
#include <cstdlib>
#include <math.h>
#include <stdio.h>
#include <cuda.h>
#include <chrono>

#include "gpu.h"

using namespace std::chrono;

#include <SDL2/SDL.h>

// Global Declarations
int image_width = 1000;
int image_height = 1000; 

// long double output_start = 0.2f;
// long double output_end = 0.5f;

// long double output_start = 0.35f;
// long double output_end = 0.36f;

int n_max = 64; // 4096
int s_max = 8; // prefer to be a power of 2



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
void writePPM(Color* pixelColors)
{
    printf("P3\n");
    printf("%d %d",image_width,image_height);
    printf("\n255\n");

    for(int j = image_height-1; j >=0; --j)
    {
        for(int i = 0; i < image_width; ++i)
        {
            Color c = pixelColors[i*image_height + j];
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
    SDL_Init(SDL_INIT_EVERYTHING);
    SDL_Window* window;
    SDL_Renderer* renderer;
    SDL_Event event;
    SDL_CreateWindowAndRenderer(image_width, image_height, 0, &window, &renderer);
    SDL_RenderSetLogicalSize(renderer, image_width, image_height);


    int N = image_width;
    int M = image_height; 

    double output_start = -2.0;
    double output_end = 2.0f;
    double factor = 1.0f;

     int palleteSize = 5;

    int *iterationCounts = (int*) malloc(N*M*sizeof(int));
    Color *pixelColours = (Color*) malloc(N*M*sizeof(Color));

    while(1)
    {
        SDL_RenderPresent(renderer);
        SDL_Event event;
        while (SDL_PollEvent(&event))
        {
            if (event.type == SDL_QUIT) {
                SDL_Quit();
                return 0;
            }
        }
    

        // Number of bytes to allocate for N doubles
        // size_t iterationBytes = N*M*sizeof(int);
        size_t pixelBytes = N*M*sizeof(Color);
        size_t palleteBytes = palleteSize*sizeof(Color);

        // Allocate memory for arrays A, B, and C on host
        // int *A = (int*)malloc(iterationBytes); // iterationCounts
        Color *B = (Color*)malloc(pixelBytes); // pixelColours
        Color *P = (Color*)malloc(palleteBytes); // colorPalete
        double* output_start_host = (double*)malloc(sizeof(double));
        double* output_end_host = (double*)malloc(sizeof(double));
        double* factor_host = (double*)malloc(sizeof(double));


        // Allocate memory for arrays d_A, d_B, and d_C on device
        Color* d_B;
        Color* d_P;
        double* d_output_start;
        double* d_output_end;

        

        // cudaMalloc(&d_A, iterationBytes);
        // cudaMalloc(&d_B, pixelBytes);
        // cudaMalloc(&d_P, palleteBytes);
        // cudaMalloc(&d_output_start, sizeof(double));
        // cudaMalloc(&d_output_end, sizeof(double));

        d_B = gpuAllocColor(N, M, pixelBytes);
        d_P = gpuAllocColor(N, M, palleteBytes);
        d_output_start = gpuAllocDouble(N, M, sizeof(double));
        d_output_end = gpuAllocDouble(N, M, sizeof(double));

        // Filldata structures
        
        Color color;
        color.R = 0;color.G = 7; color.B = 100;
        P[0] = color;

        color.R = 32;color.G = 107; color.B = 204;
        P[1] = color;

        color.R = 237;color.G = 255; color.B = 255;
        P[2] = color;

        color.R = 255;color.G = 170; color.B = 0;
        P[3] = color;

        color.R = 0;color.G = 2; color.B = 0;
        P[4] = color;

        // Fill host arrays data structures
        output_start_host[0] = output_start;
        output_end_host[0] = output_end;


        // Copy data from host arrays A and B to device arrays d_A and d_B
        // cudaMemcpy(d_A, A, iterationBytes, cudaMemcpyHostToDevice);
        // cudaMemcpy(d_B, B, pixelBytes, cudaMemcpyHostToDevice);
        // cudaMemcpy(d_P, P, palleteBytes, cudaMemcpyHostToDevice);
        // cudaMemcpy(d_output_start, output_start_host, sizeof(double), cudaMemcpyHostToDevice);
        // cudaMemcpy(d_output_end, output_end_host, sizeof(double), cudaMemcpyHostToDevice);


        gpuCopyToDevice(N,M,palleteSize,d_B,d_P,B,P,d_output_start,d_output_end,output_start_host,output_end_host);
        

        // Set execution configuration parameters
        //      thr_per_blk: number of CUDA threads per grid block
        //      blk_in_grid: number of blocks in grid
        // int thr_per_blk = 16;
        // int blk_in_grid = (N + thr_per_blk -1 )/ thr_per_blk ;

        // dim3 threads(thr_per_blk, thr_per_blk);
        // dim3 blocks(blk_in_grid, blk_in_grid);

        // auto start = high_resolution_clock::now();

        // mandelbortKernel<<< blocks, threads >>>(d_B, d_P, palleteSize, N, N, d_output_start, d_output_end, n_max, s_max);

        // cudaThreadSynchronize();

        // auto stop = high_resolution_clock::now();
        // auto duration = duration_cast<microseconds>(stop - start);

        gpuRender( d_B, d_P,palleteSize,  N,  M,  d_output_start, d_output_end, n_max, s_max);

        // std::cerr << "Mandel frame computed in: " << duration.count() << "us" << "\n";

        // Copy data from device array d_C to host array C
        // cudaMemcpy(B, d_B, pixelBytes, cudaMemcpyDeviceToHost);

        gpuCopyFromDevice(N, M, palleteSize, d_B, d_P,  B, P);

        // writePPM(B);


        // Write to scren

        for(int i = 0; i < image_width; i++)
        {
            for (int j = 0; j < image_height; j++)
            {

                Color c = B[i*image_height + j];

                 // SDL Draw
                SDL_SetRenderDrawColor(renderer, c.R, c.G, c.B, 255);
                SDL_RenderDrawPoint(renderer, i, j);
            }
        }
       


        // Free CPU memory
        // free(A);
        free(B);
        free(P);
        free(output_start_host);
        free(output_end_host);

        // Free GPU memory
        // cudaFree(d_A);
        // cudaFree(d_B);
        // cudaFree(d_P);
        // cudaFree(d_output_start);
        // cudaFree(d_output_end);


        gpuFree(d_B, d_P, d_output_start, d_output_end);

        // printf("\n---------------------------\n");
        // printf("__SUCCESS__\n");
        // printf("---------------------------\n");
        // printf("N                 = %d\n", K);
        // printf("Threads Per Block = %d\n", thr_per_blk);
        // printf("Blocks In Grid    = %d\n", blk_in_grid);
        // printf("---------------------------\n\n");


        // printf("end! \n");

        // return 0;



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
        output_start+=0.15*factor;
        output_end-=0.1*factor;
        factor *= 0.9349;
        n_max+=5;



    }

    free(iterationCounts);
    free(pixelColours);

    return 0;
}