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
void writePPM(Color* pixelColors, int image_width, int image_height)
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


// inline std::vector<Color> generateColorPalete()
// {
//     std::vector<Color> colorPalete;

//     Color color1;
//     color1.R = 0xFF;
//     color1.G = 0xFF;
//     color1.B = 0xFF;

//     Color color2;
//     color2.R = 0x00;
//     color2.G = 0x00;
//     color2.B = 0x00;

//     Color black;
//     black.R = 0x00;
//     black.G = 0x00;
//     black.B = 0x00;

//     long double sumR = 0;
//     long double sumG = 0;
//     long double sumB = 0;



//     for(int i = 0; i < n_max; i++)
//     {
//         if(i == 0)
//         {
//             colorPalete.push_back(color1);
//         }
//         else if(i == n_max-1)
//         {
//             colorPalete.push_back(color2);
//         }
//         else if(i == n_max)
//         {
//             colorPalete.push_back(black);
//         }
//         else
//         {
//             Color tempColor;

//             // Calculate log step for current iteration
//             double long rStep = log2(1 + (color1.R * log2(i)) / 10) / 8;
//             double long gStep = log2(1 + (color1.G * log2(i)) / 10) / 8;
//             double long bStep = log2(1 + (color1.B * log2(i)) / 10) / 8;

//             // 
//             sumR += rStep;
//             sumG += gStep;
//             sumB += bStep;

//             tempColor.R = sumR;
//             tempColor.G = sumG;
//             tempColor.B = sumB;

//             colorPalete.push_back(tempColor);
//         }
//     }

//     return colorPalete;
// }


int main(int argc, char* argv[])
{

    int image_width = 1080;
    int image_height = 1080; 

    int n_max = 64; // 4096
    int s_max = 8; // prefer to be a power of 2


    int N = image_width;
    int M = image_height; 

    double output_start = -2.0;
    double output_end = 2.0f;

    // long double output_start = 0.2f;
    // long double output_end = 0.5f;

    // long double output_start = 0.35f;
    // long double output_end = 0.36f;


    double factor = 1.0f;
    int palleteSize = 5;


    // SDL Stuff
    SDL_Init(SDL_INIT_EVERYTHING);
    SDL_Window* window;
    SDL_Renderer* renderer;
    SDL_Event event;
    SDL_CreateWindowAndRenderer(image_width, image_height, 0, &window, &renderer);
    SDL_RenderSetLogicalSize(renderer, image_width, image_height);


    // Number of bytes to allocate for N doubles
    // size_t iterationBytes = N*M*sizeof(int);
    size_t pixelBytes = N*M*sizeof(Color);
    size_t palleteBytes = palleteSize*sizeof(Color);

    // Allocate memory for arrays B, P on host
    Color *B = (Color*)malloc(pixelBytes); // pixelColours
    Color *P = (Color*)malloc(palleteBytes); // colorPalete
    double* output_start_host = (double*)malloc(sizeof(double));
    double* output_end_host = (double*)malloc(sizeof(double));


    // to allocate memory for arrays d_A, d_B, and d_C on device
    Color* d_B;
    Color* d_P;
    double* d_output_start;
    double* d_output_end;

    // Actually allocate GPU memory
    d_B = gpuAllocColor(N, M, pixelBytes);
    d_P = gpuAllocColor(N, M, palleteBytes);
    d_output_start = gpuAllocDouble(N, M, sizeof(double));
    d_output_end = gpuAllocDouble(N, M, sizeof(double));

    // Fill host data structures
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

    // Copy host data to the device
    gpuCopyToDevice(N,M,palleteSize,d_B,d_P,B,P,d_output_start,d_output_end,output_start_host,output_end_host);


    // Main render loop
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


        // Fill host arrays data structures
        output_start_host[0] = output_start;
        output_end_host[0] = output_end;


        gpuUpdateBounds(N,M, palleteSize,d_output_start, d_output_end,output_start_host,output_end_host);

        
        auto start = high_resolution_clock::now();
        gpuRender( d_B, d_P,palleteSize,  N,  M,  d_output_start, d_output_end, n_max, s_max);
        auto stop = high_resolution_clock::now();
        auto gpuRenderDuration = duration_cast<microseconds>(stop - start);


        start = high_resolution_clock::now();
        gpuCopyFromDevice(N, M, palleteSize, d_B, d_P,  B, P);
        stop = high_resolution_clock::now();
        auto gpuCopyFromDeviceDuration = duration_cast<microseconds>(stop - start);
        
        // For file output
        // writePPM(B, N, M);

        // Write to scren
        start = high_resolution_clock::now();
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
        stop = high_resolution_clock::now();
        auto SDLRenderDuration = duration_cast<microseconds>(stop - start);
        std::cerr << "\rgpuRender time : " << gpuRenderDuration.count()/1000 << " ms " << "gpuCopyFromDevice time : " << gpuCopyFromDeviceDuration.count()/1000 << " ms " << "SDLRender time : " << SDLRenderDuration.count()/1000 << " ms " << std::flush;
       

        // Zoom in code by https://www.youtube.com/watch?v=KnCNfBb2ODQ
        output_start+=0.15*factor;
        output_end-=0.1*factor;
        factor *= 0.9349;
        n_max+=1;
    }

    // Free CPU memory
    free(B);
    free(P);
    free(output_start_host);
    free(output_end_host); 

    // Free GPU memory
    gpuFree(d_B, d_P, d_output_start, d_output_end);

    // Try and look pretty
    std::cerr << std::endl;
    std::cerr << std::endl;

    return 0;
}