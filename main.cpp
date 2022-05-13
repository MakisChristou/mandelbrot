#include <iostream>
#include <vector>
#include <cstdlib>
#include <math.h>
#include <stdio.h>
#include <cuda.h>
#include <chrono>
#include <string>

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

//Timer Class for performance evaluation
class Timer{
	public:
	//Start timer
	std::chrono::steady_clock::time_point start, end;
	std::chrono::duration<float> duration;

	//Constructor
	Timer()
	{
		start = std::chrono::steady_clock::now();
	}

	//Destructor
	~Timer()
	{
		end = std::chrono::steady_clock::now();
		duration = end-start; //duration in seconds
		float ms = duration.count() * 1000.0f;
	}
};

int main(int argc, char* argv[])
{

    int image_width = 800;
    int image_height = 800; 

    int n_max = 64; // 4096
    int s_max = 8; // prefer to be a power of 2


    int N = image_width;
    int M = image_height; 

    double output_start = -3.0;
    double output_end = 2.66f;

    // long double output_start = 0.2f;
    // long double output_end = 0.5f;

    // double output_start = 0.35f;
    // double output_end = 0.36f;


    double factor = 1.0f;
    int palleteSize = 5;

    int count = 0; // frame counter


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
        Timer timer;

        
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

        // Segfault here ?

        // SDL_Surface* sshot = SDL_GetWindowSurface(window);
        // printf("Test1");
        // SDL_RenderReadPixels(renderer, NULL, SDL_PIXELFORMAT_BGRA8888, sshot->pixels, sshot->pitch);
        // std::string fileName = std::to_string(count) + ".bmp";
        // SDL_SaveBMP(sshot, fileName.c_str());
        // SDL_FreeSurface(sshot);


        stop = high_resolution_clock::now();
        auto SDLRenderDuration = duration_cast<microseconds>(stop - start);
        // std::cerr << "\rgpuRender time : " << gpuRenderDuration.count()/1000 << " ms " << "gpuCopyFromDevice time : " << gpuCopyFromDeviceDuration.count()/1000 << " ms " << "SDLRender time : " << SDLRenderDuration.count()/1000 << " ms " << std::flush;
       

        // Zoom in code by https://www.youtube.com/watch?v=KnCNfBb2ODQ
        
        output_start+=0.15*factor;
        output_end-=0.1*factor;
        factor *= 0.9550;


        // std::chrono::duration<float> duration = timer.duration;

        //Duration of 1 input pattern in ms
		float iterationTime = (float) std::chrono::duration_cast<std::chrono::microseconds>(timer.duration).count()/1000;
        float fps = (float)1000.0/iterationTime;
        int x = (int)(1000/iterationTime);

        // Update every 1 second
        if((count == 1) || (count % x) == 0)
            std::cerr <<"\r" << output_start << " " << output_end << " " << factor  << " " << n_max <<" fps: " << fps << std::flush;


        n_max+=1;




        count++;

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