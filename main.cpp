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
#include <GLFW/glfw3.h>

// Convert pixel coords to complex coords
inline long double map(long double input, long double output_start, long double output_end, long double input_start, long double input_end)
{
    return (output_start + ((output_end - output_start) / (input_end - input_start)) * (input - input_start));
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

// Timer Class for performance evaluation
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
    int image_width = 2080;
    int image_height = 2080; 

    int n_max = 64; // 4096
    int s_max = 8; // prefer to be a power of 2


    int N = image_width;
    int M = image_height; 

    long double output_start_x = -3.0f;
    long double output_end_x = 2.66f;
    double output_start_y = -3.0f;
    double output_end_y = 2.66f;


    double factor = 1.0f;
    int palleteSize = 5;

    // Used for zooming
    long double scaleX = 1.00f;
    long double scaleY = 1.00f;

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

    // New Coordinates
    double* output_start_host_x = (double*)malloc(sizeof(double));
    double* output_end_host_x = (double*)malloc(sizeof(double));
    double* output_start_host_y = (double*)malloc(sizeof(double));
    double* output_end_host_y = (double*)malloc(sizeof(double));

    // to allocate memory for arrays d_A, d_B, and d_C on device
    Color* d_B;
    Color* d_P;

    // New Coordinates
    double* d_output_start_x;
    double* d_output_end_x;
    double* d_output_start_y;
    double* d_output_end_y;

    // Actually allocate GPU memory
    d_B = gpuAllocColor(N, M, pixelBytes);
    d_P = gpuAllocColor(N, M, palleteBytes);


    d_output_start_x = gpuAllocDouble(N, M, sizeof(double));
    d_output_end_x = gpuAllocDouble(N, M, sizeof(double));
    d_output_start_y = gpuAllocDouble(N, M, sizeof(double));
    d_output_end_y = gpuAllocDouble(N, M, sizeof(double));


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
    gpuCopyToDevice(N, M, palleteSize, d_B, d_P, B, P, 
    d_output_start_x, d_output_end_x, output_start_host_x, output_end_host_x,
    d_output_start_y, d_output_end_y, output_start_host_y, output_end_host_y);

    // SDL2 Texture for faster rendering
    SDL_Texture* theTexture = SDL_CreateTexture( renderer, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STREAMING, image_width, image_height );
    uint32_t *textureBuffer = new uint32_t[ image_width * image_height ];


    // Main render loop
    while(1)
    {
        Timer timer;

        
        SDL_RenderPresent(renderer);
        SDL_Event event;
        while (SDL_PollEvent(&event))
        {
            // User clicks the exit button
            if (event.type == SDL_QUIT) {
                SDL_Quit();
                return 0;
            }
            // Key press
            else if(event.type == SDL_KEYDOWN)
            {   
                long double distance_x = output_end_x - output_start_x;
                long double distance_y = output_end_y - output_start_y;


                // Arrow Keys for panning
                if(event.key.keysym.sym == SDLK_RIGHT)
                {   
                    output_start_x+=distance_x/(1*100);
                    output_end_x+=distance_x/(1*100);
                }
                else if(event.key.keysym.sym == SDLK_LEFT)
                {
                    output_start_x-=distance_x/(1*100);
                    output_end_x-=distance_x/(1*100);
                }
                else if(event.key.keysym.sym == SDLK_DOWN)
                {
                    output_start_y+=distance_y/(1*100);
                    output_end_y+=distance_y/(1*100);
                }
                else if(event.key.keysym.sym == SDLK_UP)
                {
                    output_start_y-=distance_y/(1*100);
                    output_end_y-=distance_y/(1*100);
                }


                if(event.key.keysym.sym == SDLK_w)
                {
                    n_max++;
                }
                else if(event.key.keysym.sym == SDLK_s)
                {   
                    if(n_max > 1)
                        n_max--;
                }

            }

            // Scroll Wheel
            else if(event.type == SDL_MOUSEWHEEL)
            {
                // Get mouse position
                int mouse_x;
                int mouse_y;
                SDL_GetMouseState(&mouse_x, &mouse_y);

                // Convert mouse position to complex plane
                long double mouse_complex_i = map(mouse_x, output_start_x, output_end_x, 0, image_width);
                long double mouse_complex_j = map(mouse_y, output_start_y, output_end_y, 0, image_height);

                long double distance_x = output_end_x - output_start_x;
                long double distance_y = output_end_y - output_start_y;
                
                if(event.wheel.y > 0) // scroll up
                {
                    output_start_x+=distance_x/(1*100);
                    output_end_x-=distance_x/(1*100);
                    output_start_y+=distance_y/(1*100);
                    output_end_y-=distance_y/(1*100);
                }
                else if(event.wheel.y < 0) // scroll down
                {
                    output_start_x-=distance_x/(1*100);
                    output_end_x+=distance_x/(1*100);
                    output_start_y-=distance_y/(1*100);
                    output_end_y+=distance_y/(1*100);
                }

                
            }
        }


        // Update render coordinates
        output_start_host_x[0] = output_start_x;
        output_end_host_x[0] = output_end_x;
        output_start_host_y[0] = output_start_y;
        output_end_host_y[0] = output_end_y;

        // Copy render coordinates to GPU memory
        gpuUpdateBounds(N,M, palleteSize, d_output_start_x, d_output_end_x, output_start_host_x, output_end_host_x, d_output_start_y, d_output_end_y, output_start_host_y, output_end_host_y);

        
        auto start = high_resolution_clock::now();
        gpuRender( d_B, d_P,palleteSize,  N,  M, d_output_start_x, d_output_end_x, d_output_start_y, d_output_end_y, n_max, s_max);
        auto stop = high_resolution_clock::now();
        auto gpuRenderDuration = duration_cast<microseconds>(stop - start);


        start = high_resolution_clock::now();
        gpuCopyFromDevice(N, M, palleteSize, d_B, d_P,  B, P);
        stop = high_resolution_clock::now();
        auto gpuCopyFromDeviceDuration = duration_cast<microseconds>(stop - start);
        

        


        start = high_resolution_clock::now();
        
        // Write to screen
        for(int i = 0; i < image_width; i++)
        {
            for (int j = 0; j < image_height; j++)
            {
                Color c = B[i*image_height + j];
                uint32_t ir = c.R;
                uint32_t ig = c.G;
                uint32_t ib = c.B;
                textureBuffer[j * image_height + i] = 0xFF000000 | (ir<<16) | (ig<<8) | ib;
            }
        }
        
        SDL_Rect texture_rect; //create a rect
        texture_rect.x = 0;  //controls the rect's x coordinate 
        texture_rect.y = 0; // controls the rect's y coordinte
        texture_rect.w = image_width; // controls the width of the rect
        texture_rect.h = image_height; // controls the height of the rect

        SDL_UpdateTexture( theTexture , NULL, textureBuffer, image_width * sizeof (uint32_t));
        SDL_RenderCopy(renderer, theTexture, NULL, &texture_rect); 
        stop = high_resolution_clock::now();
        auto SDLRenderDuration = duration_cast<microseconds>(stop - start);

        // Segfault here ?
        // SDL_Surface* sshot = SDL_GetWindowSurface(window);
        // printf("Test1");
        // SDL_RenderReadPixels(renderer, NULL, SDL_PIXELFORMAT_BGRA8888, sshot->pixels, sshot->pitch);
        // std::string fileName = std::to_string(count) + ".bmp";
        // SDL_SaveBMP(sshot, fileName.c_str());
        // SDL_FreeSurface(sshot);
       

        // Zoom in code by https://www.youtube.com/watch?v=KnCNfBb2ODQ

        // Automatic Zoom 
        // output_start+=0.15*factor;
        // output_end-=0.1*factor;
        // factor *= 0.9550;

        //Duration of 1 input pattern in ms
		float iterationTime = (float) std::chrono::duration_cast<std::chrono::microseconds>(timer.duration).count()/1000;
        float fps = (float)1000.0/iterationTime;
        int x = (int)(1000/iterationTime);
        // Avoid division by 0
        if(x == 0)
			x = 1;	
			
        // Update every 1 second    
        if((count == 1) || (count % x) == 0)
            std::cerr << "\r" << " Iterations: " << n_max <<" fps: " << fps << " gpuRender time : " << gpuRenderDuration.count()/1000 << " ms " << "gpuCopyFromDevice time : " << gpuCopyFromDeviceDuration.count()/1000 << " ms " << "SDLRender time : " << SDLRenderDuration.count()/1000 << " ms " << std::flush;


        

        count++;

    }

    // Free CPU memory
    free(B);
    free(P);
    free(output_start_host_x);
    free(output_end_host_x);
    free(output_start_host_y);
    free(output_end_host_y);


    // Free GPU memory
    gpuFree(d_B, d_P, d_output_start_x, d_output_end_x, d_output_start_y, d_output_end_y);

    // Try and look pretty
    std::cerr << std::endl;
    std::cerr << std::endl;

    return 0;
}