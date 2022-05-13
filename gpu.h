typedef struct 
{
    int R;
    int G;
    int B;
}Color;

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


Color* gpuAllocColor(int N, int M, int bytes);
double* gpuAllocDouble(int N, int M, int bytes);
void gpuFree(Color* d_B, Color* d_P, double* d_output_start, double* d_output_end);
void gpuCopyToDevice(int N, int M, int palleteSize, Color* d_B, Color* d_P, Color* B, Color* P, double* d_output_start, double* d_output_end, double* output_start_host, double* output_end_host);
void gpuCopyFromDevice(int N, int M, int palleteSize, Color* d_B, Color* d_P, Color* B, Color* P);
void gpuRender(Color* d_B, Color* d_P, int palleteSize, int N, int M, double* d_output_start, double* d_output_end, int n_max, int s_max);
void gpuUpdateBounds(int N, int M, int palleteSize, double* d_output_start, double* d_output_end, double* output_start_host, double* output_end_host);
