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
void gpuFree(Color* d_B, Color* d_P, double* d_output_start_x, double* d_output_end_x, double* d_output_start_y, double* d_output_end_y);
void gpuCopyToDevice(int N, int M, int palleteSize, Color* d_B, Color* d_P, Color* B, Color* P, double* d_output_start_x, double* d_output_end_x, double* output_start_host_x, double* output_end_host_x, double* d_output_start_y, double* d_output_end_y, double* output_start_host_y, double* output_end_host_y);
void gpuCopyFromDevice(int N, int M, int palleteSize, Color* d_B, Color* d_P, Color* B, Color* P);
void gpuRender(Color* d_B, Color* d_P, int palleteSize, int N, int M, double* d_output_start_x, double* d_output_end_x, double* d_output_start_y, double* d_output_end_y, int n_max, int s_max);
void gpuUpdateBounds(int N, int M, int palleteSize, double* d_output_start_x, double* d_output_end_x, double* output_start_host_x, double* output_end_host_x, double* d_output_start_y, double* d_output_end_y, double* output_start_host_y, double* output_end_host_y);
