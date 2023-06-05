use clap::Parser;

/// Simple mandelbrot fractal renderer
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Args{
    /// Output image path
    #[arg(short, long, default_value_t = String::from("fractal.png"))]
    pub file_path: String,
    
    /// Output image width
    #[arg(long, default_value_t = 1000)]
    pub width: usize,

    /// Output image width
    #[arg(long, default_value_t = 1000)]
    pub height: usize,

    /// Start/End coordinate on the complex plane (both x and y)
    #[arg(long, default_value_t = String::from("-2.0, 2.0"))]
    pub bounds: String,

    /// Number of iterations to do before assuming point is in the set
    #[arg(short, long, default_value_t = 64)]
    pub n_max: u32,

    /// Anti-aliasing
    #[arg(short, long, default_value_t = 4)]
    pub s_max: u32,
    
}