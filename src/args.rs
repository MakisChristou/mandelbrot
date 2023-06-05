use std::{str::FromStr, fmt};
use clap::Parser;

#[derive(Debug, Clone)]
pub struct Bounds {
    pub output_start: f64,
    pub output_end: f64,
}

// Automatically convert from a String input to Bounds
impl FromStr for Bounds {
    type Err = &'static str;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let parts: Vec<&str> = s.split(',').collect();

        if parts.len() != 2 {
            return Err("Invalid Bounds");
        }

        let output_start = parts[0].parse::<f64>().map_err(|_| "Cannot parse start bound")?;
        let output_end = parts[1].parse::<f64>().map_err(|_| "Cannot parse end bound")?;

        Ok(
            Bounds{
                output_start,
                output_end,
            }
        )
    }
}

// Required for default_value_t
impl fmt::Display for Bounds {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{},{}", self.output_start, self.output_end)
    }
}

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
    #[arg(long, default_value_t = Bounds{output_start: -2.0, output_end: 1.5})]
    pub bounds: Bounds,

    /// Number of iterations to do before assuming point is in the set
    #[arg(short, long, default_value_t = 64)]
    pub n_max: u32,

    /// Anti-aliasing
    #[arg(short, long, default_value_t = 4)]
    pub s_max: u32,
    
}