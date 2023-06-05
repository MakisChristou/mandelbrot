use clap::Parser;
use num::Complex;
use std::{fmt, str::FromStr};

#[derive(Debug, Clone)]
pub struct Bounds {
    pub upper_left: Complex<f64>,
    pub lower_right: Complex<f64>,
}

// Automatically convert from a String input to Bounds
impl FromStr for Bounds {
    type Err = &'static str;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let parts: Vec<&str> = s.split(',').collect();

        if parts.len() != 4 {
            return Err("Invalid Bounds");
        }

        let upper_left_re = parts[0]
            .parse::<f64>()
            .map_err(|_| "Cannot parse start bound")?;
        let upper_left_im = parts[1]
            .parse::<f64>()
            .map_err(|_| "Cannot parse end bound")?;

        let lower_right_re = parts[2]
            .parse::<f64>()
            .map_err(|_| "Cannot parse start bound")?;
        let lower_right_im = parts[3]
            .parse::<f64>()
            .map_err(|_| "Cannot parse end bound")?;

        Ok(Bounds {
            upper_left: Complex {
                re: upper_left_re,
                im: upper_left_im,
            },
            lower_right: Complex {
                re: lower_right_re,
                im: lower_right_im,
            },
        })
    }
}

// Required for default_value_t
impl fmt::Display for Bounds {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{},{},{},{}",
            self.upper_left.re, self.upper_left.im, self.lower_right.re, self.lower_right.im
        )
    }
}

/// Simple mandelbrot fractal renderer
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Args {
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
    #[arg(long, default_value_t = Bounds{upper_left: Complex{re: -2.0, im: 1.5}, lower_right: Complex{re: 1.5, im: -2.0}})]
    pub bounds: Bounds,

    /// Number of iterations to do before assuming point is in the set
    #[arg(short, long, default_value_t = 64)]
    pub n_max: u32,

    /// Anti-aliasing
    #[arg(short, long, default_value_t = 4)]
    pub s_max: u32,
}
