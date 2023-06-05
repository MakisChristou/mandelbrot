pub mod args;
pub mod color;
pub mod config;
pub mod mandelbrot;

use crate::mandelbrot::Mandelbrot;
use crate::mandelbrot::Renderable;
use args::Args;
use clap::Parser;
use color::Color;
use config::Config;
use core::panic;

fn main() {
    let args = Args::parse();

    let blue_gold_palette = vec![
        Color::new(0, 0, 64),      // Dark Blue
        Color::new(0, 0, 128),     // Medium Blue
        Color::new(0, 0, 192),     // Lighter Blue
        Color::new(64, 64, 224),   // Blue with a bit of purple
        Color::new(128, 128, 255), // Very light blue
        Color::new(192, 192, 192), // Light Grey
        Color::new(224, 224, 128), // Yellowish Grey
        Color::new(255, 255, 64),  // Light Yellow
        Color::new(255, 224, 0),   // Gold
    ];

    let mandelbrot = Mandelbrot::new(Config {
        width: args.width,
        height: args.height,
        bounds: args.bounds,
        factor: 1.0,
        n_max: args.n_max,
        s_max: args.s_max,
        color_pallete: Some(blue_gold_palette),
    });

    match mandelbrot {
        Ok(mut mandelbrot) => {
            mandelbrot.render();
            match mandelbrot.save_image(&args.file_path) {
                Ok(()) => {}
                Err(e) => {
                    panic!("Image failed to be saved {:?}", e)
                }
            }
        }
        Err(e) => {
            panic!("Mandelbrot failed to start with error {:?}", e);
        }
    }
}
