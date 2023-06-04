pub mod color;
pub mod mandelbrot;
use color::Color;

use crate::mandelbrot::Mandelbrot;
use crate::mandelbrot::Renderable;

fn main() {
    
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

    let mandelbrot = Mandelbrot::new(1000, 1000, -2.0, 2.0, 1.0, 64, 4, Some(blue_gold_palette));

    match mandelbrot {
        Ok(mut mandelbrot) => {
            mandelbrot.render();
            match mandelbrot.save_image("fractal.png") {
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
