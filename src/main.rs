pub mod color;
pub mod mandelbrot;
use crate::mandelbrot::Mandelbrot;
use crate::mandelbrot::Renderable;

fn main() {
    let mandelbrot = Mandelbrot::new(1000, 1000, -2.0, 2.0, 1.0, 64, 4);

    match mandelbrot {
        Ok(mut mandelbrot) => {
            // mandelbrot.parallel_render(Some(10));
            mandelbrot.render();
            // mandelbrot.write_ppm();
            match mandelbrot.save_image("fractal.png") {
                Ok(()) => {},
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
