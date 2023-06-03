pub mod color;
pub mod mandelbrot;
use crate::mandelbrot::Mandelbrot;
use crate::mandelbrot::Renderable;

fn main() {
    let mut mandelbrot = Mandelbrot::new(1000, 1000, -2.0, 2.0, 1.0, 64, 4);

    match mandelbrot {
        Ok(mut mandelbrot) => {
            mandelbrot.render();
            mandelbrot.write_ppm();
        }
        Err(e) => {
            println!("Mandelbrot failed to start with error {:?}", e);
        }
    }
}
