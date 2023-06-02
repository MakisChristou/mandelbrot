pub mod color;
pub mod mandelbrot;

use crate::mandelbrot::Mandelbrot;

fn main() {
    let mut mandelbrot = Mandelbrot::new(1000, 1000, -2.0, 2.0, 1.0, 64, 4);
    mandelbrot.render();
    mandelbrot.write_ppm();
}
