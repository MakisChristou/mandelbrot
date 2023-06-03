use crate::color::Color;
use itertools::Itertools;

pub struct Mandelbrot {
    width: u32,
    height: u32,
    output_start: f64,
    output_end: f64,
    factor: f64,
    n_max: u32,
    s_max: u32,
    pixel_colours: Vec<Vec<Color>>,
    color_pallete: Vec<Color>,
}

#[derive(Debug)]

pub enum MandelbrotError {
    InvalidDimentions,
    InvalidRenderRange,
    InvalidIterations,
    InvalidAntiAliasing(&'static str),
    DivideByZero,
}

pub trait Renderable {
    fn render(&mut self);
}

impl Renderable for Mandelbrot {
    fn render(&mut self) {
        self.pixel_colours.clear();

        for _i in 0..self.width {
            self.pixel_colours.push(Vec::new());
        }

        for (i, j) in (0..self.width).cartesian_product(0..self.height) {
            let mut n = 0;
            let mut sum = 0;
            let mut k: f64 = 0.0;

            while k < 1.0 {
                let ii = i as f64 + k;
                let jj = j as f64 + k;

                let coordinates = self.pixels_to_coordinates(ii, jj);
                let citerations = self.iterate_mandelbrot(coordinates.0, coordinates.1);

                n = citerations.2;
                sum += n;
                k += 1.0 / (self.s_max as f64);
            }

            sum = sum / self.s_max;
            n = sum;

            let color = self.get_color(n);
            self.pixel_colours[i as usize].push(color);
        }
    }
}

impl Mandelbrot {
    pub fn new(
        width: u32,
        height: u32,
        output_start: f64,
        output_end: f64,
        factor: f64,
        n_max: u32,
        s_max: u32,
    ) -> Result<Self, MandelbrotError> {
        if width != height {
            return Err(MandelbrotError::InvalidDimentions);
        }

        if output_start >= output_end {
            return Err(MandelbrotError::InvalidRenderRange);
        }

        if n_max <= 0 {
            return Err(MandelbrotError::InvalidIterations);
        }

        if s_max <= 0 || !s_max.is_power_of_two() {
            return Err(MandelbrotError::InvalidAntiAliasing("Must be a power of 2"));
        }

        let pixel_colours: Vec<Vec<Color>> = vec![vec![]];

        // Default color pallete
        let color_pallete = vec![
            Color::new(0, 7, 100),
            Color::new(32, 107, 203),
            Color::new(237, 255, 255),
            Color::new(255, 170, 0),
            Color::new(0, 2, 0),
        ];

        Ok(Mandelbrot {
            width,
            height,
            output_start,
            output_end,
            factor,
            n_max,
            s_max,

            pixel_colours,
            color_pallete,
        })
    }

    fn map(
        input: f64,
        output_start: f64,
        output_end: f64,
        input_start: f64,
        input_end: f64,
    ) -> Result<f64, MandelbrotError> {
        if input_end - input_start == 0.0 {
            return Err(MandelbrotError::DivideByZero);
        }
        Ok(output_start
            + ((output_end - output_start) / (input_end - input_start)) * (input - input_start))
    }

    fn linear_interpolation(&self, v: &Color, u: &Color, a: f64) -> Color {
        let b: f64 = 1.0 - a;

        let R = (b * (v.R as f64) + a * (u.R as f64)) as u32;
        let G = (b * (v.G as f64) + a * (u.G as f64)) as u32;
        let B = (b * (v.B as f64) + a * (u.B as f64)) as u32;

        Color::new(R, G, B)
    }

    fn get_color(&self, iter: u32) -> Color {
        // Stolen Code from https://github.com/sevity/mandelbrot
        let MAX_COLOR: usize = self.color_pallete.len() - 1;
        let mu = (iter as f64) / self.n_max as f64;
        //scale mu to be in the range of colors
        let mu = mu * MAX_COLOR as f64;
        let i_mu = mu as usize;
        let color1 = &self.color_pallete[i_mu];
        let color2 = &self.color_pallete[std::cmp::min(i_mu + 1, MAX_COLOR)];
        let mut c = self.linear_interpolation(color1, color2, mu - i_mu as f64);

        if iter == self.n_max {
            c.R = 0;
            c.G = 0;
            c.B = 0;
        }

        c
    }

    fn pixels_to_coordinates(&self, i: f64, j: f64) -> (f64, f64) {
        let complex_i = Mandelbrot::map(
            i,
            self.output_start,
            self.output_end,
            0 as f64,
            self.width as f64,
        )
        .expect("Divide by 0 occured");

        let complex_j = Mandelbrot::map(
            j,
            self.output_start,
            self.output_end,
            0 as f64,
            self.height as f64,
        )
        .expect("Divide by 0 occured");

        (complex_i, complex_j)
    }

    fn iterate_mandelbrot(&self, complex_i: f64, complex_j: f64) -> (f64, f64, u32) {
        let mut x0 = 0.0;
        let mut y0 = 0.0;

        let mut x: f64 = 0.0;
        let mut y: f64 = 0.0;

        let mut n: u32 = 0;

        while x * x + y * y <= 4.0 && n < self.n_max {
            x = x0 * x0 - y0 * y0 + complex_i;
            y = 2.0 * x0 * y0 + complex_j;

            x0 = x;
            y0 = y;
            n += 1;
        }

        return (x, y, n);
    }

    pub fn write_ppm(&self) {
        print!("P3\n{} {}\n255\n", self.width, self.height);
        for (i, j) in (0..self.width).cartesian_product(0..self.height) {
            let c = &self.pixel_colours[j as usize][i as usize];
            print!("{} {} {}\n", c.R, c.G, c.B);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{Mandelbrot, Renderable};
    use crate::mandelbrot::MandelbrotError;

    #[test]
    fn should_instantiate_mandelbrot() {
        let mandelbrot = Mandelbrot::new(1000, 1000, -2.0, 2.0, 1.0, 64, 4);
        assert!(mandelbrot.is_ok());
    }

    #[test]
    fn width_equal_to_height() {
        let mandelbrot = Mandelbrot::new(640, 480, -2.0, 2.0, 1.0, 64, 4);
        assert!(matches!(
            mandelbrot,
            Err(MandelbrotError::InvalidDimentions)
        ));
    }

    #[test]
    fn start_before_end() {
        let mandelbrot = Mandelbrot::new(1000, 1000, 2.0, -2.0, 1.0, 64, 4);
        assert!(matches!(
            mandelbrot,
            Err(MandelbrotError::InvalidRenderRange)
        ));
    }

    #[test]
    fn valid_iterations() {
        let mandelbrot = Mandelbrot::new(1000, 1000, -2.0, 2.0, 1.0, 0, 4);
        assert!(matches!(
            mandelbrot,
            Err(MandelbrotError::InvalidIterations)
        ));
    }

    #[test]
    fn valid_anti_aliasing() {
        let mandelbrot = Mandelbrot::new(1000, 1000, -2.0, 2.0, 1.0, 64, 5);
        assert!(matches!(
            mandelbrot,
            Err(MandelbrotError::InvalidAntiAliasing("Must be a power of 2"))
        ));
    }

    #[test]
    fn does_not_diverge_to_infinity_at_zero() {
        let n_max = 64;
        let mandelbrot = Mandelbrot::new(1000, 1000, -2.0, 2.0, 1.0, n_max, 4);

        match mandelbrot {
            Ok(mut mandelbrot) => {
                let iterations = mandelbrot.iterate_mandelbrot(0.0, 0.0);
                assert_eq!(iterations.2, n_max)
            }
            Err(_) => {
                assert!(false)
            }
        }
    }

    #[test]
    fn diverges_to_infinity_at_known_point() {
        let mandelbrot = Mandelbrot::new(1000, 1000, -2.0, 2.0, 1.0, 64, 4);

        match mandelbrot {
            Ok(mut mandelbrot) => {
                let iterations = mandelbrot.iterate_mandelbrot(-100.0, 0.0);
                assert_eq!(iterations.2, 1)
            }
            Err(_) => {
                assert!(false)
            }
        }
    }
}
