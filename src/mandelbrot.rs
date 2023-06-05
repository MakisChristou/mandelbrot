use crate::args::Bounds;
use crate::color::Color;
use crate::config::Config;
use image::ImageError;
use itertools::Itertools;
use num::Complex;
use std::thread;
use std::usize;

#[derive(Clone)]
pub struct Mandelbrot {
    width: usize,
    height: usize,
    bounds: Bounds,
    factor: f64,
    n_max: u32,
    s_max: u32,
    pixel_colours: Vec<Color>,
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
    fn parallel_render(&mut self, threads: Option<usize>);
}

impl Renderable for Mandelbrot {
    fn render(&mut self) {
        self.pixel_colours = vec![Color::new(0, 0, 0); self.height * self.width];

        for (i, j) in (0..self.width).cartesian_product(0..self.height) {
            self.pixel_colours[j * self.height + i] = self.colorize_pixel(i, j);
        }
    }

    fn parallel_render(&mut self, threads: Option<usize>) {
        let mut cpus = num_cpus::get();

        if let Some(threads) = threads {
            cpus = threads;
        }

        let mut handlers = Vec::new();
        let ranges = self.split_pixels(cpus);

        // Split pixel_colors into ranges
        for (l, h) in ranges {
            // Clone whole mandelbrot for now
            let cloned_self = self.clone();

            let handler = thread::spawn(move || cloned_self.render_band(l, h));
            handlers.push(handler);
        }

        self.pixel_colours = vec![];

        for handler in handlers {
            let mut render_slice = handler.join().unwrap();
            self.pixel_colours.append(&mut render_slice);
        }
    }
}

impl Mandelbrot {
    pub fn new(config: Config) -> Result<Self, MandelbrotError> {
        if config.width != config.height {
            return Err(MandelbrotError::InvalidDimentions);
        }

        if config.bounds.upper_left.re >= config.bounds.lower_right.re
            || config.bounds.upper_left.im <= config.bounds.lower_right.im
        {
            return Err(MandelbrotError::InvalidRenderRange);
        }

        if config.n_max == 0 {
            return Err(MandelbrotError::InvalidIterations);
        }

        if config.s_max == 0 || !config.s_max.is_power_of_two() {
            return Err(MandelbrotError::InvalidAntiAliasing("Must be a power of 2"));
        }

        let pixel_colours: Vec<Color> = vec![];

        // Default color pallete
        let default_color_pallete = vec![
            Color::new(0, 7, 100),
            Color::new(32, 107, 203),
            Color::new(237, 255, 255),
            Color::new(255, 170, 0),
            Color::new(0, 2, 0),
        ];

        if let Some(color_pallete) = config.color_pallete {
            return Ok(Mandelbrot {
                width: config.width,
                height: config.height,
                bounds: config.bounds,
                factor: config.factor,
                n_max: config.n_max,
                s_max: config.s_max,
                pixel_colours,
                color_pallete,
            });
        }

        Ok(Mandelbrot {
            width: config.width,
            height: config.height,
            bounds: config.bounds,
            factor: config.factor,
            n_max: config.n_max,
            s_max: config.s_max,
            pixel_colours,
            color_pallete: default_color_pallete,
        })
    }

    fn pixel_to_point(
        bounds: (usize, usize),
        pixel: (f64, f64),
        upper_left: Complex<f64>,
        lower_right: Complex<f64>,
    ) -> Complex<f64> {
        let (width, height) = (
            lower_right.re - upper_left.re,
            upper_left.im - lower_right.im,
        );
        Complex {
            re: upper_left.re + pixel.0 * width / bounds.0 as f64,
            im: upper_left.im - pixel.1 * height / bounds.1 as f64, // Why subtraction here? pixel.1 increases as we go down,
                                                                    // but the imaginary component increases as we go up.
        }
    }

    fn render_band(&self, start: usize, end: usize) -> Vec<Color> {
        let mut local_render: Vec<Color> = Vec::new();

        for k in start..end {
            let i = k % self.width;
            let j = k / self.width;

            local_render.push(self.colorize_pixel(i, j));
        }

        local_render
    }

    fn split_pixels(&self, threads: usize) -> Vec<(usize, usize)> {
        let mut ranges = Vec::new();

        for i in 0..threads {
            let lower_bound = i * self.width * self.height / threads;
            let upper_bound = (i + 1) * self.width * self.height / threads;

            ranges.push((lower_bound, upper_bound));
        }

        ranges
    }

    fn linear_interpolation(&self, v: &Color, u: &Color, a: f64) -> Color {
        let b: f64 = 1.0 - a;

        let r = (b * (v.r as f64) + a * (u.r as f64)) as u8;
        let g = (b * (v.g as f64) + a * (u.g as f64)) as u8;
        let b = (b * (v.b as f64) + a * (u.b as f64)) as u8;

        Color::new(r, g, b)
    }

    fn colorize_pixel(&self, i: usize, j: usize) -> Color {
        let mut n;
        let mut sum = 0;
        let mut k: f64 = 0.0;

        while k < 1.0 {
            let ii = i as f64 + k;
            let jj = j as f64 + k;

            let coordinates = Self::pixel_to_point(
                (self.width, self.height),
                (ii, jj),
                self.bounds.upper_left,
                self.bounds.lower_right,
            );
            let citerations = self.iterate_mandelbrot(coordinates);

            n = citerations.2;
            sum += n;
            k += 1.0 / (self.s_max as f64);
        }

        sum /= self.s_max;
        n = sum;

        self.get_color(n)
    }

    fn get_color(&self, iter: u32) -> Color {
        // Stolen Code from https://github.com/sevity/mandelbrot
        let max_color: usize = self.color_pallete.len() - 1;
        let mu = (iter as f64) / self.n_max as f64;
        //scale mu to be in the range of colors
        let mu = mu * max_color as f64;
        let i_mu = mu as usize;
        let color1 = &self.color_pallete[i_mu];
        let color2 = &self.color_pallete[std::cmp::min(i_mu + 1, max_color)];
        let mut c = self.linear_interpolation(color1, color2, mu - i_mu as f64);

        if iter == self.n_max {
            c.r = 0;
            c.g = 0;
            c.b = 0;
        }
        c
    }

    fn iterate_mandelbrot(&self, complex: Complex<f64>) -> (f64, f64, u32) {
        let mut x0 = 0.0;
        let mut y0 = 0.0;

        let mut x: f64 = 0.0;
        let mut y: f64 = 0.0;

        let mut n: u32 = 0;

        while x * x + y * y <= 4.0 && n < self.n_max {
            x = x0 * x0 - y0 * y0 + complex.re;
            y = 2.0 * x0 * y0 + complex.im;

            x0 = x;
            y0 = y;
            n += 1;
        }

        (x, y, n)
    }

    pub fn save_image(&self, file_path: &str) -> Result<(), ImageError> {
        let mut imgbuf = image::ImageBuffer::new(self.width as u32, self.height as u32);

        // Iterate over the coordinates and pixels of the image
        for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
            let c = &self.pixel_colours[y as usize * self.height + x as usize];
            *pixel = image::Rgb([c.r, c.g, c.b]);
        }

        imgbuf.save(file_path)
    }

    pub fn write_ppm(&self) {
        print!("P3\n{} {}\n255\n", self.width, self.height);
        for (i, j) in (0..self.width).cartesian_product(0..self.height) {
            let c = &self.pixel_colours[i * self.height + j];
            println!("{} {} {}", c.r, c.g, c.b);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::Mandelbrot;
    use crate::{args::Bounds, color::Color, config::Config, mandelbrot::MandelbrotError};
    use num::Complex;

    #[test]
    fn should_instantiate_mandelbrot() {
        let mandelbrot = Mandelbrot::new(Config {
            width: 1000,
            height: 1000,
            bounds: Bounds {
                upper_left: Complex { re: -2.0, im: 2.0 },
                lower_right: Complex { re: 2.0, im: -2.0 },
            },
            factor: 1.0,
            n_max: 64,
            s_max: 4,
            color_pallete: None,
        });
        assert!(mandelbrot.is_ok());
    }

    #[test]
    fn width_equal_to_height() {
        let mandelbrot = Mandelbrot::new(Config {
            width: 640,
            height: 480,
            bounds: Bounds {
                upper_left: Complex { re: -2.0, im: 2.0 },
                lower_right: Complex { re: 2.0, im: -2.0 },
            },
            factor: 1.0,
            n_max: 64,
            s_max: 4,
            color_pallete: None,
        });
        assert!(matches!(
            mandelbrot,
            Err(MandelbrotError::InvalidDimentions)
        ));
    }

    #[test]
    fn start_before_end() {
        let mandelbrot = Mandelbrot::new(Config {
            width: 1000,
            height: 1000,
            bounds: Bounds {
                upper_left: Complex { re: 2.0, im: -2.0 },
                lower_right: Complex { re: -2.0, im: 2.0 },
            },
            factor: 1.0,
            n_max: 64,
            s_max: 4,
            color_pallete: None,
        });
        assert!(matches!(
            mandelbrot,
            Err(MandelbrotError::InvalidRenderRange)
        ));
    }

    #[test]
    fn valid_iterations() {
        let mandelbrot = Mandelbrot::new(Config {
            width: 1000,
            height: 1000,
            bounds: Bounds {
                upper_left: Complex { re: -2.0, im: 2.0 },
                lower_right: Complex { re: 2.0, im: -2.0 },
            },
            factor: 1.0,
            n_max: 0,
            s_max: 4,
            color_pallete: None,
        });
        assert!(matches!(
            mandelbrot,
            Err(MandelbrotError::InvalidIterations)
        ));
    }

    #[test]
    fn valid_anti_aliasing() {
        let mandelbrot = Mandelbrot::new(Config {
            width: 1000,
            height: 1000,
            bounds: Bounds {
                upper_left: Complex { re: -2.0, im: 2.0 },
                lower_right: Complex { re: 2.0, im: -2.0 },
            },
            factor: 1.0,
            n_max: 64,
            s_max: 5,
            color_pallete: None,
        });
        assert!(matches!(
            mandelbrot,
            Err(MandelbrotError::InvalidAntiAliasing("Must be a power of 2"))
        ));
    }

    #[test]
    fn does_not_diverge_to_infinity_at_zero() {
        let n_max = 64;
        let mandelbrot = Mandelbrot::new(Config {
            width: 1000,
            height: 1000,
            bounds: Bounds {
                upper_left: Complex { re: -2.0, im: 2.0 },
                lower_right: Complex { re: 2.0, im: -2.0 },
            },
            factor: 1.0,
            n_max: n_max,
            s_max: 4,
            color_pallete: None,
        });

        match mandelbrot {
            Ok(mandelbrot) => {
                let iterations = mandelbrot.iterate_mandelbrot(Complex { re: 0.0, im: 0.0 });
                assert_eq!(iterations.2, n_max)
            }
            Err(_) => {
                assert!(false)
            }
        }
    }

    #[test]
    fn diverges_to_infinity_at_known_point() {
        let mandelbrot = Mandelbrot::new(Config {
            width: 1000,
            height: 1000,
            bounds: Bounds {
                upper_left: Complex { re: -2.0, im: 2.0 },
                lower_right: Complex { re: 2.0, im: -2.0 },
            },
            factor: 1.0,
            n_max: 64,
            s_max: 4,
            color_pallete: None,
        });

        match mandelbrot {
            Ok(mandelbrot) => {
                let iterations = mandelbrot.iterate_mandelbrot(Complex {
                    re: -100.0,
                    im: 0.0,
                });
                assert_eq!(iterations.2, 1)
            }
            Err(_) => {
                assert!(false)
            }
        }
    }

    #[test]
    fn custom_and_default_color_palette() {
        let mandelbrot = Mandelbrot::new(Config {
            width: 1000,
            height: 1000,
            bounds: Bounds {
                upper_left: Complex { re: -2.0, im: 2.0 },
                lower_right: Complex { re: 2.0, im: -2.0 },
            },
            factor: 1.0,
            n_max: 64,
            s_max: 4,
            color_pallete: None,
        });

        let vibrant_color_palette = vec![
            Color::new(66, 30, 15),    // Dark brown
            Color::new(25, 7, 26),     // Dark purple
            Color::new(9, 1, 47),      // Deep blue
            Color::new(4, 4, 73),      // Blue
            Color::new(0, 7, 100),     // Lighter blue
            Color::new(12, 44, 138),   // Blue-Indigo
            Color::new(24, 82, 177),   // Sky Blue
            Color::new(57, 125, 209),  // Light Blue
            Color::new(134, 181, 229), // Lighter Blue
            Color::new(211, 236, 248), // Very light blue
            Color::new(241, 233, 191), // Cream
            Color::new(248, 201, 95),  // Yellow
            Color::new(255, 170, 0),   // Orange
            Color::new(204, 128, 0),   // Dark Orange
            Color::new(153, 87, 0),    // Darker Orange
            Color::new(106, 52, 3),    // Even darker orange
        ];

        let mandelbrot = Mandelbrot::new(Config {
            width: 1000,
            height: 1000,
            bounds: Bounds {
                upper_left: Complex { re: -2.0, im: 2.0 },
                lower_right: Complex { re: 2.0, im: -2.0 },
            },
            factor: 1.0,
            n_max: 64,
            s_max: 4,
            color_pallete: Some(vibrant_color_palette),
        });
    }
}
