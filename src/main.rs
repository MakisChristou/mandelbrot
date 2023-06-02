use std::f64;

struct Color {
    R: u32,
    G: u32,
    B: u32,
}

impl Color {
    fn new(R: u32, G: u32, B: u32) -> Self {
        Color { R, G, B }
    }
}

struct Mandelbrot {
    width: u32,
    height: u32,

    output_start: f64,
    output_end: f64,

    factor: f64,

    n_max: u32,
    s_max: u32,

    iteration_counts: Vec<Vec<u32>>,
    iteration_values: Vec<Vec<(f64, f64)>>,
    pixel_colours: Vec<Vec<Color>>,
    color_pallete: Vec<Color>,
}

impl Mandelbrot {
    fn new(
        width: u32,
        height: u32,
        output_start: f64,
        output_end: f64,
        factor: f64,
        n_max: u32,
        s_max: u32,
    ) -> Self {
        let iteration_counts: Vec<Vec<u32>> = vec![vec![]];
        let iteration_values: Vec<Vec<(f64, f64)>> = vec![vec![]];

        let pixel_colours: Vec<Vec<Color>> = vec![vec![]];

        // Default color pallete
        let color_pallete = vec![
            Color::new(0, 7, 100),
            Color::new(32, 107, 203),
            Color::new(237, 255, 255),
            Color::new(255, 170, 0),
            Color::new(0, 2, 0),
        ];

        Mandelbrot {
            width,
            height,
            output_start,
            output_end,
            factor,
            n_max,
            s_max,

            iteration_counts,
            iteration_values,
            pixel_colours,
            color_pallete,
        }
    }

    fn map(
        input: f64,
        output_start: f64,
        output_end: f64,
        input_start: f64,
        input_end: f64,
    ) -> f64 {
        output_start
            + ((output_end - output_start) / (input_end - input_start)) * (input - input_start)
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

    fn iterate_mandelbrot(&self, i: f64, j: f64) -> (f64, f64, u32) {
        let complex_i = Mandelbrot::map(
            i,
            self.output_start,
            self.output_end,
            0 as f64,
            self.width as f64,
        );
        let complex_j = Mandelbrot::map(
            j,
            self.output_start,
            self.output_end,
            0 as f64,
            self.height as f64,
        );

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

    fn render(&mut self) {
        let mut i = 0;

        self.iteration_counts.clear();
        self.iteration_values.clear();
        self.pixel_colours.clear();

        while i < self.width {
            self.iteration_counts.push(Vec::new());
            self.iteration_values.push(Vec::new());
            self.pixel_colours.push(Vec::new());

            let mut j = 0;

            while j < self.height {
                let mut c: (f64, f64) = (0.0, 0.0);

                let mut n = 0;
                let mut sum = 0;

                let mut k: f64 = 0.0;

                while k < 1.0 {
                    let ii = i as f64 + k;
                    let jj = j as f64 + k;

                    let citerations = self.iterate_mandelbrot(ii, jj);

                    n = citerations.2;
                    c = (citerations.0, citerations.1);

                    sum += n;

                    k += 1.0 / (self.s_max as f64);
                }

                sum = sum / self.s_max;
                n = sum;

                let color = self.get_color(n);

                self.iteration_counts[i as usize].push(n);
                self.iteration_values[i as usize].push(c);
                self.pixel_colours[i as usize].push(color);

                j += 1;
            }

            i += 1;
        }
    }

    fn write_ppm(&self) {
        print!("P3\n{} {}\n255\n", self.width, self.height);

        let mut i: u32 = 0;
        let mut j: u32 = 0;

        let mut count = 0;

        while j < self.height as u32 {
            i = 0;
            while i < self.width as u32 {
                count += 1;

                let c = &self.pixel_colours[i as usize][j as usize];
                print!("{} {} {}\n", c.R, c.G, c.B);
                i += 1;
            }
            j += 1;
        }
    }
}

fn main() {
    // let complex_integer = num::complex::Complex::new(10, 20);
    // let complex_float = num::complex::Complex::new(10.1, 20.1);

    // println!("Complex integer: {}", complex_integer);
    // println!("Complex float: {}", complex_float);

    let width = 1000;
    let height = 1000;

    let output_start = -2.0;
    let output_end = 2.0;

    let factor = 1.0;

    let n_max = 64;
    let s_max = 4;

    let mut mandelbrot = Mandelbrot::new(
        width,
        height,
        output_start,
        output_end,
        factor,
        n_max,
        s_max,
    );

    mandelbrot.render();
    mandelbrot.write_ppm();
}
