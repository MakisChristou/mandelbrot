use crate::{args::Bounds, color::Color};

pub struct Config {
    pub width: usize,
    pub height: usize,
    pub bounds: Bounds,
    pub factor: f64,
    pub n_max: u32,
    pub s_max: u32,
    pub color_pallete: Option<Vec<Color>>,
}
