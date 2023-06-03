#[derive(Clone)]
pub struct Color {
    pub R: u32,
    pub G: u32,
    pub B: u32,
}

impl Color {
    pub fn new(R: u32, G: u32, B: u32) -> Self {
        Color { R, G, B }
    }
}
