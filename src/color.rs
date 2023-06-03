#[derive(Clone)]
pub struct Color {
    pub R: u8,
    pub G: u8,
    pub B: u8,
}

impl Color {
    pub fn new(R: u8, G: u8, B: u8) -> Self {
        Color { R, G, B }
    }
}
