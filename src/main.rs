use lazy_static::lazy_static;
use std::fs::File;
use std::io::Write;

mod maze;
#[cfg(test)]
mod test;

lazy_static! {
    static ref TAN_30: f32 = 1.0f32 / 3.0f32.sqrt();
    static ref SEC_30: f32 = 2.0f32 / 3.0f32.sqrt();
}

fn main() {
    match 2 {
        2 => {
            let _ = draw_maze("/tmp/x.svg");
        }
        _ => {
            let _ = svg_check();
        }
    }
}

fn draw_maze(fname: &str) -> Result<(), std::io::Error> {
    let generator = maze::MazeGenerator::new();
    let edges = generator.generate_edges(HexCellAddress::new(0, 0));

    let mut f = File::create(fname)?;
    write!(&mut f, "<svg>\n")?;

    write!(&mut f, "<path d=\"")?;
    for edge in edges {
        let (x1, y1) = edge.0.coords_2d();
        let (x2, y2) = edge.1.coords_2d();
        write!(&mut f, "M {},{} {},{} ", x1, y1, x2, y2)?;
    }
    writeln!(
        &mut f,
        "\" style=\"stroke:#00f; fill:none; stroke-width:0.1px\"/>"
    )?;
    writeln!(&mut f, "</svg>")?;
    println!("wrote {}", fname);

    Ok(())
}

fn svg_check() -> Result<(), std::io::Error> {
    let fname = "/tmp/x.svg";
    let mut f = File::create(fname)?;
    write!(&mut f, "<svg>\n")?;
    for u in 0..10 {
        for v in 0..10 {
            let cell = HexCellAddress::new(u, v);
            let (x, y) = cell.coords_2d();
            writeln!(
                &mut f,
                "<circle cx=\"{}\" cy=\"{}\" r=\"0.1\" style=\"fill:#000\"/>",
                x, y
            )?;

            for n in cell.neighbors() {
                let (x9, y9) = n.coords_2d();
                let dx = x9 - x;
                let dy = y9 - y;

                let side = 0.05;
                let x1 = lerp(x, 0.1, x9) + side * dy;
                let y1 = lerp(y, 0.1, y9) - side * dx;
                let x8 = lerp(x, 0.9, x9) + side * dy;
                let y8 = lerp(y, 0.9, y9) - side * dx;

                writeln!(&mut f, "<path d=\"M {},{} {},{}\" style=\"fill:none; stroke:#00f; stroke-width: 0.1px\"/>",
                         x1, y1, x8, y8)?;
            }
        }
    }
    writeln!(&mut f, "</svg>")?;
    println!("wrote {}", fname);

    Ok(())
}

fn lerp(a: f32, t: f32, b: f32) -> f32 {
    a * (1.0 - t) + b * t
}

#[derive(Eq, Hash, PartialEq, Debug, Copy, Clone)]
pub struct HexCellAddress {
    pub u: i32,
    pub v: i32,
}

impl HexCellAddress {
    pub fn new(u: i32, v: i32) -> Self {
        HexCellAddress { u, v }
    }

    pub fn coords_2d(&self) -> (f32, f32) {
        let x = self.u as f32;
        let y = self.u as f32 * *TAN_30 + self.v as f32 * *SEC_30;
        (x, y)
    }

    pub fn neighbors(&self) -> [HexCellAddress; 6] {
        [
            HexCellAddress::new(self.u, self.v + 1),
            HexCellAddress::new(self.u, self.v - 1),
            HexCellAddress::new(self.u + 1, self.v - 1),
            HexCellAddress::new(self.u + 1, self.v),
            HexCellAddress::new(self.u - 1, self.v),
            HexCellAddress::new(self.u - 1, self.v + 1),
        ]
    }
}
