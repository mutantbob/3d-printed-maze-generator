use crate::blender_geometry::{BlenderGeometry, Point3D};
use lazy_static::lazy_static;
use std::fs::File;
use std::io::Write;

mod blender_geometry;
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
        3 => {
            let _ = check_blender_math("/tmp/x.py");
        }
        _ => {
            let _ = svg_check();
        }
    }
}

pub fn with_z(xy: (f32, f32), z: f32) -> Point3D {
    [xy.0, xy.1, z]
}

struct MazeTopology1 {
    after_max_u: i32,
}

impl MazeTopology1 {
    pub(crate) fn wrap(&self, pre: HexCellAddress) -> HexCellAddress {
        if true {
            return pre;
        }

        let mut u = pre.u;
        let mut v = pre.v;
        while u < 0 {
            u += self.after_max_u;
            v -= self.after_max_u / 2;
        }
        while u >= self.after_max_u {
            u -= self.after_max_u;
            v += self.after_max_u / 2;
        }
        return HexCellAddress::new(u, v);
    }
}

impl MazeTopology1 {
    fn in_bounds(&self, cell: &HexCellAddress) -> bool {
        let (_, y) = cell.coords_2d();
        cell.u >= 0 && cell.u < self.after_max_u && y >= 0. && y < 20.0
    }
}

fn draw_maze(fname: &str) -> Result<(), std::io::Error> {
    let generator = maze::MazeGenerator::new();
    let topology1 = MazeTopology1 { after_max_u: 20 };
    let edges: Vec<_> = generator
        .generate_edges(HexCellAddress::new(0, 0), |cell| {
            cell.neighbors()
                .into_iter()
                .map(|n| topology1.wrap(n))
                .filter(|n| topology1.in_bounds(n))
        })
        .into_iter()
        .map(|(c1, c2)| HexMazeEdge(c1, c2))
        .collect();

    save_edges_svg(fname, &edges)?;

    write_blender_python("/tmp/geom.py", &edges)?;

    Ok(())
}

fn write_blender_python(fname: &str, edges: &[HexMazeEdge]) -> Result<(), std::io::Error> {
    let mut blender = BlenderGeometry::new();

    for edge in edges {
        add_edge_flat(&mut blender, edge, 0.3);
    }

    let mut f = File::create(fname)?;
    blender.emit(&mut f)?;

    println!("wrote {}", fname);
    Ok(())
}

fn save_edges_svg(fname: &str, edges: &[HexMazeEdge]) -> Result<(), std::io::Error> {
    let mut f = File::create(fname)?;
    writeln!(&mut f, "<svg>")?;

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

//

fn check_blender_math(fname: &str) -> Result<(), std::io::Error> {
    let mut f = File::create(fname)?;

    let mut blender = BlenderGeometry::new();

    let edge = HexMazeEdge(HexCellAddress::new(0, 0), HexCellAddress::new(1, 0));

    add_edge_flat(&mut blender, &edge, 1.0);

    blender.emit(&mut f)?;

    Ok(())
}

fn add_edge_flat(blender: &mut BlenderGeometry, edge: &HexMazeEdge, high_z: f32) {
    let v0 = with_z(edge.0.coords_2d(), 0.0);
    let v1 = with_z(edge.1.coords_2d(), 0.0);

    let v2 = with_z(edge.coord_left(), high_z);
    let v3 = with_z(edge.coord_right(), high_z);

    blender.add_face(&[v0, v1, v2]);
    blender.add_face(&[v1, v0, v3]);
}

//

fn svg_check() -> Result<(), std::io::Error> {
    let fname = "/tmp/x.svg";
    let mut f = File::create(fname)?;
    writeln!(&mut f, "<svg>")?;
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

//

fn lerp(a: f32, t: f32, b: f32) -> f32 {
    a * (1.0 - t) + b * t
}

//

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

//

struct HexMazeEdge(HexCellAddress, HexCellAddress);

impl HexMazeEdge {
    pub fn coord_left(&self) -> (f32, f32) {
        let v1 = self.0.coords_2d();
        let v2 = self.1.coords_2d();
        let frac = 0.5 * 0.5 / 0.75_f32.sqrt();

        let dx = v2.0 - v1.0;
        let dy = v2.1 - v1.1;
        let x3 = v1.0 + dx * 0.5 - dy * frac;
        let y3 = v1.1 + dy * 0.5 + dx * frac;
        (x3, y3)
    }

    pub fn coord_right(&self) -> (f32, f32) {
        let v1 = self.0.coords_2d();
        let v2 = self.1.coords_2d();
        let frac = 0.5 * 0.5 / 0.75_f32.sqrt();

        let dx = v2.0 - v1.0;
        let dy = v2.1 - v1.1;
        let x3 = v1.0 + dx * 0.5 + dy * frac;
        let y3 = v1.1 + dy * 0.5 - dx * frac;
        (x3, y3)
    }
}
