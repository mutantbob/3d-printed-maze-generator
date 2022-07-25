use crate::blender_geometry::{BlenderGeometry, Point3D};
use lazy_static::lazy_static;
use std::cmp::Ordering;
use std::fs::File;
use std::io::Write;

mod blender_geometry;
mod experiments;
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
            let _ = experiments::check_blender_math("/tmp/x.py");
        }
        _ => {
            let _ = experiments::svg_check();
        }
    }
}

pub fn with_z(xy: (f32, f32), z: f32) -> Point3D {
    [xy.0, xy.1, z]
}

pub struct MazeTopology1 {
    after_max_u: i32,
}

impl MazeTopology1 {
    pub(crate) fn neighbors<'a>(
        &'a self,
        anchor: &HexCellAddress,
    ) -> impl Iterator<Item = HexCellAddress> + 'a {
        anchor
            .neighbors()
            .into_iter()
            .map(|n| self.wrap(n))
            .filter(|n| self.in_bounds(n))
    }
}

impl MazeTopology1 {
    pub(crate) fn all_cells<'a>(&'a self) -> impl Iterator<Item = HexCellAddress> + 'a {
        let mut min_v = 0;
        (0..self.after_max_u).flat_map(move |u| {
            if self.in_bounds(&HexCellAddress::new(u, min_v - 1)) {
                min_v -= 1;
            }
            (min_v..9999)
                .map(move |v| HexCellAddress::new(u, v))
                .take_while(|cell| self.in_bounds(cell))
        })
    }

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
        HexCellAddress::new(u, v)
    }

    fn in_bounds(&self, cell: &HexCellAddress) -> bool {
        let (_, y) = cell.coords_2d();
        cell.u >= 0 && cell.u < self.after_max_u && y >= 0. && y < 20.0
    }
}

pub fn draw_maze(fname: &str) -> Result<(), std::io::Error> {
    let generator = maze::MazeGenerator::new();
    let topology1 = MazeTopology1 { after_max_u: 20 };
    let edges: Vec<_> = generator
        .generate_edges(HexCellAddress::new(0, 0), |cell| topology1.neighbors(cell))
        .into_iter()
        .map(|(c1, c2)| HexMazeEdge(c1, c2))
        .collect();

    save_edges_svg(fname, &edges)?;

    let walls = compute_walls(&topology1, &edges);

    println!("{} edges; {} walls", edges.len(), walls.len());

    write_blender_python("/tmp/geom.py", &edges, &walls)?;

    Ok(())
}

fn compute_walls(topology: &MazeTopology1, corridors: &[HexMazeEdge]) -> Vec<HexMazeEdge> {
    let cells: Vec<_> = topology.all_cells().collect();
    println!("{} cells", cells.len());
    let mut rval = vec![];
    for cell in cells {
        let neighbors: Vec<_> = topology.neighbors(&cell).collect();
        println!("{:?} has {} neighbors", &cell, neighbors.len());
        for n in neighbors {
            if cell.cmp(&n) == Ordering::Less {
                let edge = HexMazeEdge(cell, n);
                if !corridors.iter().any(|old| *old == edge) {
                    //that edge does not have a corridor
                    rval.push(edge)
                }
            }
        }
    }

    rval
}

pub fn write_blender_python(
    fname: &str,
    edges: &[HexMazeEdge],
    walls: &[HexMazeEdge],
) -> Result<(), std::io::Error> {
    let mut blender = BlenderGeometry::new();

    for edge in edges {
        add_edge_flat(&mut blender, edge, 0.3);
    }
    for wall in walls {
        add_wall_flat(&mut blender, wall, 0.3);
    }

    let mut f = File::create(fname)?;
    blender.emit(&mut f)?;

    println!("wrote {}", fname);
    Ok(())
}

pub fn save_edges_svg(fname: &str, edges: &[HexMazeEdge]) -> Result<(), std::io::Error> {
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

pub fn add_edge_flat(blender: &mut BlenderGeometry, edge: &HexMazeEdge, high_z: f32) {
    let v0 = with_z(edge.0.coords_2d(), 0.0);
    let v1 = with_z(edge.1.coords_2d(), 0.0);

    let v2 = with_z(edge.coord_left(), high_z);
    let v3 = with_z(edge.coord_right(), high_z);

    blender.add_face(&[v0, v1, v2]);
    blender.add_face(&[v1, v0, v3]);
}

pub fn add_wall_flat(blender: &mut BlenderGeometry, edge: &HexMazeEdge, high_z: f32) {
    let v0 = with_z(edge.0.coords_2d(), 0.0);
    let v1 = with_z(edge.1.coords_2d(), 0.0);

    let v2 = with_z(edge.coord_left(), high_z);
    let v3 = with_z(edge.coord_right(), high_z);

    blender.add_face(&[v0, v3, v2]);
    blender.add_face(&[v1, v2, v3]);
}

//

//

pub fn lerp(a: f32, t: f32, b: f32) -> f32 {
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

impl PartialOrd<Self> for HexCellAddress {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for HexCellAddress {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.u < other.u {
            Ordering::Less
        } else if self.u > other.u {
            Ordering::Greater
        } else if self.v < other.v {
            Ordering::Less
        } else if self.v > other.v {
            Ordering::Greater
        } else {
            Ordering::Equal
        }
    }
}

//

pub struct HexMazeEdge(HexCellAddress, HexCellAddress);

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

impl PartialEq<Self> for HexMazeEdge {
    fn eq(&self, other: &Self) -> bool {
        if self.0 == other.0 {
            self.1 == other.1
        } else if self.0 == other.1 {
            self.1 == other.0
        } else {
            false
        }
    }
}

impl Eq for HexMazeEdge {}
