use crate::blender_geometry::{BlenderGeometry, Point3D};
use lazy_static::lazy_static;
use std::cmp::Ordering;
use std::f32::consts::TAU;
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

    pub(crate) fn wall_neighbors<'a>(
        &'a self,
        anchor: &HexCellAddress,
    ) -> impl Iterator<Item = HexCellAddress> + 'a {
        anchor
            .neighbors()
            .into_iter()
            .map(|n| self.wrap(n))
            .filter(|n| self.wall_bounds(n))
    }

    #[allow(clippy::needless_lifetimes)] // clippy is wrong.  removing the lifetime triggers an error in my version of rust
    pub(crate) fn all_cells<'a>(&'a self) -> impl Iterator<Item = HexCellAddress> + 'a {
        let mut min_v = 0;
        (0..self.after_max_u).flat_map(move |u| {
            if self.wall_bounds(&HexCellAddress::new(u, min_v - 1)) {
                min_v -= 1;
            }
            (min_v..9999)
                .map(move |v| HexCellAddress::new(u, v))
                .take_while(|cell| self.wall_bounds(cell))
        })
    }

    pub(crate) fn wrap(&self, pre: HexCellAddress) -> HexCellAddress {
        if false {
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

    fn wall_bounds(&self, cell: &HexCellAddress) -> bool {
        let (_, y) = cell.coords_2d();
        cell.u >= 0 && cell.u < self.after_max_u && y >= -1. && y < 21.0
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

    write_blender_python("/tmp/geom.py", &edges, &walls, topology1.after_max_u as f32)?;

    Ok(())
}

fn compute_walls(topology: &MazeTopology1, corridors: &[HexMazeEdge]) -> Vec<HexMazeWall> {
    let cells: Vec<_> = topology.all_cells().collect();
    println!("{} cells", cells.len());
    let mut rval = vec![];
    for cell in cells {
        let neighbors: Vec<_> = topology.wall_neighbors(&cell).collect();
        // println!("{:?} has {} neighbors", &cell, neighbors.len());
        let mut directions = vec![];
        let mut is_wall = vec![];
        for n in neighbors {
            let (edge, wallness) = {
                let edge = HexMazeEdge(cell, n);
                if !corridors.iter().any(|old| *old == edge) {
                    //that edge does not have a corridor
                    (Some(edge), true)
                } else {
                    (None, false)
                }
            };
            directions.push(edge);
            is_wall.push(wallness);
        }

        let n = directions.len();
        for (i, wall) in directions.into_iter().enumerate() {
            if let Some(wall) = wall {
                let wall_ccw = is_wall[(i + n - 1) % n];
                let wall_cw = is_wall[(i + 1) % n];
                rval.push(HexMazeWall::new(wall.0, wall.1, wall_ccw, wall_cw));
            }
        }
    }

    rval
}

pub fn write_blender_python(
    fname: &str,
    edges: &[HexMazeEdge],
    walls: &[HexMazeWall],
    max_rho: f32,
) -> Result<(), std::io::Error> {
    let mut blender = BlenderGeometry::new();

    let cylindrical = CylindricalSpace { r0: 2.0, max_rho };
    for edge in edges {
        add_edge_flat(
            &mut blender,
            edge,
            0.3,
            |xyz| cylindrical.to_blender(xyz),
            &cylindrical,
        );
    }
    for wall in walls {
        add_wall_flat(&mut blender, wall, 0.3, &cylindrical);
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

pub fn add_edge_flat<F>(
    blender: &mut BlenderGeometry,
    edge: &HexMazeEdge,
    high_z: f32,
    mapping: F,
    space: &dyn Space<(f32, f32)>,
) where
    F: Fn(Point3D) -> Point3D,
{
    let v0 = mapping(with_z(edge.0.coords_2d(), 0.0));
    let v1 = mapping(with_z(edge.1.coords_2d(), 0.0));

    let v2 = mapping(with_z(edge.coord_left(space), high_z));
    let v3 = mapping(with_z(edge.coord_right(space), high_z));

    blender.add_face(&[v0, v1, v2]);
    blender.add_face(&[v1, v0, v3]);
}

pub fn add_wall_flat(
    blender: &mut BlenderGeometry,
    wall: &HexMazeWall,
    high_z: f32,
    space: &CylindricalSpace,
) {
    let low_z = 0.0;

    let xy0 = wall.a.coords_2d();
    let v0 = with_z(xy0, low_z);
    let xy2 = wall.coord_left(space);
    let v2 = with_z(xy2, high_z);
    let xy3 = wall.coord_right(space);
    let v3 = with_z(xy3, high_z);

    match (wall.wall_ccw, wall.wall_cw) {
        (false, false) => {
            let xy8 = space.midpoint3(xy0, xy2, xy3);
            let v8 = with_z(xy8, high_z);

            let v0 = space.to_blender(v0);
            let v2 = space.to_blender(v2);
            let v3 = space.to_blender(v3);
            let v8 = space.to_blender(v8);

            blender.add_face(&[v0, v3, v8]);
            blender.add_face(&[v3, v2, v8]);
            blender.add_face(&[v2, v0, v8]);
        }
        (true, false) => {
            let v4 = with_z(space.midpoint(xy0, xy2), high_z);

            let v0 = space.to_blender(v0);
            let v2 = space.to_blender(v2);
            let v3 = space.to_blender(v3);
            let v4 = space.to_blender(v4);

            if (v0[0] - v4[0]).abs() > 2.0 {
                println!("problem {:?}", &wall);
            }

            blender.add_face(&[v0, v3, v4]);
            blender.add_face(&[v4, v3, v2]);
        }
        (false, true) => {
            let v5 = with_z(space.midpoint(xy0, xy3), high_z);

            let v0 = space.to_blender(v0);
            let v2 = space.to_blender(v2);
            let v3 = space.to_blender(v3);
            let v5 = space.to_blender(v5);

            blender.add_face(&[v0, v5, v2]);
            blender.add_face(&[v5, v3, v2]);
        }
        (true, true) => {
            let v4 = with_z(space.midpoint(xy0, xy2), high_z);
            let v5 = with_z(space.midpoint(xy0, xy3), high_z);

            let v0 = space.to_blender(v0);
            let v2 = space.to_blender(v2);
            let v3 = space.to_blender(v3);
            let v4 = space.to_blender(v4);
            let v5 = space.to_blender(v5);

            blender.add_face(&[v0, v5, v4]);
            blender.add_face(&[v4, v5, v3, v2]);
        }
    }
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
            HexCellAddress::new(self.u + 1, self.v),
            HexCellAddress::new(self.u + 1, self.v - 1),
            HexCellAddress::new(self.u, self.v - 1),
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
    pub fn coord_left(&self, space: &dyn Space<(f32, f32)>) -> (f32, f32) {
        let v1 = self.0.coords_2d();
        let v2 = self.1.coords_2d();
        let frac = 0.5 * 0.5 / 0.75_f32.sqrt();

        let (dx, dy) = space.subtract(v2, v1);

        let x3 = v1.0 + dx * 0.5 - dy * frac;
        let y3 = v1.1 + dy * 0.5 + dx * frac;
        (x3, y3)
    }

    pub fn coord_right(&self, space: &dyn Space<(f32, f32)>) -> (f32, f32) {
        let v1 = self.0.coords_2d();
        let v2 = self.1.coords_2d();
        let frac = 0.5 * 0.5 / 0.75_f32.sqrt();

        let (dx, dy) = space.subtract(v2, v1);
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

//

#[derive(Debug)]
pub struct HexMazeWall {
    pub a: HexCellAddress,
    pub b: HexCellAddress,
    pub wall_ccw: bool,
    pub wall_cw: bool,
}

impl HexMazeWall {
    pub fn new(a: HexCellAddress, b: HexCellAddress, wall_ccw: bool, wall_cw: bool) -> Self {
        HexMazeWall {
            a,
            b,
            wall_ccw,
            wall_cw,
        }
    }

    fn edge(&self) -> HexMazeEdge {
        HexMazeEdge(self.a, self.b)
    }

    pub(crate) fn coord_left(&self, space: &dyn Space<(f32, f32)>) -> (f32, f32) {
        self.edge().coord_left(space)
    }

    pub(crate) fn coord_right(&self, space: &dyn Space<(f32, f32)>) -> (f32, f32) {
        self.edge().coord_right(space)
    }
}

//

pub trait Space<C> {
    fn midpoint(&self, a: C, b: C) -> C;
    fn midpoint3(&self, a: C, b: C, c: C) -> C;
    // fn to_blender(&self, p: C) -> Point3D;

    fn subtract(&self, p1: C, p2: C) -> C;
}

pub struct CylindricalSpace {
    pub r0: f32,
    pub max_rho: f32,
}

impl Space<(f32, f32)> for CylindricalSpace {
    fn midpoint(&self, (rho1, z1): (f32, f32), (mut rho2, z2): (f32, f32)) -> (f32, f32) {
        let old = rho2;
        self.harmonize_angle(rho1, &mut rho2);
        if old != rho2 {
            println!("harmonized {} to {}", old, rho2);
        }

        (0.5 * (rho1 + rho2), 0.5 * (z1 + z2))
    }

    fn midpoint3(
        &self,
        (rho1, z1): (f32, f32),
        (mut rho2, z2): (f32, f32),
        (mut rho3, z3): (f32, f32),
    ) -> (f32, f32) {
        self.harmonize_angle(rho1, &mut rho2);
        self.harmonize_angle(rho1, &mut rho3);
        ((rho1 + rho2 + rho3) / 3.0, (z1 + z2 + z3) / 3.0)
    }

    fn subtract(&self, p1: (f32, f32), p2: (f32, f32)) -> (f32, f32) {
        let rho1 = p1.0;
        let mut rho2 = p2.0;
        self.harmonize_angle(rho1, &mut rho2);

        (rho1 - rho2, p1.1 - p2.1)
    }
}

impl CylindricalSpace {
    fn harmonize_angle(&self, rho1: f32, rho2: &mut f32) {
        let theta1 = rho1 / self.max_rho;
        let theta2 = *rho2 / self.max_rho;
        if theta1 - theta2 > 0.5 {
            *rho2 += self.max_rho;
        } else if theta2 - theta1 > 0.5 {
            *rho2 -= self.max_rho;
        }
    }

    fn to_blender(&self, [rho1, z, r]: Point3D) -> Point3D {
        let r0 = self.r0;
        let max_rho = self.max_rho;

        let theta = rho1 / max_rho * TAU;
        let x = theta.cos() * (r + r0);
        let y = theta.sin() * (r + r0);

        [x, y, z]
    }
}
