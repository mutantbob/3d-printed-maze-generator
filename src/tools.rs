use crate::Point3D;
use std::f32::consts::{PI, TAU};
use std::hash::Hash;

pub trait Topology<'a, CA> {
    type IterNeighbors: Iterator<Item = CA>;
    type IterAll: Iterator<Item = CA>;
    type IterWall: Iterator<Item = CA>;

    fn maximum_x(&self) -> f32;
    fn maximum_y(&self) -> f32;

    fn neighbors(&'a self, anchor: &CA) -> Self::IterNeighbors;

    fn wall_neighbors(&'a self, anchor: &CA) -> Self::IterWall;

    fn all_cells(&'a self) -> Self::IterAll;
}

pub trait Wall {}

pub trait CellAddress: Copy {
    fn coords_2d(&self) -> (f32, f32);
}

pub trait Space<C> {
    fn lerp(&self, a: &C, t: f32, b: &C) -> C;
    fn midpoint(&self, a: C, b: C) -> C;
    fn midpoint3(&self, a: C, b: C, c: C) -> C;
    // fn to_blender(&self, p: C) -> Point3D;

    fn subtract(&self, p1: C, p2: C) -> C;
}

pub trait BlenderMapping<COORD> {
    fn to_blender(&self, cc: COORD) -> Point3D;
}

pub trait EdgeCornerMapping<CA> {
    fn coord_left(&self, space: &dyn Space<(f32, f32)>) -> (f32, f32);
    fn coord_right(&self, space: &dyn Space<(f32, f32)>) -> (f32, f32);
}

//

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CylindricalCoodinate {
    pub rho: f32,
    pub r: f32,
    pub z: f32,
}

impl CylindricalCoodinate {
    pub fn new(rho: f32, r: f32, z: f32) -> Self {
        CylindricalCoodinate { rho, r, z }
    }
}

#[cfg(test)]
impl approx::AbsDiffEq for CylindricalCoodinate {
    type Epsilon = f32;

    fn default_epsilon() -> Self::Epsilon {
        1.0e-5
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.rho.abs_diff_eq(&other.rho, epsilon)
            && self.r.abs_diff_eq(&other.r, epsilon)
            && self.z.abs_diff_eq(&other.z, epsilon)
    }
}

//

pub struct Edge<CA>(pub CA, pub CA);

impl<CA: PartialEq<CA>> PartialEq<Self> for Edge<CA> {
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

impl<CA: PartialEq<CA>> Eq for Edge<CA> {}

//

#[derive(Debug)]
pub struct MazeWall<CA: CellAddress> {
    pub a: CA,
    pub b: CA,
    pub wall_ccw: bool,
    pub wall_cw: bool,
    pub wall_all: bool,
}

impl<CA: CellAddress> MazeWall<CA> {
    pub fn new(a: CA, b: CA, wall_ccw: bool, wall_cw: bool, wall_all: bool) -> Self {
        MazeWall {
            a,
            b,
            wall_ccw,
            wall_cw,
            wall_all,
        }
    }

    fn edge(&self) -> Edge<CA> {
        Edge(self.a, self.b)
    }
}

impl<CA: CellAddress> MazeWall<CA>
where
    Edge<CA>: EdgeCornerMapping<CA>,
{
    pub fn coord_left(&self, space: &dyn Space<(f32, f32)>) -> (f32, f32) {
        self.edge().coord_left(space)
    }

    pub fn coord_right(&self, space: &dyn Space<(f32, f32)>) -> (f32, f32) {
        self.edge().coord_right(space)
    }
}

//

//

pub struct CylindricalSpace {
    pub r0: f32,
    pub max_rho: f32,
}

impl Space<(f32, f32)> for CylindricalSpace {
    fn lerp(&self, a: &(f32, f32), t: f32, b: &(f32, f32)) -> (f32, f32) {
        let rho1 = a.0;
        let mut rho2 = b.0;
        let old = rho2;
        self.harmonize_angle(rho1, &mut rho2);
        if old != rho2 {
            println!("harmonized {} to {}", old, rho2);
        }
        let rho = lerp(rho1, t, rho2);
        let z = lerp(a.1, t, b.1);
        (rho, z)
    }

    fn midpoint(&self, (rho1, z1): (f32, f32), (mut rho2, z2): (f32, f32)) -> (f32, f32) {
        if true {
            return self.lerp(&(rho1, z1), 0.5, &(rho2, z2));
        }
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
    pub fn harmonize_angle(&self, rho1: f32, rho2: &mut f32) {
        let theta1 = rho1 / self.max_rho;
        let theta2 = *rho2 / self.max_rho;
        if theta1 - theta2 > 0.5 {
            *rho2 += self.max_rho;
        } else if theta2 - theta1 > 0.5 {
            *rho2 -= self.max_rho;
        }
    }

    pub fn scale_z(&self, z: f32) -> f32 {
        z * 2.0 * PI * self.r0 / self.max_rho
    }
}

impl BlenderMapping<CylindricalCoodinate> for CylindricalSpace {
    fn to_blender(&self, cc: CylindricalCoodinate) -> Point3D {
        let r0 = self.r0;
        let max_rho = self.max_rho;

        let theta = cc.rho / max_rho * TAU;
        let x = theta.cos() * (cc.r + r0);
        let y = theta.sin() * (cc.r + r0);

        let z = self.scale_z(cc.z);

        [x, y, z]
    }
}

//

#[derive(Debug, PartialEq, Eq, Hash)]
pub struct BidirectionalEdge {
    pub idx1: usize,
    pub idx2: usize,
}

impl BidirectionalEdge {
    pub fn new(idx1: usize, idx2: usize) -> Self {
        if idx1 <= idx2 {
            BidirectionalEdge { idx1, idx2 }
        } else {
            BidirectionalEdge {
                idx1: idx2,
                idx2: idx1,
            }
        }
    }
}

//
//
//

pub fn with_r(xz: (f32, f32), r: f32) -> CylindricalCoodinate {
    CylindricalCoodinate {
        rho: xz.0,
        z: xz.1,
        r,
    }
}

pub fn lerp(a: f32, t: f32, b: f32) -> f32 {
    a * (1.0 - t) + b * t
}

pub fn coord_left(
    v1: (f32, f32),
    v2: (f32, f32),
    frac: f32,
    space: &dyn Space<(f32, f32)>,
) -> (f32, f32) {
    let (dx, dy) = space.subtract(v2, v1);

    let x3 = v1.0 + dx * 0.5 - dy * frac;
    let y3 = v1.1 + dy * 0.5 + dx * frac;
    (x3, y3)
}
