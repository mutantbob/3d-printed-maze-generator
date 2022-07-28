use crate::{Point3D, SEC_30, TAN_30};
use std::cmp::Ordering;
use std::f32::consts::{PI, TAU};

pub fn with_z(xy: (f32, f32), z: f32) -> Point3D {
    [xy.0, xy.1, z]
}

pub struct MazeTopology1 {
    pub after_max_u: i32,
    pub max_y: f32,
}

impl MazeTopology1 {
    pub fn new(u_count: u32, v_count: u32) -> Self {
        MazeTopology1 {
            after_max_u: u_count as i32,
            max_y: 0.1 + (v_count as f32) * *SEC_30,
        }
    }

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
        cell.u >= 0 && cell.u < self.after_max_u && y >= 0. && y < self.max_y
    }

    fn wall_bounds(&self, cell: &HexCellAddress) -> bool {
        let (_, y) = cell.coords_2d();
        cell.u >= 0
            && cell.u < self.after_max_u
            && y >= -0.1 - *SEC_30
            && y < self.max_y + *SEC_30 + 0.1
    }
}

pub fn lerp(a: f32, t: f32, b: f32) -> f32 {
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

pub struct HexMazeEdge(pub HexCellAddress, pub HexCellAddress);

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

#[derive(Debug)]
pub struct HexMazeWall {
    pub a: HexCellAddress,
    pub b: HexCellAddress,
    pub wall_ccw: bool,
    pub wall_cw: bool,
    pub wall_all: bool,
}

impl HexMazeWall {
    pub fn new(
        a: HexCellAddress,
        b: HexCellAddress,
        wall_ccw: bool,
        wall_cw: bool,
        wall_all: bool,
    ) -> Self {
        HexMazeWall {
            a,
            b,
            wall_ccw,
            wall_cw,
            wall_all,
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

    pub fn to_blender(&self, [rho1, z, r]: Point3D) -> Point3D {
        let r0 = self.r0;
        let max_rho = self.max_rho;

        let theta = rho1 / max_rho * TAU;
        let x = theta.cos() * (r + r0);
        let y = theta.sin() * (r + r0);

        let z = self.scale_z(z);

        [x, y, z]
    }
}

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

enum RingCommand {
    Append(Point3D),
    Unshift(Point3D),
    Close,
    Nope,
}

#[derive(Default)]
pub struct RingAccumulator {
    pub open_strings: Vec<Vec<Point3D>>,
    pub closed_strings: Vec<Vec<Point3D>>,
}

impl RingAccumulator {
    pub(crate) fn absorb(&mut self, p1: Point3D, p2: Point3D) {
        for idx in 0..self.open_strings.len() {
            let string = &mut self.open_strings[idx];
            let cmd = Self::command_for(&p1, &p2, string);
            match cmd {
                RingCommand::Append(p) => {
                    string.push(p);
                    self.consolidate_open_strings();
                    return;
                }
                RingCommand::Unshift(p) => {
                    string.insert(0, p);
                    self.consolidate_open_strings();
                    return;
                }
                RingCommand::Close => {
                    self.closed_strings.push(self.open_strings.remove(idx));
                    return;
                }
                RingCommand::Nope => {}
            }
        }

        self.open_strings.push(vec![p1, p2]);
    }

    fn consolidate_open_strings(&mut self) {
        for i in 0..self.open_strings.len() {
            let a = self.open_strings[i].first().unwrap();
            let b = self.open_strings[i].last().unwrap();
            for j in (i + 1)..self.open_strings.len() {
                let n = self.open_strings[j].first().unwrap();
                let o = self.open_strings[j].last().unwrap();

                if Self::close_enough(a, n) {
                    let mut doomed = self.open_strings.remove(j);
                    doomed.remove(0);
                    self.open_strings[i].splice(0..0, doomed.into_iter().rev());
                    self.consolidate_open_strings();
                    return;
                } else if Self::close_enough(b, n) {
                    let mut doomed = self.open_strings.remove(j);
                    doomed.remove(0);
                    self.open_strings[i].extend(doomed);
                    self.consolidate_open_strings();
                    return;
                } else if Self::close_enough(a, o) {
                    let mut doomed = self.open_strings.remove(j);
                    doomed.pop();
                    self.open_strings[i].splice(0..0, doomed);
                    self.consolidate_open_strings();
                    return;
                } else if Self::close_enough(b, o) {
                    let mut doomed = self.open_strings.remove(j);
                    doomed.pop();
                    self.open_strings[i].extend(doomed.into_iter().rev());
                    self.consolidate_open_strings();
                    return;
                }
            }
        }
    }

    fn command_for(p1: &Point3D, p2: &Point3D, string: &[Point3D]) -> RingCommand {
        let begin = string.first().unwrap();
        let end = string.last().unwrap();
        if Self::close_enough(begin, p1) {
            if Self::close_enough(end, p2) {
                RingCommand::Close
            } else {
                RingCommand::Unshift(*p2)
            }
        } else if Self::close_enough(end, p1) {
            if Self::close_enough(begin, p2) {
                RingCommand::Close
            } else {
                RingCommand::Append(*p2)
            }
        } else if Self::close_enough(begin, p2) {
            if Self::close_enough(end, p1) {
                RingCommand::Close
            } else {
                RingCommand::Unshift(*p1)
            }
        } else if Self::close_enough(end, p2) {
            if Self::close_enough(begin, p1) {
                RingCommand::Close
            } else {
                RingCommand::Append(*p1)
            }
        } else {
            RingCommand::Nope
        }
    }

    fn close_enough(p1: &Point3D, p2: &Point3D) -> bool {
        let delta = p1
            .iter()
            .zip(p2.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, |a, b| if a < b { b } else { a });
        delta < 1.0e-6
    }
}
