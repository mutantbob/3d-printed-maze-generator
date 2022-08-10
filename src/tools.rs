use crate::{Point3D, SEC_30, TAN_30};
use std::cmp::Ordering;
use std::f32::consts::{PI, TAU};
use std::hash::Hash;

pub fn with_r(xz: (f32, f32), r: f32) -> CylindricalCoodinate {
    CylindricalCoodinate {
        rho: xz.0,
        z: xz.1,
        r,
    }
}

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

//

pub struct Edge<CA>(pub CA, pub CA);

pub struct CylindricalCoodinate {
    pub rho: f32,
    pub r: f32,
    pub z: f32,
}

pub struct HexMazeTopology {
    pub after_max_u: i32,
    pub max_y: f32,
}

impl HexMazeTopology {
    pub fn new(u_count: u32, v_count: u32) -> Self {
        HexMazeTopology {
            after_max_u: u_count as i32,
            max_y: 0.1 + (v_count as f32) * *SEC_30,
        }
    }

    pub fn wrap(&self, pre: HexCellAddress) -> HexCellAddress {
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
            && y >= -0.1 - *SEC_30 * 2.0
            && y < self.max_y + *SEC_30 * 1.0 + 0.1
    }
}

impl<'a> Topology<'a, HexCellAddress> for HexMazeTopology {
    type IterNeighbors = Box<dyn Iterator<Item = HexCellAddress> + 'a>;
    type IterAll = Box<dyn Iterator<Item = HexCellAddress> + 'a>;
    type IterWall = Box<dyn Iterator<Item = HexCellAddress> + 'a>;

    fn maximum_x(&self) -> f32 {
        self.after_max_u as f32
    }

    fn maximum_y(&self) -> f32 {
        self.max_y
    }

    fn neighbors(&'a self, anchor: &HexCellAddress) -> Self::IterNeighbors {
        Box::new(
            anchor
                .neighbors()
                .into_iter()
                .map(|n| self.wrap(n))
                .filter(|n| self.in_bounds(n)),
        )
    }

    fn wall_neighbors(&'a self, anchor: &HexCellAddress) -> Self::IterWall {
        Box::new(
            anchor
                .neighbors()
                .into_iter()
                .map(|n| self.wrap(n))
                .filter(|n| self.wall_bounds(n)),
        )
    }

    #[allow(clippy::needless_lifetimes)] // clippy is wrong.  removing the lifetime triggers an error in my version of rust
    fn all_cells(&'a self) -> Self::IterAll {
        let mut min_v = 0;
        Box::new((0..self.after_max_u).flat_map(move |u| {
            if self.wall_bounds(&HexCellAddress::new(u, min_v - 1)) {
                min_v -= 1;
            }
            (min_v..9999)
                .map(move |v| HexCellAddress::new(u, v))
                .take_while(|cell| self.wall_bounds(cell))
        }))
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

impl CellAddress for HexCellAddress {
    fn coords_2d(&self) -> (f32, f32) {
        let x = self.u as f32;
        let y = self.u as f32 * *TAN_30 + self.v as f32 * *SEC_30;
        (x, y)
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

pub type HexMazeEdge = Edge<HexCellAddress>;

pub trait EdgeCornerMapping<CA> {
    fn coord_left(&self, space: &dyn Space<(f32, f32)>) -> (f32, f32);
    fn coord_right(&self, space: &dyn Space<(f32, f32)>) -> (f32, f32);
}

impl EdgeCornerMapping<HexCellAddress> for HexMazeEdge {
    fn coord_left(&self, space: &dyn Space<(f32, f32)>) -> (f32, f32) {
        let frac = 0.5 * 0.5 / 0.75_f32.sqrt();
        let v1 = self.0.coords_2d();
        let v2 = self.1.coords_2d();

        coord_left(v1, v2, frac, space)
    }

    fn coord_right(&self, space: &dyn Space<(f32, f32)>) -> (f32, f32) {
        let frac = 0.5 * 0.5 / 0.75_f32.sqrt();
        let v1 = self.0.coords_2d();
        let v2 = self.1.coords_2d();
        if false {
            return coord_left(v1, v2, -frac, space);
        }
        let (dx, dy) = space.subtract(v2, v1);
        let x3 = v1.0 + dx * 0.5 + dy * frac;
        let y3 = v1.1 + dy * 0.5 - dx * frac;
        (x3, y3)
    }
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

pub trait Space<C> {
    fn midpoint(&self, a: C, b: C) -> C;
    fn midpoint3(&self, a: C, b: C, c: C) -> C;
    // fn to_blender(&self, p: C) -> Point3D;

    fn subtract(&self, p1: C, p2: C) -> C;
}

pub trait BlenderMapping<COORD> {
    fn to_blender(&self, cc: COORD) -> Point3D;
}

//

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
    pub fn absorb(&mut self, p1: Point3D, p2: Point3D) {
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

//

//

#[derive(Copy, Clone, PartialEq, Eq, Debug, Hash)]
pub struct SqCellAddress {
    pub u: i32,
    pub v: i32,
}

impl SqCellAddress {
    pub fn new(u: i32, v: i32) -> Self {
        SqCellAddress { u, v }
    }

    pub fn neighbors(&self) -> [SqCellAddress; 4] {
        let SqCellAddress { u, v } = *self;
        [
            SqCellAddress::new(u - 1, v),
            SqCellAddress::new(u, v - 1),
            SqCellAddress::new(u + 1, v),
            SqCellAddress::new(u, v + 1),
        ]
    }
}

impl CellAddress for SqCellAddress {
    fn coords_2d(&self) -> (f32, f32) {
        (self.u as f32, self.v as f32)
    }
}

impl EdgeCornerMapping<SqCellAddress> for Edge<SqCellAddress> {
    fn coord_left(&self, space: &dyn Space<(f32, f32)>) -> (f32, f32) {
        coord_left(self.0.coords_2d(), self.1.coords_2d(), 0.5, space)
    }

    fn coord_right(&self, space: &dyn Space<(f32, f32)>) -> (f32, f32) {
        coord_left(self.0.coords_2d(), self.1.coords_2d(), -0.5, space)
    }
}

//

pub struct SquareMazeTopology {
    pub u_count: u32,
    pub v_count: u32,
}

impl SquareMazeTopology {
    pub fn new(u_count: u32, v_count: u32) -> Self {
        SquareMazeTopology { u_count, v_count }
    }

    pub fn wrap(&self, pre: SqCellAddress) -> SqCellAddress {
        if false {
            return pre;
        }

        let mut u = pre.u;
        let v = pre.v;
        while u < 0 {
            u += self.u_count as i32;
        }
        while u >= self.u_count as i32 {
            u -= self.u_count as i32;
        }
        SqCellAddress::new(u, v)
    }

    fn in_bounds(&self, cell: &SqCellAddress) -> bool {
        cell.u >= 0 && cell.u < self.u_count as i32 && cell.v >= 0 && cell.v < self.v_count as i32
    }

    fn wall_bounds(&self, cell: &SqCellAddress) -> bool {
        cell.u >= 0
            && cell.u < self.u_count as i32
            && cell.v >= -1
            && cell.v < 1 + self.v_count as i32
    }
}

impl<'a> Topology<'a, SqCellAddress> for SquareMazeTopology {
    type IterNeighbors = Box<dyn Iterator<Item = SqCellAddress> + 'a>;
    type IterAll = Box<dyn Iterator<Item = SqCellAddress> + 'a>;
    type IterWall = Box<dyn Iterator<Item = SqCellAddress> + 'a>;

    fn maximum_x(&self) -> f32 {
        self.u_count as f32
    }

    fn maximum_y(&self) -> f32 {
        self.v_count as f32
    }

    fn neighbors(&'a self, anchor: &SqCellAddress) -> Self::IterNeighbors {
        Box::new(
            anchor
                .neighbors()
                .into_iter()
                .map(|n| self.wrap(n))
                .filter(|n| self.in_bounds(n)),
        )
    }

    fn wall_neighbors(&'a self, anchor: &SqCellAddress) -> Self::IterWall {
        Box::new(
            anchor
                .neighbors()
                .into_iter()
                .map(|n| self.wrap(n))
                .filter(|n| self.wall_bounds(n)),
        )
    }

    #[allow(clippy::needless_lifetimes)] // clippy is wrong.  removing the lifetime triggers an error in my version of rust
    fn all_cells(&'a self) -> Self::IterAll {
        Box::new((0..self.u_count).flat_map(move |u| {
            (0..self.v_count).map(move |v| SqCellAddress::new(u as i32, v as i32))
        }))
    }
}
/*
struct SMTNeighbors<'a> {
    base: Box<dyn Iterator<Item = SqCellAddress> + 'a>,
}

impl<'a> SMTNeighbors<'a> {
    pub fn new(topology: &'a SquareMazeTopology, anchor: &'a SqCellAddress) -> Self {
        SMTNeighbors {
            base: Box::new(
                anchor
                    .neighbors()
                    .into_iter()
                    .map(|n| topology.wrap(n))
                    .filter(|n| topology.in_bounds(n)),
            ),
        }
    }
}

impl<'a> Iterator for SMTNeighbors<'a> {
    type Item = SqCellAddress;

    fn next(&mut self) -> Option<Self::Item> {
        self.base.next()
    }
}
*/
