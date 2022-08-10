use crate::tools::coord_left;
use crate::{CellAddress, Edge, EdgeCornerMapping, Space, Topology};
use lazy_static::lazy_static;
use std::cmp::Ordering;

lazy_static! {
    static ref TAN_30: f32 = 1.0f32 / 3.0f32.sqrt();
    static ref SEC_30: f32 = 2.0f32 / 3.0f32.sqrt();
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
            return coord_left(v1, v2, -frac, space); // try this, later
        }
        let (dx, dy) = space.subtract(v2, v1);
        let x3 = v1.0 + dx * 0.5 + dy * frac;
        let y3 = v1.1 + dy * 0.5 - dx * frac;
        (x3, y3)
    }
}

//

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

//

pub type HexMazeEdge = Edge<HexCellAddress>;
