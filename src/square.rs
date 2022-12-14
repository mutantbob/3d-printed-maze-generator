use crate::tools::coord_left;
use crate::{CellAddress, Edge, EdgeCornerMapping, Space, Topology};

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

#[derive(Clone)]
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
            && cell.v >= -2
            && cell.v < 1 + self.v_count as i32
    }
}

impl Topology<SqCellAddress> for SquareMazeTopology {
    fn maximum_x(&self) -> f32 {
        self.u_count as f32
    }

    fn maximum_y(&self) -> f32 {
        self.v_count as f32
    }

    fn neighbors(&self, anchor: &SqCellAddress) -> Box<dyn Iterator<Item = SqCellAddress>> {
        let clone1 = (*self).clone();
        let clone2 = (*self).clone();
        Box::new(
            anchor
                .neighbors()
                .into_iter()
                .map(move |n| clone1.wrap(n))
                .filter(move |n| clone2.in_bounds(n)),
        )
    }

    fn wall_neighbors(&self, anchor: &SqCellAddress) -> Box<dyn Iterator<Item = SqCellAddress>> {
        let clone1 = (*self).clone();
        let clone2 = (*self).clone();
        Box::new(
            anchor
                .neighbors()
                .into_iter()
                .map(move |n| clone1.wrap(n))
                .filter(move |n| clone2.wall_bounds(n)),
        )
    }

    #[allow(clippy::needless_lifetimes)] // clippy is wrong.  removing the lifetime triggers an error in my version of rust
    fn all_cells(&self) -> Box<dyn Iterator<Item = SqCellAddress>> {
        let u_count = self.u_count;
        let v_count = self.v_count;
        Box::new((0..u_count).flat_map(move |u| {
            (-1..(v_count as i32)).map(move |v| SqCellAddress::new(u as i32, v as i32))
        }))
    }
}
