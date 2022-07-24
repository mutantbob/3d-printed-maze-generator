use crate::HexCellAddress;
use rand::Rng;
use std::collections::HashMap;

pub struct MazeGenerator {
    visited: HashMap<HexCellAddress, bool>,
    boundary: Vec<HexCellAddress>,
}

impl MazeGenerator {
    pub fn new() -> Self {
        MazeGenerator {
            visited: Default::default(),
            boundary: Default::default(),
        }
    }

    pub fn generate_edges(
        mut self,
        start: HexCellAddress,
    ) -> Vec<(HexCellAddress, HexCellAddress)> {
        self.boundary.push(start);

        let mut rng = rand::thread_rng();

        let mut rval = vec![];
        loop {
            if self.boundary.is_empty() {
                break;
            }
            let idx = rng.gen_range(0..self.boundary.len());
            let item = self.boundary.remove(idx);
            let mut candidates = vec![];
            let neighbors: Vec<_> = item
                .neighbors()
                .into_iter()
                .filter(|cell| self.in_bounds(cell))
                .collect();
            for n in neighbors {
                if self.visited.contains_key(&n) {
                    println!("candidate {:?}", n);
                    candidates.push(n);
                } else if !self.boundary.iter().find(|x| **x == n).is_some() {
                    println!("boundary {:?}", n);
                    self.boundary.push(n);
                }
            }

            if !candidates.is_empty() {
                let old = candidates.remove(rng.gen_range(0..candidates.len()));
                let edge = (old, item);
                rval.push(edge);
            }
            self.visited.insert(item, true);
        }

        rval
    }

    pub fn in_bounds(&self, cell: &HexCellAddress) -> bool {
        let (x, y) = cell.coords_2d();
        x >= 0. && x < 10.0 && y >= 0. && y < 20.0
    }
}
