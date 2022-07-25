use rand::Rng;
use std::collections::HashSet;
use std::fmt::Debug;
use std::hash::Hash;

pub struct MazeGenerator<T> {
    visited: HashSet<T>,
    boundary: Vec<T>,
}

impl<T> MazeGenerator<T>
where
    T: Eq + Hash + Debug + Clone,
{
    pub fn new() -> Self {
        MazeGenerator {
            visited: Default::default(),
            boundary: Default::default(),
        }
    }

    pub fn generate_edges<FN, I>(mut self, start: T, neighbors: FN) -> Vec<(T, T)>
    where
        FN: Fn(&T) -> I,
        I: IntoIterator<Item = T>,
    {
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
            let neighbors: Vec<_> = neighbors(&item)
                .into_iter()
                //.filter(|cell| topology.in_bounds(cell))
                .collect();
            for n in neighbors {
                if self.visited.contains(&n) {
                    println!("candidate {:?}", &n);
                    candidates.push(n);
                } else if !self.boundary.iter().find(|x| **x == n).is_some() {
                    println!("boundary {:?}", n.clone());
                    self.boundary.push(n);
                }
            }

            if !candidates.is_empty() {
                let old = candidates.remove(rng.gen_range(0..candidates.len()));
                let edge = (old, item.clone());
                rval.push(edge);
            }
            self.visited.insert(item);
        }

        rval
    }
}
