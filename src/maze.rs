use crate::Edge;
use rand::Rng;
use std::cmp::Ordering;
use std::collections::hash_map::Entry;
use std::collections::{HashMap, HashSet};
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

    pub fn generate_edges<FN, F2, F3, I>(
        mut self,
        rng: &mut impl Rng,
        start: T,
        neighbors: FN,
        finisher: Option<F2>,
        eligible_to_finish: F3,
    ) -> Vec<Edge<T>>
    where
        FN: Fn(&T) -> I,
        I: IntoIterator<Item = T>,
        T: Clone,
        F2: Fn(&T) -> T,
        F3: Fn(&T) -> bool,
    {
        self.boundary.push(start.clone());

        let mut rval: Vec<Edge<T>> = vec![];
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
            // println!("picked {:?}; has {} neighbors", &item, neighbors.len());
            for n in neighbors {
                if self.visited.contains(&n) {
                    // println!("candidate {:?}", &n);
                    candidates.push(n);
                } else if !self.boundary.iter().any(|x| *x == n) {
                    // is not yet in the boundary
                    // println!("boundary {:?}", n.clone());
                    self.boundary.push(n);
                }
            }

            if !candidates.is_empty() {
                let old = candidates.remove(rng.gen_range(0..candidates.len()));
                let edge = Edge(old, item.clone());
                rval.push(edge);
            }
            self.visited.insert(item);
        }

        if let Some(finisher) = finisher {
            self.augment(&mut rval, start, rng, finisher, eligible_to_finish);
        }

        rval
    }

    pub(crate) fn augment<F, F2, RNG>(
        &self,
        edges: &mut Vec<Edge<T>>,
        start: T,
        rng: &mut RNG,
        finisher: F2,
        eligible_to_finish: F,
    ) where
        F: Fn(&T) -> bool,
        F2: Fn(&T) -> T,
        RNG: Rng,
    {
        let mut distances = HashMap::new();

        let mut boundary = HashMap::new();
        boundary.insert(start, 0);
        while !boundary.is_empty() {
            let candidates = Self::pick_low_cost(&boundary);
            let chosen = &candidates[
                0//rng.gen_range(0..candidates.len())
                ];
            let distance = boundary.remove(chosen).unwrap();
            distances.insert(chosen.clone(), distance);
            for n in Self::connected(chosen, edges) {
                if distances.contains_key(&n) {
                    continue;
                }
                let entry = boundary.entry(n);
                match entry {
                    Entry::Occupied(x) => {
                        if *x.get() > distance + 1 {
                            *x.into_mut() = distance + 1;
                        }
                    }
                    Entry::Vacant(x) => {
                        x.insert(distance + 1);
                    }
                }
            }
        }

        let tmp = distances
            .into_iter()
            .filter(|(cell, _)| eligible_to_finish(cell));
        let tmp: Vec<_> = tmp.collect();
        println!(" exit ring has {} ", tmp.len());
        let far = Self::filter_min(tmp, |cost, old_cost| old_cost.cmp(&cost));
        let winner = far[rng.gen_range(0..far.len())].clone();

        let tail = finisher(&winner);
        edges.push(Edge(winner, tail));
    }

    pub fn connected(node: &T, edges: &[Edge<T>]) -> Vec<T> {
        let mut rval = Vec::new();

        for edge in edges {
            if edge.0 == *node {
                rval.push(edge.1.clone())
            } else if edge.1 == *node {
                rval.push(edge.0.clone())
            }
        }

        rval
    }

    pub fn filter_min<I, F, T2, C>(map: I, cmp: F) -> Vec<T2>
    where
        I: IntoIterator<Item = (T2, C)>,
        F: Fn(C, C) -> Ordering,
        T2: Clone,
        C: Copy,
    {
        let mut low_cost = None;
        let mut candidates = Vec::new();

        for (node, cost) in map.into_iter() {
            match low_cost {
                None => {
                    low_cost = Some(cost);
                    candidates.push(node.clone());
                }
                Some(old_cost) => match cmp(cost, old_cost) {
                    Ordering::Less => {
                        low_cost = Some(cost);
                        candidates.clear();
                        candidates.push(node.clone());
                    }
                    Ordering::Equal => {
                        candidates.push(node.clone());
                    }
                    Ordering::Greater => {
                        // don't care
                    }
                },
            }
        }

        candidates
    }

    pub fn pick_low_cost(boundary: &HashMap<T, u32>) -> Vec<T> {
        Self::filter_min(boundary, |cost, old_cost| cost.cmp(old_cost))
            .into_iter()
            .map(Clone::clone)
            .collect()
    }
}
