use crate::tools::lerp;
use crate::{Point3Ds, Vector3Ds};

/// if (v-v0).dot(normal) >= 0.0, then the point is inside the half-plane
pub struct HalfSpace {
    pub v0: Point3Ds,
    pub normal: Vector3Ds,
}

impl HalfSpace {
    pub fn inside(&self, other: &Point3Ds) -> bool {
        self.score(other) >= 0.0
    }

    pub fn score(&self, other: &Point3Ds) -> f32 {
        (*other - self.v0).dot(self.normal)
    }

    pub fn crossing(&self, v1: &Point3Ds, v2: &Point3Ds) -> Point3Ds {
        let s1 = self.score(v1);
        let s2 = self.score(v2);
        let t = (0.0 - s1) / (s2 - s1);

        Point3Ds {
            x: lerp(v1.x, t, v2.x),
            y: lerp(v1.y, t, v2.y),
            z: lerp(v1.z, t, v2.z),
            _unit: Default::default(),
        }
    }

    pub fn intersect_polygon(&self, polygon: &[Point3Ds]) -> Vec<Point3Ds> {
        let inside: Vec<_> = polygon.iter().map(|v1| self.inside(v1)).collect();
        let mut rval: Vec<Point3Ds> = vec![];
        let n_vertices = polygon.len();
        for (idx1, (v1, inside1)) in polygon.iter().zip(&inside).enumerate() {
            let idx2 = (idx1 + 1) % n_vertices;
            let v2 = &polygon[idx2];
            let inside2 = inside[idx2];

            if *inside1 {
                rval.push(*v1);
                if !inside2 {
                    let v3 = self.crossing(v1, v2);
                    rval.push(v3);
                }
            } else if inside2 {
                let v3 = self.crossing(v1, v2);
                rval.push(v3);
            }
        }
        rval
    }
}

//
//
