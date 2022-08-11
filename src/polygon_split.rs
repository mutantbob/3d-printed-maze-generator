use crate::tools::lerp;
use crate::CylindricalCoodinate;

pub fn subdivide_faces<'a, I>(faces: I, theta_resolution: f32) -> Vec<Vec<CylindricalCoodinate>>
where
    I: IntoIterator<Item = &'a [CylindricalCoodinate]>,
{
    let mut rval = vec![];
    for face in faces {
        rval.extend(subdivide_face(face, theta_resolution))
    }
    rval
}

pub fn subdivide_face(
    face: &[CylindricalCoodinate],
    theta_resolution: f32,
) -> Vec<Vec<CylindricalCoodinate>> {
    let mut rval = vec![];

    // XXX do something about coordinate wrapping

    let mut slicer = PolygonSlicer::new(face);

    let mut i = 0;
    loop {
        let q = (slicer.min_theta / theta_resolution).floor() + i as f32;
        let theta1 = if i == 0 {
            slicer.min_theta
        } else {
            q * theta_resolution
        };
        if theta1 >= slicer.max_theta {
            break;
        }
        let theta2 = (q + 1.0) * theta_resolution;
        // println!("theta1 {}; theta2 {}", theta1, theta2);

        rval.extend(slicer.slice(theta1, theta2));
        i += 1;
    }

    rval
}

pub enum VertexOvershoot {
    None,
    Gamma,
    Delta,
}

pub struct PolygonSlicer<'a> {
    pub face: &'a [CylindricalCoodinate],
    pub alpha: usize,
    pub beta: usize,
    pub min_theta: f32,
    pub max_theta: f32,
}

impl<'a> PolygonSlicer<'a> {
    fn new(face: &'a [CylindricalCoodinate]) -> Self {
        let (min_theta, max_theta) = {
            let raw_thetas = face.iter().map(|cc| cc.rho);
            let thetas = sorted_set(raw_thetas);
            (thetas[0], *thetas.last().unwrap())
        };

        // find an idx of a coordinate matching min_theta
        let idx = face
            .iter()
            .enumerate()
            .filter_map(|(i, theta)| {
                if theta.rho == min_theta {
                    Some(i)
                } else {
                    None
                }
            })
            .next()
            .unwrap();

        let theta0 = min_theta;
        let mut beta = idx;
        let mut alpha = (idx + 1) % face.len();
        // println!("rho {}; theta0 {}", face[alpha].rho, theta0);
        if face[alpha].rho > theta0 {
            (alpha, beta) = (beta, (idx + face.len() - 1) % face.len());
            // println!("rho {}; theta0 {}", face[beta].rho, theta0);
            if face[beta].rho > theta0 {
                beta = alpha;
            }
        }
        // println!("theta0 = {}", theta0);
        // println!("alpha {}; beta {:?}", alpha, beta);

        PolygonSlicer {
            face,
            alpha,
            beta,
            min_theta,
            max_theta,
        }
    }

    pub fn slice(&mut self, theta1: f32, theta2: f32) -> Vec<Vec<CylindricalCoodinate>> {
        let mut rval = vec![];

        let mut gamma = (self.alpha + 1) % self.face.len();
        let mut delta = (self.beta + self.face.len() - 1) % (self.face.len());

        let mut v1 = project(&self.face[self.alpha], &self.face[gamma], theta1);
        let mut v2 = project(&self.face[self.beta], &self.face[delta], theta1);
        let mut v3 = project(&self.face[self.alpha], &self.face[gamma], theta2);
        let mut v4 = project(&self.face[self.beta], &self.face[delta], theta2);

        for _ in 0..100 {
            let degenerate_start = v1 == v2;

            let vertex_overshoot = if gamma == delta || self.face[gamma].rho == self.face[delta].rho
            {
                VertexOvershoot::None
            } else if self.face[gamma].rho < theta2 {
                if self.face[delta].rho < theta2 {
                    if self.face[gamma].rho < self.face[delta].rho {
                        VertexOvershoot::Gamma
                    } else {
                        VertexOvershoot::Delta
                    }
                } else {
                    VertexOvershoot::Gamma
                }
            } else if self.face[delta].rho < theta2 {
                VertexOvershoot::Delta
            } else {
                VertexOvershoot::None
            };

            match vertex_overshoot {
                VertexOvershoot::Gamma => {
                    let mut new_face = vec![];
                    if !degenerate_start {
                        new_face.push(v2);
                    }
                    new_face.extend([v1, self.face[gamma], v4]);
                    rval.push(new_face);

                    self.alpha = gamma;
                    gamma = (self.alpha + 1) % self.face.len();
                    v1 = self.face[self.alpha];
                    v2 = v4;
                    v3 = project(&self.face[self.alpha], &self.face[gamma], theta2);
                    /*println!(
                        "vertex overshoot gamma; alpha :={}; gamma :={}",
                        self.alpha, gamma
                    );*/
                }
                VertexOvershoot::Delta => {
                    let mut new_face = vec![];
                    if !degenerate_start {
                        new_face.push(v2);
                    }
                    new_face.extend([v1, v3, self.face[delta]]);
                    rval.push(new_face);

                    self.beta = delta;
                    delta = (self.beta + self.face.len() - 1) % (self.face.len());
                    v2 = self.face[self.beta];
                    v1 = v3;
                    v4 = project(&self.face[self.beta], &self.face[delta], theta2);
                    /* println!(
                        "vertex overshoot delta; beta :={}; delta :={}",
                        self.beta, delta
                    );*/
                }
                VertexOvershoot::None => {
                    // println!("no overshoot");
                    break;
                }
            }
        }

        let degenerate_start = v1 == v2;

        let mut face = if degenerate_start {
            vec![v1]
        } else {
            vec![v2, v1]
        };

        if v3 != v1 && v3 != v2 {
            face.push(v3);
        }
        if v4 != v2 && v4 != v3 {
            face.push(v4);
        }

        if face.len() > 2 {
            rval.push(face);
        }

        /* println!(
            "theta2 {}; gamma.rho {}; delta.rho {}",
            theta2, self.face[gamma].rho, self.face[delta].rho
        );*/

        if theta2 >= self.face[gamma].rho {
            self.alpha = gamma;
        }
        if theta2 >= self.face[delta].rho {
            self.beta = delta;
        }

        // println!("theta2 = {}", theta2);
        // println!("alpha {}; beta {:?}", self.alpha, self.beta);

        for face in &rval {
            if degenerate_face(face) {
                panic!("degenerate polygon\n{:?}\nfrom\n{:#?}", face, self.face);
            }
        }

        rval
    }
}

fn degenerate_face(face: &[CylindricalCoodinate]) -> bool {
    face.len() > 2 && face[0] == *face.last().unwrap()
}

fn project(a: &CylindricalCoodinate, b: &CylindricalCoodinate, rho: f32) -> CylindricalCoodinate {
    let t = (rho - a.rho) / (b.rho - a.rho);
    let t = t.clamp(0.0, 1.0);
    let rho = lerp(a.rho, t, b.rho);
    let r = lerp(a.r, t, b.r);
    let z = lerp(a.z, t, b.z);

    CylindricalCoodinate { rho, r, z }
}

fn sorted_set<I>(values: I) -> Vec<f32>
where
    I: IntoIterator<Item = f32>,
{
    let mut thetas: Vec<f32> = vec![];
    for theta in values {
        let idx = thetas
            .iter()
            .enumerate()
            .filter_map(|(i, t)| if *t >= theta { Some(i) } else { None })
            .next();
        match idx {
            Some(idx) => {
                if theta != thetas[idx] {
                    thetas.insert(idx, theta)
                }
            }
            None => thetas.push(theta),
        }
    }
    thetas
}

#[cfg(test)]
mod test {
    use crate::polygon_split::{sorted_set, subdivide_face, subdivide_faces, PolygonSlicer};
    use crate::CylindricalCoodinate;
    use approx::assert_abs_diff_eq;

    #[test]
    pub fn test_sorted_set() {
        let sorted = sorted_set([4.0, 0.0, 3.0, 2.0, 4.0]);

        assert_eq!(sorted, vec![0.0, 2.0, 3.0, 4.0])
    }

    #[test]
    pub fn test_subdivide_1() {
        let faces = subdivide_face(
            &[
                CylindricalCoodinate::new(0.0, 0.0, 0.0),
                CylindricalCoodinate::new(1.0, 1.0, 0.0),
                CylindricalCoodinate::new(1.0, -1.0, 2.0),
            ],
            0.5,
        );

        assert_eq!(2, faces.len());
        assert_eq!(
            faces[0],
            vec![
                CylindricalCoodinate {
                    rho: 0.0,
                    r: 0.0,
                    z: 0.0
                },
                CylindricalCoodinate {
                    rho: 0.5,
                    r: 0.5,
                    z: 0.0
                },
                CylindricalCoodinate {
                    rho: 0.5,
                    r: -0.5,
                    z: 1.0
                }
            ],
        );

        assert_eq!(
            faces[1],
            vec![
                CylindricalCoodinate {
                    rho: 0.5,
                    r: -0.5,
                    z: 1.0
                },
                CylindricalCoodinate {
                    rho: 0.5,
                    r: 0.5,
                    z: 0.0
                },
                CylindricalCoodinate::new(1.0, 1.0, 0.0),
                CylindricalCoodinate::new(1.0, -1.0, 2.0),
            ],
        );
    }

    #[test]
    pub fn test_subdivide_2() {
        let faces = subdivide_face(
            &[
                CylindricalCoodinate::new(0.0, 0.0, 0.0),
                CylindricalCoodinate::new(1.0, 1.0, 0.0),
                CylindricalCoodinate::new(2.0, -1.0, 2.0),
            ],
            1.0,
        );

        assert_eq!(2, faces.len());
        assert_eq!(
            faces[0],
            vec![
                CylindricalCoodinate {
                    rho: 0.0,
                    r: 0.0,
                    z: 0.0
                },
                CylindricalCoodinate {
                    rho: 1.0,
                    r: 1.0,
                    z: 0.0
                },
                CylindricalCoodinate {
                    rho: 1.0,
                    r: -0.5,
                    z: 1.0
                }
            ],
        );

        assert_eq!(
            faces[1],
            vec![
                CylindricalCoodinate {
                    rho: 1.0,
                    r: -0.5,
                    z: 1.0
                },
                CylindricalCoodinate {
                    rho: 1.0,
                    r: 1.0,
                    z: 0.0
                },
                CylindricalCoodinate::new(2.0, -1.0, 2.0),
            ],
        );
    }

    #[test]
    pub fn test_subdivide_3() {
        let faces = subdivide_face(
            &[
                CylindricalCoodinate::new(-1.0, 1.0, 0.0),
                CylindricalCoodinate::new(0.0, 0.0, 0.0),
                CylindricalCoodinate::new(-1.0, -1.0, 2.0),
            ],
            0.5,
        );

        assert_eq!(2, faces.len());
        assert_eq!(
            faces[0],
            vec![
                CylindricalCoodinate {
                    rho: -1.0,
                    r: -1.0,
                    z: 2.0
                },
                CylindricalCoodinate {
                    rho: -1.0,
                    r: 1.0,
                    z: 0.0
                },
                CylindricalCoodinate {
                    rho: -0.5,
                    r: 0.5,
                    z: 0.0
                },
                CylindricalCoodinate {
                    rho: -0.5,
                    r: -0.5,
                    z: 1.0
                },
            ],
        );

        assert_eq!(
            faces[1],
            vec![
                CylindricalCoodinate {
                    rho: -0.5,
                    r: -0.5,
                    z: 1.0
                },
                CylindricalCoodinate {
                    rho: -0.5,
                    r: 0.5,
                    z: 0.0
                },
                CylindricalCoodinate::new(0.0, 0.0, 0.0),
            ],
        );
    }

    #[test]
    pub fn test_subdivide_4() {
        let faces = subdivide_face(
            &[
                CylindricalCoodinate::new(-0.75, 1.0, 0.0),
                CylindricalCoodinate::new(0.0, 0.0, 0.0),
                CylindricalCoodinate::new(-1.0, -1.0, 2.0),
            ],
            0.5,
        );

        assert_eq!(3, faces.len());

        assert_eq!(
            faces[0],
            vec![
                CylindricalCoodinate::new(-1.0, -1.0, 2.0),
                CylindricalCoodinate::new(-0.75, 1.0, 0.0),
                CylindricalCoodinate::new(-0.5, -0.5, 1.0),
            ]
        );

        // println!("faces[1] = {:?}", faces[1]);
        assert_eq!(3, faces[1].len());
        assert_abs_diff_eq!(faces[1][0], CylindricalCoodinate::new(-0.5, -0.5, 1.0));
        assert_abs_diff_eq!(faces[1][1], CylindricalCoodinate::new(-0.75, 1.0, 0.0));
        assert_abs_diff_eq!(faces[1][2], CylindricalCoodinate::new(-0.5, 2.0 / 3.0, 0.0),);

        assert_eq!(3, faces[2].len());
        assert_abs_diff_eq!(faces[2][0], CylindricalCoodinate::new(-0.5, -0.5, 1.0));
        assert_abs_diff_eq!(faces[2][1], CylindricalCoodinate::new(-0.5, 2.0 / 3.0, 0.0),);
        assert_abs_diff_eq!(faces[2][2], CylindricalCoodinate::new(0.0, 0.0, 0.0));
    }

    #[test]
    pub fn test_subdivide_5() {
        let face1 = [
            CylindricalCoodinate {
                rho: 0.0,
                r: 0.0,
                z: -1.0,
            },
            CylindricalCoodinate {
                rho: -0.375,
                r: 2.0,
                z: -1.375,
            },
            CylindricalCoodinate {
                rho: 0.375,
                r: 2.0,
                z: -1.375,
            },
        ];
        let face2 = [
            CylindricalCoodinate {
                rho: -0.375,
                r: 2.0,
                z: -1.375,
            },
            CylindricalCoodinate {
                rho: -0.5,
                r: 2.0,
                z: -1.5,
            },
            CylindricalCoodinate {
                rho: 0.5,
                r: 2.0,
                z: -1.5,
            },
            CylindricalCoodinate {
                rho: 0.375,
                r: 2.0,
                z: -1.375,
            },
        ];
        let faces = [face1.as_slice(), face2.as_slice()];

        let faces = subdivide_faces(faces, 0.25);

        println!("ts5 {:#?}", &faces);
        for face in faces {
            assert!(face.len() > 2)
        }
    }

    #[test]
    pub fn test_subdivide_6() {
        let face_pre = [
            CylindricalCoodinate::new(0.0, 0.0, 0.75),
            CylindricalCoodinate::new(0.0, 0.0, 0.0),
            CylindricalCoodinate::new(3.0, 0.75, 0.0),
        ];
        let faces = subdivide_face(&face_pre, 2.0);

        assert_eq!(2, faces.len());

        assert_eq!(4, faces[0].len());
        assert_abs_diff_eq!(faces[0][0], CylindricalCoodinate::new(0.0, 0.0, 0.75),);
        assert_abs_diff_eq!(faces[0][1], CylindricalCoodinate::new(0.0, 0.0, 0.0),);
        assert_abs_diff_eq!(faces[0][2], CylindricalCoodinate::new(2.0, 0.50, 0.0),);
        assert_abs_diff_eq!(faces[0][3], CylindricalCoodinate::new(2.0, 0.50, 0.25),);

        assert_eq!(3, faces[1].len());
    }

    #[test]
    pub fn test_slice_1() {
        let face_pre = [
            CylindricalCoodinate::new(0.0, 0.0, 0.75),
            CylindricalCoodinate::new(0.0, 0.0, 0.0),
            CylindricalCoodinate::new(3.0, 0.75, 0.0),
        ];
        let mut slicer = PolygonSlicer::new(&face_pre);

        let faces = slicer.slice(2.0, 4.0);
        assert_eq!(1, faces.len());
        assert_eq!(3, faces[0].len());
    }

    #[test]
    pub fn test_subdivide_7() {
        let face_pre = [
            CylindricalCoodinate {
                rho: 0.375,
                r: 2.0,
                z: 2.375,
            },
            CylindricalCoodinate {
                rho: 0.5,
                r: 2.0,
                z: 2.5,
            },
            CylindricalCoodinate {
                rho: -0.5,
                r: 2.0,
                z: 2.5,
            },
            CylindricalCoodinate {
                rho: -0.375,
                r: 2.0,
                z: 2.375,
            },
        ];

        let _faces = subdivide_face(&face_pre, 0.25);
    }
}
