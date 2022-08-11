use crate::tools::lerp;
use crate::{
    with_r, BlenderGeometry, BlenderMapping, CellAddress, CylindricalCoodinate, Edge,
    EdgeCornerMapping, HexCellAddress, MazeWall, Point3D, Space, SqCellAddress, Topology,
};

pub fn compute_walls<'a, TOP: Topology<'a, CA>, CA: CellAddress + PartialEq>(
    topology: &'a TOP,
    corridors: &[Edge<CA>],
) -> Vec<MazeWall<CA>>
where
    Edge<CA>: EdgeCornerMapping<CA>,
{
    let cells: Vec<_> = topology.all_cells().collect();
    println!("{} cells", cells.len());
    let mut rval: Vec<MazeWall<CA>> = vec![];
    for cell in cells {
        let neighbors: Vec<_> = topology.wall_neighbors(&cell).collect();
        // println!("{:?} has {} neighbors", &cell, neighbors.len());
        let mut directions = vec![];
        let mut is_wall = vec![];
        for n in neighbors {
            let (edge, wallness) = {
                let edge = Edge(cell, n);
                if !corridors.iter().any(|old| *old == edge) {
                    //that edge does not have a corridor
                    (Some(edge), true)
                } else {
                    (None, false)
                }
            };
            directions.push(edge);
            is_wall.push(wallness);
        }
        let wall_all = !is_wall.iter().any(|b| !*b);

        let n = directions.len();
        for (i, wall) in directions.into_iter().enumerate() {
            if let Some(wall) = wall {
                let wall_ccw = is_wall[(i + n - 1) % n];
                let wall_cw = is_wall[(i + 1) % n];
                rval.push(MazeWall::new(wall.0, wall.1, wall_ccw, wall_cw, wall_all));
            }
        }
    }

    rval
}

pub fn add_edge_flat<CA, SPACE>(
    blender: &mut BlenderGeometry,
    edge: &Edge<CA>,
    delta_r: f32,
    space: &SPACE,
) where
    CA: CellAddress,
    SPACE: Space<(f32, f32)> + BlenderMapping<CylindricalCoodinate>,
    Edge<CA>: EdgeCornerMapping<CA> + CorridorPolygons<SPACE>,
{
    for face in edge.calculate_path_polygons(delta_r, space) {
        blender.add_face(face.as_slice());
    }
}

pub fn add_wall<CA: CellAddress, SPACE: Space<(f32, f32)> + BlenderMapping<CylindricalCoodinate>>(
    blender: &mut BlenderGeometry,
    wall: &MazeWall<CA>,
    delta_r: f32,
    space: &SPACE,
) where
    Edge<CA>: EdgeCornerMapping<CA>,
    MazeWall<CA>: WallPolygons<SPACE>,
    SPACE: Space<(f32, f32)> + BlenderMapping<CylindricalCoodinate>,
{
    for face in wall.calculate_wall_polygons(delta_r, space) {
        blender.add_face(face.as_slice());
    }
}

//

pub trait WallPolygons<SPACE>
where
    SPACE: Space<(f32, f32)> + BlenderMapping<CylindricalCoodinate>,
{
    fn calculate_wall_polygons(&self, delta_r: f32, space: &SPACE) -> Vec<Vec<Point3D>>;
}

impl<SPACE> WallPolygons<SPACE> for MazeWall<HexCellAddress>
where
    SPACE: Space<(f32, f32)> + BlenderMapping<CylindricalCoodinate>,
{
    fn calculate_wall_polygons(&self, delta_r: f32, space: &SPACE) -> Vec<Vec<Point3D>> {
        let wall = self;

        let xy0 = wall.a.coords_2d();
        let v0 = with_r(xy0, 0.0);
        let xy2 = wall.coord_left(space);
        let v2 = with_r(xy2, delta_r);
        let xy3 = wall.coord_right(space);
        let v3 = with_r(xy3, delta_r);

        if wall.wall_all {
            let v0 = with_r(xy0, delta_r);

            let v0 = space.to_blender(v0);
            let v2 = space.to_blender(v2);
            let v3 = space.to_blender(v3);

            vec![vec![v0, v3, v2]]
        } else {
            match (wall.wall_ccw, wall.wall_cw) {
                (false, false) => {
                    let xy8 = space.midpoint3(xy0, xy2, xy3);
                    let v8 = with_r(xy8, delta_r);

                    let v0 = space.to_blender(v0);
                    let v2 = space.to_blender(v2);
                    let v3 = space.to_blender(v3);
                    let v8 = space.to_blender(v8);

                    vec![vec![v0, v3, v8], vec![v3, v2, v8], vec![v2, v0, v8]]
                }
                (true, false) => {
                    let v4 = with_r(space.midpoint(xy0, xy2), delta_r);

                    let v0 = space.to_blender(v0);
                    let v2 = space.to_blender(v2);
                    let v3 = space.to_blender(v3);
                    let v4 = space.to_blender(v4);

                    vec![vec![v0, v3, v4], vec![v4, v3, v2]]
                }
                (false, true) => {
                    let v5 = with_r(space.midpoint(xy0, xy3), delta_r);

                    let v0 = space.to_blender(v0);
                    let v2 = space.to_blender(v2);
                    let v3 = space.to_blender(v3);
                    let v5 = space.to_blender(v5);

                    vec![vec![v0, v5, v2], vec![v5, v3, v2]]
                }
                (true, true) => {
                    let v4 = with_r(space.midpoint(xy0, xy2), delta_r);
                    let v5 = with_r(space.midpoint(xy0, xy3), delta_r);

                    let v0 = space.to_blender(v0);
                    let v2 = space.to_blender(v2);
                    let v3 = space.to_blender(v3);
                    let v4 = space.to_blender(v4);
                    let v5 = space.to_blender(v5);

                    vec![vec![v0, v5, v4], vec![v4, v5, v3, v2]]
                }
            }
        }
    }
}

impl<SPACE> WallPolygons<SPACE> for MazeWall<SqCellAddress>
where
    SPACE: Space<(f32, f32)> + BlenderMapping<CylindricalCoodinate>,
{
    fn calculate_wall_polygons(&self, delta_r: f32, space: &SPACE) -> Vec<Vec<Point3D>> {
        let wall = self;

        let xy0 = wall.a.coords_2d();
        let v0 = with_r(xy0, if wall.wall_all { delta_r } else { 0.0 });
        let xy2 = wall.coord_left(space);
        let v2 = with_r(xy2, delta_r);
        let xy3 = wall.coord_right(space);
        let v3 = with_r(xy3, delta_r);

        let frac = 0.75;
        let xy4 = space.lerp(&xy0, frac, &xy2);
        let v4 = with_r(xy4, delta_r);
        let xy5 = space.lerp(&xy0, frac, &xy3);
        let v5 = with_r(xy5, delta_r);

        if true {
            let face1 = [v0, v5, v4];
            let face2 = [v5, v3, v2, v4];
            let faces = [face1.as_slice(), face2.as_slice()];

            println!("debug wall polygons; subdivide {:?}", faces);

            let faces = subdivide_faces(faces, 0.25);

            multi_face_to_blender(space, faces)
        } else {
            let v0 = space.to_blender(v0);
            let v2 = space.to_blender(v2);
            let v3 = space.to_blender(v3);
            let v4 = space.to_blender(v4);
            let v5 = space.to_blender(v5);

            vec![vec![v0, v5, v4], vec![v5, v3, v2, v4]]
        }
    }
}

fn multi_face_to_blender<'a, SPACE, I>(
    space: &SPACE,
    faces: impl IntoIterator<Item = I>,
) -> Vec<Vec<Point3D>>
where
    SPACE: BlenderMapping<CylindricalCoodinate>,
    I: IntoIterator<Item = CylindricalCoodinate>,
{
    faces
        .into_iter()
        .map(|face| face_to_blender(space, face.into_iter()))
        .collect()
}

fn face_to_blender<'a, SPACE>(
    space: &SPACE,
    face: impl IntoIterator<Item = CylindricalCoodinate>,
) -> Vec<Point3D>
where
    SPACE: BlenderMapping<CylindricalCoodinate>,
{
    face.into_iter().map(|cc| space.to_blender(cc)).collect()
}

//

pub trait CorridorPolygons<SPACE>
where
    SPACE: Space<(f32, f32)> + BlenderMapping<CylindricalCoodinate>,
{
    fn calculate_path_polygons(&self, delta_r: f32, space: &SPACE) -> Vec<Vec<Point3D>>;
}

impl<SPACE> CorridorPolygons<SPACE> for Edge<HexCellAddress>
where
    SPACE: Space<(f32, f32)> + BlenderMapping<CylindricalCoodinate>,
{
    fn calculate_path_polygons(&self, delta_r: f32, space: &SPACE) -> Vec<Vec<Point3D>> {
        let edge = self;
        let tz0 = edge.0.coords_2d();
        let mut tz1 = edge.1.coords_2d();
        space.maybe_wrap(&mut tz1, &tz0);
        let tzr0 = with_r(tz0, 0.0);
        // let v0 = space.to_blender(tzr0);
        let tzr1 = with_r(tz1, 0.0);
        // let v1 = space.to_blender(tzr1);

        let tzr2 = with_r(edge.coord_left(space), delta_r);
        // let v2 = space.to_blender(tzr2);
        let tzr3 = with_r(edge.coord_right(space), delta_r);
        // let v3 = space.to_blender(tzr3);

        let face1 = [tzr0, tzr1, tzr2];
        let face2 = [tzr1, tzr0, tzr3];
        let faces = [face1.as_slice(), face2.as_slice()];
        let faces = subdivide_faces(faces, 0.25);

        //vec![vec![v0, v1, v2], vec![v1, v0, v3]]
        multi_face_to_blender(space, faces)
    }
}

fn subdivide_faces<'a, I>(faces: I, theta_resolution: f32) -> Vec<Vec<CylindricalCoodinate>>
where
    I: IntoIterator<Item = &'a [CylindricalCoodinate]>,
{
    let mut rval = vec![];
    for face in faces {
        rval.extend(subdivide_face(face, theta_resolution))
    }
    rval
}

/* This only works for polygons that are convex on the rho/x dimension */
fn subdivide_face(
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
        println!("theta1 {}; theta2 {}", theta1, theta2);

        rval.extend(slicer.slice(theta1, theta2));
        i += 1;
    }

    rval
}

//

enum VertexOvershoot {
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
        println!("rho {}; theta0 {}", face[alpha].rho, theta0);
        if face[alpha].rho > theta0 {
            (alpha, beta) = (beta, (idx + face.len() - 1) % face.len());
            println!("rho {}; theta0 {}", face[beta].rho, theta0);
            if face[beta].rho > theta0 {
                beta = alpha;
            }
        }
        println!("theta0 = {}", theta0);
        println!("alpha {}; beta {:?}", alpha, beta);

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
                    println!(
                        "vertex overshoot gamma; alpha :={}; gamma :={}",
                        self.alpha, gamma
                    );
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
                    println!(
                        "vertex overshoot delta; beta :={}; delta :={}",
                        self.beta, delta
                    );
                }
                VertexOvershoot::None => break,
            }
        }

        let degenerate_start = v1 == v2;

        let mut face = if degenerate_start {
            vec![v1]
        } else {
            println!(
                "v2 {:?} = project({:?}, {:?}, {})",
                &v2, &self.face[self.beta], &self.face[delta], theta1
            );
            vec![v2, v1]
        };

        if v3 != v1 {
            face.push(v3);
        }
        if v4 != v2 && v4 != v3 {
            face.push(v4);
        }

        if face.len() > 2 {
            rval.push(face);
        } else {
            println!("degenerate polygon {:#?}", face);
        }

        println!(
            "theta2 {}; gamma.rho {}; delta.rho {}",
            theta2, self.face[gamma].rho, self.face[delta].rho
        );

        if theta2 >= self.face[gamma].rho {
            self.alpha = gamma;
        }
        if theta2 >= self.face[delta].rho {
            self.beta = delta;
        }

        println!("theta2 = {}", theta2);
        println!("alpha {}; beta {:?}", self.alpha, self.beta);

        rval
    }
}

//

fn project(a: &CylindricalCoodinate, b: &CylindricalCoodinate, rho: f32) -> CylindricalCoodinate {
    let t = (rho - a.rho) / (b.rho - a.rho);
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

impl<SPACE> CorridorPolygons<SPACE> for Edge<SqCellAddress>
where
    SPACE: Space<(f32, f32)> + BlenderMapping<CylindricalCoodinate>,
{
    fn calculate_path_polygons(&self, delta_r: f32, space: &SPACE) -> Vec<Vec<Point3D>> {
        let edge = self;
        let xy0 = edge.0.coords_2d();
        let v0 = with_r(xy0, 0.0);
        let mut xy1 = edge.1.coords_2d();
        space.maybe_wrap(&mut xy1, &xy0);
        let v1 = with_r(xy1, 0.0);

        let xy2 = edge.coord_left(space);
        let v2 = with_r(xy2, delta_r);
        let xy3 = edge.coord_right(space);
        let v3 = with_r(xy3, delta_r);

        let frac = 0.75;

        let xy4 = space.lerp(&xy0, frac, &xy2);
        let v4 = with_r(xy4, delta_r);
        let xy5 = space.lerp(&xy0, frac, &xy3);
        let v5 = with_r(xy5, delta_r);

        let xy6 = space.lerp(&xy1, frac, &xy2);
        let v6 = with_r(xy6, delta_r);
        let xy7 = space.lerp(&xy1, frac, &xy3);
        let v7 = with_r(xy7, delta_r);

        if true {
            let faces = subdivide_faces(
                vec![
                    [v0, v1, v6, v4].as_slice(),
                    &[v4, v6, v2],
                    &[v0, v5, v7, v1],
                    &[v7, v5, v3],
                ],
                0.25,
            );

            multi_face_to_blender(space, faces)
        } else {
            let v0 = space.to_blender(v0);
            let v1 = space.to_blender(v1);
            let v2 = space.to_blender(v2);
            let v3 = space.to_blender(v3);
            let v4 = space.to_blender(v4);
            let v5 = space.to_blender(v5);
            let v6 = space.to_blender(v6);
            let v7 = space.to_blender(v7);

            vec![
                vec![v0, v1, v6, v4],
                vec![v4, v6, v2],
                vec![v0, v5, v7, v1],
                vec![v7, v5, v3],
            ]
        }
    }
}

#[cfg(test)]
mod test {
    use crate::walls::{sorted_set, subdivide_face, subdivide_faces, PolygonSlicer};
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
}
