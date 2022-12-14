use crate::polygon_split::subdivide_faces;
use crate::{
    with_r, BlenderGeometry, BlenderMapping, CellAddress, CylindricalCoodinate, Edge,
    EdgeCornerMapping, HexCellAddress, MazeWall, Point3Ds, Space, SqCellAddress, Topology,
};

pub fn compute_walls<CA: CellAddress + PartialEq>(
    topology: &dyn Topology<CA>,
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
    space: &SPACE,
    radius_high: f32,
    radius_groove: f32,
) where
    CA: CellAddress,
    SPACE: Space<(f32, f32)> + BlenderMapping<CylindricalCoodinate>,
    Edge<CA>: EdgeCornerMapping<CA> + CorridorPolygons<SPACE>,
{
    for face in edge.calculate_path_polygons(space, radius_high, radius_groove) {
        blender.add_face(face.as_slice());
    }
}

pub fn add_wall<CA: CellAddress, SPACE: Space<(f32, f32)> + BlenderMapping<CylindricalCoodinate>>(
    blender: &mut BlenderGeometry,
    wall: &MazeWall<CA>,
    space: &SPACE,
    radius_high: f32,
    radius_groove: f32,
) where
    Edge<CA>: EdgeCornerMapping<CA>,
    MazeWall<CA>: WallPolygons<SPACE>,
    SPACE: Space<(f32, f32)> + BlenderMapping<CylindricalCoodinate>,
{
    for face in wall.calculate_wall_polygons(space, radius_high, radius_groove) {
        blender.add_face(face.as_slice());
    }
}

//

pub trait WallPolygons<SPACE>
where
    SPACE: Space<(f32, f32)> + BlenderMapping<CylindricalCoodinate>,
{
    fn calculate_wall_polygons(
        &self,
        space: &SPACE,
        radius_high: f32,
        radius_groove: f32,
    ) -> Vec<Vec<Point3Ds>>;
}

impl<SPACE> WallPolygons<SPACE> for MazeWall<HexCellAddress>
where
    SPACE: Space<(f32, f32)> + BlenderMapping<CylindricalCoodinate>,
{
    fn calculate_wall_polygons(
        &self,
        space: &SPACE,
        radius_high: f32,
        radius_groove: f32,
    ) -> Vec<Vec<Point3Ds>> {
        let wall = self;

        let xy0 = wall.a.coords_2d();
        let v0 = with_r(xy0, radius_groove);
        let xy2 = wall.coord_left(space);
        let v2 = with_r(xy2, radius_high);
        let xy3 = wall.coord_right(space);
        let v3 = with_r(xy3, radius_high);

        let faces = if wall.wall_all {
            let v0 = with_r(xy0, radius_high);

            vec![vec![v0, v3, v2]]
        } else {
            match (wall.wall_ccw, wall.wall_cw) {
                (false, false) => {
                    let xy8 = space.midpoint3(xy0, xy2, xy3);
                    let v8 = with_r(xy8, radius_high);

                    vec![vec![v0, v3, v8], vec![v3, v2, v8], vec![v2, v0, v8]]
                }
                (true, false) => {
                    let v4 = with_r(space.midpoint(xy0, xy2), radius_high);

                    vec![vec![v0, v3, v4], vec![v4, v3, v2]]
                }
                (false, true) => {
                    let v5 = with_r(space.midpoint(xy0, xy3), radius_high);

                    vec![vec![v0, v5, v2], vec![v5, v3, v2]]
                }
                (true, true) => {
                    let v4 = with_r(space.midpoint(xy0, xy2), radius_high);
                    let v5 = with_r(space.midpoint(xy0, xy3), radius_high);

                    vec![vec![v0, v5, v4], vec![v4, v5, v3, v2]]
                }
            }
        };

        let faces = subdivide_faces(faces.iter().map(|face| face.as_slice()), 0.25);

        multi_face_to_blender(space, faces)
    }
}

impl<SPACE> WallPolygons<SPACE> for MazeWall<SqCellAddress>
where
    SPACE: Space<(f32, f32)> + BlenderMapping<CylindricalCoodinate>,
{
    fn calculate_wall_polygons(
        &self,
        space: &SPACE,
        radius_high: f32,
        radius_groove: f32,
    ) -> Vec<Vec<Point3Ds>> {
        let wall = self;

        let xy0 = wall.a.coords_2d();
        let v0 = with_r(
            xy0,
            if wall.wall_all {
                radius_high
            } else {
                radius_groove
            },
        );
        let xy2 = wall.coord_left(space);
        let v2 = with_r(xy2, radius_high);
        let xy3 = wall.coord_right(space);
        let v3 = with_r(xy3, radius_high);

        let frac = 0.75;
        let xy4 = space.lerp(&xy0, frac, &xy2);
        let v4 = with_r(xy4, radius_high);
        let xy5 = space.lerp(&xy0, frac, &xy3);
        let v5 = with_r(xy5, radius_high);

        let face1 = [v0, v5, v4];
        let face2 = [v5, v3, v2, v4];
        let faces = [face1.as_slice(), face2.as_slice()];

        //println!("debug wall polygons; subdivide {:?}", faces);

        let faces = subdivide_faces(faces, 0.25);

        multi_face_to_blender(space, faces)
    }
}

pub fn multi_face_to_blender<SPACE, I>(
    space: &SPACE,
    faces: impl IntoIterator<Item = I>,
) -> Vec<Vec<Point3Ds>>
where
    SPACE: BlenderMapping<CylindricalCoodinate>,
    I: IntoIterator<Item = CylindricalCoodinate>,
{
    faces
        .into_iter()
        .map(|face| face_to_blender(space, face.into_iter()))
        .collect()
}

fn face_to_blender<SPACE>(
    space: &SPACE,
    face: impl IntoIterator<Item = CylindricalCoodinate>,
) -> Vec<Point3Ds>
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
    fn calculate_path_polygons(
        &self,
        space: &SPACE,
        radius_high: f32,
        radius_groove: f32,
    ) -> Vec<Vec<Point3Ds>>;
}

impl<SPACE> CorridorPolygons<SPACE> for Edge<HexCellAddress>
where
    SPACE: Space<(f32, f32)> + BlenderMapping<CylindricalCoodinate>,
{
    fn calculate_path_polygons(
        &self,
        space: &SPACE,
        radius_high: f32,
        radius_groove: f32,
    ) -> Vec<Vec<Point3Ds>> {
        let edge = self;
        let tz0 = edge.0.coords_2d();
        let mut tz1 = edge.1.coords_2d();
        space.maybe_wrap(&mut tz1, &tz0);
        let tzr0 = with_r(tz0, radius_groove);
        // let v0 = space.to_blender(tzr0);
        let tzr1 = with_r(tz1, radius_groove);
        // let v1 = space.to_blender(tzr1);

        let tzr2 = with_r(edge.coord_left(space), radius_high);
        // let v2 = space.to_blender(tzr2);
        let tzr3 = with_r(edge.coord_right(space), radius_high);
        // let v3 = space.to_blender(tzr3);

        let face1 = [tzr0, tzr1, tzr2];
        let face2 = [tzr1, tzr0, tzr3];
        let faces = [face1.as_slice(), face2.as_slice()];
        let faces = subdivide_faces(faces, 0.25);

        //vec![vec![v0, v1, v2], vec![v1, v0, v3]]
        multi_face_to_blender(space, faces)
    }
}

impl<SPACE> CorridorPolygons<SPACE> for Edge<SqCellAddress>
where
    SPACE: Space<(f32, f32)> + BlenderMapping<CylindricalCoodinate>,
{
    fn calculate_path_polygons(
        &self,
        space: &SPACE,
        radius_high: f32,
        radius_groove: f32,
    ) -> Vec<Vec<Point3Ds>> {
        let edge = self;
        let xy0 = edge.0.coords_2d();
        let v0 = with_r(xy0, radius_groove);
        let mut xy1 = edge.1.coords_2d();
        space.maybe_wrap(&mut xy1, &xy0);
        let v1 = with_r(xy1, radius_groove);

        let xy2 = edge.coord_left(space);
        let v2 = with_r(xy2, radius_high);
        let xy3 = edge.coord_right(space);
        let v3 = with_r(xy3, radius_high);

        let frac = 0.75;

        let xy4 = space.lerp(&xy0, frac, &xy2);
        let v4 = with_r(xy4, radius_high);
        let xy5 = space.lerp(&xy0, frac, &xy3);
        let v5 = with_r(xy5, radius_high);

        let xy6 = space.lerp(&xy1, frac, &xy2);
        let v6 = with_r(xy6, radius_high);
        let xy7 = space.lerp(&xy1, frac, &xy3);
        let v7 = with_r(xy7, radius_high);

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
    }
}
