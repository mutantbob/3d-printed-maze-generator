extern crate core;

use crate::cut::HalfSpace;
use crate::cylinder_shell::ShellDimensions;
use crate::walls::{add_edge_flat, add_wall, compute_walls};
use blender_geometry::BlenderGeometry;
use euclid::Vector3D;
use hexagonal::{HexCellAddress, HexMazeTopology};
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use ring::RingAccumulator;
use square::{SqCellAddress, SquareMazeTopology};
use std::cmp::Ordering;
use std::collections::HashMap;
use std::f32::consts::{PI, TAU};
use std::fs::File;
use std::io::Write;
use tools::{
    with_r, BidirectionalEdge, BlenderMapping, CellAddress, CylindricalCoodinate, CylindricalSpace,
    Edge, EdgeCornerMapping, MazeWall, Space, Topology,
};
use walls::{CorridorPolygons, WallPolygons};

mod blender_geometry;
mod cut;
mod cylinder_shell;
mod experiments;
mod hexagonal;
mod maze;
mod polygon_split;
mod ring;
mod square;
#[cfg(test)]
mod test;
mod tools;
mod walls;

pub type Point3Ds = euclid::Point3D<f32, ()>;
pub type Vector3Ds = Vector3D<f32, ()>;

fn main() {
    match 5 {
        2 => {
            let _ = craft_hex_maze_1("/tmp/x.svg", "/tmp/geom-hex.py");
        }
        3 => {
            let _ = experiments::check_blender_math("/tmp/x.py");
        }
        4 => {
            let _ = craft_square_maze_1("/tmp/square.svg", "/tmp/geom-sq.py");
        }
        5 => {
            let _ = craft_shells();
        }
        6 => {
            check_pin_mesh().unwrap();
        }
        _ => {
            let _ = experiments::svg_check();
        }
    }
}

fn check_pin_mesh() -> Result<(), std::io::Error> {
    let mesh = cylinder_shell::pin_slices((14.0, 12.0), 4.5, 0.4, 0.1);

    let fname = "/tmp/pin.py";
    let mut f = File::create(fname)?;
    // f.write_fmt()
    std::io::Write::write(&mut f, mesh.generate_script_for_3_1("pin1").as_bytes())?;
    println!("wrote {}", fname);

    Ok(())
}

pub fn craft_shells() -> Result<(), std::io::Error> {
    let cap_thickness = 2.0;
    let overall_length = 61.0 + 11.0 + cap_thickness;
    let mesh = cylinder_shell::make_cylinder_shell(&ShellDimensions {
        angular_resolution: 12 * 4,
        outer_radius: 13.0,
        inner_radius: 11.0,
        overall_length,
        cap_thickness,
        pin_length: 2.0,
        pin_tip_z: overall_length - 4.0,
        pin_slope: 0.5,
    });

    let fname = "/tmp/shell1.py";
    let mut f = File::create(fname)?;
    // f.write_fmt()
    std::io::Write::write(&mut f, mesh.generate_script_for_3_1("shell1").as_bytes())?;
    println!("wrote {}", fname);
    Ok(())
}

pub struct MazeDimensions {
    /// the inner void where you hide the prize
    inner_radius: f32,
    /// the groove of the maze path
    groove_radius: f32,
    /// the main radius of the maze cylinder
    maze_outer_radius: f32,
    /// z coordinate of the "bottom" of the maze (before we augment it with a rest groove and grip cap)
    bottom_z: f32,
    /// Z coordinate of the "top" of the maze
    top_z: f32,
    /// Z coordinate of the bottom of the inside of the cylinder
    pocket_z: f32,
    /// z coordinate of the top of the maze grip
    grip_top: f32,
    ///
    cap_outer_radius: f32,
}

pub fn craft_hex_maze_1(svg_fname: &str, blender_fname_py: &str) -> Result<(), std::io::Error> {
    let topology = HexMazeTopology::new(14, 14);

    let shell_r = 10.0;
    let groove_r = 8.0;

    let cylindrical = CylindricalSpace {
        r0: groove_r,
        max_rho: topology.maximum_x(),
    };

    // let prescale_top_z = topology.maximum_y() + 1.5;
    let bottom_z = -11.0;
    let maze_dimensions = MazeDimensions {
        inner_radius: 7.5,
        groove_radius: groove_r,
        maze_outer_radius: shell_r,
        bottom_z,
        top_z: 61.0,
        pocket_z: -9.0,
        grip_top: -7.0,
        cap_outer_radius: 12.0,
    };

    craft_hex_maze(
        svg_fname,
        blender_fname_py,
        topology,
        ChaCha8Rng::from_seed([7; 32]),
        cylindrical,
        maze_dimensions,
    )
}

pub fn craft_hex_maze(
    fname: &str,
    blender_fname_py: &str,
    topology1: HexMazeTopology,
    mut rng: ChaCha8Rng,
    cylindrical: CylindricalSpace,
    maze_dimensions: MazeDimensions,
) -> Result<(), std::io::Error> {
    let generator = maze::MazeGenerator::new();

    let start = HexCellAddress::new(0, 0);
    let finisher = Some(|end: &HexCellAddress| HexCellAddress::new(end.u, end.v + 1));
    let eligible_to_finish = |cell: &HexCellAddress| {
        let (_, y) = cell.coords_2d();
        y + 0.2 >= topology1.max_y
    };
    let mut edges: Vec<_> = generator
        .generate_edges(
            &mut rng,
            start,
            |cell| topology1.neighbors(cell),
            finisher,
            eligible_to_finish,
        )
        .into_iter()
        // .map(|(c1, c2)| Edge::<HexCellAddress>(c1, c2))
        .collect();

    edges.push(Edge::<HexCellAddress>(
        start,
        HexCellAddress::new(start.u, start.v - 1),
    ));

    save_edges_svg(fname, &edges)?;

    let walls = compute_walls(&topology1, &edges);

    println!("{} edges; {} walls", edges.len(), walls.len());

    write_blender_python(
        blender_fname_py,
        &edges,
        &walls,
        &cylindrical,
        &maze_dimensions,
    )?;

    Ok(())
}

pub fn craft_square_maze_1(fname: &str, blender_fname_py: &str) -> Result<(), std::io::Error> {
    let topology = SquareMazeTopology::new(14, 11);
    let max_rho = topology.maximum_x();

    let cylindrical = CylindricalSpace { r0: 14.0, max_rho };

    let bottom_z = -14.0; //cylindrical.scale_z(-2.1);
    let maze_dimensions = MazeDimensions {
        inner_radius: 12.5,
        groove_radius: 14.0,
        maze_outer_radius: 16.0,
        bottom_z,
        top_z: if true {
            68.05
        } else {
            cylindrical.scale_z(topology.maximum_y() + 2.5)
        },
        pocket_z: bottom_z + 2.0,
        grip_top: -10.0,
        cap_outer_radius: 19.0,
    };

    draw_square_maze(
        fname,
        blender_fname_py,
        topology,
        ChaCha8Rng::from_seed([7; 32]),
        cylindrical,
        maze_dimensions,
    )
}

pub fn draw_square_maze(
    fname: &str,
    blender_fname_py: &str,
    topology: SquareMazeTopology,
    mut rng: ChaCha8Rng,
    cylindrical: CylindricalSpace,
    maze_dimensions: MazeDimensions,
) -> Result<(), std::io::Error> {
    let generator = maze::MazeGenerator::new();

    let start = SqCellAddress::new(0, 0);
    let mut edges: Vec<_> = generator
        .generate_edges(
            &mut rng,
            start,
            |cell| topology.neighbors(cell),
            Some(|end: &SqCellAddress| SqCellAddress::new(end.u, end.v + 1)),
            |cell| cell.v + 1 >= topology.v_count as i32,
        )
        .into_iter()
        // .map(|(c1, c2)| Edge(c1, c2))
        .collect();

    edges.push(Edge(start, SqCellAddress::new(start.u, start.v - 1)));

    save_edges_svg(fname, &edges)?;

    let walls = compute_walls(&topology, &edges);

    println!("{} edges; {} walls", edges.len(), walls.len());

    write_blender_python(
        blender_fname_py,
        edges.as_slice(),
        &walls,
        &cylindrical,
        &maze_dimensions,
    )?;

    Ok(())
}

pub fn maze_to_mesh<CA: CellAddress>(
    edges: &[Edge<CA>],
    walls: &[MazeWall<CA>],
    cylindrical: &CylindricalSpace,
    maze_dimensions: &MazeDimensions,
) -> BlenderGeometry
where
    Edge<CA>: EdgeCornerMapping<CA> + CorridorPolygons<CylindricalSpace>,
    MazeWall<CA>: WallPolygons<CylindricalSpace>,
{
    let mut blender = BlenderGeometry::new();

    for edge in edges {
        add_edge_flat(
            &mut blender,
            edge,
            cylindrical,
            maze_dimensions.maze_outer_radius,
            maze_dimensions.groove_radius,
        );
    }
    for wall in walls {
        add_wall(
            &mut blender,
            wall,
            cylindrical,
            maze_dimensions.maze_outer_radius,
            maze_dimensions.groove_radius,
        );
    }

    if true {
        blender = finish_cylinder(blender, maze_dimensions);
    }
    blender
}

pub fn write_blender_python<CA: CellAddress>(
    fname: &str,
    edges: &[Edge<CA>],
    walls: &[MazeWall<CA>],
    cylindrical: &CylindricalSpace,
    maze_dimensions: &MazeDimensions,
) -> Result<(), std::io::Error>
where
    Edge<CA>: EdgeCornerMapping<CA> + CorridorPolygons<CylindricalSpace>,
    MazeWall<CA>: WallPolygons<CylindricalSpace>,
{
    let blender = maze_to_mesh(edges, walls, cylindrical, maze_dimensions);

    // println!("epsilon = {:?}", blender.epsilon);

    let mut f = File::create(fname)?;
    f.write_all(blender.generate_script_for_3_1("mesh").as_bytes())?;

    println!("wrote {}", fname);
    Ok(())
}

fn finish_cylinder(mesh: BlenderGeometry, maze_dimensions: &MazeDimensions) -> BlenderGeometry {
    let bottom_z = maze_dimensions.bottom_z;
    let top_z = maze_dimensions.top_z;

    let mesh = clip(
        &mesh,
        &HalfSpace {
            v0: Point3Ds::new(0.0, 0.0, top_z),
            normal: Vector3Ds::new(0.0, 0.0, -1.0),
        },
    );

    let mut mesh = clip(
        &mesh,
        &HalfSpace {
            v0: Point3Ds::new(0.0, 0.0, maze_dimensions.grip_top),
            normal: [0.0, 0.0, 1.0].into(),
        },
    );

    let mid_z = 0.5 * (bottom_z + top_z);

    let edge_counts = count_edges(&mesh);

    let mut accum = RingAccumulator::default();

    for (edge, count) in edge_counts.into_iter() {
        match 2.cmp(&count) {
            Ordering::Less => {
                panic!("count = {} for {:?}", count, edge);
            }
            Ordering::Equal => {}
            Ordering::Greater => {
                let v1 = *mesh.get_vertex(edge.idx1);
                let v2 = *mesh.get_vertex(edge.idx2);

                if (v1.x - v2.x).abs().max((v1.y - v2.y).abs()) < 1.0e-6 {
                    println!("alarmingly vertical edge, this could cause problems")
                }

                let low = v1.z < mid_z;
                extrude_edge(
                    &v1,
                    &v2,
                    &mut mesh,
                    &mut accum,
                    low,
                    if low { maze_dimensions.grip_top } else { top_z },
                );
            }
        }
    }

    println!(
        "ring accumulator {} open; {} closed",
        accum.open_strings.len(),
        accum.closed_strings.len()
    );

    match 2.cmp(&accum.closed_strings.len()) {
        Ordering::Less => {
            println!("oh no, too few edges");

            for ring in &accum.closed_strings {
                mesh.add_face(ring);
            }
        }
        Ordering::Equal => {
            let string2 = accum.closed_strings.remove(1);
            let string1 = accum.closed_strings.remove(0);
            let (bottom, top) = {
                if string1.first().unwrap().z < string2.first().unwrap().z {
                    (string1, string2)
                } else {
                    (string2, string1)
                }
            };
            // mesh.add_face(&bottom);
            finish_cylinder_top(
                &mut mesh,
                &top,
                maze_dimensions.inner_radius,
                maze_dimensions.top_z,
                maze_dimensions.pocket_z,
            );
            finish_cylinder_bottom(
                &mut mesh,
                &bottom,
                maze_dimensions.cap_outer_radius,
                maze_dimensions.grip_top,
                maze_dimensions.bottom_z,
            )
        }
        Ordering::Greater => {
            println!("how did this happen?");
        }
    }
    if accum.closed_strings.len() > 2 {
        println!("how did this happen?");
    } else {
        for ring in &accum.closed_strings {
            mesh.add_face(ring);
        }
    }

    mesh
}

fn extrude_edge(
    v1: &Point3Ds,
    v2: &Point3Ds,
    mesh: &mut BlenderGeometry,
    accum: &mut RingAccumulator,
    hippo: bool,
    z0: f32,
) {
    let Point3Ds {
        x: x1,
        y: y1,
        // z: z1,
        ..
    } = *v1;

    let Point3Ds {
        x: x2,
        y: y2,
        // z: z2,
        ..
    } = *v2;

    let clockwise = is_clockwise((x1, y1), (x2, y2));

    let v3 = Point3Ds::new(x1, y1, z0);
    let v4 = Point3Ds::new(x2, y2, z0);

    let (v1, v2, v3, v4) = if hippo != clockwise {
        (v1, v2, v3, v4)
    } else {
        (v2, v1, v4, v3)
    };

    accum.absorb(v4, v3); // this affects the order of the final ring

    let mut new_face = vec![*v2, *v1];
    if !mesh.close_enough(&v3, v1) {
        new_face.push(v3)
    }
    if !mesh.close_enough(&v4, v2) && !mesh.close_enough(&v3, &v4) {
        new_face.push(v4)
    }

    if new_face.len() > 2 {
        mesh.add_face(new_face.as_slice())
    }
}

pub fn clip(mesh: &BlenderGeometry, half_space: &HalfSpace) -> BlenderGeometry {
    let mut rval = BlenderGeometry::new();

    for face in mesh.face_iter() {
        let polygon: Vec<_> = face.iter().map(|&idx| *mesh.get_vertex(idx)).collect();
        let p2 = half_space.intersect_polygon(polygon.as_slice());
        if p2.len() > 2 {
            rval.add_face(&p2);
        }
    }

    println!("{} faces after clip", rval.face_count());

    rval
}

fn count_edges(mesh: &BlenderGeometry) -> HashMap<BidirectionalEdge, i32> {
    let mut edge_counts = HashMap::new();

    for face in mesh.face_iter() {
        for (i, v1) in face.iter().enumerate() {
            let v2 = face[(i + 1) % face.len()];
            let edge = BidirectionalEdge::new(*v1, v2);
            let entry = edge_counts.entry(edge);
            let count = entry.or_insert(0);
            *count += 1;
        }
    }
    edge_counts
}

fn is_clockwise((x1, y1): (f32, f32), (x2, y2): (f32, f32)) -> bool {
    let cross = x1 * y2 - y1 * x2;
    cross < 0.0
}

fn finish_cylinder_top(
    mesh: &mut BlenderGeometry,
    outer_ring: &[Point3Ds],
    core_radius: f32,
    high_z: f32,
    pocket_z: f32,
) {
    let interior_face_count = 12 * 4;

    let ring: Vec<_> = (0..interior_face_count)
        .map(|i| {
            let theta = TAU * i as f32 / interior_face_count as f32;
            let x = core_radius * theta.cos();
            let y = core_radius * theta.sin();
            (x, y)
        })
        .collect();

    let top_ring: Vec<_> = ring
        .iter()
        .map(|(x, y)| Point3Ds::new(*x, *y, high_z))
        .collect();
    let bottom_ring: Vec<_> = ring
        .into_iter()
        .map(|(x, y)| Point3Ds::new(x, y, pocket_z))
        .collect();

    mesh.add_face(&bottom_ring);

    for (i, (v1, v3)) in top_ring.iter().zip(bottom_ring.iter()).enumerate() {
        let i2 = (i + 1) % top_ring.len();
        let v2 = top_ring[i2];
        let v4 = bottom_ring[i2];

        mesh.add_face([*v1, v2, v4, *v3].as_slice())
    }

    for face in bridge_rings(outer_ring, &top_ring) {
        mesh.add_face(&face)
    }

    // mesh.add_face(&top);
}

pub fn finish_cylinder_bottom(
    mesh: &mut BlenderGeometry,
    bottom: &[Point3Ds],
    outer_radius: f32,
    grip_top: f32,
    bottom_z: f32,
) {
    println!("grip top {}; bottom_z {}", grip_top, bottom_z);

    let interior_face_count = 12 * 4;

    let ring: Vec<_> = (0..interior_face_count)
        .map(|i| {
            let theta = -TAU * i as f32 / interior_face_count as f32;
            let x = outer_radius * theta.cos();
            let y = outer_radius * theta.sin();
            (x, y)
        })
        .collect();

    let top_ring: Vec<_> = ring
        .iter()
        .map(|(x, y)| Point3Ds::new(*x, *y, grip_top))
        .collect();
    let bottom_ring: Vec<_> = ring
        .into_iter()
        .map(|(x, y)| Point3Ds::new(x, y, bottom_z))
        .collect();

    mesh.add_face(&bottom_ring);

    for (i, (v1, v3)) in top_ring.iter().zip(bottom_ring.iter()).enumerate() {
        let i2 = (i + 1) % top_ring.len();
        let v2 = top_ring[i2];
        let v4 = bottom_ring[i2];

        mesh.add_face([*v1, v2, v4, *v3].as_slice())
    }

    for face in bridge_rings(bottom, &top_ring) {
        mesh.add_face(&face)
    }
}

pub fn bridge_rings(ring1: &[Point3Ds], ring2: &[Point3Ds]) -> impl Iterator<Item = Vec<Point3Ds>> {
    let thetas1: Vec<_> = ring1.iter().map(|v1| f32::atan2(v1.y, v1.x)).collect();
    let thetas2: Vec<_> = ring2.iter().map(|v2| f32::atan2(v2.y, v2.x)).collect();

    let clockwise = ring_is_clockwise(&thetas1);

    let theta0 = thetas1[0];

    let base2 = thetas2
        .iter()
        .enumerate()
        .fold((None, 0.0), |(old, min), (curr, rot)| {
            let delta = radians_wrap(rot - theta0, -PI).abs();
            match old {
                Some(idx) => {
                    if delta < min {
                        (Some(curr), delta)
                    } else {
                        (Some(idx), min)
                    }
                }
                None => (Some(curr), delta),
            }
        })
        .0
        .unwrap();

    let mut cursor1 = 0;
    let mut cursor2 = base2;
    let mut rval = vec![];
    while cursor1 <
        // 20
        ring1.len()
        || cursor2 < base2 + ring2.len()
    {
        let i1b = (cursor1 + 1) % ring1.len();
        let i2b = (cursor2 + 1) % ring2.len();
        let c1 = cursor1 % ring1.len();
        let v1 = ring1[c1];
        let c2 = cursor2 % ring2.len();
        let v2 = ring2[c2];
        let v3 = &ring1[i1b];
        let v4 = &ring2[i2b];
        let theta1 = radians_middle(thetas1[c1], thetas1[i1b]);
        let theta2 = radians_middle(thetas2[c2], thetas2[i2b]);
        /* println!(
            "{}<{}<{};\t{}<{}<{}",
            thetas1[c1], theta2, thetas2[i2b], thetas2[c2], theta1, thetas1[i1b],
        );*/

        if clockwise == a_cw_of_b(theta1, theta2) {
            rval.push(vec![v1, *v3, v2]);
            cursor1 += 1;
        } else {
            rval.push(vec![v1, *v4, v2]);
            cursor2 += 1;
        }
        // println!("{}, {} ; {}, {}", cursor1, cursor2, theta1, theta2);

        if cursor2 > 1000 {
            panic!("bugger");
        }
    }

    rval.into_iter()
}

fn ring_is_clockwise(ring: &[f32]) -> bool {
    let mut theta = ring[0];
    for th2 in ring {
        theta = radians_wrap(*th2, theta - PI);
    }
    theta > 0.0
}

pub fn radians_wrap(theta: f32, min: f32) -> f32 {
    if theta < min {
        theta + TAU
    } else if theta - TAU >= min {
        theta - TAU
    } else {
        theta
    }
}

pub fn radians_middle(t1: f32, t2: f32) -> f32 {
    let t2 = radians_wrap(t2, t1 - PI);
    0.5 * (t1 + t2)
}

pub fn a_cw_of_b(a: f32, b: f32) -> bool {
    radians_wrap(b - a, 0.0) < PI
}

//

pub fn save_edges_svg<CA: CellAddress>(
    fname: &str,
    edges: &[Edge<CA>],
) -> Result<(), std::io::Error> {
    let mut f = File::create(fname)?;
    writeln!(&mut f, "<svg>")?;

    write!(&mut f, "<path d=\"")?;
    for edge in edges {
        let (x1, y1) = edge.0.coords_2d();
        let (x2, y2) = edge.1.coords_2d();
        write!(&mut f, "M {},{} {},{} ", x1, y1, x2, y2)?;
    }
    writeln!(
        &mut f,
        "\" style=\"stroke:#00f; fill:none; stroke-width:0.1px\"/>"
    )?;
    writeln!(&mut f, "</svg>")?;
    println!("wrote {}", fname);

    Ok(())
}

//
