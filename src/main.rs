extern crate core;

use crate::walls::{add_edge_flat, add_wall, compute_walls};
use blender_geometry::{BlenderGeometry, Point3D};
use hexagonal::{HexCellAddress, HexMazeTopology};
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use ring::RingAccumulator;
use square::{SqCellAddress, SquareMazeTopology};
use std::cmp::Ordering;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use tools::{
    with_r, BidirectionalEdge, BlenderMapping, CellAddress, CylindricalCoodinate, CylindricalSpace,
    Edge, EdgeCornerMapping, MazeWall, Space, Topology,
};
use walls::{CorridorPolygons, WallPolygons};

mod blender_geometry;
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

fn main() {
    match 4 {
        2 => {
            let _ = craft_hex_maze_1("/tmp/x.svg", "/tmp/geom-hex.py");
        }
        3 => {
            let _ = experiments::check_blender_math("/tmp/x.py");
        }
        4 => {
            let _ = craft_square_maze_1("/tmp/square.svg", "/tmp/geom-sq.py");
        }
        _ => {
            let _ = experiments::svg_check();
        }
    }
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
}

pub fn craft_hex_maze_1(svg_fname: &str, blender_fname_py: &str) -> Result<(), std::io::Error> {
    let topology = HexMazeTopology::new(14, 14);

    let shell_r = 10.0;
    let groove_r = 8.0;

    let cylindrical = CylindricalSpace {
        r0: groove_r,
        max_rho: topology.maximum_x(),
    };

    let maze_dimensions = MazeDimensions {
        inner_radius: 7.5,
        groove_radius: groove_r,
        maze_outer_radius: shell_r,
        bottom_z: cylindrical.scale_z(-3.0),
        top_z: cylindrical.scale_z(topology.maximum_y() + 2.5),
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

    let maze_dimensions = MazeDimensions {
        inner_radius: 12.5,
        groove_radius: 14.0,
        maze_outer_radius: 16.0,
        bottom_z: cylindrical.scale_z(-2.1),
        top_z: cylindrical.scale_z(topology.maximum_y() + 2.5),
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
        finish_cylinder(
            &mut blender,
            maze_dimensions.bottom_z,
            maze_dimensions.top_z,
        );
    }

    println!("epsilon = {:?}", blender.epsilon);

    let mut f = File::create(fname)?;
    blender.emit(&mut f)?;

    println!("wrote {}", fname);
    Ok(())
}

fn finish_cylinder(mesh: &mut BlenderGeometry, bottom_z: f32, top_z: f32) {
    let edge_counts = count_edges(mesh);

    let mut accum = RingAccumulator::default();

    for (edge, count) in edge_counts.into_iter() {
        match 2.cmp(&count) {
            Ordering::Less => {
                panic!("count = {} for {:?}", count, edge);
            }
            Ordering::Equal => {}
            Ordering::Greater => {
                let &[x1, y1, z1] = mesh.get_vertex(edge.idx1);
                let &[x2, y2, z2] = mesh.get_vertex(edge.idx2);

                if (x1 - x2).abs().max((y1 - y2).abs()) < 1.0e-6 {
                    println!("alarmingly vertical edge, this could cause problems")
                }

                let low = z1 < 0.5 * (bottom_z + top_z);
                let clockwise = is_clockwise((x1, y1), (x2, y2));
                let z0 = if low { bottom_z } else { top_z };

                let v3 = [x1, y1, z0];
                let v4 = [x2, y2, z0];
                if low != clockwise {
                    accum.absorb(v4, v3); // this affects the order of the final ring
                    mesh.add_face(&[[x1, y1, z1], v3, v4, [x2, y2, z2]])
                } else {
                    accum.absorb(v3, v4);
                    mesh.add_face(&[[x1, y1, z1], [x2, y2, z2], v4, v3])
                }
            }
        }
    }

    println!(
        "ring accumulator {} open; {} closed",
        accum.open_strings.len(),
        accum.closed_strings.len()
    );

    if accum.closed_strings.len() > 2 {
        println!("how did this happen?");
    } else {
        for ring in &accum.closed_strings {
            mesh.add_face(ring);
        }
    }
}

fn count_edges(mesh: &mut BlenderGeometry) -> HashMap<BidirectionalEdge, i32> {
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
