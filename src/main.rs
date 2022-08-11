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
mod ring;
mod square;
#[cfg(test)]
mod test;
mod tools;
mod walls;

fn main() {
    match 4 {
        2 => {
            let _ = draw_hex_maze("/tmp/x.svg");
        }
        3 => {
            let _ = experiments::check_blender_math("/tmp/x.py");
        }
        4 => {
            let _ = draw_square_maze("/tmp/square.svg");
        }
        _ => {
            let _ = experiments::svg_check();
        }
    }
}

pub fn draw_hex_maze(fname: &str) -> Result<(), std::io::Error> {
    let generator = maze::MazeGenerator::new();
    let topology1 = HexMazeTopology::new(14, 14);

    let mut rng = ChaCha8Rng::from_seed([7; 32]);

    let mut edges: Vec<_> = generator
        .generate_edges(
            &mut rng,
            HexCellAddress::new(0, 0),
            |cell| topology1.neighbors(cell),
            Some(|end: &HexCellAddress| HexCellAddress::new(end.u, end.v + 1)),
            |cell| {
                let (_, y) = cell.coords_2d();
                y + 0.2 >= topology1.max_y
            },
        )
        .into_iter()
        .map(|(c1, c2)| Edge::<HexCellAddress>(c1, c2))
        .collect();

    edges.push(Edge::<HexCellAddress>(
        HexCellAddress::new(0, 0),
        HexCellAddress::new(0, -1),
    ));

    save_edges_svg(fname, &edges)?;

    let walls = compute_walls(&topology1, &edges);

    println!("{} edges; {} walls", edges.len(), walls.len());

    let max_rho = topology1.maximum_x();
    let groove_depth = 2.0;
    let cylindrical = CylindricalSpace {
        r0: 10.0 - groove_depth,
        max_rho,
    };
    write_blender_python(
        "/tmp/geom-hex.py",
        &edges,
        &walls,
        &cylindrical,
        groove_depth,
        -3.0,
        topology1.maximum_y() + 2.5,
    )?;

    Ok(())
}

pub fn draw_square_maze(fname: &str) -> Result<(), std::io::Error> {
    let generator = maze::MazeGenerator::new();

    let topology = SquareMazeTopology::new(14, 12);

    let mut rng = ChaCha8Rng::from_seed([7; 32]);

    let mut edges: Vec<_> = generator
        .generate_edges(
            &mut rng,
            SqCellAddress::new(0, 0),
            |cell| topology.neighbors(cell),
            Some(|end: &SqCellAddress| SqCellAddress::new(end.u, end.v + 1)),
            |cell| cell.v + 1 >= topology.v_count as i32,
        )
        .into_iter()
        .map(|(c1, c2)| Edge(c1, c2))
        .collect();

    edges.push(Edge(SqCellAddress::new(0, 0), SqCellAddress::new(0, -1)));

    save_edges_svg(fname, &edges)?;

    let walls = compute_walls(&topology, &edges);

    println!("{} edges; {} walls", edges.len(), walls.len());

    let max_rho = topology.maximum_x();
    let groove_depth = 2.0;
    let cylindrical = CylindricalSpace {
        r0: 15.0 - groove_depth,
        max_rho,
    };

    write_blender_python(
        "/tmp/geom-sq.py",
        edges.as_slice(),
        &walls,
        &cylindrical,
        groove_depth,
        -2.1,
        topology.maximum_y() + 2.5,
    )?;

    Ok(())
}

pub fn write_blender_python<CA: CellAddress>(
    fname: &str,
    edges: &[Edge<CA>],
    walls: &[MazeWall<CA>],
    cylindrical: &CylindricalSpace,
    groove_depth: f32,
    prescale_bottom_z: f32,
    prescale_top_z: f32,
) -> Result<(), std::io::Error>
where
    Edge<CA>: EdgeCornerMapping<CA> + CorridorPolygons<CylindricalSpace>,
    MazeWall<CA>: WallPolygons<CylindricalSpace>,
{
    let mut blender = BlenderGeometry::new();

    for edge in edges {
        add_edge_flat(&mut blender, edge, groove_depth, cylindrical);
    }
    for wall in walls {
        add_wall(&mut blender, wall, groove_depth, cylindrical);
    }

    if true {
        finish_cylinder(
            &mut blender,
            cylindrical.scale_z(prescale_bottom_z),
            cylindrical.scale_z(prescale_top_z),
        );
    }

    println!("epsilon = {:?}", blender.epsilon);

    let mut f = File::create(fname)?;
    blender.emit(&mut f)?;

    println!("wrote {}", fname);
    Ok(())
}

fn finish_cylinder(mesh: &mut BlenderGeometry, bottom_z: f32, top_z: f32) {
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
