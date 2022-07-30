use crate::blender_geometry::{BlenderGeometry, Point3D};
use crate::tools::with_z;
use lazy_static::lazy_static;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use tools::{
    BidirectionalEdge, CylindricalSpace, HexCellAddress, HexMazeEdge, HexMazeWall, MazeTopology1,
    RingAccumulator, Space,
};

mod blender_geometry;
mod experiments;
mod maze;
#[cfg(test)]
mod test;
mod tools;

lazy_static! {
    static ref TAN_30: f32 = 1.0f32 / 3.0f32.sqrt();
    static ref SEC_30: f32 = 2.0f32 / 3.0f32.sqrt();
}

fn main() {
    match 2 {
        2 => {
            let _ = draw_maze("/tmp/x.svg");
        }
        3 => {
            let _ = experiments::check_blender_math("/tmp/x.py");
        }
        _ => {
            let _ = experiments::svg_check();
        }
    }
}

pub fn draw_maze(fname: &str) -> Result<(), std::io::Error> {
    let generator = maze::MazeGenerator::new();
    let topology1 = MazeTopology1::new(14, 14);
    let mut edges: Vec<_> = generator
        .generate_edges(
            HexCellAddress::new(0, 0),
            |cell| topology1.neighbors(cell),
            Some(|end: &HexCellAddress| HexCellAddress::new(end.u, end.v + 1)),
            |cell| {
                let (_, y) = cell.coords_2d();
                y + 0.2 >= topology1.max_y
            },
        )
        .into_iter()
        .map(|(c1, c2)| HexMazeEdge(c1, c2))
        .collect();

    edges.push(HexMazeEdge(
        HexCellAddress::new(0, 0),
        HexCellAddress::new(0, -1),
    ));

    save_edges_svg(fname, &edges)?;

    let walls = compute_walls(&topology1, &edges);

    println!("{} edges; {} walls", edges.len(), walls.len());

    write_blender_python("/tmp/geom.py", &edges, &walls, &topology1)?;

    Ok(())
}

fn compute_walls(topology: &MazeTopology1, corridors: &[HexMazeEdge]) -> Vec<HexMazeWall> {
    let cells: Vec<_> = topology.all_cells().collect();
    println!("{} cells", cells.len());
    let mut rval = vec![];
    for cell in cells {
        let neighbors: Vec<_> = topology.wall_neighbors(&cell).collect();
        // println!("{:?} has {} neighbors", &cell, neighbors.len());
        let mut directions = vec![];
        let mut is_wall = vec![];
        for n in neighbors {
            let (edge, wallness) = {
                let edge = HexMazeEdge(cell, n);
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
                rval.push(HexMazeWall::new(
                    wall.0, wall.1, wall_ccw, wall_cw, wall_all,
                ));
            }
        }
    }

    rval
}

pub fn write_blender_python(
    fname: &str,
    edges: &[HexMazeEdge],
    walls: &[HexMazeWall],
    topology: &MazeTopology1,
) -> Result<(), std::io::Error> {
    let max_rho = topology.after_max_u as f32;
    let mut blender = BlenderGeometry::new();

    let high_z = 2.0;
    let cylindrical = CylindricalSpace {
        r0: 10.0 - high_z,
        max_rho,
    };
    for edge in edges {
        add_edge_flat(
            &mut blender,
            edge,
            high_z,
            |xyz| cylindrical.to_blender(xyz),
            &cylindrical,
        );
    }
    for wall in walls {
        add_wall(&mut blender, wall, high_z, &cylindrical);
    }

    if true {
        finish_cylinder(
            &mut blender,
            cylindrical.scale_z(-3.0),
            cylindrical.scale_z(topology.max_y + 2.5),
        );
    }

    println!("epsilon = {:?}", blender.epsilon);

    let mut f = File::create(fname)?;
    blender.emit(&mut f)?;

    println!("wrote {}", fname);
    Ok(())
}

fn finish_cylinder(mesh: &mut BlenderGeometry, low_z: f32, high_z: f32) {
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

                let z0 = if z1 < 0.5 * (low_z + high_z) {
                    low_z
                } else {
                    high_z
                };

                let v3 = [x1, y1, z0];
                let v4 = [x2, y2, z0];
                accum.absorb(v3, v4);
                mesh.add_face(&[[x1, y1, z1], v3, v4, [x2, y2, z2]])
            }
        }
    }

    println!(
        "ring accumulator {} open; {} closed",
        accum.open_strings.len(),
        accum.closed_strings.len()
    );

    if accum.closed_strings.len() > 2 {
        panic!("how did this happen?");
    }

    for ring in &accum.closed_strings {
        mesh.add_face(ring);
    }
}

pub fn save_edges_svg(fname: &str, edges: &[HexMazeEdge]) -> Result<(), std::io::Error> {
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

pub fn add_edge_flat<F>(
    blender: &mut BlenderGeometry,
    edge: &HexMazeEdge,
    high_z: f32,
    mapping: F,
    space: &dyn Space<(f32, f32)>,
) where
    F: Fn(Point3D) -> Point3D,
{
    let v0 = mapping(with_z(edge.0.coords_2d(), 0.0));
    let v1 = mapping(with_z(edge.1.coords_2d(), 0.0));

    let v2 = mapping(with_z(edge.coord_left(space), high_z));
    let v3 = mapping(with_z(edge.coord_right(space), high_z));

    blender.add_face(&[v0, v1, v2]);
    blender.add_face(&[v1, v0, v3]);
}

pub fn add_wall(
    blender: &mut BlenderGeometry,
    wall: &HexMazeWall,
    high_z: f32,
    space: &CylindricalSpace,
) {
    let low_z = 0.0;

    let xy0 = wall.a.coords_2d();
    let v0 = with_z(xy0, low_z);
    let xy2 = wall.coord_left(space);
    let v2 = with_z(xy2, high_z);
    let xy3 = wall.coord_right(space);
    let v3 = with_z(xy3, high_z);

    if wall.wall_all {
        let v0 = with_z(xy0, high_z);

        let v0 = space.to_blender(v0);
        let v2 = space.to_blender(v2);
        let v3 = space.to_blender(v3);

        blender.add_face(&[v0, v3, v2]);
    } else {
        match (wall.wall_ccw, wall.wall_cw) {
            (false, false) => {
                let xy8 = space.midpoint3(xy0, xy2, xy3);
                let v8 = with_z(xy8, high_z);

                let v0 = space.to_blender(v0);
                let v2 = space.to_blender(v2);
                let v3 = space.to_blender(v3);
                let v8 = space.to_blender(v8);

                blender.add_face(&[v0, v3, v8]);
                blender.add_face(&[v3, v2, v8]);
                blender.add_face(&[v2, v0, v8]);
            }
            (true, false) => {
                let v4 = with_z(space.midpoint(xy0, xy2), high_z);

                let v0 = space.to_blender(v0);
                let v2 = space.to_blender(v2);
                let v3 = space.to_blender(v3);
                let v4 = space.to_blender(v4);

                /*if (v0[0] - v4[0]).abs() > 2.0 {
                    println!("problem {:?}", &wall);
                }*/

                blender.add_face(&[v0, v3, v4]);
                blender.add_face(&[v4, v3, v2]);
            }
            (false, true) => {
                let v5 = with_z(space.midpoint(xy0, xy3), high_z);

                let v0 = space.to_blender(v0);
                let v2 = space.to_blender(v2);
                let v3 = space.to_blender(v3);
                let v5 = space.to_blender(v5);

                blender.add_face(&[v0, v5, v2]);
                blender.add_face(&[v5, v3, v2]);
            }
            (true, true) => {
                let v4 = with_z(space.midpoint(xy0, xy2), high_z);
                let v5 = with_z(space.midpoint(xy0, xy3), high_z);

                let v0 = space.to_blender(v0);
                let v2 = space.to_blender(v2);
                let v3 = space.to_blender(v3);
                let v4 = space.to_blender(v4);
                let v5 = space.to_blender(v5);

                blender.add_face(&[v0, v5, v4]);
                blender.add_face(&[v4, v5, v3, v2]);
            }
        }
    }
}
