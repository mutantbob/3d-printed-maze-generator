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

mod blender_geometry;
mod experiments;
mod hexagonal;
mod maze;
mod ring;
mod square;
#[cfg(test)]
mod test;
mod tools;

fn main() {
    match 2 {
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

    write_blender_python("/tmp/geom-hex.py", &edges, &walls, &topology1)?;

    Ok(())
}

pub fn draw_square_maze(fname: &str) -> Result<(), std::io::Error> {
    let generator = maze::MazeGenerator::new();

    let topology = SquareMazeTopology::new(20, 14);

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

    write_blender_python("/tmp/geom-sq.py", edges.as_slice(), &walls, &topology)?;

    Ok(())
}

fn compute_walls<'a, TOP: Topology<'a, CA>, CA: CellAddress + PartialEq>(
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

pub fn write_blender_python<'a, CA: CellAddress>(
    fname: &str,
    edges: &[Edge<CA>],
    walls: &[MazeWall<CA>],
    topology: &impl Topology<'a, CA>,
) -> Result<(), std::io::Error>
where
    Edge<CA>: EdgeCornerMapping<CA>,
{
    let max_rho = topology.maximum_x() as f32;
    let mut blender = BlenderGeometry::new();

    let groove_depth = 2.0;
    let cylindrical = CylindricalSpace {
        r0: 10.0 - groove_depth,
        max_rho,
    };
    for edge in edges {
        add_edge_flat::<CA, _>(
            &mut blender,
            edge,
            groove_depth,
            |xyz| cylindrical.to_blender(xyz),
            &cylindrical,
        );
    }
    for wall in walls {
        add_wall(&mut blender, wall, groove_depth, &cylindrical);
    }

    if true {
        finish_cylinder(
            &mut blender,
            cylindrical.scale_z(-3.0),
            cylindrical.scale_z(topology.maximum_y() + 2.5),
        );
    }

    println!("epsilon = {:?}", blender.epsilon);

    let mut f = File::create(fname)?;
    blender.emit(&mut f)?;

    println!("wrote {}", fname);
    Ok(())
}

fn finish_cylinder(mesh: &mut BlenderGeometry, groove_r: f32, cylinder_r: f32) {
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

                let z0 = if z1 < 0.5 * (groove_r + cylinder_r) {
                    groove_r
                } else {
                    cylinder_r
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
        println!("how did this happen?");
    } else {
        for ring in &accum.closed_strings {
            mesh.add_face(ring);
        }
    }
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

pub fn add_edge_flat<CA, F>(
    blender: &mut BlenderGeometry,
    edge: &Edge<CA>,
    delta_r: f32,
    mapping: F,
    space: &dyn Space<(f32, f32)>,
) where
    CA: CellAddress,
    F: Fn(CylindricalCoodinate) -> Point3D,
    Edge<CA>: EdgeCornerMapping<CA>,
{
    let v0 = mapping(with_r(edge.0.coords_2d(), 0.0));
    let v1 = mapping(with_r(edge.1.coords_2d(), 0.0));

    let v2 = mapping(with_r(edge.coord_left(space), delta_r));
    let v3 = mapping(with_r(edge.coord_right(space), delta_r));

    blender.add_face(&[v0, v1, v2]);
    blender.add_face(&[v1, v0, v3]);
}

pub fn add_wall<CA: CellAddress, SPACE: Space<(f32, f32)> + BlenderMapping<CylindricalCoodinate>>(
    blender: &mut BlenderGeometry,
    wall: &MazeWall<CA>,
    delta_r: f32,
    space: &SPACE,
) where
    Edge<CA>: EdgeCornerMapping<CA>,
{
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

        blender.add_face(&[v0, v3, v2]);
    } else {
        match (wall.wall_ccw, wall.wall_cw) {
            (false, false) => {
                let xy8 = space.midpoint3(xy0, xy2, xy3);
                let v8 = with_r(xy8, delta_r);

                let v0 = space.to_blender(v0);
                let v2 = space.to_blender(v2);
                let v3 = space.to_blender(v3);
                let v8 = space.to_blender(v8);

                blender.add_face(&[v0, v3, v8]);
                blender.add_face(&[v3, v2, v8]);
                blender.add_face(&[v2, v0, v8]);
            }
            (true, false) => {
                let v4 = with_r(space.midpoint(xy0, xy2), delta_r);

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
                let v5 = with_r(space.midpoint(xy0, xy3), delta_r);

                let v0 = space.to_blender(v0);
                let v2 = space.to_blender(v2);
                let v3 = space.to_blender(v3);
                let v5 = space.to_blender(v5);

                blender.add_face(&[v0, v5, v2]);
                blender.add_face(&[v5, v3, v2]);
            }
            (true, true) => {
                let v4 = with_r(space.midpoint(xy0, xy2), delta_r);
                let v5 = with_r(space.midpoint(xy0, xy3), delta_r);

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
