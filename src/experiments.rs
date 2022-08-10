use crate::hexagonal::HexCellAddress;
use crate::tools::CylindricalSpace;
use crate::{BlenderGeometry, CellAddress, Edge};
use std::fs::File;
use std::io::Write;

pub fn check_blender_math(fname: &str) -> Result<(), std::io::Error> {
    let mut f = File::create(fname)?;

    let mut blender = BlenderGeometry::new();

    let edge = Edge::<HexCellAddress>(HexCellAddress::new(0, 0), HexCellAddress::new(1, 0));

    crate::walls::add_edge_flat(
        &mut blender,
        &edge,
        1.0,
        &CylindricalSpace {
            r0: 10.0,
            max_rho: 20.0,
        },
    );

    blender.emit(&mut f)?;

    Ok(())
}

pub fn svg_check() -> Result<(), std::io::Error> {
    let fname = "/tmp/x.svg";
    let mut f = File::create(fname)?;
    writeln!(&mut f, "<svg>")?;
    for u in 0..10 {
        for v in 0..10 {
            let cell = HexCellAddress::new(u, v);
            let (x, y) = cell.coords_2d();
            writeln!(
                &mut f,
                "<circle cx=\"{}\" cy=\"{}\" r=\"0.1\" style=\"fill:#000\"/>",
                x, y
            )?;

            for n in cell.neighbors() {
                let (x9, y9) = n.coords_2d();
                let dx = x9 - x;
                let dy = y9 - y;

                let side = 0.05;
                let x1 = crate::tools::lerp(x, 0.1, x9) + side * dy;
                let y1 = crate::tools::lerp(y, 0.1, y9) - side * dx;
                let x8 = crate::tools::lerp(x, 0.9, x9) + side * dy;
                let y8 = crate::tools::lerp(y, 0.9, y9) - side * dx;

                writeln!(&mut f, "<path d=\"M {},{} {},{}\" style=\"fill:none; stroke:#00f; stroke-width: 0.1px\"/>",
                         x1, y1, x8, y8)?;
            }
        }
    }
    writeln!(&mut f, "</svg>")?;
    println!("wrote {}", fname);

    Ok(())
}
