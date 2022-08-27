use crate::{BlenderGeometry, Point3Ds};
use std::f32::consts::TAU;

pub struct ShellDimensions {
    /// how many facets on a cylinder
    pub angular_resolution: u16,
    /// radius of the outer cylinder
    pub outer_radius: f32,
    /// radius of the inner empty space
    pub inner_radius: f32,
    /// total length of the cylinder
    pub overall_length: f32,
    /// thickness of the closed end of the cylinder
    pub cap_thickness: f32,
}

pub fn make_cylinder_shell(dimensions: &ShellDimensions) -> BlenderGeometry {
    let cis: Vec<_> = (0..dimensions.angular_resolution)
        .map(|idx| {
            let theta = TAU * idx as f32 / dimensions.angular_resolution as f32;
            (theta.cos(), theta.sin())
        })
        .collect();

    let cap_ring_outer: Vec<_> = cis
        .iter()
        .map(|(x, y)| {
            Point3Ds::new(
                x * dimensions.outer_radius,
                y * dimensions.outer_radius,
                0.0,
            )
        })
        .collect();
    let cap_ring_inner: Vec<_> = cis
        .iter()
        .map(|(x, y)| {
            Point3Ds::new(
                x * dimensions.inner_radius,
                y * dimensions.inner_radius,
                dimensions.cap_thickness,
            )
        })
        .collect();

    let top_ring_outer: Vec<_> = cis
        .iter()
        .map(|(x, y)| {
            Point3Ds::new(
                x * dimensions.outer_radius,
                y * dimensions.outer_radius,
                dimensions.overall_length,
            )
        })
        .collect();
    let top_ring_inner: Vec<_> = cis
        .iter()
        .map(|(x, y)| {
            Point3Ds::new(
                x * dimensions.inner_radius,
                y * dimensions.inner_radius,
                dimensions.overall_length,
            )
        })
        .collect();

    //

    let mut rval = BlenderGeometry::new();

    rval.add_face(cap_ring_outer.iter().rev());
    for face in rings_to_quads(&top_ring_outer, &cap_ring_outer) {
        rval.add_face(&face);
    }
    for face in rings_to_quads(&top_ring_inner, &top_ring_outer) {
        rval.add_face(&face);
    }
    for face in rings_to_quads(&cap_ring_inner, &top_ring_inner) {
        rval.add_face(&face);
    }
    rval.add_face(&cap_ring_inner);

    rval
}

pub fn rings_to_quads(a: &[Point3Ds], b: &[Point3Ds]) -> Vec<Vec<Point3Ds>> {
    a.iter()
        .zip(b.iter())
        .enumerate()
        .map(|(i, (v1, v2))| {
            let i2 = (i + 1) % a.len();
            let v3 = a[i2];
            let v4 = b[i2];
            vec![*v1, *v2, v4, v3]
        })
        .collect()
}

pub fn pin_slices(
    pin_tip_xz: (f32, f32),
    pin_length: f32,
    pin_dy_dx: f32,
    grid_resolution: f32,
) -> BlenderGeometry {
    let radius = pin_length * pin_dy_dx;

    let z_count_half = (radius / grid_resolution).ceil() as i32;

    let mut rval = BlenderGeometry::new();

    let pin = ConeXAxis {
        tip_x: pin_tip_xz.0,
        tip_z: pin_tip_xz.1,
        dy_dx: pin_dy_dx,
    };

    let vertices: Vec<Vec<_>> = (-z_count_half..=z_count_half)
        .map(|zq| {
            let z = zq as f32 * grid_resolution + pin.tip_z;
            (-z_count_half..=z_count_half)
                .map(|yq| {
                    let y = yq as f32 * grid_resolution;
                    let x = pin.x_for(y, z);
                    Point3Ds::new(x, y, z)
                })
                .collect()
        })
        .collect();

    for (i, row) in vertices.iter().enumerate().skip(1) {
        let i2 = i - 1;
        for (j, p4) in row.iter().enumerate().skip(1) {
            let j2 = j - 1;
            let p2 = vertices[i][j2];
            let p3 = vertices[i2][j];
            let p1 = vertices[i2][j2];
            rval.add_face(&[p1, p2, *p4, p3]);
        }
    }

    rval
}

pub struct ConeXAxis {
    pub tip_x: f32,
    pub tip_z: f32,
    pub dy_dx: f32,
}

impl ConeXAxis {
    pub fn x_for(&self, y: f32, z: f32) -> f32 {
        let dz = z - self.tip_z;
        let dx = (dz * dz + y * y).sqrt() / self.dy_dx;
        self.tip_x + dx
    }
}
