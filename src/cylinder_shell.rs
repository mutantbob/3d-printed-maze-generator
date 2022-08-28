use crate::tools::lerp;
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

    pub pin_length: f32,
    pub pin_tip_z: f32,
    pub pin_slope: f32,
}

//

#[derive(Debug)]
struct QuantizedYIterator {
    i: f32,
    y1: f32,
    y2: f32,
    resolution: f32,
    yq0: f32,
    after_last: bool,
}

impl QuantizedYIterator {
    pub fn new(y1: f32, y2: f32, resolution: f32) -> Self {
        let resolution = if y1 < y2 { resolution } else { -resolution };
        let yq0 = (y1 / resolution).floor();
        QuantizedYIterator {
            i: 0.0,
            y1,
            y2,
            resolution,
            yq0,
            after_last: false,
        }
    }

    pub fn inside(&self, y: f32) -> bool {
        if self.y1 < self.y2 {
            self.y1 <= y && y < self.y2
        } else {
            self.y2 < y && y <= self.y1
        }
    }
}

impl Iterator for QuantizedYIterator {
    type Item = f32;

    fn next(&mut self) -> Option<Self::Item> {
        let y = if self.i == 0.0 {
            self.y1
        } else {
            (self.yq0 + self.i) * self.resolution
        };
        //println!("y {};\tQYI {:?}", y, self);
        if self.after_last {
            None
        } else if self.inside(y) {
            self.i += 1.0;
            Some(y)
        } else {
            self.after_last = true;
            Some(self.y2)
        }
    }
}

//

pub struct VerticalRectangle {
    x1: f32,
    y1: f32,
    z_high: f32,
    z_low: f32,
    x2: f32,
    y2: f32,
}

impl VerticalRectangle {
    pub fn max_x(&self) -> f32 {
        self.x1.max(self.x2)
    }

    pub fn interpolate_x(&self, y: f32) -> f32 {
        let t = (y - self.y1) / (self.y2 - self.y1);
        lerp(self.x1, t, self.x2)
    }

    pub fn as_quad(&self) -> Vec<Point3Ds> {
        vec![
            Point3Ds::new(self.x1, self.y1, self.z_low),
            Point3Ds::new(self.x1, self.y1, self.z_high),
            Point3Ds::new(self.x2, self.y2, self.z_high),
            Point3Ds::new(self.x2, self.y2, self.z_low),
        ]
    }
}

//

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
    {
        let pin = ConeXAxis {
            tip_x: dimensions.inner_radius - dimensions.pin_length,
            tip_z: dimensions.pin_tip_z,
            dy_dx: dimensions.pin_slope,
        };
        let yz_radius = dimensions.pin_slope * dimensions.pin_length;
        for face in rings_to_quads(&cap_ring_inner, &top_ring_inner) {
            let vr = VerticalRectangle {
                // this only works because of how these quads are built
                x1: face[0].x,
                y1: face[0].y,
                z_high: dimensions.overall_length,
                z_low: dimensions.cap_thickness,
                x2: face[2].x,
                y2: face[2].y,
            };
            for face2 in maybe_subdivide_for_pin(vr, &pin, 0.1, yz_radius) {
                rval.add_face(&face2)
            }
        }
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

pub fn maybe_subdivide_for_pin(
    face: VerticalRectangle,
    pin: &ConeXAxis,
    xz_resolution: f32,
    yz_radius: f32,
) -> Vec<Vec<Point3Ds>> {
    let max_x = face.max_x();
    let beaver_alpha = face.y1.abs() < yz_radius;
    let beaver_dos = face.y2.abs() < yz_radius;
    if max_x > pin.tip_x && (beaver_alpha || beaver_dos) {
        // might overlap with pin

        let mut rval = vec![];

        let steps = (yz_radius / xz_resolution).ceil() as i32;
        let z1 = pin.tip_z - steps as f32 * xz_resolution;
        let z2 = pin.tip_z + steps as f32 * xz_resolution;

        // x_grid[y][z]
        let mut x_grid: Vec<Vec<_>> = QuantizedYIterator::new(face.y1, face.y2, xz_resolution)
            .map(|y| {
                (-steps..=steps)
                    .map(|zq| {
                        let z1 = pin.tip_z - zq as f32 * xz_resolution;
                        let x1 = face.interpolate_x(y);
                        let x2 = pin.x_for(y, z1);
                        Point3Ds::new(x1.min(x2), y, z1)
                    })
                    .collect()
            })
            .collect();

        if beaver_alpha {
            let mut new_face = vec![
                Point3Ds::new(face.x1, face.y1, face.z_high),
                Point3Ds::new(face.x2, face.y2, face.z_high),
                Point3Ds::new(face.x2, face.y2, face.z_low),
                Point3Ds::new(face.x1, face.y1, face.z_low),
            ];

            if true {
                for y in QuantizedYIterator::new(face.y1, face.y2, xz_resolution) {
                    let x = face.interpolate_x(y);
                    new_face.push(Point3Ds::new(x, y, z1));
                }
                new_face.pop();
            }

            if true {
                for xyz in x_grid[x_grid.len() - 2].iter().rev().skip(1) {
                    new_face.push(*xyz);
                }
                x_grid.pop(); // that last column is non-manifold with the adjacent piece
            }
            for y in QuantizedYIterator::new(face.y2, face.y1, xz_resolution).skip(2) {
                let x = face.interpolate_x(y);
                new_face.push(Point3Ds::new(x, y, z2));
            }

            rval.push(new_face);
        } else if beaver_dos {
            let mut new_face = vec![
                Point3Ds::new(face.x2, face.y2, face.z_low),
                Point3Ds::new(face.x1, face.y1, face.z_low),
                Point3Ds::new(face.x1, face.y1, face.z_high),
                Point3Ds::new(face.x2, face.y2, face.z_high),
            ];
            for y in QuantizedYIterator::new(face.y2, face.y1, xz_resolution) {
                let x = face.interpolate_x(y);
                new_face.push(Point3Ds::new(x, y, z2));
            }
            new_face.pop();

            if true {
                x_grid.remove(0); // that first column is non-manifold with the adjacent piece
                for xyz in x_grid[0].iter().skip(1) {
                    new_face.push(*xyz);
                }
            }

            if true {
                for y in QuantizedYIterator::new(face.y1, face.y2, xz_resolution).skip(2) {
                    let x = face.interpolate_x(y);
                    new_face.push(Point3Ds::new(x, y, z1));
                }
            }

            rval.push(new_face);
        } else if z1 > face.z_low {
            let mut new_face = vec![
                Point3Ds::new(face.x1, face.y1, z1),
                Point3Ds::new(face.x1, face.y1, face.z_low),
                Point3Ds::new(face.x2, face.y2, face.z_low),
                Point3Ds::new(face.x2, face.y2, z1),
            ];

            for y in QuantizedYIterator::new(face.y1, face.y2, xz_resolution) {
                let x = face.interpolate_x(y);

                new_face.push(Point3Ds::new(x, y, z1));
            }

            rval.push(new_face);
        }

        if true {
            for face in grid_array_to_quads(&x_grid) {
                rval.push(face);
            }
        }

        return rval;
    };

    vec![face.as_quad()]
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

    for face in grid_array_to_quads(&vertices) {
        rval.add_face(&face);
    }

    rval
}

pub fn grid_array_to_quads(
    vertices: &'_ [Vec<Point3Ds>],
) -> impl Iterator<Item = Vec<Point3Ds>> + '_ {
    vertices
        .iter()
        .enumerate()
        .skip(1)
        .flat_map(move |(i, row)| {
            let i2 = i - 1;
            row.iter().enumerate().skip(1).map(move |(j, p4)| {
                let j2 = j - 1;
                let p2 = vertices[i][j2];
                let p3 = vertices[i2][j];
                let p1 = vertices[i2][j2];
                vec![p1, p2, *p4, p3]
            })
        })
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
