use std::io::Write;

pub type Point3D = [f32; 3];

#[derive(Default)]
pub struct BlenderGeometry {
    vertices: Vec<Point3D>,
    faces: Vec<Vec<usize>>,
    pub epsilon: Option<f32>,
}

const MAX_EPSILON: f32 = 5.0e-5;

impl BlenderGeometry {
    pub fn new() -> Self {
        BlenderGeometry::default()
    }

    pub fn add_face(&mut self, vertices: &[Point3D]) {
        let indices = vertices
            .iter()
            .map(|xyz| self.get_or_create_vertex_index(xyz))
            .collect();
        self.faces.push(indices);
    }

    fn get_or_create_vertex_index(&mut self, xyz: &Point3D) -> usize {
        let idx = self.get_vertex_index(xyz);
        match idx {
            None => {
                let rval = self.vertices.len();
                self.vertices.push(*xyz);
                rval
            }
            Some(idx) => idx,
        }
    }

    fn get_vertex_index(&mut self, xyz: &Point3D) -> Option<usize> {
        self.vertices
            .clone()
            .iter()
            .enumerate()
            .find(|(_, old)| self.close_enough(*old, xyz))
            .map(|(idx, _)| idx)
    }

    pub(crate) fn emit(&self, sink: &mut dyn Write) -> Result<(), std::io::Error> {
        writeln!(sink, "vertices = [")?;

        for vertex in &self.vertices {
            writeln!(sink, "  [{},{},{}],", vertex[0], vertex[1], vertex[2])?;
        }
        writeln!(sink, "]")?;

        writeln!(sink, "faces = [")?;
        for face in &self.faces {
            write!(sink, "  [")?;
            for idx in face {
                write!(sink, "{}, ", idx)?;
            }
            writeln!(sink, "],")?;
        }
        writeln!(sink, "]")?;

        Ok(())
    }

    pub fn face_iter(&self) -> impl Iterator<Item = &[usize]> {
        self.faces.iter().map(|face| face.as_slice())
    }

    pub(crate) fn get_vertex(&self, idx: usize) -> &Point3D {
        &self.vertices[idx]
    }

    fn close_enough(&mut self, p0: &Point3D, p1: &Point3D) -> bool {
        let dx = (p0[0] - p1[0]).abs();
        let dy = (p0[1] - p1[1]).abs();
        let dz = (p0[2] - p1[2]).abs();
        let delta = dx.max(dy).max(dz);
        if delta <= MAX_EPSILON {
            return true;
        }
        match self.epsilon {
            None => self.epsilon = Some(delta),
            Some(old_epsilon) => {
                if delta < old_epsilon {
                    self.epsilon = Some(delta)
                }
            }
        }
        false
    }
}
