use crate::Point3Ds;
use std::io::Write;

#[derive(Default)]
pub struct BlenderGeometry {
    vertices: Vec<Point3Ds>,
    faces: Vec<Vec<usize>>,
    // pub epsilon: Option<f32>,
}

const MAX_EPSILON: f32 = 5.0e-5;

impl BlenderGeometry {
    pub fn new() -> Self {
        BlenderGeometry::default()
    }

    pub fn add_face(&mut self, vertices: &[Point3Ds]) {
        let indices = vertices
            .iter()
            .map(|xyz| self.get_or_create_vertex_index(xyz))
            .collect();
        self.faces.push(indices);
    }

    fn get_or_create_vertex_index(&mut self, xyz: &Point3Ds) -> usize {
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

    fn get_vertex_index(&mut self, xyz: &Point3Ds) -> Option<usize> {
        self.vertices
            .clone()
            .iter()
            .enumerate()
            .find(|(_, old)| self.close_enough(*old, xyz))
            .map(|(idx, _)| idx)
    }

    pub fn emit(&self, sink: &mut dyn Write) -> Result<(), std::io::Error> {
        writeln!(sink, "vertices = [")?;

        for vertex in &self.vertices {
            writeln!(sink, "  [{},{},{}],", vertex.x, vertex.y, vertex.z)?;
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

    pub fn get_vertex(&self, idx: usize) -> &Point3Ds {
        &self.vertices[idx]
    }

    pub fn close_enough(&self, p0: &Point3Ds, p1: &Point3Ds) -> bool {
        let dx = (p0.x - p1.x).abs();
        let dy = (p0.y - p1.y).abs();
        let dz = (p0.z - p1.z).abs();
        let delta = dx.max(dy).max(dz);
        if delta <= MAX_EPSILON {
            return true;
        }
        /* match self.epsilon {
            None => self.epsilon = Some(delta),
            Some(old_epsilon) => {
                if delta < old_epsilon {
                    self.epsilon = Some(delta)
                }
            }
        }*/
        false
    }

    pub fn face_count(&self) -> usize {
        self.faces.len()
    }
}
