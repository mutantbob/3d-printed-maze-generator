use std::io::Write;

pub type Point3D = [f32; 3];

#[derive(Default)]
pub struct BlenderGeometry {
    vertices: Vec<Point3D>,
    faces: Vec<Vec<usize>>,
}

impl BlenderGeometry {
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
}

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
            .iter()
            .enumerate()
            .find(|(_, old)| **old == *xyz)
            .map(|(idx, _)| idx)
    }
}
