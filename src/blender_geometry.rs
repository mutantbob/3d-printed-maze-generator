use crate::Point3Ds;

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

    pub fn add_face<'a, I>(&mut self, vertices: I)
    where
        I: IntoIterator<Item = &'a Point3Ds>,
    {
        let indices: Vec<_> = vertices
            .into_iter()
            .map(|xyz| self.get_or_create_vertex_index(xyz))
            .collect();
        let indices = suppress_duplicate_vertex(&indices);
        if BlenderGeometry::degenerate_face(&indices) {
            eprintln!("suppressing degenerate face {:?}", &indices);
        } else {
            self.faces.push(indices);
        }
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

    pub fn as_python(&self) -> String {
        let mut rval = String::new();
        self.emit_python(&mut rval).unwrap();
        rval
    }

    pub fn emit_python(&self, sink: &mut dyn std::fmt::Write) -> Result<(), std::fmt::Error> {
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

    pub fn generate_script_for_3_1(&self, mesh_name: &str) -> String {
        format!(
            r#"import bpy

{}

name = "{}"
mesh = bpy.data.meshes.new(name)
mesh.from_pydata(vertices, [] , faces)

obj = bpy.data.objects.get(name)
if obj is None:
  obj = bpy.data.objects.new(name, mesh)
else:
  obj.data = mesh

try:
  bpy.context.scene.collection.objects.link(obj)
except:
  pass
"#,
            self.as_python(),
            mesh_name
        )
    }

    pub fn degenerate_face(face: &[usize]) -> bool {
        if face.len() <= 2 {
            return true;
        }

        for (i, v1) in face.iter().enumerate() {
            for v2 in face[(i + 1)..].iter() {
                if *v1 == *v2 {
                    return true;
                }
            }
        }
        false
    }
}

fn suppress_duplicate_vertex(face: &[usize]) -> Vec<usize> {
    let mut rval: Vec<usize> = vec![];
    for (i, &vi) in face.iter().enumerate() {
        let prev = if i > 0 {
            face[i - 1]
        } else {
            *face.last().unwrap()
        };
        if vi == prev {
            eprintln!("suppressing duplicate vertex {} in face", vi)
        } else {
            rval.push(vi);
        }
    }
    rval
}
