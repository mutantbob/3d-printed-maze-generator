use crate::{MazeDimensions, SeedableRng};
use configparser::ini::Ini;
use rand::Rng;
use rand_chacha::ChaCha8Rng;

pub struct MazeConfig {
    pub config: Ini,
}

impl MazeConfig {
    pub fn from_file(fname: &str) -> Result<Self, String> {
        let mut config = Ini::new();
        config.load(fname)?;

        Ok(MazeConfig { config })
    }

    pub fn output(&self) -> OutputFileNames {
        OutputFileNames::from_config(&self.config)
    }

    pub fn maze_dimensions(&self) -> Result<MazeDimensions, String> {
        maze_dimensions_for(&self.config)
    }

    pub fn rng(&self) -> impl Rng {
        rng_from_config(&self.config)
    }
}

//

pub fn config_expect_maze_dimension(config: &Ini, key: &str) -> Result<f32, String> {
    config
        .getfloat("maze dimensions", key)
        .map(|v| v.unwrap_or_else(|| panic!("missing {}", key)) as f32)
}

pub fn maze_dimensions_for(config: &Ini) -> Result<MazeDimensions, String> {
    let bottom_z = config_expect_maze_dimension(config, "bottom z")?;

    let maze_dimensions = MazeDimensions {
        inner_radius: config_expect_maze_dimension(config, "inner radius")?,
        groove_radius: config_expect_maze_dimension(config, "groove radius")?,
        maze_outer_radius: config_expect_maze_dimension(config, "maze outer radius")?,
        bottom_z,
        top_z: config_expect_maze_dimension(config, "top z")?,
        pocket_z: config
            .getfloat("maze dimensions", "pocket z")?
            .unwrap_or(bottom_z as f64 + 2.0) as f32,
        grip_top: config_expect_maze_dimension(config, "grip top")?,
        cap_outer_radius: config_expect_maze_dimension(config, "cap outer radius")?,
    };
    Ok(maze_dimensions)
}

//

pub struct OutputFileNames {
    pub svg_fname: String,
    pub blender_fname: String,
    pub mesh_name: String,
}

impl OutputFileNames {
    pub fn from_config(config: &Ini) -> Self {
        let svg_fname = config
            .get("output", "svg fname")
            .unwrap_or_else(|| "/tmp/maze.svg".into());
        let blender_fname = config
            .get("output", "blender fname")
            .unwrap_or_else(|| "/tmp/geom-maze".into());
        let mesh_name = config
            .get("output", "mesh name")
            .unwrap_or_else(|| "maze".into());

        OutputFileNames {
            svg_fname,
            blender_fname,
            mesh_name,
        }
    }
}

pub fn rng_from_config(config: &Ini) -> ChaCha8Rng {
    ChaCha8Rng::from_seed(crate::convert_to_seed(config.get("rng", "seed")))
}
