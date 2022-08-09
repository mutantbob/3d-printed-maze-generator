use crate::robinson::{decompose, RobinsonTriangle, RobinsonTriangleA};
use std::error::Error;
use std::f32::consts::PI;
use std::fmt::Write as _;
use std::fs::File;
use std::io::Write as _;

mod robinson;

pub fn main() -> Result<(), Box<dyn Error>> {
    let seed = match 1 {
        0 => {
            let start = RobinsonTriangle::default();

            vec![start]
        }
        _ => {
            let tan18 = (PI / 10.0).tan();
            let alpha = RobinsonTriangle::A(RobinsonTriangleA::new(
                (0.0, 0.0).into(),
                (1.0, -tan18).into(),
                (1.0, tan18).into(),
            ));
            let dos = RobinsonTriangle::A(RobinsonTriangleA::new(
                (1.0, tan18).into(),
                (0.0, 2.0 * tan18).into(),
                (0.0, 0.0).into(),
            ));
            vec![alpha, dos]
        }
    };

    let levels = 4;
    let tris = decompose(levels, &seed);

    let ofname = "/tmp/x.svg";
    let mut f = File::create(ofname).unwrap();

    write!(&mut f, "{}", svg_header(800.0, 800.0))?;

    writeln!(&mut f, r#"<g transform ="scale(800)">"#)?;
    {
        let mut d = String::new();
        for tri in tris {
            writeln!(
                &mut d,
                "M {},{} {},{} {},{} Z",
                tri.a().x,
                tri.a().y,
                tri.b().x,
                tri.b().y,
                tri.c().x,
                tri.c().y
            )?;
        }
        let style = "fill:none; stroke:#000; stroke-width:0.01px";
        writeln!(&mut f, r#"<path d="{}" style="{}" />"#, d, style)?;
    }
    writeln!(&mut f, "</g>")?;
    writeln!(&mut f, "</svg>")?;
    drop(f);

    println!("wrote {}", ofname);
    Ok(())
}

pub fn svg_header(width: f32, height: f32) -> String {
    // let winW = Math.min(1600, width + 800);
    // let win_h = Math.min(1000, height + 300);
    let win_w = 1600;
    let win_h = 1000;
    return format!(
        r##"<svg   xmlns:dc="#http://purl.org/dc/elements/1.1/"
   xmlns:cc="http://creativecommons.org/ns#"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
   xmlns:svg="http://www.w3.org/2000/svg"
   xmlns="http://www.w3.org/2000/svg"
   xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
   xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
   width="{}"
   height="{}"
   id="svg2"
   version="1.1"
   sodipodi:docname="New document 1">
    <sodipodi:namedview
     borderlayer="true"
     inkscape:window-width="{}"
     inkscape:window-height="{}"
     inkscape:zoom="1"
     inkscape:cx="{}"
     inkscape:cy="{}"
     inkscape:window-x="10"
     inkscape:window-y="10"
     inkscape:window-maximized="0"
      >
 <inkscape:grid empspacing="10" type="xygrid" visible="true" enabled="true" />
</sodipodi:namedview>
"##,
        width,
        height,
        win_w,
        win_h,
        width / 2.0,
        height / 2.0
    );
}
