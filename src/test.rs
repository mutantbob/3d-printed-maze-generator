use crate::hexagonal::{HexCellAddress, HexMazeTopology};
use crate::polygon_split::subdivide_faces;
use crate::ring::RingAccumulator;
use crate::tools::{with_r, CylindricalSpace, MazeWall, Space};
use crate::walls::multi_face_to_blender;
use crate::{
    BlenderGeometry, BlenderMapping, CellAddress, CylindricalCoodinate, Edge, EdgeCornerMapping,
    Topology, WallPolygons,
};
use assert_approx_eq::assert_approx_eq;

#[test]
pub fn test1() {
    let topology1 = crate::HexMazeTopology::new(10, 10);
    let c1 = HexCellAddress::new(0, 5);
    let c2 = HexCellAddress::new(10, 0);
    let c2w = topology1.wrap(c2);
    assert_eq!(c1.coords_2d().1, c2.coords_2d().1);
    assert_eq!(c1, c2w)
}

#[test]
pub fn test_corners_1() {
    let space = cspace20();
    let c0 = HexCellAddress::new(0, 0);
    let c1 = HexCellAddress::new(1, 0);
    let edge = Edge(c0, c1);

    let (x, y) = edge.coord_left(&space);
    assert_eq!((0.0, 0.0), c0.coords_2d());
    assert_eq!((1.0, 1.0 / 3.0_f32.sqrt()), c1.coords_2d());

    assert_approx_eq!(1.0 / 3.0, x);
    assert_approx_eq!(1.0 / 3.0_f32.sqrt(), y);

    let (x, y) = edge.coord_right(&space);
    assert_approx_eq!(2.0 / 3.0, x);
    assert_approx_eq!(0.0, y);
}

const fn cspace20() -> CylindricalSpace {
    CylindricalSpace {
        r0: 10.0,
        max_rho: 20.0,
    }
}

#[test]
pub fn test_corners_2() {
    let c0 = HexCellAddress::new(4, 0);
    let c1 = HexCellAddress::new(5, 0);
    let edge = Edge(c0, c1);

    let (x0, y0) = c0.coords_2d();
    assert_eq!((4.0, 2.309401), (x0, y0));

    let (x, y) = edge.coord_left(&cspace20());

    assert_approx_eq!(x0 + 1.0 / 3.0, x);
    assert_approx_eq!(y0 + 1.0 / 3.0_f32.sqrt(), y);
}

///  _0
/// 1
#[test]
pub fn test_corners_3() {
    let c1 = HexCellAddress::new(4, 0);
    let c0 = HexCellAddress::new(5, 0);
    let edge = Edge(c0, c1);

    let (x0, y0) = c1.coords_2d();
    assert_eq!((4.0, 2.309401), (x0, y0));

    let (x, y) = edge.coord_right(&cspace20());

    assert_approx_eq!(x0 + 1.0 / 3.0, x);
    assert_approx_eq!(y0 + 1.0 / 3.0_f32.sqrt(), y);
}

/// 0
/// |
/// 1
#[test]
pub fn test_corners_4() {
    let c1 = HexCellAddress::new(4, 0);
    let c0 = HexCellAddress::new(4, 1);
    let edge = Edge(c0, c1);

    let (x0, y0) = c1.coords_2d();
    assert_eq!((4.0, 2.309401), (x0, y0));

    let space = cspace20();
    {
        let (x, y) = edge.coord_left(&space);

        assert_approx_eq!(x0 + 1.0 / 3.0, x);
        assert_approx_eq!(y0 + 1.0 / 3.0_f32.sqrt(), y);
    }

    {
        let (x, y) = edge.coord_right(&space);

        assert_approx_eq!(x0 - 1.0 / 3.0, x);
        assert_approx_eq!(y0 + 1.0 / 3.0_f32.sqrt(), y);
    }
}

#[test]
pub fn test_topology_1() {
    let topology = HexMazeTopology::new(20, 10);

    let ns: Vec<_> = topology.neighbors(&HexCellAddress::new(0, 3)).collect();

    assert_eq!(6, ns.len());
}

#[test]
pub fn test_harmonize_angle() {
    let space = CylindricalSpace {
        r0: 10.0,
        max_rho: 20.0,
    };

    let mut r2 = 21.0;
    space.harmonize_angle(0.0, &mut r2);

    assert_eq!(1.0, r2);
}

#[test]
pub fn test_wall() {
    let space = cspace20();
    let low_z = 0.0;
    let high_z = 0.3;

    let wall = MazeWall::new(
        HexCellAddress::new(19, -9),
        HexCellAddress::new(0, 1),
        true,
        false,
        false,
    );

    let xy0 = wall.a.coords_2d();
    let v0 = with_r(xy0, low_z);
    let xy2 = wall.coord_left(&space);
    // let v2 = with_z(xy2, high_z);
    // let xy3 = wall.coord_right(&space);
    // let v3 = with_z(xy3, high_z);

    let xy4 = space.midpoint(xy0, xy2);
    println!("{:?}  {:?}  {:?}", xy0, xy4, xy2);
    let v4 = with_r(xy4, high_z);

    let v0 = space.to_blender(v0);
    // let v2 = space.to_blender(v2);
    // let v3 = space.to_blender(v3);
    let v4 = space.to_blender(v4);

    if (v0.x - v4.x).abs() > 2.0 {
        panic!("problem {:?}", &wall);
    }
}

#[test]
fn test_ring_accumulator() {
    let p1 = [0.0, 0.0, 0.0].into();
    let p2 = [1.0, 0.0, 0.0].into();
    let p3 = [2.0, 0.0, 0.0].into();
    let p4 = [3.0, 1.0, 0.0].into();

    {
        let mut accum = RingAccumulator::default();

        accum.absorb(p1, p2);
        accum.absorb(p2, p3);
        accum.absorb(p3, p1);
        assert_eq!(1, accum.closed_strings.len());
        assert_eq!(0, accum.open_strings.len());
        assert_eq!(3, accum.closed_strings[0].len());
    }

    {
        let mut accum = RingAccumulator::default();

        accum.absorb(p1, p2);
        accum.absorb(p3, p2);
        accum.absorb(p3, p1);
        assert_eq!(1, accum.closed_strings.len());
        assert_eq!(0, accum.open_strings.len());
        assert_eq!(3, accum.closed_strings[0].len());
    }

    {
        let mut accum = RingAccumulator::default();

        accum.absorb(p1, p2);
        accum.absorb(p1, p3);
        accum.absorb(p3, p2);
        assert_eq!(1, accum.closed_strings.len());
        assert_eq!(0, accum.open_strings.len());
        assert_eq!(3, accum.closed_strings[0].len());
    }

    {
        let mut accum = RingAccumulator::default();

        accum.absorb(p1, p2);
        accum.absorb(p3, p1);
        accum.absorb(p3, p2);
        assert_eq!(1, accum.closed_strings.len());
        assert_eq!(0, accum.open_strings.len());
        assert_eq!(3, accum.closed_strings[0].len());
    }

    {
        let mut accum = RingAccumulator::default();

        accum.absorb(p1, p2);
        accum.absorb(p3, p4);
        assert_eq!(2, accum.open_strings.len());
        accum.absorb(p2, p3);
        assert_eq!(0, accum.closed_strings.len());
        assert_eq!(1, accum.open_strings.len());
        accum.absorb(p4, p1);
        assert_eq!(1, accum.closed_strings.len());
        assert_eq!(0, accum.open_strings.len());
        assert_eq!(4, accum.closed_strings[0].len());
    }
}

#[test]
pub fn test_angle_compare() {
    assert_eq!(6.0, crate::radians_wrap(6.0, 0.0));
    assert!(crate::a_cw_of_b(0.0, 1.0));
    assert!(!crate::a_cw_of_b(1.0, 0.0));

    assert!(!crate::a_cw_of_b(0.0, 6.0));
    assert!(crate::a_cw_of_b(0.0, 0.0));
}

#[test]
pub fn test_wall_polygons() {
    let wall = MazeWall::new(
        HexCellAddress::new(4, 2),
        HexCellAddress::new(4, 1),
        false,
        false,
        false,
    );
    let space = CylindricalSpace {
        r0: 10.0,
        max_rho: 80.0,
    };
    let radius_high = 10.0;
    let radius_groove = 9.0;
    let faces = wall.calculate_wall_polygons(&space, radius_high, radius_groove);

    let mut geom = BlenderGeometry::new();
    for face in faces {
        geom.add_face(face.as_slice());
    }
}

#[test]
pub fn test_subdivide() {
    let face1 = [
        CylindricalCoodinate {
            rho: 4.0,
            r: 9.0,
            z: 4.618802,
        },
        CylindricalCoodinate {
            rho: 3.6666667,
            r: 10.0,
            z: 4.041452,
        },
        CylindricalCoodinate {
            rho: 4.0000005,
            r: 10.0,
            z: 4.233902,
        },
    ];

    let smaller = subdivide_faces(
        [
            face1.as_slice(),
            // face2.as_slice(),
            // face3.as_slice(),
        ],
        0.25,
    );

    let space = CylindricalSpace {
        r0: 10.0,
        max_rho: 80.0,
    };
    let smaller = multi_face_to_blender(&space, smaller);

    println!("problem? {:#?}", smaller);

    let mut geom = BlenderGeometry::new();
    for face in smaller {
        geom.add_face(face.as_slice());
    }
}
