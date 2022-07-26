use crate::tools::{
    with_z, CylindricalSpace, HexCellAddress, HexMazeEdge, HexMazeWall, MazeTopology1,
    RingAccumulator, Space,
};
use assert_approx_eq::assert_approx_eq;

#[test]
pub fn test1() {
    let topology1 = crate::MazeTopology1::new(10, 10);
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
    let edge = HexMazeEdge(c0, c1);

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
    let edge = HexMazeEdge(c0, c1);

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
    let edge = HexMazeEdge(c0, c1);

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
    let edge = HexMazeEdge(c0, c1);

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
    let topology = MazeTopology1::new(20, 10);

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

    let wall = HexMazeWall::new(
        HexCellAddress::new(19, -9),
        HexCellAddress::new(0, 1),
        true,
        false,
        false,
    );

    let xy0 = wall.a.coords_2d();
    let v0 = with_z(xy0, low_z);
    let xy2 = wall.coord_left(&space);
    // let v2 = with_z(xy2, high_z);
    // let xy3 = wall.coord_right(&space);
    // let v3 = with_z(xy3, high_z);

    let xy4 = space.midpoint(xy0, xy2);
    println!("{:?}  {:?}  {:?}", xy0, xy4, xy2);
    let v4 = with_z(xy4, high_z);

    let v0 = space.to_blender(v0);
    // let v2 = space.to_blender(v2);
    // let v3 = space.to_blender(v3);
    let v4 = space.to_blender(v4);

    if (v0[0] - v4[0]).abs() > 2.0 {
        panic!("problem {:?}", &wall);
    }
}

#[test]
fn test_ring_accumulator() {
    let p1 = [0.0, 0.0, 0.0];
    let p2 = [1.0, 0.0, 0.0];
    let p3 = [2.0, 0.0, 0.0];
    let p4 = [3.0, 1.0, 0.0];

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
