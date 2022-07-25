use crate::{
    with_z, CylindricalSpace, HexCellAddress, HexMazeEdge, HexMazeWall, MazeTopology1, Space,
};
use assert_approx_eq::assert_approx_eq;

#[test]
pub fn test1() {
    let topology1 = crate::MazeTopology1 { after_max_u: 10 };
    let c1 = HexCellAddress::new(0, 5);
    let c2 = HexCellAddress::new(10, 0);
    let c2w = topology1.wrap(c2);
    assert_eq!(c1.coords_2d().1, c2.coords_2d().1);
    assert_eq!(c1, c2w)
}

#[test]
pub fn test_corners_1() {
    let c0 = HexCellAddress::new(0, 0);
    let c1 = HexCellAddress::new(1, 0);
    let edge = HexMazeEdge(c0, c1);

    let (x, y) = edge.coord_left();
    assert_eq!((0.0, 0.0), c0.coords_2d());
    assert_eq!((1.0, 1.0 / 3.0_f32.sqrt()), c1.coords_2d());

    let c2 = HexCellAddress::new(0, 1);
    // let c2 = HexCellAddress::new(1, -1);
    let avg = (
        (c0.coords_2d().0 + c1.coords_2d().0 + c2.coords_2d().0) / 3.0,
        (c0.coords_2d().1 + c1.coords_2d().1 + c2.coords_2d().1) / 3.0,
    );

    assert_approx_eq!(1.0 / 3.0, x);
    assert_approx_eq!(1.0 / 3.0_f32.sqrt(), y);

    let (x, y) = edge.coord_right();
    assert_approx_eq!(2.0 / 3.0, x);
    assert_approx_eq!(0.0, y);
}

#[test]
pub fn test_corners_2() {
    let c0 = HexCellAddress::new(4, 0);
    let c1 = HexCellAddress::new(5, 0);
    let edge = HexMazeEdge(c0, c1);

    let (x0, y0) = c0.coords_2d();
    assert_eq!((4.0, 2.309401), (x0, y0));

    let (x, y) = edge.coord_left();

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

    let (x, y) = edge.coord_right();

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

    {
        let (x, y) = edge.coord_left();

        assert_approx_eq!(x0 + 1.0 / 3.0, x);
        assert_approx_eq!(y0 + 1.0 / 3.0_f32.sqrt(), y);
    }

    {
        let (x, y) = edge.coord_right();

        assert_approx_eq!(x0 - 1.0 / 3.0, x);
        assert_approx_eq!(y0 + 1.0 / 3.0_f32.sqrt(), y);
    }
}

#[test]
pub fn test_topology_1() {
    let topology = MazeTopology1 { after_max_u: 20 };

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
    let low_z = 0.0;
    let high_z = 0.3;

    let wall = HexMazeWall::new(
        HexCellAddress::new(19, -9),
        HexCellAddress::new(0, 1),
        true,
        false,
    );
    let space = CylindricalSpace {
        r0: 10.0,
        max_rho: 20.0,
    };

    let xy0 = wall.a.coords_2d();
    let v0 = with_z(xy0, low_z);
    let xy2 = wall.coord_left();
    let v2 = with_z(xy2, high_z);
    let xy3 = wall.coord_right();
    let v3 = with_z(xy3, high_z);

    let xy4 = space.midpoint(xy0, xy2);
    println!("{:?}  {:?}  {:?}", xy0, xy4, xy2);
    let v4 = with_z(xy4, high_z);

    let v0 = space.to_blender(v0);
    let v2 = space.to_blender(v2);
    let v3 = space.to_blender(v3);
    let v4 = space.to_blender(v4);

    if (v0[0] - v4[0]).abs() > 2.0 {
        panic!("problem {:?}", &wall);
    }
}
