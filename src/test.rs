use crate::{HexCellAddress, HexMazeEdge, MazeTopology1};
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

    let ns: Vec<_> = topology.neighbors(HexCellAddress::new(0, 3)).collect();

    assert_eq!(6, ns.len());
}
