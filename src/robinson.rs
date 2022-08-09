use euclid::Point2D;
use std::f32::consts::PI;
use std::ops::{Add, Mul, Sub};

pub type Point2Ds = Point2D<f32, ()>;

pub const GR0: f32 = 0.618034; //(5.0_f32.sqrt() - 1.0) / 2.0;

pub fn interpolate2<T: Copy, U>(a: Point2D<T, U>, t: T, b: Point2D<T, U>) -> Point2D<T, U>
where
    T: Add<Output = T> + Mul<Output = T> + Sub<Output = T> + Copy,
    T: From<f32>,
{
    let x = interpolate(a.x, t, b.x);
    let y = interpolate(a.y, t, b.y);
    (x, y).into()
}

pub fn interpolate<T>(a: T, t: T, b: T) -> T
where
    T: Add<Output = T> + Mul<Output = T> + Sub<Output = T> + Copy,
    T: From<f32>,
{
    let one = <T as From<f32>>::from(1.0);
    a * (one - t) + b * t
}

#[derive(Copy, Clone, Debug)]
pub enum RobinsonTriangle {
    A(RobinsonTriangleA),
    O(RobinsonTriangleO),
}

impl RobinsonTriangle {
    pub fn a(&self) -> &Point2Ds {
        match self {
            RobinsonTriangle::A(tri) => &tri.a,
            RobinsonTriangle::O(tri) => &tri.a,
        }
    }

    pub fn b(&self) -> &Point2Ds {
        match self {
            RobinsonTriangle::A(tri) => &tri.b,
            RobinsonTriangle::O(tri) => &tri.b,
        }
    }

    pub fn c(&self) -> &Point2Ds {
        match self {
            RobinsonTriangle::A(tri) => &tri.c,
            RobinsonTriangle::O(tri) => &tri.c,
        }
    }

    pub fn decompose_alpha(&self) -> Vec<RobinsonTriangle> {
        match self {
            RobinsonTriangle::A(tri) => tri.decompose_alpha(),
            RobinsonTriangle::O(tri) => tri.decompose_alpha(),
        }
    }

    pub fn decompose_beta(&self) -> Vec<RobinsonTriangle> {
        match self {
            RobinsonTriangle::A(tri) => tri.decompose_beta(),
            RobinsonTriangle::O(tri) => tri.decompose_beta(),
        }
    }
}

impl Default for RobinsonTriangle {
    fn default() -> Self {
        RobinsonTriangle::A(RobinsonTriangleA::sample())
    }
}

//

#[derive(Copy, Clone, Debug)]
pub struct RobinsonTriangleA {
    pub a: Point2Ds,
    pub b: Point2Ds,
    pub c: Point2Ds,
}

impl RobinsonTriangleA {
    pub fn new(a: Point2Ds, b: Point2Ds, c: Point2Ds) -> Self {
        RobinsonTriangleA { a, b, c }
    }

    pub fn sample() -> Self {
        let d36 = PI / 5.0;
        RobinsonTriangleA {
            a: (0.0, 0.0).into(),
            b: (d36.cos(), d36.sin()).into(),
            c: (1.0, 0.0).into(),
        }
    }

    pub fn decompose_alpha(&self) -> Vec<RobinsonTriangle> {
        let p4 = interpolate2(self.a, GR0, self.c);

        vec![
            RobinsonTriangle::A(RobinsonTriangleA::new(self.b, self.c, p4)),
            RobinsonTriangle::O(RobinsonTriangleO::new(p4, self.a, self.b)),
        ]
    }

    pub fn decompose_beta(&self) -> Vec<RobinsonTriangle> {
        vec![RobinsonTriangle::A(*self)]
    }
}

#[derive(Copy, Clone, Debug)]
pub struct RobinsonTriangleO {
    pub a: Point2Ds,
    pub b: Point2Ds,
    pub c: Point2Ds,
}

impl RobinsonTriangleO {
    pub fn new(a: Point2Ds, b: Point2Ds, c: Point2Ds) -> Self {
        RobinsonTriangleO { a, b, c }
    }

    pub fn decompose_alpha(&self) -> Vec<RobinsonTriangle> {
        vec![RobinsonTriangle::O(*self)]
    }

    pub fn decompose_beta(&self) -> Vec<RobinsonTriangle> {
        let p4 = interpolate2(self.c, GR0, self.b);

        vec![
            RobinsonTriangle::A(RobinsonTriangleA::new(self.c, p4, self.a)),
            RobinsonTriangle::O(RobinsonTriangleO::new(p4, self.a, self.b)),
        ]
    }
}

//

pub fn decompose(levels: usize, seed: &[RobinsonTriangle]) -> Vec<RobinsonTriangle> {
    if levels == 0 {
        seed.to_vec()
    } else {
        let mut rval = vec![];
        for level in 0..levels {
            rval = decompose1(if level == 0 { seed } else { &rval }, 0 != (level & 1));
        }
        rval
    }
}

pub fn decompose1(triangles: &[RobinsonTriangle], alpha_not_beta: bool) -> Vec<RobinsonTriangle> {
    let mut rval = vec![];
    for tri in triangles {
        if alpha_not_beta {
            rval.extend(tri.decompose_alpha())
        } else {
            rval.extend(tri.decompose_beta())
        }
    }
    rval
}

#[cfg(test)]

mod test {
    use crate::robinson::GR0;

    #[test]
    fn test_gr0() {
        assert_eq!(0.0, GR0 - (1.25_f32.sqrt() - 0.5));
    }
}
