//! 3D point representation and vector arithmetic.
//!
//! This module defines [`Point`], a simple 3D coordinate type used
//! throughout the crate for positions and vectors. It implements the
//! standard arithmetic operators (`+`, `-`, `*`, `/`) so points can be
//! used interchangeably as positions or displacement vectors, and it
//! implements the [`Position`] trait so it can be used generically
//! alongside types like `Sphere` in neighbour-search code.
//!
//! # Examples
//!
//! ```
//! use garnet::geo::point::{Point, Position};
//!
//! let a = Point::new(1.0, 2.0, 3.0);
//! let b = Point::new(4.0, 5.0, 6.0);
//! let sum = a + b;
//! assert_eq!(sum.position().x, 5.0);
//! ```

use std::ops::{Add, Div, Mul, Sub};

use nalgebra::{Rotation3, Vector3};

/// A trait for types that have a 3D position
/// This is used for both Point and Sphere to guarantee
/// that we can access their x,y,z fields during neighbour searching
pub trait Position {
    fn position(&self) -> Vector3<f64>;
}

impl<T: Position> Position for &T {
    fn position(&self) -> Vector3<f64> {
        (*self).position()
    }
}

/// A 3D coordinate
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point {
    location: Vector3<f64>,
}

impl Point {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Point {
            location: Vector3::new(x, y, z),
        }
    }

    pub fn zero() -> Self {
        Point {
            location: Vector3::new(0.0, 0.0, 0.0),
        }
    }

    /// Euclidean length, also the L2 norm
    pub fn norm(&self) -> f64 {
        self.location.norm()
    }

    /// Add another position-like value to this point,
    /// used to accumulate coordinates in centroid calculation.
    pub fn add_position<P: Position>(self, other: &P) -> Point {
        Point {
            location: self.location + other.position(),
        }
    }

    /// Rotate this point using an nalgebra rotation matrix.
    pub fn rotate(&self, rotation: &Rotation3<f64>) -> Point {
        let v: Vector3<f64> = (*self).into();
        (rotation * v).into()
    }
}

impl Default for Point {
    fn default() -> Self {
        Point::zero()
    }
}

/// Operator overloading for Point <-> Point and &Point <-> &Point
/// In both cases the result is a Point
///
/// Note that the asymmetric operations are not defined
/// (i.e., Point + &Point is not defined)
impl Add for Point {
    type Output = Point;
    fn add(self, rhs: Point) -> Point {
        Point {
            location: self.position() + rhs.position(),
        }
    }
}

impl Add<&Point> for &Point {
    type Output = Point;
    fn add(self, rhs: &Point) -> Point {
        *self + *rhs
    }
}

impl Sub for Point {
    type Output = Point;
    fn sub(self, rhs: Point) -> Point {
        Point {
            location: self.position() - rhs.position(),
        }
    }
}

impl Sub<&Point> for &Point {
    type Output = Point;
    fn sub(self, rhs: &Point) -> Point {
        *self - *rhs
    }
}

impl Mul<f64> for Point {
    type Output = Point;
    fn mul(self, s: f64) -> Point {
        Point {
            location: self.position() * s,
        }
    }
}

impl Mul<f64> for &Point {
    type Output = Point;
    fn mul(self, s: f64) -> Point {
        *self * s
    }
}

impl Div<f64> for Point {
    type Output = Point;
    fn div(self, s: f64) -> Point {
        Point {
            location: self.position() / s,
        }
    }
}

impl Div<f64> for &Point {
    type Output = Point;
    fn div(self, s: f64) -> Point {
        *self / s
    }
}

/// The Position trait is used to allow efficient neighbor searching
impl Position for Point {
    fn position(&self) -> Vector3<f64> {
        self.location
    }
}

/// Point is used to represent a point in 3D space
/// To perform rotations and translations on Points using nalgebra types,
/// first convert the Point to a Vector3.
impl From<Point> for Vector3<f64> {
    fn from(p: Point) -> Self {
        p.position()
    }
}

impl From<Vector3<f64>> for Point {
    fn from(v: Vector3<f64>) -> Self {
        Point::new(v.x, v.y, v.z)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_point_fields() {
        let p = Point::new(1.0, 2.0, 3.0);
        assert_eq!(p.location.x, 1.0);
        assert_eq!(p.location.y, 2.0);
        assert_eq!(p.location.z, 3.0);
    }

    #[test]
    fn test_point_add() {
        let a = Point::new(1.0, 2.0, 3.0);
        let b = Point::new(4.0, 5.0, 6.0);
        let r = a + b;
        assert_eq!(r.location.x, 5.0);
        assert_eq!(r.location.y, 7.0);
        assert_eq!(r.location.z, 9.0);
    }

    #[test]
    fn test_point_sub() {
        let a = Point::new(4.0, 5.0, 6.0);
        let b = Point::new(1.0, 2.0, 3.0);
        let r = a - b;
        assert_eq!(r.location.x, 3.0);
        assert_eq!(r.location.y, 3.0);
        assert_eq!(r.location.z, 3.0);
    }

    #[test]
    fn test_point_scale() {
        let p = Point::new(2.0, 4.0, 6.0);
        let r = p / 2.0;
        assert_eq!(r.location.x, 1.0);
        assert_eq!(r.location.y, 2.0);
        assert_eq!(r.location.z, 3.0);
    }

    #[test]
    fn test_point_div() {
        let p = Point::new(2.0, 4.0, 6.0);
        let r = p / 2.0;
        assert_eq!(r.location.x, 1.0);
        assert_eq!(r.location.y, 2.0);
        assert_eq!(r.location.z, 3.0);
    }
    #[test]
    fn test_point_norm() {
        let p = Point::new(3.0, 4.0, 0.0);
        assert!((p.norm() - 5.0).abs() < 1e-9);
    }
}
