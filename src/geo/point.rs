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
//! use garnet::geo::point::Point;
//!
//! let a = Point::new(1.0, 2.0, 3.0);
//! let b = Point::new(4.0, 5.0, 6.0);
//! let sum = a + b;
//! assert_eq!(sum.x, 5.0);
//! ```

use std::ops::{Add, Div, Mul, Sub};

use nalgebra::{Rotation3, Vector3};

/// A trait for types that have a 3D position
/// This is used for both Point and Sphere to guarantee
/// that we can access their x,y,z fields during neighbour searching
pub trait Position {
    fn x(&self) -> f64;
    fn y(&self) -> f64;
    fn z(&self) -> f64;
}

impl<T: Position> Position for &T {
    fn x(&self) -> f64 {
        (*self).x()
    }

    fn y(&self) -> f64 {
        (*self).y()
    }

    fn z(&self) -> f64 {
        (*self).z()
    }
}

/// A 3D coordinate
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Point {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Point { x, y, z }
    }

    pub fn zero() -> Self {
        Point {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }

    /// Euclidean length
    pub fn norm(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    /// Add another position-like value to this point,
    /// used to accumulate coordinates in centroid calculation.
    pub fn add_position<P: Position>(self, other: &P) -> Point {
        Point::new(self.x + other.x(), self.y + other.y(), self.z + other.z())
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
        Point::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
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
        Point::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
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
        Point::new(self.x * s, self.y * s, self.z * s)
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
        Point::new(self.x / s, self.y / s, self.z / s)
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
    fn x(&self) -> f64 {
        self.x
    }
    fn y(&self) -> f64 {
        self.y
    }
    fn z(&self) -> f64 {
        self.z
    }
}

/// Point is used to represent a point in 3D space
/// To perform rotations and translations on Points using nalgebra types,
/// first convert the Point to a Vector3.
impl From<Point> for Vector3<f64> {
    fn from(p: Point) -> Self {
        Vector3::new(p.x, p.y, p.z)
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
        assert_eq!(p.x, 1.0);
        assert_eq!(p.y, 2.0);
        assert_eq!(p.z, 3.0);
    }

    #[test]
    fn test_point_add() {
        let a = Point::new(1.0, 2.0, 3.0);
        let b = Point::new(4.0, 5.0, 6.0);
        let r = a + b;
        assert_eq!(r.x, 5.0);
        assert_eq!(r.y, 7.0);
        assert_eq!(r.z, 9.0);
    }

    #[test]
    fn test_point_sub() {
        let a = Point::new(4.0, 5.0, 6.0);
        let b = Point::new(1.0, 2.0, 3.0);
        let r = a - b;
        assert_eq!(r.x, 3.0);
        assert_eq!(r.y, 3.0);
        assert_eq!(r.z, 3.0);
    }

    #[test]
    fn test_point_scale() {
        let p = Point::new(2.0, 4.0, 6.0);
        let r = p / 2.0;
        assert_eq!(r.x, 1.0);
        assert_eq!(r.y, 2.0);
        assert_eq!(r.z, 3.0);
    }

    #[test]
    fn test_point_div() {
        let p = Point::new(2.0, 4.0, 6.0);
        let r = p / 2.0;
        assert_eq!(r.x, 1.0);
        assert_eq!(r.y, 2.0);
        assert_eq!(r.z, 3.0);
    }
    #[test]
    fn test_point_norm() {
        let p = Point::new(3.0, 4.0, 0.0);
        assert!((p.norm() - 5.0).abs() < 1e-9);
    }
}
