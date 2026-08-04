use std::ops::{Add, Div, Mul, Sub};

/// A trait for types that have a 3D position
/// This is used for both Point and Sphere to guarantee
/// that we can access their x,y,z fields during neighbour searching
pub trait Position {
    fn x(&self) -> f64;
    fn y(&self) -> f64;
    fn z(&self) -> f64;
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
}

impl Default for Point {
    fn default() -> Self {
        Point::zero()
    }
}

/// Operator overloading for Point
impl Add for Point {
    type Output = Point;
    fn add(self, rhs: Point) -> Point {
        Point::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl Sub for Point {
    type Output = Point;
    fn sub(self, rhs: Point) -> Point {
        Point::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl Mul<f64> for Point {
    type Output = Point;
    fn mul(self, s: f64) -> Point {
        Point::new(self.x * s, self.y * s, self.z * s)
    }
}

impl Div<f64> for Point {
    type Output = Point;
    fn div(self, s: f64) -> Point {
        Point::new(self.x / s, self.y / s, self.z / s)
    }
}

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
