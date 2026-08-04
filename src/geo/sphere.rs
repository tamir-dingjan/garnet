use crate::geo::point::{Point, Position};

use std::ops::Deref;

/// A sphere, has a 3D central coordinate and a radius
/// Used as the base for Atom and Probe
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sphere {
    pub coord: Point,
    pub r: f64,
}

impl Sphere {
    pub fn new(coord: Point, r: f64) -> Self {
        Sphere { coord, r }
    }
}

impl Deref for Sphere {
    type Target = Point;
    fn deref(&self) -> &Self::Target {
        &self.coord
    }
}

impl Default for Sphere {
    fn default() -> Self {
        Sphere::new(Point::default(), 0.0)
    }
}

impl Position for &Sphere {
    fn x(&self) -> f64 {
        self.coord.x
    }
    fn y(&self) -> f64 {
        self.coord.y
    }
    fn z(&self) -> f64 {
        self.coord.z
    }
}

#[test]
fn test_sphere_fields() {
    let s = Sphere::new(Point::new(1.0, 2.0, 3.0), 1.4);
    assert_eq!(s.x, 1.0);
    assert_eq!(s.y, 2.0);
    assert_eq!(s.z, 3.0);
    assert_eq!(s.r, 1.4);
}
