//! Sphere representation and operations.
//!
//! This module defines [`Sphere`], a simple 3D sphere type used
//! throughout the crate for atom and probe representations.
//! It implements the [`Position`] trait so it can be used generically
//! alongside types like `Point` in neighbour-search code.
//!
//! # Examples
//!
//! ```
//! use garnet::geo::sphere::Sphere;
//! use garnet::geo::point::{Point, Position};
//!
//! let s = Sphere::new(Point::new(1.0, 2.0, 3.0), 1.4);
//! assert_eq!(s.position().x, 1.0);
//! ```

use nalgebra::Vector3;

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

impl Position for Sphere {
    fn position(&self) -> Vector3<f64> {
        self.coord.position()
    }
}

#[test]
fn test_sphere_fields() {
    let s = Sphere::new(Point::new(1.0, 2.0, 3.0), 1.4);
    assert_eq!(s.position().x, 1.0);
    assert_eq!(s.position().y, 2.0);
    assert_eq!(s.position().z, 3.0);
    assert_eq!(s.r, 1.4);
}
