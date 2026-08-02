use std::ops::{Add, Deref, Div, Mul, Sub};

/// XS is an "extra small" constant used as a zero guard
pub const XS: f64 = 1e-60;

/// A trait for types that have a 3D position
pub trait HasPosition {
    fn x(&self) -> f64;
    fn y(&self) -> f64;
    fn z(&self) -> f64;
}

/// A 3D coordinate
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Coor {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Coor {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Coor { x, y, z }
    }

    pub fn zero() -> Self {
        Coor {
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

impl Default for Coor {
    fn default() -> Self {
        Coor::zero()
    }
}

/// Operator overloading for Coor
impl Add for Coor {
    type Output = Coor;
    fn add(self, rhs: Coor) -> Coor {
        Coor::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl Sub for Coor {
    type Output = Coor;
    fn sub(self, rhs: Coor) -> Coor {
        Coor::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl Mul<f64> for Coor {
    type Output = Coor;
    fn mul(self, s: f64) -> Coor {
        Coor::new(self.x * s, self.y * s, self.z * s)
    }
}

impl Div<f64> for Coor {
    type Output = Coor;
    fn div(self, s: f64) -> Coor {
        Coor::new(self.x / s, self.y / s, self.z / s)
    }
}

impl HasPosition for Coor {
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

/// A sphere, has a 3D central coordinate and a radius
/// Used as the base for Atom and Probe
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sphere {
    pub coor: Coor,
    pub r: f64,
}

impl Sphere {
    pub fn new(coor: Coor, r: f64) -> Self {
        Sphere { coor, r }
    }
}

impl Deref for Sphere {
    type Target = Coor;
    fn deref(&self) -> &Self::Target {
        &self.coor
    }
}

impl Default for Sphere {
    fn default() -> Self {
        Sphere::new(Coor::default(), 0.0)
    }
}

impl HasPosition for &Sphere {
    fn x(&self) -> f64 {
        self.coor.x
    }
    fn y(&self) -> f64 {
        self.coor.y
    }
    fn z(&self) -> f64 {
        self.coor.z
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_coor_fields() {
        let c = Coor::new(1.0, 2.0, 3.0);
        assert_eq!(c.x, 1.0);
        assert_eq!(c.y, 2.0);
        assert_eq!(c.z, 3.0);
    }

    #[test]
    fn test_coor_add() {
        let a = Coor::new(1.0, 2.0, 3.0);
        let b = Coor::new(4.0, 5.0, 6.0);
        let r = a + b;
        assert_eq!(r.x, 5.0);
        assert_eq!(r.y, 7.0);
        assert_eq!(r.z, 9.0);
    }

    #[test]
    fn test_coor_sub() {
        let a = Coor::new(4.0, 5.0, 6.0);
        let b = Coor::new(1.0, 2.0, 3.0);
        let r = a - b;
        assert_eq!(r.x, 3.0);
        assert_eq!(r.y, 3.0);
        assert_eq!(r.z, 3.0);
    }

    #[test]
    fn test_coor_scale() {
        let c = Coor::new(2.0, 4.0, 6.0);
        let r = c / 2.0;
        assert_eq!(r.x, 1.0);
        assert_eq!(r.y, 2.0);
        assert_eq!(r.z, 3.0);
    }

    #[test]
    fn test_coor_div() {
        let c = Coor::new(2.0, 4.0, 6.0);
        let r = c / 2.0;
        assert_eq!(r.x, 1.0);
        assert_eq!(r.y, 2.0);
        assert_eq!(r.z, 3.0);
    }
    #[test]
    fn test_coor_norm() {
        let c = Coor::new(3.0, 4.0, 0.0);
        assert!((c.norm() - 5.0).abs() < 1e-9);
    }

    #[test]
    fn test_sphere_fields() {
        let s = Sphere::new(Coor::new(1.0, 2.0, 3.0), 1.4);
        assert_eq!(s.x, 1.0);
        assert_eq!(s.y, 2.0);
        assert_eq!(s.z, 3.0);
        assert_eq!(s.r, 1.4);
    }
}
