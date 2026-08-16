use nalgebra::Vector3;

use crate::geo::point::{Point, Position};

/// XS is an "extra small" constant used as a zero guard
pub const XS: f64 = 1e-60;

/// Euclidean distance between two points
pub fn distance(a: &impl Position, b: &impl Position) -> f64 {
    let delta: Vector3<f64> = a.position() - b.position();
    delta.norm()
}

/// Dot product of two vectors
pub fn dot(a: &Point, b: &Point) -> f64 {
    a.position().dot(&b.position())
}

/// Cross product of two vectors
pub fn cross(a: &Point, b: &Point) -> Point {
    Point::new(
        a.position().y * b.position().z - a.position().z * b.position().y,
        a.position().z * b.position().x - a.position().x * b.position().z,
        a.position().x * b.position().y - a.position().y * b.position().x,
    )
}

/// Return unit vector
/// If the normalized vector is close to zero,
/// return the zero vector instead.
pub fn normalize(v: &Point) -> Point {
    let n: f64 = v.norm();
    if n < XS { Point::zero() } else { *v / n }
}

/// Compute the angle between two vectors in radians
pub fn angle(a: &Point, b: &Point) -> f64 {
    let dot_product = dot(a, b);
    let n = a.norm() * b.norm();
    if n < XS {
        0.0
    } else {
        (dot_product / n).clamp(-1.0, 1.0)
    }
    .acos()
}

/// Compute the volume of a tetrahedron defined by four points
///
/// The volume is computed as:
/// V = ( (b-a) × (c-a) ) · (d-a) / 6
pub fn tetra_volume(a: &Point, b: &Point, c: &Point, d: &Point) -> f64 {
    let bma = b - a;
    let cma = c - a;
    let dma = d - a;

    let cross = cross(&bma, &cma);
    (cross.position().x * dma.position().x
        + cross.position().y * dma.position().y
        + cross.position().z * dma.position().z)
        / 6.0
}

pub fn centroid(points: &[impl Position]) -> Option<Point> {
    if points.is_empty() {
        return None;
    }
    let n = points.len() as f64;
    Some(
        points
            .iter()
            .fold(Point::zero(), |acc, p| acc.add_position(p))
            / n,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geo::point::Point;

    #[test]
    fn test_distance_3_4_0() {
        let a = Point::new(0.0, 0.0, 0.0);
        let b = Point::new(3.0, 4.0, 0.0);
        assert!((distance(&a, &b) - 5.0).abs() < 1e-9);
    }

    #[test]
    fn test_dot_perpendicular_is_zero() {
        let a = Point::new(1.0, 0.0, 0.0);
        let b = Point::new(0.0, 1.0, 0.0);
        assert_eq!(dot(&a, &b), 0.0);
    }

    #[test]
    fn test_dot_parallel_is_one() {
        let a = Point::new(1.0, 0.0, 0.0);
        assert!((dot(&a, &a) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_cross_x_y_gives_z() {
        let x = Point::new(1.0, 0.0, 0.0);
        let y = Point::new(0.0, 1.0, 0.0);
        let z = cross(&x, &y);
        assert!((z.position().x - 0.0).abs() < 1e-9);
        assert!((z.position().y - 0.0).abs() < 1e-9);
        assert!((z.position().z - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_normalize_produces_unit_vector() {
        let v = Point::new(3.0, 4.0, 0.0);
        let n = normalize(&v);
        assert!((n.norm() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_normalize_zero_vector_returns_zero() {
        let v = Point::new(0.0, 0.0, 0.0);
        let n = normalize(&v);
        assert_eq!(n.position().x, 0.0);
        assert_eq!(n.position().y, 0.0);
        assert_eq!(n.position().z, 0.0);
    }

    #[test]
    fn test_tetra_volume() {
        let coords = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
            Point::new(0.0, 0.0, 1.0),
        ];
        let v = tetra_volume(&coords[0], &coords[1], &coords[2], &coords[3]);
        assert!((v - 0.16666667).abs() < 1e-6);
    }
}
