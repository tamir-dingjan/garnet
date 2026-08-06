use crate::geo::point::Point;

/// XS is an "extra small" constant used as a zero guard
pub const XS: f64 = 1e-60;

/// Euclidean distance between two points
pub fn distance(a: &Point, b: &Point) -> f64 {
    let dx = a.x - b.x;
    let dy = a.y - b.y;
    let dz = a.z - b.z;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Dot product of two vectors
pub fn dot(a: &Point, b: &Point) -> f64 {
    a.x * b.x + a.y * b.y + a.z * b.z
}

/// Cross product of two vectors
pub fn cross(a: &Point, b: &Point) -> Point {
    Point::new(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
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
        assert!((z.x - 0.0).abs() < 1e-9);
        assert!((z.y - 0.0).abs() < 1e-9);
        assert!((z.z - 1.0).abs() < 1e-9);
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
        assert_eq!(n.x, 0.0);
        assert_eq!(n.y, 0.0);
        assert_eq!(n.z, 0.0);
    }
}
