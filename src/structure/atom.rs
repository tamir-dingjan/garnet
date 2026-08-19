use nalgebra::Vector3;

use crate::geo::point::Point;
use crate::geo::point::Position;
use crate::structure::residue::{Residue, ResidueName};
use std::ops::Deref;
use std::sync::{Arc, Weak};

/// Represent the pre-freeze mutable atom used
/// during parsing and surface marking.
#[derive(Debug, Clone, PartialEq)]
pub struct RawAtom {
    pub coor: Point,
    pub r: f64,
    pub name: String,
    pub resn: ResidueName,
    pub chain_id: char,
    pub resi: i32,
    pub het: bool,
    pub model: u32,
    pub is_surface: bool,
}

/// Represents an atom in a protein structure after
/// parsing and freezing.
#[derive(Debug)]
pub struct Atom {
    pub coor: Point,
    pub r: f64,                 // Van der Waals radius (A)
    pub name: String,           // Atom name, max 4 chars e.g. " CA "
    pub resn: ResidueName,      // Residue name, max 3 chars e.g. "ALA"
    pub chain_id: char,         // Chain ID, single char e.g. 'A'
    pub resi: i32,              // Residue sequence number
    pub residue: Weak<Residue>, // Use a Weak to avoid circular strong reference to Residue
    pub het: bool,              // True for HETATM records
    pub model: u32,             // Model number
    pub is_surface: bool,       // True if atom is on the surface
}

impl Deref for Atom {
    type Target = Point;
    fn deref(&self) -> &Self::Target {
        &self.coor
    }
}

impl Deref for RawAtom {
    type Target = Point;
    fn deref(&self) -> &Self::Target {
        &self.coor
    }
}

impl Default for RawAtom {
    fn default() -> Self {
        RawAtom {
            coor: Point::zero(),
            r: 1.7,
            name: "    ".to_string(),
            resn: ResidueName::ALA,
            chain_id: ' ',
            resi: 0,
            het: false,
            model: 1,
            is_surface: false,
        }
    }
}

impl Default for Atom {
    fn default() -> Self {
        Atom {
            coor: Point::zero(),
            r: 1.7,
            name: "    ".to_string(),
            resn: ResidueName::ALA,
            chain_id: ' ',
            resi: 0,
            residue: Weak::new(),
            het: false,
            model: 1,
            is_surface: false,
        }
    }
}

impl Position for RawAtom {
    fn position(&self) -> Vector3<f64> {
        self.coor.position()
    }
}

impl Position for Arc<Atom> {
    fn position(&self) -> Vector3<f64> {
        self.coor.position()
    }
}

impl Atom {
    // Access the residue that this atom belongs to.
    // We use a Weak to avoid circular references between Atom and Residue.
    pub fn residue(&self) -> Arc<Residue> {
        self.residue
            .upgrade()
            .expect("atom should always point to a residue")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geo::point::Point;

    #[test]
    fn test_atom_default_fields() {
        let a = Atom::default();
        assert_eq!(a.name, "    ");
        assert_eq!(a.resn, ResidueName::ALA);
        assert_eq!(a.chain_id, ' ');
        assert_eq!(a.resi, 0);
        assert!(!a.het);
        assert_eq!(a.model, 1);
        assert!(!a.is_surface);
    }

    #[test]
    fn test_atom_position_via_coor() {
        let a = Atom {
            coor: Point::new(1.0, 2.0, 3.0),
            r: 1.7,
            ..Atom::default()
        };
        assert_eq!(a.position().x, 1.0);
        assert_eq!(a.position().y, 2.0);
        assert_eq!(a.position().z, 3.0);
    }
}
