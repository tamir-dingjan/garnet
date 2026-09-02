use nalgebra::Vector3;

use crate::geo::point::Point;
use crate::geo::point::Position;
use crate::structure::residue::{Residue, ResidueName};
use std::ops::Deref;
use std::sync::{Arc, Weak};

const DEFAULT_ATOM_RADIUS: f64 = 0.0;

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

/// Default implementation for [`RawAtom`] used for testing
/// These values are neutral placeholders to represent an "unset" state
impl Default for RawAtom {
    fn default() -> Self {
        RawAtom {
            coor: Point::zero(), // Origin point (0,0,0)
            r: DEFAULT_ATOM_RADIUS,
            name: "".to_string(), // Empty string
            resn: ResidueName::Other(Default::default()),
            chain_id: ' ',     // Space character
            resi: 0,           // Residue numbers start at 1 in PDB file parsing, so 0 is "unset"
            het: false,        // Not a HETATOM record, as for most protein atoms
            model: 0,          // Model's are 1-indexed in PDB file parsing, so 0 is "unset"
            is_surface: false, // Not a surface atom, as for most protein atoms
        }
    }
}

impl Default for Atom {
    fn default() -> Self {
        Atom {
            coor: Point::zero(), // Origin point (0,0,0)
            r: DEFAULT_ATOM_RADIUS,
            name: "".to_string(), // Empty string
            resn: ResidueName::Other(Default::default()),
            chain_id: ' ',        // Space character
            resi: 0,              // Residue numbers start at 1 in PDB file parsing, so 0 is "unset"
            residue: Weak::new(), // Default Atom types do not have a Weak<Residue> ready to access
            het: false,           // Not a HETATOM record, as for most protein atoms
            model: 0,             // Model's are 1-indexed in PDB file parsing, so 0 is "unset"
            is_surface: false,    // Not a surface atom, as for most protein atoms
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
        assert_eq!(a.name, "");
        assert_eq!(a.resn, ResidueName::Other(Default::default()));
        assert_eq!(a.chain_id, ' ');
        assert_eq!(a.resi, 0);
        assert!(!a.het);
        assert_eq!(a.model, 0);
        assert!(!a.is_surface);
    }

    #[test]
    fn test_atom_position_via_coor() {
        let a = Atom {
            coor: Point::new(1.0, 2.0, 3.0),
            r: DEFAULT_ATOM_RADIUS,
            ..Atom::default()
        };
        assert_eq!(a.position().x, 1.0);
        assert_eq!(a.position().y, 2.0);
        assert_eq!(a.position().z, 3.0);
    }
}
