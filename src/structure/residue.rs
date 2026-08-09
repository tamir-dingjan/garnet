//! Residue represents a residue in a protein structure.

use std::cmp::Ordering;
use std::hash::{Hash, Hasher};
use std::sync::Arc;

use crate::geo::ops::centroid;
use crate::geo::point::Point;
use crate::structure::atom::Atom;

/// Represents a residue in a protein structure.
/// Identity is determined by (chain_id, resi) - residue name
/// can differ across alignments
#[derive(Debug, Clone)]
pub struct Residue {
    pub resn: String,              // Residue name, max 3 chars e.g. "ALA"
    pub chain_id: char,            // Chain identifier, e.g. 'A' for chain A
    pub resi: i32,                 // Residue sequence number
    pub atoms: Vec<Arc<Atom>>,     // Atoms in this residue
    pub calpha: Option<Arc<Atom>>, // C-alpha atom, if available
    pub cons: u32,                 // Conservation count
    pub defective: bool,           // Defective residues are missing atoms
}

impl Residue {
    /// A residue key used for comparison, consisting of (chain_id, resi)
    fn key(&self) -> (char, i32) {
        (self.chain_id, self.resi)
    }

    pub fn one_letter(&self) -> char {
        let residue = self.resn.as_str();
        three_to_one(residue)
    }

    pub fn atoms(&self) -> &[Arc<Atom>] {
        &self.atoms
    }

    pub fn calpha(&self) -> Option<&Arc<Atom>> {
        self.calpha.as_ref()
    }
}

impl Default for Residue {
    fn default() -> Self {
        Residue {
            chain_id: ' ',
            resn: "   ".to_string(),
            resi: 0,
            atoms: Vec::new(),
            calpha: None,
            cons: 0,
            defective: false,
        }
    }
}

impl PartialEq for Residue {
    fn eq(&self, other: &Self) -> bool {
        self.key() == other.key()
    }
}

impl Eq for Residue {}

impl PartialOrd for Residue {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Compares residues by their key (chain_id, resi)
/// This uses the cmp() method on the key tuple
impl Ord for Residue {
    fn cmp(&self, other: &Self) -> Ordering {
        self.key().cmp(&other.key())
    }
}

impl Hash for Residue {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.key().hash(state);
    }
}

/// Convert the three-letter residue name to the one-letter residue name
/// Returns '?' if the residue name is not recognized
pub fn three_to_one(residue_name: &str) -> char {
    match residue_name.trim() {
        "ALA" => 'A',
        "ARG" => 'R',
        "ASN" => 'N',
        "ASP" => 'D',
        "CYS" => 'C',
        "GLU" => 'E',
        "GLN" => 'Q',
        "GLY" => 'G',
        "HIS" => 'H',
        "ILE" => 'I',
        "LEU" => 'L',
        "LYS" => 'K',
        "MET" => 'M',
        "PHE" => 'F',
        "PRO" => 'P',
        "SER" => 'S',
        "THR" => 'T',
        "TRP" => 'W',
        "TYR" => 'Y',
        "VAL" => 'V',
        _ => '?',
    }
}

/// Lookup for the minimum number of heavy atoms required for a residue.
/// Defaults to 0 if the residue name is not recognized
pub fn min_heavy_atoms(residue: &str) -> usize {
    match residue.trim() {
        "ALA" => 5,
        "ARG" => 11,
        "ASN" => 8,
        "ASP" => 8,
        "CYS" => 6,
        "GLU" => 9,
        "GLN" => 9,
        "GLY" => 4,
        "HIS" => 10,
        "ILE" => 8,
        "LEU" => 8,
        "LYS" => 9,
        "MET" => 8,
        "PHE" => 11,
        "PRO" => 7,
        "SER" => 6,
        "THR" => 7,
        "TRP" => 14,
        "TYR" => 12,
        "VAL" => 7,
        _ => 0,
    }
}

/// Lookup to check if a residue contains an aromatic ring
pub fn residue_contains_aromatic_ring(resn: &str) -> bool {
    matches!(resn, "HIS" | "PHE" | "TYR" | "TRP")
}

/// Returns the aromatic ring descriptor for the given residue
/// The descriptor is placed at the centroid of the aromatic ring atoms,
/// which are selected by atom name.
pub fn aromatic_ring_centroid(atom: &Arc<Atom>) -> Option<Point> {
    // Identify the aromatic ring atoms by the residue name
    let ring_atoms = match atom.resn.as_str() {
        "HIS" => vec!["CG", "ND1", "CD2", "CE1", "NE2"],
        "PHE" => vec!["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
        "TYR" => vec!["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
        "TRP" => vec!["CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
        _ => return None,
    };

    // Collect the aromatic ring atom coordinates
    let ring_atoms = atom
        .residue()
        .atoms()
        .iter()
        .filter(|a| ring_atoms.contains(&a.name.as_str()))
        .map(|a| a.coor)
        .collect::<Vec<_>>();
    if ring_atoms.is_empty() {
        return None;
    }

    // Return the ring centroid, or a None if the centroid fails to compute
    centroid(&ring_atoms)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_residue_ordering_by_resi() {
        let r1 = Residue {
            chain_id: 'A',
            resi: 1,
            ..Residue::default()
        };
        let r2 = Residue {
            chain_id: 'A',
            resi: 2,
            ..Residue::default()
        };
        assert!(r1 < r2);
    }

    #[test]
    fn test_residue_chain_sorts_before_resi() {
        let r_a = Residue {
            chain_id: 'A',
            resi: 99,
            ..Residue::default()
        };
        let r_b = Residue {
            chain_id: 'B',
            resi: 1,
            ..Residue::default()
        };
        assert!(r_a < r_b);
    }

    #[test]
    fn test_residue_equality_by_chain_and_resi() {
        let r1 = Residue {
            chain_id: 'A',
            resn: "GLY".to_string(),
            resi: 42,
            ..Residue::default()
        };
        let r2 = Residue {
            chain_id: 'A',
            resn: "ALA".to_string(),
            resi: 42,
            ..Residue::default()
        };
        assert_eq!(r1, r2);
    }

    #[test]
    fn test_residue_hashable() {
        use std::collections::HashSet;
        let mut s = HashSet::new();
        s.insert(Residue {
            chain_id: 'A',
            resi: 1,
            ..Residue::default()
        });
        s.insert(Residue {
            chain_id: 'A',
            resi: 1,
            ..Residue::default()
        });
        s.insert(Residue {
            chain_id: 'A',
            resi: 2,
            ..Residue::default()
        });
        assert_eq!(s.len(), 2);
    }
}
