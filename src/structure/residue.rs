//! Residue represents a residue in a protein structure.

use arrayvec::ArrayString;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::hash::{Hash, Hasher};
use std::sync::Arc;

use crate::geo::ops::centroid;
use crate::geo::point::Point;
use crate::structure::atom::Atom;

/// Residue names
#[derive(Debug, Hash, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum ResidueName {
    ALA,
    ARG,
    ASN,
    ASP,
    CYS,
    GLU,
    GLN,
    GLY,
    HIS,
    ILE,
    LEU,
    LYS,
    MET,
    PHE,
    PRO,
    SER,
    THR,
    TRP,
    TYR,
    VAL,
    Other(ArrayString<4>),
}

impl std::fmt::Display for ResidueName {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

impl ResidueName {
    pub fn parse(s: &str) -> Self {
        match s {
            "ALA" => Self::ALA,
            "ARG" => Self::ARG,
            "ASN" => Self::ASN,
            "ASP" => Self::ASP,
            "CYS" => Self::CYS,
            "GLU" => Self::GLU,
            "GLN" => Self::GLN,
            "GLY" => Self::GLY,
            "HIS" => Self::HIS,
            "ILE" => Self::ILE,
            "LEU" => Self::LEU,
            "LYS" => Self::LYS,
            "MET" => Self::MET,
            "PHE" => Self::PHE,
            "PRO" => Self::PRO,
            "SER" => Self::SER,
            "THR" => Self::THR,
            "TRP" => Self::TRP,
            "TYR" => Self::TYR,
            "VAL" => Self::VAL,
            _ => Self::Other(ArrayString::from(s).unwrap_or_default()),
        }
    }
    pub fn as_str(&self) -> &str {
        match self {
            Self::ALA => "ALA",
            Self::ARG => "ARG",
            Self::ASN => "ASN",
            Self::ASP => "ASP",
            Self::CYS => "CYS",
            Self::GLU => "GLU",
            Self::GLN => "GLN",
            Self::GLY => "GLY",
            Self::HIS => "HIS",
            Self::ILE => "ILE",
            Self::LEU => "LEU",
            Self::LYS => "LYS",
            Self::MET => "MET",
            Self::PHE => "PHE",
            Self::PRO => "PRO",
            Self::SER => "SER",
            Self::THR => "THR",
            Self::TRP => "TRP",
            Self::TYR => "TYR",
            Self::VAL => "VAL",
            Self::Other(s) => s.as_str(),
        }
    }
    pub fn one_letter(&self) -> Option<char> {
        match self {
            Self::ALA => Some('A'),
            Self::ARG => Some('R'),
            Self::ASN => Some('N'),
            Self::ASP => Some('D'),
            Self::CYS => Some('C'),
            Self::GLU => Some('E'),
            Self::GLN => Some('Q'),
            Self::GLY => Some('G'),
            Self::HIS => Some('H'),
            Self::ILE => Some('I'),
            Self::LEU => Some('L'),
            Self::LYS => Some('K'),
            Self::MET => Some('M'),
            Self::PHE => Some('F'),
            Self::PRO => Some('P'),
            Self::SER => Some('S'),
            Self::THR => Some('T'),
            Self::TRP => Some('W'),
            Self::TYR => Some('Y'),
            Self::VAL => Some('V'),
            Self::Other(_) => None,
        }
    }

    pub fn min_heavy_atoms(&self) -> Option<usize> {
        match self {
            Self::ALA => Some(5),
            Self::ARG => Some(11),
            Self::ASN => Some(8),
            Self::ASP => Some(8),
            Self::CYS => Some(6),
            Self::GLU => Some(9),
            Self::GLN => Some(9),
            Self::GLY => Some(4),
            Self::HIS => Some(10),
            Self::ILE => Some(8),
            Self::LEU => Some(8),
            Self::LYS => Some(9),
            Self::MET => Some(8),
            Self::PHE => Some(11),
            Self::PRO => Some(7),
            Self::SER => Some(6),
            Self::THR => Some(7),
            Self::TRP => Some(14),
            Self::TYR => Some(12),
            Self::VAL => Some(7),
            Self::Other(_) => None,
        }
    }
}

/// Represents a residue in a protein structure.
/// Identity is determined by (chain_id, resi) - residue name
/// can differ across alignments
#[derive(Debug, Clone)]
pub struct Residue {
    pub resn: ResidueName,         // Residue name, max 3 chars e.g. "ALA"
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

    pub fn one_letter(&self) -> Option<char> {
        self.resn.one_letter()
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
            resn: ResidueName::ALA,
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

/// Lookup to check if a residue contains an aromatic ring
pub fn residue_contains_aromatic_ring(resn: &ResidueName) -> bool {
    matches!(
        resn,
        ResidueName::HIS | ResidueName::PHE | ResidueName::TYR | ResidueName::TRP
    )
}

/// Returns the aromatic ring descriptor for the given residue
/// The descriptor is placed at the centroid of the aromatic ring atoms,
/// which are selected by atom name.
pub fn aromatic_ring_centroid(atom: &Arc<Atom>) -> Option<Point> {
    // Identify the aromatic ring atoms by the residue name
    let ring_atoms = match atom.resn {
        ResidueName::HIS => vec!["CG", "ND1", "CD2", "CE1", "NE2"],
        ResidueName::PHE => vec!["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
        ResidueName::TYR => vec!["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
        ResidueName::TRP => vec!["CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
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
            resn: ResidueName::GLY,
            resi: 42,
            ..Residue::default()
        };
        let r2 = Residue {
            chain_id: 'A',
            resn: ResidueName::ALA,
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
