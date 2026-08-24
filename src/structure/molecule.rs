/// Molecule is a type that allows reading from a PDB file and storing the
/// atoms and residues. Because atoms and residues refer to eachother,
/// they are stored using 'Arc' with a cyclic reference. The ownership
/// model is:
/// - Arc<Atom> owns a Weak<Residue> reference to its residue
/// - Arc<Residue> owns the Vec<Arc<Atom>> containing its atoms
/// - Molecule owns the Vec<Arc<Residue>> containing all its residues
///
/// Molecule is parsed from a PDB file. Atom data is read into RawAtom types
/// and once surface exposure is assigned, "frozen" into Atom types.
/// This conversion wires up the Atom<->Residue cyclic references in
/// finalize_structure().
///
/// Molecule also has a function to select binding sites and mask off
/// atoms that are not part of the binding site by setting them as non-surface
/// atoms.
use crate::geo::grid::Grid;
use crate::geo::ops::distance;
use crate::geo::point::Position;
use crate::io::pdb::{ParseOptions, parse_pdb_file};
use crate::structure::atom::{Atom, RawAtom};
use crate::structure::chainset::ChainSet;
use crate::structure::residue::{Residue, ResidueName};
use anyhow::{Context, Result};

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

const ATOM_OVERLAP_TOLERANCE: f64 = 1e-3;

const HELPER_GRID_SPACING: f64 = 2.0;

/// C-alpha atom name
pub const CA_ATOM_NAME: &str = "CA";

pub struct Molecule {
    pub path: String,                // Path to the PDB file
    pub opts: ParseOptions,          // Parsing options for the PDB file
    pub name: String,                // File name without extension
    pub atoms: Vec<Arc<Atom>>,       // All atoms (includes both ATOM and HETATM records)
    pub residues: Vec<Arc<Residue>>, // Residues built from atoms
    pub grid: Option<Grid<RawAtom>>, // Grid for spatial queries over atoms, built on demand
    pub raw_atoms: Vec<RawAtom>,
}

impl Molecule {
    pub fn new(name: impl Into<String>) -> Self {
        let name = name.into();
        Molecule {
            path: String::new(),
            opts: ParseOptions::default(),
            name,
            atoms: Vec::new(),
            residues: Vec::new(),
            grid: None,
            raw_atoms: Vec::new(),
        }
    }

    /// Read a PDB file from `path` using `opts`.
    pub fn read_pdb(&mut self, path: &str, opts: ParseOptions) -> Result<()> {
        self.path = path.to_string();
        self.opts = opts;

        self.raw_atoms = parse_pdb_file(path, &self.opts)?;
        self.atoms.clear();
        self.residues.clear();
        self.grid = None;

        Ok(())
    }

    /// Collect the chain IDs of all atoms in the molecule.
    pub fn chain_ids(&self) -> ChainSet {
        self.raw_atoms.iter().map(|atom| atom.chain_id).collect()
    }

    /// Build the `grid` for spatial queries over atoms.
    pub fn build_grid(&mut self, spacing: f64) {
        let mut g = Grid::new(spacing);
        g.build(&self.raw_atoms);
        self.grid = Some(g);
    }

    pub fn raw_atoms(&self) -> &[RawAtom] {
        &self.raw_atoms
    }

    pub fn raw_atoms_mut(&mut self) -> &mut [RawAtom] {
        &mut self.raw_atoms
    }

    fn get_residue_atoms(&self, resn: &ResidueName, resi: i32, chain_id: char) -> Vec<&RawAtom> {
        self.raw_atoms
            .iter()
            .filter(|a| &a.resn == resn && a.resi == resi && a.chain_id == chain_id)
            .collect()
    }

    /// Since most ligands are present only as HETATM lines, we load a helper
    /// version of this molecule that does read HETATM lines to select the
    /// binding site residues.
    fn get_binding_site_residues(
        &self,
        resn: &ResidueName,
        resi: i32,
        chain_id: char,
        distance: f64,
    ) -> Result<HashSet<(char, i32, ResidueName)>> {
        let mut het_opts = self.opts.clone();
        het_opts.include_hetero = true;
        het_opts.chain_ids = None;

        let mut helper = Molecule::new("helper");
        helper
            .read_pdb(&self.path, het_opts)
            .with_context(|| "Failed to read PDB for binding site residues")?;

        helper.build_grid(HELPER_GRID_SPACING);

        let binding_site_center = helper.get_residue_atoms(resn, resi, chain_id);
        if binding_site_center.is_empty() {
            anyhow::bail!(
                "No atoms found for binding site residue {} {} {}",
                resn,
                resi,
                chain_id
            );
        }

        // Residue keys are stable across the whole-structure helper and the
        // chain-selected molecule, unlike atom vector indices.
        let binding_site_residues = binding_site_center
            .iter()
            .flat_map(|atom| {
                helper.grid.as_ref().unwrap().neighbors_within(
                    &helper.raw_atoms,
                    atom.coor,
                    distance,
                )
            })
            .map(|idx| {
                let atom = &helper.raw_atoms[idx];
                (atom.chain_id, atom.resi, atom.resn.clone())
            })
            .collect();

        Ok(binding_site_residues)
    }

    /// Hides atoms beyond the given distance from the binding site center.
    ///
    /// Atoms that are beyond the given distance from the binding site center will
    /// have their `is_surface` flag set to `false`. This prevents Descriptor
    /// generation for these atoms.
    pub fn hide_atoms_beyond_binding_site(
        &mut self,
        resn: ResidueName,
        resi: i32,
        chain_id: char,
        distance: f64,
    ) -> Result<()> {
        // Get the binding site residues
        let binding_site_residues =
            self.get_binding_site_residues(&resn, resi, chain_id, distance)?;

        // Hide atoms beyond the binding site by marking them as not surface-exposed
        for atom in self.raw_atoms_mut() {
            let residue_key = (atom.chain_id, atom.resi, atom.resn.clone());
            if !binding_site_residues.contains(&residue_key) {
                atom.is_surface = false;
            }
        }
        Ok(())
    }

    /// Takes the raw atoms and groups them into residues,
    /// including marking of defective residues.
    /// This is called after surface exposure has been computed.
    pub fn finalize_structure(&mut self) {
        // Group the atoms by their chain/residue/name key
        // Using the group_index HashMap lets us preserve residue ordering
        // in the grouped Vec
        let mut grouped: Vec<((char, i32, ResidueName), Vec<usize>)> = Vec::new();

        let mut group_index: HashMap<(char, i32, ResidueName), usize> = HashMap::new();

        for (idx, atom) in self.raw_atoms.iter().enumerate() {
            let key = (atom.chain_id, atom.resi, atom.resn.clone());
            if let Some(existing) = group_index.get(&key) {
                grouped[*existing].1.push(idx);
            } else {
                let slot = grouped.len();
                group_index.insert(key.clone(), slot);
                grouped.push((key, vec![idx]));
            }
        }

        // Setup for populating the frozen collections
        let defective_keys = self.find_defective_residue_keys();
        let mut frozen_atoms: Vec<Option<Arc<Atom>>> = vec![None; self.raw_atoms.len()];
        let mut frozen_residues: Vec<Arc<Residue>> = Vec::new();

        for ((chain_id, resi, resn), atom_indices) in grouped {
            let defective = defective_keys.contains(&(chain_id, resi));
            // Define the atoms and their residue together
            // The Arc<Atom> gets a Weak<Residue> reference
            // The Arc<Residue> owns the Vec<Arc<Atom>>
            // Mol owns the Vec<Arc<Residue>>
            let residue = Arc::new_cyclic(|weak_residue| {
                let atoms: Vec<Arc<Atom>> = atom_indices
                    .iter()
                    .map(|&raw_idx| {
                        let raw = &self.raw_atoms[raw_idx];
                        Arc::new(Atom {
                            coor: raw.coor,
                            r: raw.r,
                            name: raw.name.clone(),
                            resn: raw.resn.clone(),
                            chain_id: raw.chain_id,
                            resi: raw.resi,
                            residue: weak_residue.clone(),
                            het: raw.het,
                            model: raw.model,
                            is_surface: raw.is_surface,
                        })
                    })
                    .collect();
                let calpha = atoms.iter().find(|atom| atom.name == CA_ATOM_NAME).cloned();

                Residue {
                    resn: resn.clone(),
                    chain_id,
                    resi,
                    atoms,
                    calpha,
                    cons: 0,
                    defective,
                }
            });
            // Copy the residue's atoms into frozen atoms
            for (atom_slot, &raw_idx) in residue.atoms().iter().zip(atom_indices.iter()) {
                frozen_atoms[raw_idx] = Some(atom_slot.clone());
            }
            // Move the residue into frozen residues
            frozen_residues.push(residue);
        }

        // Move the frozen collections into the molecule
        self.atoms = frozen_atoms
            .into_iter()
            .map(|atom| atom.expect("every raw atom should freeze"))
            .collect();
        self.residues = frozen_residues;
    }

    fn find_defective_residue_keys(&self) -> HashSet<(char, i32)> {
        let mut grouped: HashMap<(char, i32, ResidueName), Vec<usize>> = HashMap::new();
        let mut order: Vec<(char, i32, ResidueName)> = Vec::new();
        let mut defective = HashSet::new();

        for (idx, atom) in self.raw_atoms.iter().enumerate() {
            let key = (atom.chain_id, atom.resi, atom.resn.clone());
            if !grouped.contains_key(&key) {
                order.push(key.clone());
            }
            grouped.entry(key).or_default().push(idx);
        }

        // Mark residues missing atoms as defective
        // Residues with non-normal names are not marked defective
        for key in &order {
            if let Some(num_heavy_atoms) = &key.2.min_heavy_atoms()
                && grouped[key].len() < *num_heavy_atoms
            {
                defective.insert((key.0, key.1));
            }
        }

        // Mark overlapping residues as defective
        // This check marks the first residue in each overlapping pair
        // as defective
        for pair in order.windows(2) {
            let left_atoms = &grouped[&pair[0]]
                .iter()
                .map(|&idx| &self.raw_atoms[idx])
                .collect::<Vec<_>>();
            let right_atoms = &grouped[&pair[1]]
                .iter()
                .map(|&idx| &self.raw_atoms[idx])
                .collect::<Vec<_>>();
            if self.atoms_overlap(
                left_atoms.iter(),
                right_atoms.iter(),
                ATOM_OVERLAP_TOLERANCE,
            ) {
                defective.insert((pair[1].0, pair[1].1));
            }
        }

        defective
    }

    /// Utility to check if two collections of atoms overlap spatially
    /// within the given tolerance
    /// The T: Position trait bound allows this to be used for RawAtom and Atom types
    pub fn atoms_overlap<'a, T>(
        &self,
        atoms1: impl Iterator<Item = &'a T>,
        atoms2: impl Iterator<Item = &'a T> + Clone,
        tolerance: f64,
    ) -> bool
    where
        T: Position + 'a,
    {
        for atom1 in atoms1 {
            // Clone for atoms2 resets the iterator to its original position for each loop
            for atom2 in atoms2.clone() {
                if distance(&atom1, &atom2) < tolerance {
                    return true;
                }
            }
        }
        false
    }
}

#[cfg(test)]
mod tests {
    fn pdb_atom_line(
        record: &str,
        serial: usize,
        atom_name: &str,
        residue_name: &str,
        chain_id: char,
        residue_number: i32,
        x: f64,
    ) -> String {
        format!(
            "{record:<6}{serial:>5} {atom_name:^4} {residue_name:>3} {chain_id}{residue_number:>4}    {x:>8.3}{y:>8.3}{z:>8.3}  1.00 20.00           C ",
            y = 0.0,
            z = 0.0,
        )
    }

    #[test]
    fn binding_site_on_another_chain_keeps_full_selected_chain_residue() {
        use super::*;
        use std::io::Write;
        use tempfile::Builder;

        let mut pdb = Builder::new().suffix(".pdb").tempfile().unwrap();
        let lines = [
            pdb_atom_line("ATOM", 1, "CA", "ALA", 'A', 1, 0.0),
            pdb_atom_line("ATOM", 2, "CB", "ALA", 'A', 1, 4.0),
            pdb_atom_line("ATOM", 3, "CA", "GLY", 'A', 2, 10.0),
            pdb_atom_line("HETATM", 4, "C1", "LIG", 'Z', 9, 0.5),
        ];
        writeln!(pdb, "{}", lines.join("\n")).unwrap();

        let mut molecule = Molecule::new("selected-chain");
        molecule
            .read_pdb(
                pdb.path().to_str().unwrap(),
                ParseOptions {
                    chain_ids: Some(ChainSet::single('A')),
                    ..ParseOptions::default()
                },
            )
            .unwrap();
        for atom in molecule.raw_atoms_mut() {
            atom.is_surface = true;
        }

        molecule
            .hide_atoms_beyond_binding_site(ResidueName::parse("LIG"), 9, 'Z', 2.0)
            .unwrap();

        assert!(molecule.raw_atoms()[0].is_surface);
        assert!(molecule.raw_atoms()[1].is_surface);
        assert!(!molecule.raw_atoms()[2].is_surface);
    }

    #[test]
    fn finalize_structure_links_atoms_and_residues() {
        use super::*;
        use crate::geo::point::Point;
        use crate::structure::atom::RawAtom;
        use std::sync::Arc;

        let mut mol = Molecule::new("mini");
        mol.raw_atoms = vec![
            RawAtom {
                coor: Point::new(0.0, 0.0, 0.0),
                r: 1.5,
                name: "N".into(),
                resn: ResidueName::parse("GLY"),
                chain_id: 'A',
                resi: 7,
                het: false,
                model: 1,
                is_surface: false,
            },
            RawAtom {
                coor: Point::new(1.0, 0.0, 0.0),
                r: 1.5,
                name: "CA".into(),
                resn: ResidueName::parse("GLY"),
                chain_id: 'A',
                resi: 7,
                het: false,
                model: 1,
                is_surface: false,
            },
            RawAtom {
                coor: Point::new(2.0, 0.0, 0.0),
                r: 1.5,
                name: "C".into(),
                resn: ResidueName::parse("GLY"),
                chain_id: 'A',
                resi: 7,
                het: false,
                model: 1,
                is_surface: false,
            },
        ];
        mol.finalize_structure();

        assert_eq!(mol.atoms.len(), 3);
        assert_eq!(mol.residues.len(), 1);

        let atom = mol.atoms[1].clone();
        let residue = atom.residue();

        assert_eq!(residue.resi, 7);
        assert_eq!(residue.atoms().len(), 3);
        assert!(Arc::ptr_eq(residue.calpha().unwrap(), &atom));
    }
}
