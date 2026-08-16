//! Surface detection uses the Shrake-Rupley algorithm.
//! Each atom is modeled as a sphere with a van der Waals radius `r`.
//! The probe is modeled as a sphere with a probe radius `p`.
//! Typically the probe represents a water molecule, with `p=1.4`.
//! The probe center traces out the solvent-accessible surface area
//! at a shell of distance `r+p` from each atom center.
//!
//! The algorithm starts by building a sphere around each atom center
//! with a radius of `r+p`, and choosing candidate probe positions
//! uniformly distributed on the sphere surface. These candidate positions
//! are all the points where a probe could contact the atom.
//! The uniform distribution will be placed with a Fibbonacci spiral.
//!
//! Then for each candidate position `P`, check whether a probe at this position
//! would overlap with any other atom `B`, i.e., whether
//! the distance between `P` and `B` is less than `r_B+p`.
//! Probe positions which overlap with any atom are rejected.
//!
//! Overlapping atoms are fetched from a precomputed atom overlap grid.
//! The fetch gets all atoms whose centers are within `r_A + 2p + r_B`
//! of the A atom center.
//!
//! Atoms which have at least one probe position without overlaps
//! are considered to be on the solvent-accessible surface.
//! The fraction of remaining probe positions after rejection
//! estimates the per-atom solvent-accessible surface area.
//!
//! The computed surface probe positions and labeled surface atoms are saved
//! to binary format `.srf` files with the following layout:
//!
//! [5 bytes]  b"MOL> "
//! [4 bytes u32] length of the molecule name
//! [variable bytes] molecule name
//! [4 bytes u32] number of atoms
//! For each surface atom:
//!   [4 bytes f32] x
//!   [4 bytes f32] y
//!   [4 bytes f32] z
//!   [4 bytes f32] vdW radius r
//!   [4 bytes u8 ] four-character atom name (as ASCII)
//!   [3 bytes u8 ] three-letter resn (as ASCII)
//!   [1 byte  u8 ] chain_id (as ASCII)
//!   [4 bytes i32] resi
//!   [1 byte  u8 ] HET flag
//!   [1 byte  u8 ] surface flag
//! [3 bytes]  b"P> "
//! [4 bytes u32] number of probe lists - should match the number of atoms
//! For each probe list:
//!   [4 bytes u32] number of spheres in the probe list
//!   For each sphere:
//!     [4 bytes f32] x
//!     [4 bytes f32] y
//!     [4 bytes f32] z
//!     [4 bytes f32] r
//!

use crate::geo::ops::distance;
use crate::geo::point::Point;
use crate::geo::sphere::Sphere;
use crate::structure::atom::RawAtom;
use crate::structure::molecule::Molecule;

/// Number of points to sample on the Fibonacci sphere
const N_POINTS: usize = 500;

/// Maximum van der Waals radius to use for neighbor searching
const MAX_VDW: f64 = 2.0;

/// Minimum number of probes required to mark an atom as a surface atom
const MIN_PROBE_COUNT_FOR_SURFACE_MARK: usize = 1;

/// The depth of the surface layer for marking surface atoms
const SURFACE_DEPTH: f64 = 3.0;

/// Use the Fibonacci spiral to place points on the surface of a sphere
/// with a uniform distribution.
fn uniform_sphere_points(n: usize) -> Vec<Point> {
    use std::f64::consts::PI;
    let golden_angle = PI * (1.0 + 5.0_f64.sqrt());

    (0..n)
        .map(|i| {
            let theta = (1.0 - 2.0 * (i as f64 + 0.5) / n as f64).acos();
            let phi = golden_angle * i as f64;
            Point::new(
                theta.sin() * phi.cos(),
                theta.sin() * phi.sin(),
                theta.cos(),
            )
        })
        .collect()
}

/// Generate surface probe positions around a single atom
///
/// Candidate positions are sampled uniformly on the surface of a sphere
/// of radius `atom.r + probe_radius` centered at `atom.coor`.
fn generate_probes(atom: &RawAtom, probe_radius: f64) -> Vec<Sphere> {
    let unit_pts: Vec<Point> = uniform_sphere_points(N_POINTS);
    let shell_r: f64 = atom.r + probe_radius;
    // Scale the unit sphere to the size of the shell radius
    // and translate to the atom center coordinate
    let probes: Vec<Sphere> = unit_pts
        .into_iter()
        .map(|pt| Sphere {
            coord: pt * shell_r + atom.coor,
            r: probe_radius,
        })
        .collect();
    probes
}

/// Remove probes which are within the surface accessible sphere
/// of a given atom.
///
/// Measure the distance between the probe and the atom center.
/// If this distance is less than the atom radius plus probe radius,
/// the probe is considered to be within the accessible sphere and is removed.
fn remove_overlapping_probes(mut probes: Vec<Sphere>, atom: &RawAtom) -> Vec<Sphere> {
    probes.retain(|p| distance(p, atom) > atom.r + p.r);
    probes
}

/// Get the surface probes for a molecule.
///
/// Each atom has a vector of surface accessible probe positions.
pub fn get_surface_probes(mol: &Molecule, probe_radius: f64) -> Vec<Vec<Sphere>> {
    let grid = mol
        .grid
        .as_ref()
        .expect("molecule grid must be initialised before getting surface probes");

    let mut surface_probes = Vec::new();

    for atom in &mol.raw_atoms {
        let mut atom_probes = generate_probes(atom, probe_radius);
        let search_r = atom.r + probe_radius + MAX_VDW;
        let neighbors = grid.neighbors_within(&mol.raw_atoms, atom.coor, search_r);

        for &nbor_idx in &neighbors {
            let nbor_atom = &mol.raw_atoms[nbor_idx];
            atom_probes = remove_overlapping_probes(atom_probes, nbor_atom);
        }
        // Debug
        // How many probes are left accessible?
        // log::debug!("Atom {} ({}) has {} accessible probes", atom.name, atom.resn, atom_probes.len());

        surface_probes.push(atom_probes);
    }
    surface_probes
}

/// Mark surface exposed atoms by one of two criteria:
/// 1. The atom must have at least `MIN_PROBE_COUNT_FOR_SURFACE_MARK` accessible probes
/// 2. The atom is within PROBE_RADIUS + SURFACE_DEPTH of any probe belonging to a surface-exposed atom
///
/// The second check is only run once, after the first check has already marked surface atoms
pub fn mark_surface_atoms(mol: &mut Molecule, probes: &[Vec<Sphere>]) {
    // Mark surface exposed atoms by probe counts
    for (atom, atom_probes) in mol.raw_atoms.iter_mut().zip(probes.iter()) {
        if atom_probes.len() > MIN_PROBE_COUNT_FOR_SURFACE_MARK {
            atom.is_surface = true;
        }
    }

    // Mark surface exposed atoms by proximity to other surface atoms
    let grid = mol.grid.as_ref().unwrap();
    let mut neighbors_to_mark_exposed: Vec<usize> = Vec::new();
    for (atom, atom_probes) in mol.raw_atoms.iter().zip(probes.iter()) {
        if atom.is_surface {
            // Find atoms within PROBE_RADIUS + SURFACE_DEPTH of a probe point for this exposed atom
            for probe in atom_probes {
                for neighbor in
                    grid.neighbors_within(&mol.raw_atoms, probe.coord, probe.r + SURFACE_DEPTH)
                {
                    if !mol.raw_atoms[neighbor].is_surface {
                        neighbors_to_mark_exposed.push(neighbor);
                    }
                }
            }
        }
    }
    // Mark the neighbors found as surface atoms - no need to deduplicate
    for neighbor in neighbors_to_mark_exposed {
        mol.raw_atoms[neighbor].is_surface = true;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        geo::{
            point::{Point, Position},
            sphere::Sphere,
        },
        io::pdb::ParseOptions,
    };

    /// Call to mock an Atom at position (x,y,z) with radius r
    fn atom_at(x: f64, y: f64, z: f64, r: f64) -> RawAtom {
        RawAtom {
            coor: Point::new(x, y, z),
            r,
            name: "CA".to_string(),
            resn: "ALA".to_string(),
            chain_id: 'A',
            resi: 1,
            het: false,
            model: 1,
            is_surface: false,
        }
    }

    // Test uniform placement of candidate probe positions on a sphere
    // Each point from uniform_sphere_points(n) should be satisfy the unit
    // sphere equation x^2 + y^2 + z^2 = 1.
    #[test]
    fn sphere_points_lie_on_unit_sphere() {
        let pts: Vec<Point> = uniform_sphere_points(100);
        assert_eq!(pts.len(), 100, "should return exactly 100 sphere points");
        for pt in pts {
            let r = pt.position().x * pt.position().x
                + pt.position().y * pt.position().y
                + pt.position().z * pt.position().z;
            assert!(
                (r - 1.0).abs() < 1e-6,
                "point ({}, {}, {}) should lie on the unit sphere",
                pt.position().x,
                pt.position().y,
                pt.position().z
            );
        }
    }

    // Test that points generated around an atom are located at a distance of
    // atom_radius + probe_radius from the atom center
    #[test]
    fn probes_to_atom_distance_is_sum_of_radii() {
        let atom = atom_at(0.0, 0.0, 0.0, 1.0);
        let probe_radius = 0.5;
        let probes: Vec<Sphere> = generate_probes(&atom, probe_radius);
        for probe in probes {
            let d = distance(&probe.coord, &atom.coor);
            assert!(
                (d - (atom.r + probe.r)).abs() < 1e-6,
                "probe-to-atom distance should be sum of radii"
            );
        }
    }

    // Test that probes which overlap atoms are removed
    // Atom1 is at (0,0,0), with radius 1.
    // With a probe radius of 1, the candidate probe positions
    // will be placed at a maximum X-coordinate of 2, with the probe
    // radius reaching to a maximum X-coordinate of 3.
    // This means Atom2 at (3,0,0), with radius 1 should overlap with
    // at least one of the candidate probe positions.
    // The exact number of overlapping probe positions depends on the
    // logic used to generate uniform sphere points.
    #[test]
    fn probes_overlapping_atoms_are_removed() {
        let atom1 = atom_at(0.0, 0.0, 0.0, 1.0);
        let atom2 = atom_at(3.0, 0.0, 0.0, 1.0);
        let probe_radius = 1.0;
        let all_probes: Vec<Sphere> = generate_probes(&atom1, probe_radius);
        let full_set_probes = all_probes.clone();
        let pruned_probes: Vec<Sphere> = remove_overlapping_probes(all_probes, &atom2);
        assert!(
            pruned_probes.len() < full_set_probes.len(),
            "probes should be removed when overlapping atoms are present"
        );
    }

    fn tetra_molecule() -> Molecule {
        // Four atoms arranged in a tetrahedron shape
        let r: f64 = 1.7;
        let s: f64 = 2.7;
        let atoms = vec![
            atom_at(s, s, s, r),
            atom_at(s, -s, -s, r),
            atom_at(-s, s, -s, r),
            atom_at(-s, -s, s, r),
        ];
        Molecule {
            path: String::new(),
            opts: ParseOptions::default(),
            name: "tetra".to_string(),
            raw_atoms: atoms,
            atoms: Vec::new(),
            residues: Vec::new(),
            grid: None,
        }
    }

    #[test]
    fn mark_surface_atoms_works() {
        let mut mol = tetra_molecule();
        mol.build_grid(10.0);
        let probe_radius = 1.0;
        let probes: Vec<Vec<Sphere>> = get_surface_probes(&mol, probe_radius);

        mark_surface_atoms(&mut mol, &probes);
        mol.finalize_structure();
        let surface_count = mol.atoms.iter().filter(|a| a.is_surface).count();

        assert!(
            surface_count > 0,
            "should mark at least one atom as surface"
        );
    }

    #[test]
    fn finalize_structure_preserves_surface_flags() {
        use crate::geo::point::Point;
        use crate::structure::atom::RawAtom;
        use crate::structure::molecule::Molecule;

        let mut mol = Molecule::new("surface");
        mol.raw_atoms = vec![
            RawAtom {
                coor: Point::new(0.0, 0.0, 0.0),
                r: 1.7,
                name: "CA".into(),
                resn: "ALA".into(),
                chain_id: 'A',
                resi: 1,
                het: false,
                model: 1,
                is_surface: true,
            },
            RawAtom {
                coor: Point::new(3.0, 0.0, 0.0),
                r: 1.7,
                name: "CB".into(),
                resn: "ALA".into(),
                chain_id: 'A',
                resi: 1,
                het: false,
                model: 1,
                is_surface: false,
            },
        ];

        mol.finalize_structure();

        assert!(mol.atoms[0].is_surface);
        assert!(!mol.atoms[1].is_surface);
    }
}
