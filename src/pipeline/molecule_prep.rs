//! Grouping together the repeated stages for loading
//! and preparing a molecule from a PDB file.
//!

use anyhow::{Context, Result};
use std::path::Path;

use crate::{
    geo::sphere::Sphere,
    io::pdb::ParseOptions,
    pipeline::config::Config,
    structure::{
        chainset::ChainSet,
        molecule::Molecule,
        site::BindingSiteSpec,
        surface::{get_surface_probes, mark_surface_atoms},
    },
};

pub struct PreparedMolecule {
    pub name: String,
    pub molecule: Molecule,
    pub probes: Vec<Vec<Sphere>>,
}

pub fn load_and_prepare_molecule(
    pdb_path: &str,
    chain: &ChainSet,
    site: Option<BindingSiteSpec>,
    config: &Config,
) -> Result<PreparedMolecule> {
    let mol_name = Path::new(&pdb_path)
        .file_stem()
        .unwrap_or_else(|| panic!("Could not get file name from {}", &pdb_path))
        .to_str()
        .unwrap_or_else(|| panic!("Could not parse file name from {}", &pdb_path));

    // Load molecules
    let opts = ParseOptions {
        chain_ids: Some(chain.clone()),
        ..Default::default()
    };
    log::info!("Loading {}...", mol_name);
    let mut mol = Molecule::new(mol_name);
    mol.read_pdb(pdb_path, opts)
        .with_context(|| "Failed to load molecule 1")?;

    // Surface extraction
    log::info!("Extracting surface probes");
    mol.build_grid(config.grid_spacing);

    let probes = get_surface_probes(&mol, config.probe_radius);

    mark_surface_atoms(&mut mol, &probes);

    // If the user specified a binding site, hide atoms beyond the binding site
    if let Some(site) = site {
        mol.hide_atoms_beyond_binding_site(
            site.residue_name,
            site.residue_number,
            site.chain_id,
            site.distance_angstroms,
        )?;
    }

    Ok(PreparedMolecule {
        name: mol_name.to_string(),
        molecule: mol,
        probes,
    })
}
