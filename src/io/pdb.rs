//! PDB file parsing utilities.

use crate::geo::point::Point;
use crate::structure::atom::RawAtom;
use crate::structure::chainset::ChainSet;
use anyhow::{Context, Result};

/// Van der Waals radii for each element type
///
/// Sourced from:
/// A. Bondi; van der Waals Volumes and Radii. J. Phys. Chem. 1 March 1964; 68 (3): 441–451. https://doi.org/10.1021/j100785a001
fn vdw_radius(element: &str) -> f64 {
    match element.to_uppercase().as_str() {
        "C" => 1.70,
        "N" => 1.55,
        "O" => 1.52,
        "S" => 1.80,
        "H" => 1.20,
        "P" => 1.80,
        "F" => 1.47,
        "CL" => 1.75,
        "BR" => 1.85,
        "I" => 1.98,
        _ => 1.70,
    }
}

fn element_from_line(line: &str) -> &str {
    // Columns 77-78 hold the element symbol, if present
    if line.len() >= 78 {
        let e = line[76..78].trim(); // 0-indexed
        if !e.is_empty() {
            return e;
        }
    }
    // Fallback to the first character of the atom name
    // Last resort is an empty string
    if line.len() >= 16 {
        &line[12..16].trim()[..1]
    } else {
        ""
    }
}

/// Parsing options control which records are included
#[derive(Debug, Clone)]
pub struct ParseOptions {
    /// Include HETATM records (default: false)
    pub include_hetero: bool,
    /// Chain IDs to include (default: all chains)
    /// If Some("AB"), only chains 'A' and 'B' will be included
    /// If None, all chains will be included
    pub chain_ids: Option<ChainSet>,
    /// Only include the first model (default: true)
    /// Will stop reading after the first ENDMDL record.
    pub first_model_only: bool,
}

impl Default for ParseOptions {
    fn default() -> Self {
        ParseOptions {
            include_hetero: false,
            chain_ids: None,
            first_model_only: true,
        }
    }
}

/// Parse ATOM/HETATM records from a PDB file on disk
pub fn parse_pdb_file(path: &str, opts: &ParseOptions) -> Result<Vec<RawAtom>> {
    // Check the file extension
    if !path.ends_with(".pdb") {
        return Err(anyhow::anyhow!("file must have .pdb extension"));
    }

    let content = std::fs::read_to_string(path)
        .with_context(|| format!("failed to read PDB file: {path}"))?;
    parse_pdb_str(&content, opts)
}

/// Parse ATOM/HETATM records from an in-memory string
pub fn parse_pdb_str(content: &str, opts: &ParseOptions) -> Result<Vec<RawAtom>> {
    let mut atoms: Vec<RawAtom> = Vec::new();
    let mut model_num: u32 = 1;

    for line in content.lines() {
        // Get the record type of this line
        // Guard the type check behind a length check
        if line.len() < 6 {
            continue;
        }

        let record = line[..6].trim();
        let mut is_atom: bool = false;
        let mut is_hetatm: bool = false;

        match record {
            "MODEL" => {
                model_num = line[6..].trim().parse().unwrap_or(model_num + 1);
                continue;
            }
            "ENDMDL" => {
                if opts.first_model_only {
                    break;
                }
                continue;
            }
            "ATOM" => {
                is_atom = true;
            }
            "HETATM" => {
                is_hetatm = true;
                if !opts.include_hetero {
                    continue;
                }
            }
            _ => {}
        }

        if !is_atom && !is_hetatm {
            continue; // If the line is neither ATOM nor HETATM, skip
        } else if line.len() < 54 {
            continue; // If the line is too short to contain coords, skip
        }

        // Parse the fixed-column fields
        let name = line.get(12..16).unwrap_or("    ").trim().to_string();
        let altloc = line.get(16..17).unwrap_or(" ").to_string();
        let resn = line.get(17..20).unwrap_or("   ").trim().to_string();
        let chain_id = line.chars().nth(21).unwrap_or(' ');
        let resi: i32 = line
            .get(22..26)
            .and_then(|s| s.trim().parse().ok())
            .unwrap_or(0);
        let x: f64 = line
            .get(30..38)
            .and_then(|s| s.trim().parse().ok())
            .unwrap_or(0.0);
        let y: f64 = line
            .get(38..46)
            .and_then(|s| s.trim().parse().ok())
            .unwrap_or(0.0);
        let z: f64 = line
            .get(46..54)
            .and_then(|s| s.trim().parse().ok())
            .unwrap_or(0.0);

        // Skip hydrogens ~= atoms whose name begins with "H"
        if name.starts_with("H") {
            continue;
        }

        // Check if this atom is a secondary alternate location ("B", etc)
        // If so, skip it (we only use the primary alternate location, or blank)
        if altloc != " " && altloc != "A" {
            continue;
        }

        // Check if this chain is selected for parsing
        if let Some(ref chains_to_parse) = opts.chain_ids
            && !chains_to_parse.contains(chain_id)
        {
            continue;
        }

        let r = vdw_radius(element_from_line(line));

        atoms.push(RawAtom {
            coor: Point::new(x, y, z),
            r,
            name,
            resn,
            chain_id,
            resi,
            het: is_hetatm,
            model: model_num,
            is_surface: false,
        });
    }
    Ok(atoms)
}

#[cfg(test)]
mod tests {
    use crate::geo::point::Position;

    use super::*;

    // A sample PDB file chunk
    const SAMPLE: &str = "\
ATOM      1  N   ASP H  18     159.469 114.021 101.822  1.00 45.61           N  \n\
ATOM      2  CA  ASP H  18     158.201 113.284 101.580  1.00 45.61           C  \n\
ATOM      3  C   ASP H  18     157.800 112.535 102.852  1.00 45.61           C  \n\
ATOM      4  O   ASP H  18     158.616 112.497 103.787  1.00 45.61           O  \n\
ATOM      5  CB  ASP H  18     158.315 112.337 100.384  1.00 45.61           C  \n\
ATOM      6  CG  ASP H  18     159.296 111.194 100.589  1.00 45.61           C  \n\
ATOM      7  OD1 ASP H  18     160.309 111.402 101.286  1.00 45.61           O  \n\
ATOM      8  OD2 ASP H  18     159.035 110.101 100.050  1.00 45.61           O  \n\
ATOM      9  N   VAL H  19     156.600 111.954 102.875  1.00 36.05           N  \n\
ATOM     10  CA  VAL H  19     156.098 111.255 104.095  1.00 36.05           C  \n\
ATOM     11  C   VAL H  19     156.627 109.830 104.048  1.00 36.05           C  \n\
ATOM     12  O   VAL H  19     156.447 109.189 103.008  1.00 36.05           O  \n\
ATOM     13  CB  VAL H  19     154.563 111.297 104.213  1.00 36.05           C  \n\
ATOM     14  CG1 VAL H  19     154.058 110.438 105.356  1.00 36.05           C  \n\
ATOM     15  CG2 VAL H  19     154.049 112.719 104.350  1.00 36.05           C  \n\
ATOM     16  N   GLN H  20     157.287 109.364 105.107  1.00 37.49           N  \n\
ATOM     17  CA  GLN H  20     157.842 107.987 105.123  1.00 37.49           C  \n\
ATOM     18  C   GLN H  20     157.681 107.377 106.512  1.00 37.49           C  \n\
ATOM     19  O   GLN H  20     157.828 108.118 107.495  1.00 37.49           O  \n\
ATOM     20  CB  GLN H  20     159.314 107.981 104.720  1.00 37.49           C  \n\
ATOM     21  CG  GLN H  20     159.552 108.108 103.224  1.00 37.49           C  \n\
ATOM     22  CD  GLN H  20     161.026 108.114 102.909  1.00 37.49           C  \n\
ATOM     23  OE1 GLN H  20     161.867 107.982 103.795  1.00 37.49           O  \n\
ATOM     24  NE2 GLN H  20     161.350 108.267 101.636  1.00 37.49           N  \n\
HETATM 9089  C1  CLR R 401     133.505 127.415 142.016  1.00 59.74           C  \n\
HETATM 9090  C2  CLR R 401     133.682 128.373 140.835  1.00 59.74           C  \n\
HETATM 9091  C3  CLR R 401     132.435 129.182 140.578  1.00 59.74           C  \n\
HETATM 9092  C4  CLR R 401     132.055 129.941 141.834  1.00 59.74           C  \n\
HETATM 9093  C5  CLR R 401     131.870 129.004 143.007  1.00 59.74           C  \n\
HETATM 9094  C6  CLR R 401     130.729 128.968 143.679  1.00 59.74           C  \n\
HETATM 9095  C7  CLR R 401     130.505 128.228 144.959  1.00 59.74           C  \n\
HETATM 9096  C8  CLR R 401     131.769 127.564 145.499  1.00 59.74           C  \n\
HETATM 9097  C9  CLR R 401     132.614 127.004 144.339  1.00 59.74           C  \n\
HETATM 9098  C10 CLR R 401     133.053 128.099 143.328  1.00 59.74           C  \n\
HETATM 9099  C11 CLR R 401     133.799 126.177 144.858  1.00 59.74           C  \n\
HETATM 9100  C12 CLR R 401     133.379 125.077 145.839  1.00 59.74           C  \n\
HETATM 9101  C13 CLR R 401     132.592 125.640 147.026  1.00 59.74           C  \n\
HETATM 9102  C14 CLR R 401     131.401 126.423 146.437  1.00 59.74           C  \n\
HETATM 9103  C15 CLR R 401     130.512 126.711 147.643  1.00 59.74           C  \n\
HETATM 9104  C16 CLR R 401     130.670 125.455 148.522  1.00 59.74           C  \n\
HETATM 9105  C17 CLR R 401     131.819 124.613 147.906  1.00 59.74           C  \n\
HETATM 9106  C18 CLR R 401     133.484 126.527 147.909  1.00 59.74           C  \n\
HETATM 9107  C19 CLR R 401     134.196 128.970 143.881  1.00 59.74           C  \n\
HETATM 9108  C20 CLR R 401     132.555 123.771 148.965  1.00 59.74           C  \n\
HETATM 9109  C21 CLR R 401     133.514 122.745 148.368  1.00 59.74           C  \n\
HETATM 9110  C22 CLR R 401     131.537 123.068 149.872  1.00 59.74           C  \n\
HETATM 9111  C23 CLR R 401     132.131 122.074 150.834  1.00 59.74           C  \n\
HETATM 9112  C24 CLR R 401     131.105 121.186 151.485  1.00 59.74           C  \n\
HETATM 9113  C25 CLR R 401     131.677 120.011 152.274  1.00 59.74           C  \n\
HETATM 9114  C26 CLR R 401     132.774 119.311 151.488  1.00 59.74           C  \n\
HETATM 9115  C27 CLR R 401     130.588 119.020 152.656  1.00 59.74           C  \n\
HETATM 9116  O1  CLR R 401     132.654 130.098 139.496  1.00 59.74           O  \n";

    const TWO_CHAINS: &str = "\
    ATOM   8977  N   HIS R 318     121.052 146.327 147.256  1.00 82.23           N  \n\
    ATOM   8978  CA  HIS R 318     120.475 145.797 145.991  1.00 82.23           C  \n\
    ATOM   8979  C   HIS R 318     119.685 144.510 146.256  1.00 82.23           C  \n\
    ATOM   8980  O   HIS R 318     118.624 144.347 145.626  1.00 82.23           O  \n\
    ATOM   8981  CB  HIS R 318     121.572 145.630 144.931  1.00 82.23           C  \n\
    ATOM   8982  CG  HIS R 318     121.114 144.980 143.670  1.00 82.23           C  \n\
    ATOM   8983  ND1 HIS R 318     121.407 143.664 143.371  1.00 82.23           N  \n\
    ATOM   8984  CD2 HIS R 318     120.392 145.453 142.632  1.00 82.23           C  \n\
    ATOM   8985  CE1 HIS R 318     120.884 143.354 142.203  1.00 82.23           C  \n\
    ATOM   8986  NE2 HIS R 318     120.256 144.435 141.729  1.00 82.23           N  \n\
    ATOM   8987  N   ALA R 319     120.168 143.629 147.141  1.00 79.57           N  \n\
    ATOM   8988  CA  ALA R 319     119.357 142.437 147.463  1.00 79.57           C  \n\
    ATOM   8989  C   ALA R 319     118.282 142.850 148.470  1.00 79.57           C  \n\
    ATOM   8990  O   ALA R 319     118.653 143.438 149.502  1.00 79.57           O  \n\
    ATOM   8991  CB  ALA R 319     120.231 141.355 148.028  1.00 79.57           C  \n\
    TER    8992      ALA R 319                                                      \n\
    ATOM   8993  N   ARG P   1     134.005 108.227 180.952  1.00120.35           N  \n\
    ATOM   8994  CA  ARG P   1     133.406 108.053 182.303  1.00120.35           C  \n\
    ATOM   8995  C   ARG P   1     132.314 109.110 182.499  1.00120.35           C  \n\
    ATOM   8996  O   ARG P   1     132.600 110.148 183.129  1.00120.35           O  \n\
    ATOM   8997  CB  ARG P   1     134.495 108.156 183.375  1.00120.35           C  \n\
    ATOM   8998  CG  ARG P   1     134.109 107.561 184.721  1.00120.35           C  \n\
    ATOM   8999  CD  ARG P   1     135.252 107.649 185.713  1.00120.35           C  \n\
    ATOM   9000  NE  ARG P   1     134.938 107.005 186.980  1.00120.35           N  \n\
    ATOM   9001  CZ  ARG P   1     135.758 106.952 188.024  1.00120.35           C  \n\
    ATOM   9002  NH1 ARG P   1     136.955 107.509 187.957  1.00120.35           N  \n\
    ATOM   9003  NH2 ARG P   1     135.377 106.342 189.132  1.00120.35           N  \n\
    ATOM   9004  N   PRO P   2     131.081 108.899 181.986  1.00115.18           N  \n\
    ATOM   9005  CA  PRO P   2     130.023 109.908 182.080  1.00115.18           C  \n\
    ATOM   9006  C   PRO P   2     129.620 110.166 183.537  1.00115.18           C  \n\
    ATOM   9007  O   PRO P   2     129.178 109.240 184.194  1.00115.18           O  \n\
    ATOM   9008  CB  PRO P   2     128.849 109.281 181.315  1.00115.18           C  \n\
    ATOM   9009  CG  PRO P   2     129.485 108.212 180.456  1.00115.18           C  \n\
    ATOM   9010  CD  PRO P   2     130.641 107.694 181.284  1.00115.18           C  \n";

    const TWO_MODELS: &str = "\
MODEL        1\n\
ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 10.00           C  \n\
ENDMDL\n\
MODEL        2\n\
ATOM      2  CA  ALA A   1       4.000   5.000   6.000  1.00 10.00           C  \n\
ENDMDL\n";

    #[test]
    fn test_count_atoms_no_hetero() {
        let atoms = parse_pdb_str(SAMPLE, &ParseOptions::default()).unwrap();
        assert_eq!(atoms.len(), 24);
    }

    #[test]
    fn test_count_atoms_with_hetero() {
        let opts = ParseOptions {
            include_hetero: true,
            ..ParseOptions::default()
        };
        let atoms = parse_pdb_str(SAMPLE, &opts).unwrap();
        assert_eq!(atoms.len(), 52);
    }

    #[test]
    fn test_hetatm_flag_set() {
        let opts = ParseOptions {
            include_hetero: true,
            ..ParseOptions::default()
        };
        let atoms = parse_pdb_str(SAMPLE, &opts).unwrap();
        let hets: Vec<_> = atoms.iter().filter(|a| a.het).collect();
        assert_eq!(hets.len(), 28);
        assert_eq!(hets[0].resn, "CLR");
    }

    #[test]
    fn test_coordinates_parsed() {
        let atoms = parse_pdb_str(SAMPLE, &ParseOptions::default()).unwrap();
        let n = atoms
            .iter()
            .find(|a| a.name == "N" && a.resi == 18)
            .unwrap();
        assert!((n.position().x - 159.469).abs() < 0.001);
        assert!((n.position().y - 114.021).abs() < 0.001);
        assert!((n.position().z - 101.822).abs() < 0.001);
    }

    #[test]
    fn test_residue_fields() {
        let atoms = parse_pdb_str(SAMPLE, &ParseOptions::default()).unwrap();
        let ca = atoms.iter().find(|a| a.name == "CA").unwrap();
        assert_eq!(ca.resn, "ASP");
        assert_eq!(ca.chain_id, 'H');
        assert_eq!(ca.resi, 18);
        assert_eq!(ca.het, false);
    }

    #[test]
    fn test_chain_filter() {
        let opts = ParseOptions {
            chain_ids: Some(ChainSet::single('R')),
            ..ParseOptions::default()
        };
        let atoms = parse_pdb_str(TWO_CHAINS, &opts).unwrap();
        assert_eq!(atoms.len(), 15);
        assert!(atoms.iter().all(|a| a.chain_id == 'R'))
    }

    #[test]
    fn test_first_model_only() {
        let atoms = parse_pdb_str(TWO_MODELS, &ParseOptions::default()).unwrap();
        assert_eq!(atoms.len(), 1);
        assert!((atoms[0].position().x - 1.0).abs() < 0.001);
        assert!((atoms[0].position().y - 2.0).abs() < 0.001);
        assert!((atoms[0].position().z - 3.0).abs() < 0.001);
    }

    #[test]
    fn test_all_models() {
        let opts = ParseOptions {
            first_model_only: false,
            ..ParseOptions::default()
        };
        let atoms = parse_pdb_str(TWO_MODELS, &opts).unwrap();
        assert_eq!(atoms.len(), 2);
        assert!(atoms[0].model == 1);
        assert!(atoms[1].model == 2);
        assert!((atoms[1].position().x - 4.0).abs() < 0.001);
        assert!((atoms[1].position().y - 5.0).abs() < 0.001);
        assert!((atoms[1].position().z - 6.0).abs() < 0.001);
    }
}
