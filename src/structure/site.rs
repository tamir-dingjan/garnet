use std::fmt;
use std::str::FromStr;

use anyhow::{Context, Result, bail};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BindingSiteSpec {
    pub residue_name: String,
    pub residue_number: i32,
    pub chain_id: char,
    pub distance_angstroms: f64,
}

impl FromStr for BindingSiteSpec {
    type Err = anyhow::Error;

    fn from_str(input: &str) -> Result<Self> {
        let fields: Vec<&str> = input.split(':').collect();
        if fields.len() != 4 {
            bail!(
                "binding site must have exactly four fields (RESN:RESI:CHAIN:DISTANCE), got {}",
                fields.len()
            );
        }

        let residue_name = fields[0].trim();
        if residue_name.is_empty() {
            bail!("binding site residue name cannot be empty");
        }

        let residue_number = fields[1]
            .parse::<i32>()
            .with_context(|| format!("invalid binding site residue number: {:?}", fields[1]))?;

        let mut chain_chars = fields[2].chars();
        let chain_id = chain_chars
            .next()
            .ok_or_else(|| anyhow::anyhow!("binding site chain ID cannot be empty"))?;
        if chain_chars.next().is_some() {
            bail!("binding site chain ID must be exactly one character");
        }

        let distance_angstroms = fields[3]
            .parse::<f64>()
            .with_context(|| format!("invalid binding site distance: {:?}", fields[3]))?;
        if !distance_angstroms.is_finite() {
            bail!("binding site distance must be finite");
        }
        if distance_angstroms <= 0.0 {
            bail!("binding site distance must be positive");
        }

        Ok(Self {
            residue_name: residue_name.to_string(),
            residue_number,
            chain_id,
            distance_angstroms,
        })
    }
}

impl fmt::Display for BindingSiteSpec {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "{}:{}:{}:{}",
            self.residue_name, self.residue_number, self.chain_id, self.distance_angstroms
        )
    }
}

#[cfg(test)]
mod tests {
    use super::BindingSiteSpec;

    #[test]
    fn parses_and_formats_binding_site() {
        let site: BindingSiteSpec = "ATP:501:A:6.0".parse().expect("valid binding site");

        assert_eq!(
            site,
            BindingSiteSpec {
                residue_name: "ATP".to_string(),
                residue_number: 501,
                chain_id: 'A',
                distance_angstroms: 6.0,
            }
        );
        assert_eq!(site.to_string(), "ATP:501:A:6");
    }

    #[test]
    fn trims_surrounding_residue_name_whitespace() {
        let site: BindingSiteSpec = "  ATP  :501:A:6"
            .parse()
            .expect("surrounding residue-name whitespace should be accepted");

        assert_eq!(site.residue_name, "ATP");
        assert_eq!(site.to_string(), "ATP:501:A:6");
    }

    #[test]
    fn rejects_invalid_binding_sites() {
        for input in [
            "ATP:501:A",
            ":501:A:6",
            "   :501:A:6",
            "ATP:x:A:6",
            "ATP:501::6",
            "ATP:501:AB:6",
            "ATP:501:A:bad",
            "ATP:501:A:0",
            "ATP:501:A:-1",
            "ATP:501:A:NaN",
            "ATP:501:A:inf",
        ] {
            assert!(
                input.parse::<BindingSiteSpec>().is_err(),
                "{input:?} should be rejected"
            );
        }
    }
}
