//! ChainSet represents a set of chain identifiers (e.g. `A`, `B`, `C`),
//! with parsing utilities to convert from a comma-separated string.

use std::collections::BTreeSet;
use std::fmt;
use std::str::FromStr;

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct ChainSet(BTreeSet<char>);

// Define error types for ChainSet parsing
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ChainSetParseError {
    Empty,
    EmptyToken(String),
    MultiCharToken(String),
}

impl fmt::Display for ChainSetParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ChainSetParseError::Empty => write!(f, "chain selection cannot be empty"),
            ChainSetParseError::EmptyToken(s) => {
                write!(f, "chain selection {s:?} contains an empty chain token")
            }
            ChainSetParseError::MultiCharToken(s) => {
                write!(f, "chain token {s:?} is not a single character")
            }
        }
    }
}

impl std::error::Error for ChainSetParseError {}

impl ChainSet {
    pub fn single(c: char) -> Self {
        ChainSet(BTreeSet::from([c]))
    }

    pub fn parse(s: &str) -> Result<Self, ChainSetParseError> {
        s.parse()
    }

    pub fn contains(&self, c: char) -> bool {
        self.0.contains(&c)
    }

    pub fn iter(&self) -> impl Iterator<Item = &char> + '_ {
        self.0.iter()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// If the set contains a single character, return it; otherwise, return `None`.
    pub fn as_single(&self) -> Option<char> {
        if self.0.len() == 1 {
            self.0.iter().next().copied()
        } else {
            None
        }
    }
}

impl FromStr for ChainSet {
    type Err = ChainSetParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.trim().is_empty() {
            return Err(ChainSetParseError::Empty);
        }
        let mut set = BTreeSet::new();
        for token in s.split(',') {
            let token = token.trim();
            if token.is_empty() {
                // An empty token string after trimming, e.g., ",  ,"
                return Err(ChainSetParseError::EmptyToken(s.to_string()));
            }
            let mut chars = token.chars();

            if chars.clone().count() > 1 {
                return Err(ChainSetParseError::MultiCharToken(token.to_string()));
            } else {
                set.insert(chars.next().unwrap());
            }
        }
        Ok(ChainSet(set))
    }
}

impl fmt::Display for ChainSet {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let parts: Vec<String> = self.0.iter().map(|c| c.to_string()).collect();
        write!(f, "{}", parts.join(","))
    }
}

impl FromIterator<char> for ChainSet {
    fn from_iter<I: IntoIterator<Item = char>>(iter: I) -> Self {
        ChainSet(iter.into_iter().collect())
    }
}

impl From<char> for ChainSet {
    fn from(c: char) -> Self {
        ChainSet::single(c)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_single_chain() {
        let set: ChainSet = "A".parse().unwrap();
        assert_eq!(set.as_single(), Some('A'));
    }

    #[test]
    fn parses_multiple_chains_any_order() {
        let set1: ChainSet = "A,B".parse().unwrap();
        let set2: ChainSet = "B,A".parse().unwrap();
        assert_eq!(set1, set2);
        assert_eq!(set1.to_string(), "A,B");
    }

    #[test]
    fn rejects_empty_string() {
        assert_eq!("".parse::<ChainSet>(), Err(ChainSetParseError::Empty));
    }

    #[test]
    fn rejects_empty_token() {
        assert!(matches!(
            "A,,B".parse::<ChainSet>(),
            Err(ChainSetParseError::EmptyToken(_))
        ));
    }

    #[test]
    fn rejects_multichar_token() {
        assert!(matches!(
            "AB".parse::<ChainSet>(),
            Err(ChainSetParseError::MultiCharToken(_))
        ));
    }

    #[test]
    fn contains_checks_membership() {
        let set: ChainSet = "A,B".parse().unwrap();
        assert!(set.contains('A'));
        assert!(!set.contains('C'));
    }
}
