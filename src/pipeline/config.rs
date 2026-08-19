//! Config struct for the pipeline command parameters.

#[derive(Debug, Clone, Copy)]
pub struct Config {
    /// Probe radius for surface detection
    pub probe_radius: f64,

    /// Grid spacing for tetrahedron construction
    pub grid_spacing: f64,

    /// Tetrahedron construction distance cutoff
    pub thedron_cutoff: f64,

    /// Radius to search for compatible neighboring descriptor pairs
    /// when seeding a subgraph from a matched tetrahedron pair.
    pub graph_neighbor_search_radius: f64,

    /// Maximum distance between a transformed neighbor and its candidate
    /// match for the neighbor to be added to the subgraph.
    pub graph_neighbor_add_tolerance: f64,

    /// Maximum difference between corresponding inter-descriptor distances
    /// when adding extra edges between subgraph nodes.
    pub graph_extend_edges_tolerance: f64,

    /// Maximum average inter-descriptor distance when adding extra edges
    /// between subgraph nodes.
    pub graph_extend_edges_proximity: f64,

    /// Number of shared descriptor pairs required to merge two alignments
    /// into one cluster.
    pub cluster_descpair_threshold: i32,

    /// Distance tolerance for extending an alignment to new residue pairs.
    pub extend_dist: f64,

    /// BLOSUM score threshold for a candidate triple-residue while extending.
    pub extend_blosum_threshold: i32,

    /// Backbone vector angle tolerance while extending.
    pub extend_angle_threshold: f64,
}

/// Default configuration for the pipeline.
/// Most of these values are based on the default parameters of 
/// the ProBiS alignment algorithm.
impl Default for Config {
    fn default() -> Self {
        Self {
            probe_radius: 1.4, // Radius of a water molecule for surface detection
            grid_spacing: 5.0, 
            thedron_cutoff: 9.0,
            graph_neighbor_search_radius: 15.0,
            graph_neighbor_add_tolerance: 2.0,
            graph_extend_edges_tolerance: 2.0,
            graph_extend_edges_proximity: 15.0,
            cluster_descpair_threshold: 5,
            extend_dist: 10.0,
            extend_blosum_threshold: 23,
            // 45 deg in radians = pi/4
            extend_angle_threshold: std::f64::consts::PI / 4.0,
        }
    }
}
