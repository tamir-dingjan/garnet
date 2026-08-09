//! This module owns the grid data structure used for quick neighbour searching.

use crate::geo::point::{Point, Position};
use std::collections::HashMap;

/// The grid of cells is used for quick adjacency lookups
/// The cells are hashed by their grid coordinates,
/// and each contains a Vec of indices into the collection of
/// <T> used to populate the grid.
pub struct Grid<T: Position> {
    pub spacing: f64,
    cells: HashMap<(i32, i32, i32), Vec<usize>>,
    _marker: std::marker::PhantomData<T>, // Reference to <T> to apply the same drop rules and lifetime
}

impl<T: Position> Grid<T> {
    /// Creates a new Grid with the given cell spacing in angstroms.
    pub fn new(spacing: f64) -> Self {
        Grid {
            spacing,
            cells: HashMap::new(),
            _marker: std::marker::PhantomData,
        }
    }

    /// Builds the grid from a slice of items
    /// Items are hashed into grid cells based on the cell position
    /// and stored by their index in the provided &[T].
    /// This means that the grid does not own the items,
    /// and the items must remain valid for scope where the grid methods are called.
    pub fn build(&mut self, items: &[T]) {
        self.cells.clear();
        for (idx, item) in items.iter().enumerate() {
            let key = self.cell_key(item.x(), item.y(), item.z());
            self.cells.entry(key).or_default().push(idx);
        }
    }

    /// Returns the cell key for a given position.
    /// The cell key is a tuple of the grid cell indices: (i32, i32, i32).
    fn cell_key(&self, x: f64, y: f64, z: f64) -> (i32, i32, i32) {
        (
            (x / self.spacing).floor() as i32,
            (y / self.spacing).floor() as i32,
            (z / self.spacing).floor() as i32,
        )
    }

    /// Return the indices of all items within a search radius of
    /// the given position (qx, qy, qz)
    /// This method uses a grid cell search to find neighbors efficiently.
    /// It first gets all the grid cells which cover the search radius,
    /// and then iterates over them to find neighbors within the search radius.
    pub fn neighbors_within(&self, items: &[T], position: Point, radius: f64) -> Vec<usize> {
        let mut neighbors = Vec::new();

        // How many cells are included in the search radius?
        // Add 1 cell to cover items on the edge of a cell.
        let cell_radius = (radius / self.spacing).ceil() as i32 + 1;

        // Get the origin cell position for the search query
        let (oi, oj, ok) = self.cell_key(position.x(), position.y(), position.z());

        // iterate over cells within the search radius of the origin cell
        for i in -cell_radius..=cell_radius {
            for j in -cell_radius..=cell_radius {
                for k in -cell_radius..=cell_radius {
                    let search_cell = (oi + i, oj + j, ok + k);
                    // test each of the indices in the search cell for distance from the query
                    if let Some(indices) = self.cells.get(&search_cell) {
                        for &idx in indices {
                            let item = &items[idx];
                            // compare the euclidean distance to the search radius
                            let dx = item.x() - position.x();
                            let dy = item.y() - position.y();
                            let dz = item.z() - position.z();
                            let distance = (dx * dx + dy * dy + dz * dz).sqrt();
                            if distance <= radius {
                                // the index is within the search radius
                                neighbors.push(idx);
                            }
                        }
                    }
                }
            }
        }
        neighbors
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geo::point::Point;

    fn make_points(coords: &[(f64, f64, f64)]) -> Vec<Point> {
        coords
            .iter()
            .map(|&(x, y, z)| Point::new(x, y, z))
            .collect()
    }

    // Test that a point at (1,1,1) is found when querying (0,0,0) with radius 2.0
    #[test]
    fn grid_finds_nearby_point() {
        let pts = make_points(&[(1.0, 1.0, 1.0), (10.0, 10.0, 10.0)]);
        let mut grid = Grid::new(3.0);
        grid.build(&pts);
        let search_position = Point::new(0.0, 0.0, 0.0);
        let hits = grid.neighbors_within(&pts, search_position, 2.0);
        assert!(
            hits.contains(&0),
            "should find index 0 for point near the origin"
        );
        assert!(
            !hits.contains(&1),
            "should NOT find index 1 for point far from the origin"
        );
    }

    // Test that a query far from all points is empty
    #[test]
    fn grid_empty_when_querying_far_from_all_points() {
        let pts = make_points(&[(10.0, 10.0, 10.0)]);
        let mut grid = Grid::new(3.0);
        grid.build(&pts);
        let search_position = Point::new(0.0, 0.0, 0.0);
        let hits = grid.neighbors_within(&pts, search_position, 2.0);
        assert!(
            hits.is_empty(),
            "should be empty when querying far from all points"
        );
    }

    // Test that searching around a point finds itself
    #[test]
    fn grid_finds_point_when_searching_around_itself() {
        let pts = make_points(&[(1.0, 2.0, 3.0)]);
        let mut grid = Grid::new(3.0);
        grid.build(&pts);
        let search_position = Point::new(1.0, 2.0, 3.0);
        let hits = grid.neighbors_within(&pts, search_position, 3.0);
        assert!(
            hits.contains(&0),
            "should find index 0 for point when searching around itself"
        );
    }

    // Test that a point exactly at the search boundary is found
    // Note that the search boundary is not the same as the grid cell spacing
    // Test point at (2,0,0), query radius 2.0 from (0,0,0) in grid with spacing 3.0
    #[test]
    fn grid_finds_point_at_search_boundary() {
        let pts = make_points(&[(2.0, 0.0, 0.0)]);
        let mut grid = Grid::new(3.0);
        grid.build(&pts);
        let search_position = Point::new(0.0, 0.0, 0.0);
        let hits = grid.neighbors_within(&pts, search_position, 2.0);
        assert!(
            hits.contains(&0),
            "should find index 0 for point exactly at the search radius boundary"
        );
    }

    // Test that a point exactly at the grid boundary is found
    // Test point at (2,0,0), query radius 2.0 from (0,0,0) in grid with spacing 2.0
    #[test]
    fn grid_find_point_on_grid_boundary() {
        let pts = make_points(&[(2.0, 0.0, 0.0)]);
        let mut grid = Grid::new(2.0);
        grid.build(&pts);
        let search_position = Point::new(0.0, 0.0, 0.0);
        let hits = grid.neighbors_within(&pts, search_position, 2.0);
        assert!(
            hits.contains(&0),
            "should find index 0 for point exactly at the grid boundary"
        );
    }

    // Test that rebuilding clears the old contents
    #[test]
    fn grid_rebuild_clears_old_contents() {
        // Set up grid with old points
        let pts_old = make_points(&[(1.0, 1.0, 1.0)]);
        let mut grid = Grid::new(3.0);
        grid.build(&pts_old);

        // Rebuild with new points
        let pts_new = make_points(&[(20.0, 20.0, 20.0)]);
        grid.build(&pts_new);

        // Querying near the old point location should find nothing
        let search_position = Point::new(0.0, 0.0, 0.0);
        let hits = grid.neighbors_within(&pts_new, search_position, 3.0);
        assert!(hits.is_empty(), "old contents should be gone after rebuild");
    }
}
