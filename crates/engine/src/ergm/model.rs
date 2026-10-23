//! A model: a network state, its terms, and the dyads that may change.

use crate::error::{Error, Result};
use crate::matrix::Dense;

use super::state::{Constraints, NetState};
use super::terms::Terms;

/// An exponential random graph model attached to an observed network.
#[derive(Debug, Clone)]
pub struct Model {
    state: NetState,
    terms: Terms,
    constraints: Option<Constraints>,
}

impl Model {
    /// Attaches `terms` to the observed network `state`. `constraints`, if
    /// given, must have been built for the same number of nodes and
    /// direction.
    pub fn new(state: NetState, terms: Terms, constraints: Option<Constraints>) -> Result<Model> {
        if let Some(c) = &constraints {
            if c.node_count() != state.node_count() || c.is_directed() != state.is_directed() {
                return Err(Error::InvalidArgument(
                    "fixed dyads were built for a different network".into(),
                ));
            }
        }
        Ok(Model {
            state,
            terms,
            constraints,
        })
    }

    /// The observed network.
    #[inline]
    pub fn state(&self) -> &NetState {
        &self.state
    }

    /// The terms.
    #[inline]
    pub fn terms(&self) -> &Terms {
        &self.terms
    }

    /// The fixed-dyad constraints, if any.
    #[inline]
    pub fn constraints(&self) -> Option<&Constraints> {
        self.constraints.as_ref()
    }

    /// Number of statistics.
    #[inline]
    pub fn width(&self) -> usize {
        self.terms.width()
    }

    /// Whether dyad (i, j) may be toggled.
    #[inline]
    pub fn is_free(&self, i: usize, j: usize) -> bool {
        match &self.constraints {
            Some(c) => !c.is_fixed(i, j),
            None => true,
        }
    }

    /// The free dyads in canonical order.
    pub fn free_dyads(&self) -> Vec<(u32, u32)> {
        match &self.constraints {
            Some(c) => c.free_dyads().to_vec(),
            None => self.state.dyads(),
        }
    }

    /// Number of free dyads.
    pub fn free_dyad_count(&self) -> usize {
        match &self.constraints {
            Some(c) => c.free_dyads().len(),
            None => self.state.dyad_count(),
        }
    }

    /// The observed statistics.
    pub fn stats(&self) -> Vec<f64> {
        self.terms.stats(&self.state)
    }

    /// The change statistics of every free dyad of the observed network
    /// (one row per dyad, in canonical order), the dyads' observed states,
    /// and the dyads themselves. This is the design of the
    /// pseudolikelihood and of conditional tie probabilities.
    pub fn change_matrix(&self) -> (Dense, Vec<f64>, Vec<(u32, u32)>) {
        let dyads = self.free_dyads();
        let p = self.width();
        let mut x = Dense::zeros(dyads.len(), p);
        let mut y = Vec::with_capacity(dyads.len());
        let mut row = vec![0.0; p];
        for (r, &(i, j)) in dyads.iter().enumerate() {
            self.terms
                .change(&self.state, i as usize, j as usize, &mut row);
            for (c, v) in row.iter().enumerate() {
                x.set(r, c, *v);
            }
            y.push(self.state.has(i as usize, j as usize) as u8 as f64);
        }
        (x, y, dyads)
    }
}
