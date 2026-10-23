//! Network panel data: waves, covariates, and their centering.
//!
//! Centering follows the RSiena manual (Ripley et al. 2026,
//! `ripley2026manual`, section 5.3.2): actor covariates are centered by the
//! mean over all actors and, for changing covariates, over all periods;
//! dyadic covariates by the mean of their off-diagonal cells; the
//! similarity sim_ij = 1 - |v_i - v_j| / range by its mean over all ordered
//! pairs of distinct actors (and periods). Here `range` is the difference
//! between the largest and smallest value over all actors and periods.
//!
//! A changing covariate has one value per actor and period, the value in
//! effect during period m (from wave m to wave m + 1), so M waves give
//! M - 1 columns. Values for the last wave may be supplied and are then
//! ignored, as RSiena ignores them, and a changing covariate needs at least
//! three waves. RSiena centers the similarity of changing covariates by a
//! mean computed differently (not documented in the manual), so the `simX`
//! effect of a changing covariate can differ from RSiena's; see
//! validation/manifest.d/saom.yaml.

use super::network::Network;
use crate::error::{Error, Result};

/// An actor covariate, constant or changing over waves.
#[derive(Debug, Clone)]
pub struct ActorCovariate {
    /// Raw values: `n` for a constant covariate, `n * waves` (wave-major:
    /// `values[w * n + i]`) for a changing one.
    values: Vec<f64>,
    /// Centered values, laid out like `values`.
    centered: Vec<f64>,
    waves: usize,
    n: usize,
    mean: f64,
    range: f64,
    sim_mean: f64,
}

impl ActorCovariate {
    /// Builds a covariate from raw values: `n` values for a constant
    /// covariate, `n * waves` for a changing one (wave-major).
    pub fn new(n: usize, waves: usize, values: &[f64]) -> Result<ActorCovariate> {
        if waves == 0 || values.len() != n * waves {
            return Err(Error::LengthMismatch {
                what: "actor covariate values",
                expected: n * waves,
                found: values.len(),
            });
        }
        if let Some(k) = values.iter().position(|v| !v.is_finite()) {
            return Err(Error::InvalidArgument(format!(
                "actor covariate value {k} is missing or not finite; missing covariate \
                 values are not supported"
            )));
        }
        let count = values.len() as f64;
        let mean = values.iter().sum::<f64>() / count;
        let min = values.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let range = max - min;
        let centered: Vec<f64> = values.iter().map(|v| v - mean).collect();
        // Similarity mean over ordered pairs of distinct actors and waves.
        let mut sim_mean = 0.0;
        if range > 0.0 && n > 1 {
            let mut total = 0.0;
            for w in 0..waves {
                let v = &values[w * n..(w + 1) * n];
                for i in 0..n {
                    for j in 0..n {
                        if i != j {
                            total += 1.0 - (v[i] - v[j]).abs() / range;
                        }
                    }
                }
            }
            sim_mean = total / (waves as f64 * (n * (n - 1)) as f64);
        }
        Ok(ActorCovariate {
            values: values.to_vec(),
            centered,
            waves,
            n,
            mean,
            range,
            sim_mean,
        })
    }

    /// Whether the covariate changes over waves.
    pub fn is_changing(&self) -> bool {
        self.waves > 1
    }

    /// Mean subtracted for centering.
    pub fn mean(&self) -> f64 {
        self.mean
    }

    /// Observed range (largest minus smallest value).
    pub fn range(&self) -> f64 {
        self.range
    }

    /// Mean similarity over ordered pairs of distinct actors (and waves).
    pub fn sim_mean(&self) -> f64 {
        self.sim_mean
    }

    /// Centered values in effect during `period` (0-based): wave `period`
    /// for a changing covariate.
    pub fn centered(&self, period: usize) -> &[f64] {
        let w = if self.waves > 1 {
            period.min(self.waves - 1)
        } else {
            0
        };
        &self.centered[w * self.n..(w + 1) * self.n]
    }

    /// Raw values in effect during `period`.
    pub fn raw(&self, period: usize) -> &[f64] {
        let w = if self.waves > 1 {
            period.min(self.waves - 1)
        } else {
            0
        };
        &self.values[w * self.n..(w + 1) * self.n]
    }

    /// Centered similarity of actors i and j in `period`.
    pub fn similarity(&self, period: usize, i: usize, j: usize) -> f64 {
        let v = self.raw(period);
        1.0 - (v[i] - v[j]).abs() / self.range - self.sim_mean
    }
}

/// A constant dyadic covariate w_ij.
#[derive(Debug, Clone)]
pub struct DyadicCovariate {
    /// Row-major values, `values[i * n + j]`.
    values: Vec<f64>,
    n: usize,
    mean: f64,
}

impl DyadicCovariate {
    /// Builds a dyadic covariate from a row-major n-by-n matrix. The
    /// diagonal is ignored.
    pub fn new(n: usize, values: &[f64]) -> Result<DyadicCovariate> {
        if values.len() != n * n {
            return Err(Error::LengthMismatch {
                what: "dyadic covariate values",
                expected: n * n,
                found: values.len(),
            });
        }
        let mut total = 0.0;
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    continue;
                }
                let v = values[i * n + j];
                if !v.is_finite() {
                    return Err(Error::InvalidArgument(format!(
                        "dyadic covariate value ({i}, {j}) is missing or not finite; missing \
                         covariate values are not supported"
                    )));
                }
                total += v;
            }
        }
        let mean = if n > 1 {
            total / (n * (n - 1)) as f64
        } else {
            0.0
        };
        Ok(DyadicCovariate {
            values: values.to_vec(),
            n,
            mean,
        })
    }

    /// Mean of the off-diagonal cells.
    pub fn mean(&self) -> f64 {
        self.mean
    }

    /// Centered value w_ij - mean.
    #[inline]
    pub fn centered(&self, i: usize, j: usize) -> f64 {
        self.values[i * self.n + j] - self.mean
    }
}

/// A panel of networks on the same actors with their covariates.
#[derive(Debug, Clone)]
pub struct Panel {
    n: usize,
    waves: Vec<Network>,
    actor: Vec<ActorCovariate>,
    dyadic: Vec<DyadicCovariate>,
}

impl Panel {
    /// Builds a panel from row-major 0/1 adjacency matrices, one per wave.
    /// At least two waves are needed.
    pub fn from_adjacency(n: usize, waves: &[Vec<u8>]) -> Result<Panel> {
        if waves.len() < 2 {
            return Err(Error::InvalidArgument(
                "a network panel needs at least two waves".into(),
            ));
        }
        let mut nets = Vec::with_capacity(waves.len());
        for (w, adj) in waves.iter().enumerate() {
            if adj.len() != n * n {
                return Err(Error::LengthMismatch {
                    what: "adjacency matrix",
                    expected: n * n,
                    found: adj.len(),
                });
            }
            let _ = w;
            nets.push(Network::from_adjacency(n, adj));
        }
        Ok(Panel {
            n,
            waves: nets,
            actor: Vec::new(),
            dyadic: Vec::new(),
        })
    }

    /// Builds a panel from edge lists of 0-based indices, one per wave.
    pub fn from_edge_lists(n: usize, waves: &[Vec<(usize, usize)>]) -> Result<Panel> {
        let mut adj = Vec::with_capacity(waves.len());
        for edges in waves {
            let mut a = vec![0u8; n * n];
            for &(i, j) in edges {
                if i >= n || j >= n {
                    return Err(Error::NodeOutOfRange {
                        edge: 0,
                        node: i.max(j) as i64,
                        n,
                    });
                }
                if i != j {
                    a[i * n + j] = 1;
                }
            }
            adj.push(a);
        }
        Panel::from_adjacency(n, &adj)
    }

    /// Adds an actor covariate (constant: `n` values; changing: `n * (M - 1)`
    /// values, one column per period, or `n * M`, whose last column is
    /// ignored; wave-major, where M is the number of waves). Returns its
    /// index.
    pub fn add_actor_covariate(&mut self, values: &[f64]) -> Result<usize> {
        let (n, m) = (self.n, self.waves.len());
        let periods = m.saturating_sub(1);
        let values = if values.len() == n {
            values
        } else if m < 3 && (values.len() == n * m || values.len() == n * periods) && n > 0 {
            return Err(Error::InvalidArgument(
                "a changing covariate needs at least three waves (with two waves, use a \
                 constant covariate)"
                    .into(),
            ));
        } else if values.len() == n * m {
            &values[..n * periods]
        } else {
            values
        };
        let waves = if values.len() == n { 1 } else { periods };
        self.actor.push(ActorCovariate::new(n, waves, values)?);
        Ok(self.actor.len() - 1)
    }

    /// Adds a constant dyadic covariate (row-major n-by-n). Returns its
    /// index.
    pub fn add_dyadic_covariate(&mut self, values: &[f64]) -> Result<usize> {
        self.dyadic.push(DyadicCovariate::new(self.n, values)?);
        Ok(self.dyadic.len() - 1)
    }

    /// Number of actors.
    pub fn n(&self) -> usize {
        self.n
    }

    /// Number of waves M.
    pub fn wave_count(&self) -> usize {
        self.waves.len()
    }

    /// Number of periods M - 1.
    pub fn period_count(&self) -> usize {
        self.waves.len() - 1
    }

    /// The network observed at wave `w` (0-based).
    pub fn wave(&self, w: usize) -> &Network {
        &self.waves[w]
    }

    /// Actor covariates.
    pub fn actor_covariates(&self) -> &[ActorCovariate] {
        &self.actor
    }

    /// Dyadic covariates.
    pub fn dyadic_covariates(&self) -> &[DyadicCovariate] {
        &self.dyadic
    }

    /// Number of tie variables that differ between the start and end of
    /// each period: the target statistic of the rate parameters and the
    /// stopping criterion of conditional simulation.
    pub fn distances(&self) -> Vec<f64> {
        (0..self.period_count())
            .map(|m| self.waves[m].distance(&self.waves[m + 1]) as f64)
            .collect()
    }

    /// Ties created and dissolved in each period.
    pub fn creations_and_dissolutions(&self) -> Vec<(usize, usize)> {
        let n = self.n;
        (0..self.period_count())
            .map(|m| {
                let (a, b) = (&self.waves[m], &self.waves[m + 1]);
                let mut up = 0;
                let mut down = 0;
                for i in 0..n {
                    for j in 0..n {
                        match (a.has(i, j), b.has(i, j)) {
                            (false, true) => up += 1,
                            (true, false) => down += 1,
                            _ => {}
                        }
                    }
                }
                (up, down)
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn covariates_are_centered_as_rsiena_centers_them() {
        let c = ActorCovariate::new(3, 1, &[1.0, 2.0, 6.0]).unwrap();
        assert_eq!(c.mean(), 3.0);
        assert_eq!(c.range(), 5.0);
        assert_eq!(c.centered(0), &[-2.0, -1.0, 3.0]);
        // Similarities: (1,2) 0.8, (1,6) 0, (2,6) 0.2; each twice.
        assert!((c.sim_mean() - (2.0 * (0.8 + 0.0 + 0.2) / 6.0)).abs() < 1e-12);
        assert!((c.similarity(0, 0, 1) - (0.8 - c.sim_mean())).abs() < 1e-12);

        let changing = ActorCovariate::new(2, 2, &[0.0, 1.0, 3.0, 4.0]).unwrap();
        assert_eq!(changing.mean(), 2.0);
        assert_eq!(changing.centered(1), &[1.0, 2.0]);
        assert_eq!(changing.centered(5), &[1.0, 2.0]);
        assert_eq!(changing.range(), 4.0);
        // Every pair differs by 1, so sim = 0.75 everywhere.
        assert!((changing.sim_mean() - 0.75).abs() < 1e-12);

        let d = DyadicCovariate::new(2, &[9.0, 1.0, 3.0, 9.0]).unwrap();
        assert_eq!(d.mean(), 2.0);
        assert_eq!(d.centered(0, 1), -1.0);
        assert!(ActorCovariate::new(2, 1, &[1.0, f64::NAN]).is_err());
    }

    #[test]
    fn distances_count_changed_ties() {
        let p = Panel::from_edge_lists(3, &[vec![(0, 1), (1, 2)], vec![(0, 1), (2, 1), (1, 1)]])
            .unwrap();
        assert_eq!(p.distances(), vec![2.0]);
        assert_eq!(p.creations_and_dissolutions(), vec![(1, 1)]);
        assert!(Panel::from_edge_lists(3, &[vec![]]).is_err());
    }
}
