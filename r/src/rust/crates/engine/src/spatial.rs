//! Geospatial networks: distances between located nodes, proximity
//! networks, edge lengths, group centroids, and dispersion.
//!
//! Nodes carry coordinates either as longitude and latitude in degrees on
//! a sphere ([`Crs::LonLat`]) or as planar x and y ([`Crs::Planar`]).
//! Spatial networks in general are reviewed by Barthelemy (2011,
//! `barthelemy2011spatial`).
//!
//! Sources:
//! - Great-circle distances use the haversine formula (Sinnott 1984,
//!   `sinnott1984virtues`), d = 2 R asin(sqrt(h)) with
//!   h = sin^2(dphi/2) + cos(phi_1) cos(phi_2) sin^2(dlambda/2), which is
//!   well conditioned for small distances; h is clamped to [0, 1] so
//!   rounding cannot push nearly antipodal points out of the domain of
//!   asin. Near antipodal points the formula is ill conditioned and loses
//!   about half of the floating-point digits (an error of up to about
//!   0.2 m on the Earth); elsewhere its rounding error is far below a
//!   millimetre. The default radius is the mean Earth radius
//!   R_1 = 6371008.7714 m of the Geodetic Reference System 1980 (Moritz
//!   2000, `moritz2000geodetic`), so distances are in metres; any other
//!   radius (for example 6371.0088 for kilometres) may be given. A sphere
//!   is an approximation: great-circle distances on it differ from
//!   ellipsoidal (geodesic) distances by up to about 0.5%.
//! - Planar distances are Euclidean, in the units of the coordinates.
//! - Proximity networks join each node to every node within a distance
//!   (a fixed-radius or "disk" network) or to its k nearest neighbours.
//!   Neighbours are found with a k-d tree (Bentley 1975, `bentley1975multidimensional`;
//!   Friedman, Bentley, and Finkel 1977, `friedman1977algorithm`) built over
//!   the planar points, or over the points mapped to the unit sphere in
//!   three dimensions for longitude and latitude: the chord between two
//!   points on the sphere, 2 sin(d / 2R), increases with their great-circle
//!   distance d, so the tree's Euclidean searches find the same
//!   neighbours. Every candidate's distance is then computed exactly with
//!   the requested metric, and the result is decided by those exact
//!   distances, so the tree changes the speed, never the answer.
//! - The centroid of a group is the mean of its planar coordinates, or for
//!   longitude and latitude the normalized mean of the members' unit
//!   vectors converted back to longitude and latitude (the geographic
//!   midpoint, which minimizes the sum of squared chord distances). It is
//!   undefined (NaN) when the mean vector vanishes, as for two antipodal
//!   points.
//! - Dispersion of a group is reported as the mean distance of its members
//!   from the centroid and as the standard distance, the root mean squared
//!   distance from the centroid (Bachi 1963, `bachi1963standard`).
//!
//! Conventions:
//! - Coordinates must be finite. Latitudes must lie in [-90, 90];
//!   longitudes may be any finite number of degrees.
//! - Within-distance networks contain one undirected edge {i, j} for every
//!   pair of distinct nodes at distance at most r (so coincident points are
//!   joined), listed with i < j in increasing order of (i, j).
//! - k-nearest-neighbour networks contain k directed edges from every node
//!   (fewer if there are fewer other nodes), to its nearest other nodes in
//!   order of increasing distance; ties are broken by node index, so the
//!   result is unique. The language bindings derive the symmetric "either"
//!   and "mutual" versions from these arcs.
//! - Edge lengths are the distance between an edge's endpoints in input
//!   order; a self-loop has length 0.
//! - Group labels are any integers; groups are reported in order of their
//!   first member.

use std::collections::HashMap;

use crate::error::{Error, Result};
use crate::graph::Graph;
use crate::parallel;

/// Mean Earth radius R_1 of the Geodetic Reference System 1980, in metres
/// (Moritz 2000).
pub const EARTH_RADIUS_M: f64 = 6_371_008.771_4;

/// How node coordinates are expressed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Crs {
    /// Longitude and latitude in degrees on a sphere of the given radius.
    LonLat {
        /// Sphere radius; distances are in its units.
        radius: f64,
    },
    /// Planar coordinates; distances are Euclidean.
    Planar,
}

impl Crs {
    /// Decodes the C interface codes: 0 = longitude and latitude (with
    /// `radius`), 1 = planar.
    pub fn from_code(code: i32, radius: f64) -> Result<Crs> {
        match code {
            0 => {
                if !(radius.is_finite() && radius > 0.0) {
                    return Err(Error::InvalidArgument(format!(
                        "the sphere radius must be a positive number, not {radius}"
                    )));
                }
                Ok(Crs::LonLat { radius })
            }
            1 => Ok(Crs::Planar),
            _ => Err(Error::InvalidArgument(format!(
                "unknown coordinate system code {code}"
            ))),
        }
    }
}

/// Great-circle distance between two points given in degrees, on a sphere
/// of radius `radius`, by the haversine formula.
pub fn haversine(lon1: f64, lat1: f64, lon2: f64, lat2: f64, radius: f64) -> f64 {
    let (p1, p2) = (lat1.to_radians(), lat2.to_radians());
    let dp = p2 - p1;
    let dl = (lon2 - lon1).to_radians();
    let s1 = (dp / 2.0).sin();
    let s2 = (dl / 2.0).sin();
    let h = (s1 * s1 + p1.cos() * p2.cos() * s2 * s2).clamp(0.0, 1.0);
    2.0 * radius * h.sqrt().asin()
}

/// Located points.
#[derive(Debug, Clone)]
pub struct Points<'a> {
    x: &'a [f64],
    y: &'a [f64],
    crs: Crs,
    /// Largest coordinate magnitude (at least 1), for rounding slack.
    scale: f64,
}

impl<'a> Points<'a> {
    /// Checks and wraps coordinates: x is longitude and y latitude for
    /// [`Crs::LonLat`].
    pub fn new(x: &'a [f64], y: &'a [f64], crs: Crs) -> Result<Points<'a>> {
        if x.len() != y.len() {
            return Err(Error::LengthMismatch {
                what: "y coordinates",
                expected: x.len(),
                found: y.len(),
            });
        }
        if let Some(i) = (0..x.len()).find(|&i| !x[i].is_finite() || !y[i].is_finite()) {
            return Err(Error::InvalidArgument(format!(
                "node {i} has a missing or non-finite coordinate"
            )));
        }
        if let Crs::LonLat { .. } = crs {
            if let Some(i) = y.iter().position(|&v| !(-90.0..=90.0).contains(&v)) {
                return Err(Error::InvalidArgument(format!(
                    "node {i} has latitude {}, outside [-90, 90]",
                    y[i]
                )));
            }
        }
        let scale = match crs {
            Crs::LonLat { .. } => 1.0,
            Crs::Planar => x.iter().chain(y).fold(1.0f64, |m, v| m.max(v.abs())),
        };
        Ok(Points { x, y, crs, scale })
    }

    /// Number of points.
    pub fn len(&self) -> usize {
        self.x.len()
    }

    /// Whether there are no points.
    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }

    /// Distance between points i and j.
    #[inline]
    pub fn distance(&self, i: usize, j: usize) -> f64 {
        match self.crs {
            Crs::LonLat { radius } => haversine(self.x[i], self.y[i], self.x[j], self.y[j], radius),
            Crs::Planar => (self.x[i] - self.x[j]).hypot(self.y[i] - self.y[j]),
        }
    }

    /// The point as a vector in the space the tree searches: the plane, or
    /// the unit sphere in three dimensions.
    fn embed(&self, i: usize) -> [f64; 3] {
        match self.crs {
            Crs::LonLat { .. } => unit_vector(self.x[i], self.y[i]),
            Crs::Planar => [self.x[i], self.y[i], 0.0],
        }
    }

    /// The tree-space (chord or planar) distance corresponding to a metric
    /// distance, slightly enlarged so that rounding never excludes a point
    /// that the exact metric includes.
    fn search_radius(&self, r: f64) -> f64 {
        let base = match self.crs {
            Crs::LonLat { radius } => {
                let theta = (r / radius).min(std::f64::consts::PI);
                2.0 * (theta / 2.0).sin()
            }
            Crs::Planar => r,
        };
        base * (1.0 + 1e-9) + 1e-12 * self.scale()
    }

    fn scale(&self) -> f64 {
        self.scale
    }
}

fn unit_vector(lon: f64, lat: f64) -> [f64; 3] {
    let (l, p) = (lon.to_radians(), lat.to_radians());
    [p.cos() * l.cos(), p.cos() * l.sin(), p.sin()]
}

/// Length of every edge of `g` (in input order) between its endpoints'
/// locations.
pub fn edge_lengths(g: &Graph, points: &Points) -> Result<Vec<f64>> {
    if points.len() != g.node_count() {
        return Err(Error::LengthMismatch {
            what: "coordinates",
            expected: g.node_count(),
            found: points.len(),
        });
    }
    Ok(g.edge_sources()
        .iter()
        .zip(g.edge_targets())
        .map(|(&s, &t)| points.distance(s as usize, t as usize))
        .collect())
}

/// Edges of a proximity network: endpoints and their distance.
#[derive(Debug, Clone, PartialEq)]
pub struct ProximityEdges {
    /// Edge sources.
    pub from: Vec<u32>,
    /// Edge targets.
    pub to: Vec<u32>,
    /// Distances between the endpoints.
    pub distance: Vec<f64>,
}

/// Every pair of distinct points at distance at most `r`, as undirected
/// edges with the smaller index first, in increasing order.
pub fn within_distance(points: &Points, r: f64, threads: usize) -> Result<ProximityEdges> {
    if r.is_nan() || r < 0.0 {
        return Err(Error::InvalidArgument(format!(
            "the distance must be a non-negative number, not {r}"
        )));
    }
    let n = points.len();
    let tree = KdTree::new(points);
    let reach = points.search_radius(r);
    let parts = parallel::for_each_item(
        n,
        threads,
        || (Vec::new(), Vec::new()),
        |(found, scratch): &mut (Vec<(u32, u32, f64)>, Vec<u32>), i| {
            scratch.clear();
            tree.within(&points.embed(i), reach, scratch);
            for &j in scratch.iter() {
                let j = j as usize;
                if j > i {
                    let d = points.distance(i, j);
                    if d <= r {
                        found.push((i as u32, j as u32, d));
                    }
                }
            }
        },
    );
    let mut all: Vec<(u32, u32, f64)> = parts.into_iter().flat_map(|p| p.0).collect();
    all.sort_by_key(|e| (e.0, e.1));
    Ok(unzip(all))
}

/// The `k` nearest other points of every point, as directed edges in
/// order of source, then increasing distance, ties broken by index.
pub fn nearest_neighbors(points: &Points, k: usize, threads: usize) -> Result<ProximityEdges> {
    let n = points.len();
    let k = k.min(n.saturating_sub(1));
    if k == 0 {
        return Ok(unzip(Vec::new()));
    }
    let tree = KdTree::new(points);
    let parts = parallel::for_each_item(
        n,
        threads,
        || (Vec::new(), Vec::new(), Vec::new()),
        |(found, scratch, cand): &mut KnnScratch, i| {
            let q = points.embed(i);
            // The k-th nearest in tree space bounds the exact search: every
            // point whose exact distance ties or beats the k-th must lie
            // within that tree distance (up to rounding, covered by the
            // slack in search_radius).
            let kth = tree.kth_nearest(&q, i, k);
            let reach = kth * (1.0 + 1e-9) + 1e-12 * points.scale();
            scratch.clear();
            tree.within(&q, reach, scratch);
            cand.clear();
            for &j in scratch.iter() {
                if j as usize != i {
                    cand.push((points.distance(i, j as usize), j));
                }
            }
            cand.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
            for &(d, j) in cand.iter().take(k) {
                found.push((i as u32, j, d));
            }
        },
    );
    let mut all: Vec<(u32, u32, f64)> = Vec::with_capacity(n * k);
    for p in parts {
        all.extend(p.0);
    }
    // Stable: keeps each source's neighbours in distance order.
    all.sort_by_key(|e| e.0);
    Ok(unzip(all))
}

type KnnScratch = (Vec<(u32, u32, f64)>, Vec<u32>, Vec<(f64, u32)>);

fn unzip(all: Vec<(u32, u32, f64)>) -> ProximityEdges {
    ProximityEdges {
        from: all.iter().map(|e| e.0).collect(),
        to: all.iter().map(|e| e.1).collect(),
        distance: all.iter().map(|e| e.2).collect(),
    }
}

/// A static k-d tree over embedded points: `order` is a permutation of the
/// points arranged so that each range's median splits it along the axis of
/// largest spread.
struct KdTree {
    pts: Vec<[f64; 3]>,
    order: Vec<u32>,
    axis: Vec<u8>,
}

impl KdTree {
    fn new(points: &Points) -> KdTree {
        let n = points.len();
        let pts: Vec<[f64; 3]> = (0..n).map(|i| points.embed(i)).collect();
        let mut order: Vec<u32> = (0..n as u32).collect();
        let mut axis = vec![0u8; n];
        // Explicit stack of ranges, so deep trees cannot overflow.
        let mut stack = vec![(0usize, n)];
        while let Some((lo, hi)) = stack.pop() {
            if hi - lo <= 1 {
                continue;
            }
            let mut best = (0usize, -1.0f64);
            for a in [0usize, 1, 2] {
                let (mut mn, mut mx) = (f64::INFINITY, f64::NEG_INFINITY);
                for &i in &order[lo..hi] {
                    let v = pts[i as usize][a];
                    mn = mn.min(v);
                    mx = mx.max(v);
                }
                if mx - mn > best.1 {
                    best = (a, mx - mn);
                }
            }
            let a = best.0;
            let mid = (lo + hi) / 2;
            order[lo..hi].select_nth_unstable_by(mid - lo, |&p, &q| {
                pts[p as usize][a]
                    .total_cmp(&pts[q as usize][a])
                    .then(p.cmp(&q))
            });
            axis[mid] = a as u8;
            stack.push((lo, mid));
            stack.push((mid + 1, hi));
        }
        KdTree { pts, order, axis }
    }

    fn dist2(&self, a: &[f64; 3], i: u32) -> f64 {
        let p = &self.pts[i as usize];
        (a[0] - p[0]).powi(2) + (a[1] - p[1]).powi(2) + (a[2] - p[2]).powi(2)
    }

    /// Appends every point within tree distance `r` of `q` (self included).
    fn within(&self, q: &[f64; 3], r: f64, out: &mut Vec<u32>) {
        let r2 = r * r;
        let mut stack = vec![(0usize, self.order.len())];
        while let Some((lo, hi)) = stack.pop() {
            if lo >= hi {
                continue;
            }
            let mid = (lo + hi) / 2;
            let p = self.order[mid];
            if self.dist2(q, p) <= r2 {
                out.push(p);
            }
            if hi - lo == 1 {
                continue;
            }
            let a = self.axis[mid] as usize;
            let diff = q[a] - self.pts[p as usize][a];
            if diff <= r {
                stack.push((lo, mid));
            }
            if diff >= -r {
                stack.push((mid + 1, hi));
            }
        }
    }

    /// Tree distance from `q` to its k-th nearest point other than `skip`.
    fn kth_nearest(&self, q: &[f64; 3], skip: usize, k: usize) -> f64 {
        // Max-heap of the k best squared distances found so far.
        let mut best: std::collections::BinaryHeap<Ordered> = std::collections::BinaryHeap::new();
        // Ranges with a lower bound on the squared distance to their points.
        let mut stack = vec![(0usize, self.order.len(), 0.0f64)];
        while let Some((lo, hi, floor)) = stack.pop() {
            let bound = if best.len() < k {
                f64::INFINITY
            } else {
                best.peek().map_or(f64::INFINITY, |b| b.0)
            };
            if lo >= hi || floor > bound {
                continue;
            }
            let mid = (lo + hi) / 2;
            let p = self.order[mid];
            if p as usize != skip {
                let d = self.dist2(q, p);
                if best.len() < k {
                    best.push(Ordered(d));
                } else if d < bound {
                    best.pop();
                    best.push(Ordered(d));
                }
            }
            if hi - lo == 1 {
                continue;
            }
            let a = self.axis[mid] as usize;
            let diff = q[a] - self.pts[p as usize][a];
            let (near, far) = if diff <= 0.0 {
                ((lo, mid), (mid + 1, hi))
            } else {
                ((mid + 1, hi), (lo, mid))
            };
            // The far side is pushed first so the near side is searched
            // first; its bound is checked again when it is popped.
            stack.push((far.0, far.1, floor.max(diff * diff)));
            stack.push((near.0, near.1, floor));
        }
        best.peek().map_or(f64::INFINITY, |b| b.0).sqrt()
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct Ordered(f64);
impl Eq for Ordered {}
impl PartialOrd for Ordered {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Ordered {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.total_cmp(&other.0)
    }
}

/// Centroid and dispersion of each group of points.
#[derive(Debug, Clone, PartialEq)]
pub struct GroupGeometry {
    /// The group labels, in order of first member.
    pub groups: Vec<i64>,
    /// Number of members of each group.
    pub sizes: Vec<usize>,
    /// Centroid x (longitude) of each group.
    pub x: Vec<f64>,
    /// Centroid y (latitude) of each group.
    pub y: Vec<f64>,
    /// Mean distance of the members from the centroid.
    pub mean_distance: Vec<f64>,
    /// Root mean squared distance from the centroid (standard distance).
    pub standard_distance: Vec<f64>,
}

/// Centroids and dispersion of the groups given by `labels` (one per
/// point). See the module documentation.
pub fn group_geometry(points: &Points, labels: &[i64]) -> Result<GroupGeometry> {
    let n = points.len();
    if labels.len() != n {
        return Err(Error::LengthMismatch {
            what: "groups",
            expected: n,
            found: labels.len(),
        });
    }
    let mut index: HashMap<i64, usize> = HashMap::new();
    let mut groups = Vec::new();
    let member: Vec<usize> = labels
        .iter()
        .map(|&l| {
            *index.entry(l).or_insert_with(|| {
                groups.push(l);
                groups.len() - 1
            })
        })
        .collect();
    let k = groups.len();
    let mut sizes = vec![0usize; k];
    let mut acc = vec![[0.0f64; 3]; k];
    for (i, &g) in member.iter().enumerate() {
        sizes[g] += 1;
        let v = match points.crs {
            Crs::LonLat { .. } => unit_vector(points.x[i], points.y[i]),
            Crs::Planar => [points.x[i], points.y[i], 0.0],
        };
        for (slot, x) in acc[g].iter_mut().zip(v) {
            *slot += x;
        }
    }
    let mut cx = vec![f64::NAN; k];
    let mut cy = vec![f64::NAN; k];
    for g in 0..k {
        let m = sizes[g] as f64;
        let v = [acc[g][0] / m, acc[g][1] / m, acc[g][2] / m];
        match points.crs {
            Crs::Planar => {
                cx[g] = v[0];
                cy[g] = v[1];
            }
            Crs::LonLat { .. } => {
                let norm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
                if norm > 1e-12 {
                    cx[g] = v[1].atan2(v[0]).to_degrees();
                    cy[g] = (v[2] / norm).clamp(-1.0, 1.0).asin().to_degrees();
                }
            }
        }
    }
    let mut sum = vec![0.0; k];
    let mut sum2 = vec![0.0; k];
    for (i, &g) in member.iter().enumerate() {
        let d = match points.crs {
            Crs::LonLat { radius } => haversine(points.x[i], points.y[i], cx[g], cy[g], radius),
            Crs::Planar => (points.x[i] - cx[g]).hypot(points.y[i] - cy[g]),
        };
        sum[g] += d;
        sum2[g] += d * d;
    }
    let mean_distance = (0..k).map(|g| sum[g] / sizes[g] as f64).collect();
    let standard_distance = (0..k).map(|g| (sum2[g] / sizes[g] as f64).sqrt()).collect();
    Ok(GroupGeometry {
        groups,
        sizes,
        x: cx,
        y: cy,
        mean_distance,
        standard_distance,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Stream;

    #[test]
    fn haversine_known_values() {
        // A quarter of a great circle along the equator, and pole to
        // equator.
        let r = EARTH_RADIUS_M;
        let quarter = std::f64::consts::FRAC_PI_2 * r;
        assert!((haversine(0.0, 0.0, 90.0, 0.0, r) - quarter).abs() < 1e-6);
        assert!((haversine(10.0, 90.0, -40.0, 0.0, r) - quarter).abs() < 1e-6);
        // Antipodes: half a circumference, without NaN.
        let half = std::f64::consts::PI * r;
        assert!((haversine(0.0, 0.0, 180.0, 0.0, r) - half).abs() < 1e-6);
        // Near antipodes the formula loses about half the digits: the
        // documented error bound is under a metre on the Earth.
        assert!((haversine(-120.0, 35.0, 60.0, -35.0, r) - half).abs() < 1.0);
        // Longitudes are periodic.
        assert!(
            (haversine(179.0, 10.0, -179.0, 10.0, 1.0) - haversine(0.0, 10.0, 2.0, 10.0, 1.0))
                .abs()
                < 1e-14
        );
        assert_eq!(haversine(5.0, 5.0, 5.0, 5.0, r), 0.0);
    }

    fn brute_within(p: &Points, r: f64) -> Vec<(u32, u32)> {
        let mut out = Vec::new();
        for i in 0..p.len() {
            for j in i + 1..p.len() {
                if p.distance(i, j) <= r {
                    out.push((i as u32, j as u32));
                }
            }
        }
        out
    }

    fn brute_knn(p: &Points, k: usize) -> Vec<(u32, u32)> {
        let mut out = Vec::new();
        for i in 0..p.len() {
            let mut c: Vec<(f64, u32)> = (0..p.len())
                .filter(|&j| j != i)
                .map(|j| (p.distance(i, j), j as u32))
                .collect();
            c.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
            out.extend(c.iter().take(k).map(|&(_, j)| (i as u32, j)));
        }
        out
    }

    fn random_points(n: usize, seed: u64, lonlat: bool) -> (Vec<f64>, Vec<f64>) {
        let mut rng = Stream::new(seed, 0);
        let mut x = Vec::new();
        let mut y = Vec::new();
        for i in 0..n {
            // Some exact duplicates and grid ties.
            if i % 17 == 3 && i > 0 {
                x.push(x[i - 1]);
                y.push(y[i - 1]);
                continue;
            }
            if lonlat {
                x.push(rng.uniform() * 360.0 - 180.0);
                y.push((rng.uniform() * 2.0 - 1.0).asin().to_degrees());
            } else {
                x.push((rng.uniform() * 20.0).round() / 2.0);
                y.push((rng.uniform() * 20.0).round() / 2.0);
            }
        }
        (x, y)
    }

    #[test]
    fn tree_searches_match_brute_force() {
        for &(lonlat, r) in &[(false, 1.5), (true, 2.0e6)] {
            let (x, y) = random_points(400, 7, lonlat);
            let crs = if lonlat {
                Crs::LonLat {
                    radius: EARTH_RADIUS_M,
                }
            } else {
                Crs::Planar
            };
            let p = Points::new(&x, &y, crs).unwrap();
            let w = within_distance(&p, r, 3).unwrap();
            let got: Vec<(u32, u32)> = w.from.iter().copied().zip(w.to.iter().copied()).collect();
            assert_eq!(got, brute_within(&p, r));
            for k in [1, 4, 9] {
                let nn = nearest_neighbors(&p, k, 2).unwrap();
                let got: Vec<(u32, u32)> =
                    nn.from.iter().copied().zip(nn.to.iter().copied()).collect();
                assert_eq!(got, brute_knn(&p, k), "k = {k}");
            }
        }
    }

    #[test]
    #[ignore = "timing check for large inputs; run with --ignored --release"]
    fn hundred_thousand_points() {
        for lonlat in [false, true] {
            let mut rng = Stream::new(11, 0);
            let (x, y): (Vec<f64>, Vec<f64>) = (0..100_000)
                .map(|_| {
                    if lonlat {
                        (
                            rng.uniform() * 360.0 - 180.0,
                            (rng.uniform() * 2.0 - 1.0).asin().to_degrees(),
                        )
                    } else {
                        (rng.uniform() * 100.0, rng.uniform() * 100.0)
                    }
                })
                .unzip();
            let crs = if lonlat {
                Crs::LonLat {
                    radius: EARTH_RADIUS_M,
                }
            } else {
                Crs::Planar
            };
            let p = Points::new(&x, &y, crs).unwrap();
            let start = std::time::Instant::now();
            let nn = nearest_neighbors(&p, 10, 0).unwrap();
            assert_eq!(nn.from.len(), 1_000_000);
            let r = if lonlat { 20_000.0 } else { 0.5 };
            let w = within_distance(&p, r, 0).unwrap();
            eprintln!(
                "lonlat {lonlat}: {} within-distance edges, {:?}",
                w.from.len(),
                start.elapsed()
            );
        }
    }

    #[test]
    fn small_and_degenerate_inputs() {
        let x = [0.0, 3.0, 0.0];
        let y = [0.0, 4.0, 0.0];
        let p = Points::new(&x, &y, Crs::Planar).unwrap();
        let w = within_distance(&p, 5.0, 1).unwrap();
        assert_eq!(w.from, vec![0, 0, 1]);
        assert_eq!(w.to, vec![1, 2, 2]);
        assert_eq!(w.distance, vec![5.0, 0.0, 5.0]);
        let nn = nearest_neighbors(&p, 5, 1).unwrap();
        assert_eq!(nn.from, vec![0, 0, 1, 1, 2, 2]);
        assert_eq!(nn.to, vec![2, 1, 0, 2, 0, 1]);
        assert!(within_distance(&p, -1.0, 1).is_err());
        assert!(Points::new(&[0.0], &[91.0], Crs::LonLat { radius: 1.0 }).is_err());
        assert!(Points::new(&[f64::NAN], &[0.0], Crs::Planar).is_err());
        let one = Points::new(&[1.0], &[1.0], Crs::Planar).unwrap();
        assert!(nearest_neighbors(&one, 3, 1).unwrap().from.is_empty());
    }

    #[test]
    fn edge_lengths_and_groups() {
        let g = Graph::from_edges(3, &[0, 1, 2], &[1, 2, 2], None, false).unwrap();
        let x = [0.0, 3.0, 3.0];
        let y = [0.0, 0.0, 4.0];
        let p = Points::new(&x, &y, Crs::Planar).unwrap();
        assert_eq!(edge_lengths(&g, &p).unwrap(), vec![3.0, 4.0, 0.0]);
        let gg = group_geometry(&p, &[5, 5, 9]).unwrap();
        assert_eq!(gg.groups, vec![5, 9]);
        assert_eq!(gg.x, vec![1.5, 3.0]);
        assert_eq!(gg.y, vec![0.0, 4.0]);
        assert_eq!(gg.mean_distance, vec![1.5, 0.0]);
        assert_eq!(gg.standard_distance, vec![1.5, 0.0]);
        // Geographic midpoint of two points on the equator; antipodes are
        // undefined.
        let lon = [10.0, 30.0, 0.0, 180.0];
        let lat = [0.0, 0.0, 0.0, 0.0];
        let s = Points::new(&lon, &lat, Crs::LonLat { radius: 1.0 }).unwrap();
        let gg = group_geometry(&s, &[1, 1, 2, 2]).unwrap();
        assert!((gg.x[0] - 20.0).abs() < 1e-12 && gg.y[0].abs() < 1e-12);
        assert!((gg.mean_distance[0] - 10f64.to_radians()).abs() < 1e-12);
        assert!(gg.x[1].is_nan() && gg.mean_distance[1].is_nan());
    }
}
