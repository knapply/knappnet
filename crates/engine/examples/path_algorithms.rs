//! Times betweenness under each shortest-path algorithm on one random
//! network with small integer weights.
//!
//!     cargo run --release --example path_algorithms

use std::time::Instant;

use knappnet::{betweenness, BetweennessOptions, Graph, PathAlgorithm};

fn main() {
    let n = 4000usize;
    let m = 32_000usize;
    let mut state = 7u64;
    let mut next = |k: u64| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state % k
    };
    let src: Vec<i64> = (0..m).map(|_| next(n as u64) as i64).collect();
    let tgt: Vec<i64> = (0..m).map(|_| next(n as u64) as i64).collect();
    let w: Vec<f64> = (0..m).map(|_| (1 + next(8)) as f64).collect();
    let g = Graph::from_edges(n, &src, &tgt, Some(&w), true).unwrap();
    println!("n = {n}, m = {m}, integer weights 1..8");
    for (label, weighted, algorithm) in [
        ("unweighted (BFS)", false, PathAlgorithm::Auto),
        ("weighted, Dial", true, PathAlgorithm::Dial),
        ("weighted, Dijkstra", true, PathAlgorithm::Dijkstra),
    ] {
        let opts = BetweennessOptions {
            weighted,
            algorithm,
            threads: 1,
            ..Default::default()
        };
        let t = Instant::now();
        betweenness(&g, opts).unwrap();
        println!("{label:<22}{:>8.3} s", t.elapsed().as_secs_f64());
    }
}
