//! Timing of frame validation and construction on large random networks
//! (1 and 10 million edges, integer and text identifiers, with and
//! without the optional roles). Ignored by default; run with
//!
//! ```text
//! cargo test --release --test frame_perf -- --ignored --nocapture
//! ```
//!
//! Prints the time of validation alone, of validation plus the graph
//! build, and of the summary, for one thread and for all cores.

extern crate knappnet as engine;

use std::time::Instant;

use engine::frame::{build, summarize, Column, EdgeRole, Frame, NodeRole, Strings, WeightType};
use engine::rng::Stream;

fn frame(n: usize, m: usize, text: bool, rich: bool) -> Frame {
    let mut rng = Stream::new(7, 0);
    let ids: Vec<i64> = (0..n as i64).map(|i| 1000 + i).collect();
    let mut src = Vec::with_capacity(m);
    let mut tgt = Vec::with_capacity(m);
    for _ in 0..m {
        src.push(ids[rng.below(n as u64) as usize]);
        tgt.push(ids[rng.below(n as u64) as usize]);
    }
    let column = |v: &[i64]| {
        if text {
            let s: Vec<Option<String>> = v.iter().map(|x| Some(format!("node{x}"))).collect();
            Column::Str(Strings::from_strs(&s))
        } else {
            Column::Int(v.to_vec(), None)
        }
    };
    let mut f = Frame::new(true);
    f.add_node_column("id", column(&ids)).unwrap();
    f.set_node_role(NodeRole::Id, "id").unwrap();
    f.add_edge_column("from", column(&src)).unwrap();
    f.add_edge_column("to", column(&tgt)).unwrap();
    f.set_edge_role(EdgeRole::Source, "from").unwrap();
    f.set_edge_role(EdgeRole::Target, "to").unwrap();
    if rich {
        let w: Vec<f64> = (0..m).map(|_| 1.0 + rng.below(9) as f64).collect();
        let rel: Vec<i64> = (0..m).map(|_| rng.below(3) as i64).collect();
        let start: Vec<f64> = (0..m).map(|_| rng.below(1000) as f64).collect();
        let end: Vec<f64> = start.iter().map(|s| s + 5.0).collect();
        let modes: Vec<i64> = (0..n).map(|_| rng.below(2) as i64).collect();
        f.add_edge_column("w", Column::Float(w)).unwrap();
        f.add_edge_column("rel", Column::Int(rel, None)).unwrap();
        f.add_edge_column("s", Column::Float(start)).unwrap();
        f.add_edge_column("e", Column::Float(end)).unwrap();
        f.add_node_column("mode", Column::Int(modes, None)).unwrap();
        f.set_edge_role(EdgeRole::Weight, "w").unwrap();
        f.set_weight_type(Some(WeightType::Strength));
        f.set_edge_role(EdgeRole::Relation, "rel").unwrap();
        f.set_edge_role(EdgeRole::Start, "s").unwrap();
        f.set_edge_role(EdgeRole::End, "e").unwrap();
        f.set_node_role(NodeRole::Mode, "mode").unwrap();
    }
    f
}

fn time<T>(f: impl FnOnce() -> T) -> (T, f64) {
    let t = Instant::now();
    let r = f();
    (r, t.elapsed().as_secs_f64())
}

#[test]
#[ignore]
fn frame_timings() {
    for &(n, m) in &[(100_000usize, 1_000_000usize), (1_000_000, 10_000_000)] {
        for &(text, rich) in &[(false, false), (true, false), (true, true)] {
            let f = frame(n, m, text, rich);
            for &threads in &[1usize, 0] {
                let (v, t_validate) = time(|| f.validate(threads));
                assert!(v.report.is_valid());
                let (b, t_build) = time(|| build(&f, threads).unwrap());
                let (_, t_summary) = time(|| summarize(&f, &b));
                println!(
                    "{m:>9} edges, {} ids, {} roles, {:>3} thread(s): validate {t_validate:.2}s, validate+build {t_build:.2}s, summary {t_summary:.2}s",
                    if text { "text" } else { "int " },
                    if rich { "all " } else { "none" },
                    if threads == 0 { "all".to_string() } else { threads.to_string() },
                );
            }
        }
    }
}
