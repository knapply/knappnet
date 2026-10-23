//! Timing of file reading and writing on a network the size of the
//! benchmark in `bench/bench.R` (10,000 nodes, 50,000 edges). Ignored by
//! default; run with
//!
//! ```text
//! cargo test --release --test io_perf -- --ignored --nocapture
//! ```
//!
//! Prints the median of several runs for each format, with and without
//! attributes.

extern crate knappnet as engine;

use std::time::Instant;

use engine::io::{read_str, write_string, Format, Network, NetworkBuilder, ReadOptions, Value};
use engine::rng::Stream;

fn network(attributes: bool) -> Network {
    let (n, m) = (10_000u64, 50_000);
    let mut rng = Stream::new(1, 0);
    let mut b = NetworkBuilder::new(false);
    for v in 1..=n {
        let attrs = if attributes {
            vec![
                ("group".to_string(), Value::Int((v % 7) as i64)),
                ("score".to_string(), Value::Float(rng.uniform())),
            ]
        } else {
            Vec::new()
        };
        b.add_node(&v.to_string(), attrs);
    }
    for _ in 0..m {
        let a = 1 + rng.below(n);
        let c = 1 + rng.below(n);
        let attrs = if attributes {
            vec![(
                "weight".to_string(),
                Value::Float(1.0 + rng.below(5) as f64),
            )]
        } else {
            Vec::new()
        };
        b.add_edge(&a.to_string(), &c.to_string(), attrs);
    }
    b.finish().unwrap()
}

fn median_ms(mut f: impl FnMut()) -> f64 {
    let mut times: Vec<f64> = (0..7)
        .map(|_| {
            let t = Instant::now();
            f();
            t.elapsed().as_secs_f64() * 1e3
        })
        .collect();
    times.sort_by(|a, b| a.partial_cmp(b).unwrap());
    times[times.len() / 2]
}

#[test]
#[ignore = "timing; run with --release --ignored --nocapture"]
fn read_and_write_times() {
    for attributes in [false, true] {
        let net = network(attributes);
        for format in [
            Format::GraphMl,
            Format::Gexf,
            Format::Gml,
            Format::Pajek,
            Format::Dl,
            Format::NodeLink,
        ] {
            let text = write_string(format, &net, &Default::default()).unwrap();
            let write = median_ms(|| {
                write_string(format, &net, &Default::default()).unwrap();
            });
            let read = median_ms(|| {
                read_str(format, &text, &ReadOptions::default()).unwrap();
            });
            println!(
                "{:<9} attributes {:<5} {:>6.1} MB  read {:>7.1} ms  write {:>7.1} ms",
                format.name(),
                attributes,
                text.len() as f64 / 1e6,
                read,
                write
            );
        }
    }
}

#[test]
#[ignore]
fn tmp_profile_graphml() {
    let net = network(false);
    let text = write_string(Format::GraphMl, &net, &Default::default()).unwrap();
    for _ in 0..3 {
        read_str(Format::GraphMl, &text, &ReadOptions::default()).unwrap();
    }
}
