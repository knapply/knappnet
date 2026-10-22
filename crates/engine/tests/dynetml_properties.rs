//! Randomized round trips through DyNetML: for random multimode networks
//! (several node classes, networks, periods, directed and undirected
//! networks, self-loops, parallel links, and awkward text), reading what was
//! written gives back the same tables.

extern crate knappnet as engine;

use engine::io::{
    read_str, write_string, Column, Format, Network, NetworkBuilder, ReadOptions, Value,
    WriteOptions,
};
use engine::rng::Stream;

const CLASSES: [&str; 4] = ["Agent", "Organization", "Location", "Task"];
const TEXT: [&str; 5] = ["plain", "a & b", "<tag>", "quote \"q\"", "line\nbreak"];

/// Builds a network as the DyNetML reader would produce it, so the round
/// trip must be exact.
fn random_network(seed: u64) -> Network {
    let mut rng = Stream::new(seed, 0);
    let periods = 1 + rng.below(3) as usize;
    let multi = periods > 1;
    let n = 1 + rng.below(12) as usize;
    let shared_ids = rng.below(2) == 0;
    // Each network links one ordered pair of classes in one direction.
    let mut nets: Vec<(usize, usize, bool)> = Vec::new();
    for _ in 0..1 + rng.below(3) {
        let (a, b) = (rng.below(4) as usize, rng.below(4) as usize);
        nets.push((a, b, rng.below(2) == 0));
    }
    let any_dir = nets.iter().any(|x| x.2);
    let any_undir = nets.iter().any(|x| !x.2);
    let directed = any_dir || !any_undir;
    let mixed = any_dir && any_undir;
    let weighted = rng.below(2) == 0;

    let mut b = NetworkBuilder::new(directed);
    // Nodes per period and class, in class order as the writer groups them.
    let mut members: Vec<Vec<Vec<String>>> = vec![vec![Vec::new(); 4]; periods];
    let mut rows = Vec::new();
    for p in 0..periods {
        for v in 0..n {
            let c = rng.below(4) as usize;
            let raw = if shared_ids {
                format!("n{}", rng.below(4))
            } else {
                format!("n{v}")
            };
            rows.push((p, c, raw));
        }
    }
    // Unique (period, class, id); order by period then class.
    rows.sort();
    rows.dedup();
    let collide = rows
        .iter()
        .any(|r| rows.iter().filter(|s| s.0 == r.0 && s.2 == r.2).count() > 1);
    let qualified = multi || collide;
    for (p, c, raw) in &rows {
        let mut parts = Vec::new();
        if multi {
            parts.push(format!("t{p}"));
        }
        if collide {
            parts.push(CLASSES[*c].to_string());
        }
        parts.push(raw.clone());
        let id = parts.join("::");
        let mut attrs = Vec::new();
        if qualified {
            attrs.push(("node_id".to_string(), Value::Str(raw.clone())));
        }
        attrs.push(("class".into(), Value::Str(CLASSES[*c].into())));
        attrs.push(("type".into(), Value::Str(CLASSES[*c].into())));
        if multi {
            attrs.push(("period".into(), Value::Str(format!("t{p}"))));
        }
        if rng.below(2) == 0 {
            attrs.push((
                "title".into(),
                Value::Str(TEXT[rng.below(5) as usize].into()),
            ));
        }
        if rng.below(2) == 0 {
            attrs.push(("score".into(), Value::Float(rng.uniform() * 100.0 - 50.0)));
        }
        if rng.below(3) == 0 {
            attrs.push((
                "note".into(),
                Value::Str(TEXT[rng.below(5) as usize].into()),
            ));
        }
        b.add_node(&id, attrs);
        members[*p][*c].push(id);
    }
    b.node_kind("title", engine::io::Kind::Str);
    b.node_kind("score", engine::io::Kind::Float);
    b.node_kind("note", engine::io::Kind::Str);

    b.edge_kind("network", engine::io::Kind::Str);
    if multi {
        b.edge_kind("period", engine::io::Kind::Str);
    }
    if mixed {
        b.edge_kind("directed", engine::io::Kind::Bool);
    }
    if weighted {
        b.edge_kind("weight", engine::io::Kind::Float);
    }
    // Links grouped by network, as the reader produces them.
    for (p, period) in members.iter().enumerate() {
        for (k, &(ca, cb, dir)) in nets.iter().enumerate() {
            let (sa, sb) = (&period[ca], &period[cb]);
            if sa.is_empty() || sb.is_empty() {
                continue;
            }
            for _ in 0..rng.below(6) {
                let f = sa[rng.below(sa.len() as u64) as usize].clone();
                let t = sb[rng.below(sb.len() as u64) as usize].clone();
                let mut attrs = vec![("network".to_string(), Value::Str(format!("net{k}")))];
                if multi {
                    attrs.push(("period".into(), Value::Str(format!("t{p}"))));
                }
                if mixed {
                    attrs.push(("directed".into(), Value::Bool(dir)));
                }
                if weighted {
                    // Integers and binary fractions print and parse exactly.
                    attrs.push(("weight".into(), Value::Float(rng.below(20) as f64 / 4.0)));
                }
                b.add_edge(&f, &t, attrs.clone());
                if mixed && !dir && f != t {
                    b.add_edge(&t, &f, attrs);
                }
            }
        }
    }
    b.finish().unwrap()
}

#[test]
fn dynetml_round_trips_random_multimode_networks() {
    for seed in 0..500 {
        let net = random_network(seed);
        let text = write_string(Format::DynetMl, &net, &WriteOptions::default()).unwrap();
        let back = read_str(Format::DynetMl, &text, &ReadOptions::default())
            .unwrap_or_else(|e| panic!("seed {seed}: {e}\n{text}"));
        // A mixed table whose directed networks have no links is written as
        // undirected networks only, and reads back undirected with each
        // link once.
        if let Some(Column::Bool(d)) = net.edges.column("directed") {
            if !d.is_empty() && d.iter().all(|x| *x == Some(false)) {
                assert!(!back.directed, "seed {seed}");
                let from = net.edges.strings("from").unwrap();
                let to = net.edges.strings("to").unwrap();
                let loops = from.iter().zip(&to).filter(|(f, t)| f == t).count();
                assert_eq!(back.edges.nrow(), loops + (from.len() - loops) / 2);
                continue;
            }
        }
        assert_eq!(back.directed, net.directed, "seed {seed}");
        // Columns with no values at all may be dropped or typed differently;
        // every other column must come back unchanged, and nothing else may
        // appear.
        let filled = |t: &engine::io::Table| -> Vec<String> {
            t.names
                .iter()
                .zip(&t.columns)
                .filter(|(_, c)| (0..c.len()).any(|i| c.get(i) != Value::Null))
                .map(|(n, _)| n.clone())
                .collect::<std::collections::BTreeSet<_>>()
                .into_iter()
                .collect()
        };
        assert_eq!(
            filled(&back.nodes),
            filled(&net.nodes),
            "seed {seed}\n{text}"
        );
        for name in filled(&net.nodes) {
            assert_eq!(
                back.nodes.column(&name),
                net.nodes.column(&name),
                "seed {seed} {name}\n{text}"
            );
        }
        // A file without links cannot say whether links would be weighted.
        assert_eq!(back.edges.nrow(), net.edges.nrow(), "seed {seed}");
        if net.edges.nrow() > 0 {
            let mut expected = net.edges.clone();
            // Only networks with links are written, so a `directed` column
            // that is true everywhere (no undirected network has links)
            // cannot be told apart from a directed file.
            if let Some(Column::Bool(d)) = expected.column("directed") {
                if d.iter().all(|x| *x == Some(true)) {
                    let pos = expected.names.iter().position(|n| n == "directed").unwrap();
                    expected.names.remove(pos);
                    expected.columns.remove(pos);
                }
            }
            assert_eq!(back.edges, expected, "seed {seed}\n{text}");
        }
        assert!(matches!(back.nodes.column("class"), Some(Column::Str(_))));
    }
}
