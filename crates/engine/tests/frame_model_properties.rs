//! Randomized property tests for the parts of the data model beyond one
//! node table and one edge table: multilayer edges (classification,
//! per-layer counts and subgraphs, the layer-aware parallel-edge rule),
//! unobserved pairs (excluded from the built graph, kept and counted),
//! collections of ego networks (per-ego counts and subgraphs), and the
//! node spells and node changes tables. Each property is checked against
//! a direct brute-force reading of the specification in
//! `docs/DATA_MODEL.md`, and every validation is run with one and with
//! four threads to check that the diagnostics do not depend on the
//! thread count.

extern crate knappnet as engine;

use std::collections::HashMap;

use engine::frame::diagnostics::ROW_CAP;
use engine::frame::multiplicity::multiplicity;
use engine::frame::{
    build, ego_subgraph, layer_subgraph, missing_pairs, summarize, ChangeRole, Code, Column,
    EdgeRole, Frame, LayerRole, NodeRole, Severity, SpellRole, Strings, Validation, Which,
};

/// SplitMix64 (Steele, Lea, and Flood 2014, `steele2014fast`).
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
    fn flag(&mut self) -> bool {
        self.below(2) == 0
    }
    fn shuffle<T>(&mut self, v: &mut [T]) {
        for i in (1..v.len()).rev() {
            let j = self.below(i as u64 + 1) as usize;
            v.swap(i, j);
        }
    }
}

fn text(values: &[String]) -> Column {
    let v: Vec<Option<&str>> = values.iter().map(|s| Some(s.as_str())).collect();
    Column::Str(Strings::from_strs(&v))
}

/// Validates with one and four threads and checks the two agree.
fn validate(f: &Frame) -> Validation {
    let v1 = f.validate(1);
    let v4 = f.validate(4);
    assert_eq!(v1.report, v4.report, "diagnostics depend on threads");
    v1
}

fn error_codes(v: &Validation) -> Vec<Code> {
    v.report
        .diagnostics
        .iter()
        .filter(|d| d.severity() == Severity::Error)
        .map(|d| d.code)
        .collect()
}

/// Checks that `code` is the only error and lists exactly `rows`.
fn only_error(v: &Validation, code: Code, rows: &[usize]) {
    assert_eq!(error_codes(v), vec![code], "{:?}", v.report);
    let d = v.report.find(code).unwrap();
    assert_eq!(d.count, rows.len());
    assert_eq!(d.rows, rows[..rows.len().min(ROW_CAP)].to_vec());
}

/// A node table of integer identifiers `10, 11, ...` and an edge table.
fn base(directed: bool, n: usize, src: &[usize], tgt: &[usize]) -> Frame {
    let mut f = Frame::new(directed);
    let ids: Vec<i64> = (0..n as i64).map(|i| 10 + i).collect();
    f.add_node_column("id", Column::Int(ids, None)).unwrap();
    f.set_node_role(NodeRole::Id, "id").unwrap();
    let col = |v: &[usize]| Column::Int(v.iter().map(|&i| 10 + i as i64).collect(), None);
    f.add_edge_column("from", col(src)).unwrap();
    f.add_edge_column("to", col(tgt)).unwrap();
    f.set_edge_role(EdgeRole::Source, "from").unwrap();
    f.set_edge_role(EdgeRole::Target, "to").unwrap();
    f
}

/// Random endpoints for `m` edges among `n` nodes.
fn endpoints(rng: &mut Rng, n: usize, m: usize) -> (Vec<usize>, Vec<usize>) {
    let src = (0..m).map(|_| rng.below(n as u64) as usize).collect();
    let tgt = (0..m).map(|_| rng.below(n as u64) as usize).collect();
    (src, tgt)
}

/// The class of an edge from its actors and layers: 0 intralayer,
/// 1 coupling, 2 interlayer.
fn class(a: usize, b: usize, la: usize, lb: usize) -> usize {
    if la == lb {
        0
    } else if a == b {
        1
    } else {
        2
    }
}

#[test]
fn multilayer_edges_are_classified_and_split_by_layer() {
    let mut rng = Rng(21);
    for case in 0..300 {
        let directed = rng.flag();
        let n = 1 + rng.below(12) as usize;
        let m = rng.below(40) as usize;
        let nl = 1 + rng.below(4) as usize;
        let (src, tgt) = endpoints(&mut rng, n, m);
        // Mostly intralayer edges, some couplings and interlayer edges.
        let ls: Vec<usize> = (0..m).map(|_| rng.below(nl as u64) as usize).collect();
        let lt: Vec<usize> = (0..m)
            .map(|e| {
                if rng.below(3) == 0 {
                    rng.below(nl as u64) as usize
                } else {
                    ls[e]
                }
            })
            .collect();
        let name = |k: usize| format!("L{k}");
        let mut f = base(directed, n, &src, &tgt);
        f.add_edge_column("ls", text(&ls.iter().map(|&k| name(k)).collect::<Vec<_>>()))
            .unwrap();
        f.add_edge_column("lt", text(&lt.iter().map(|&k| name(k)).collect::<Vec<_>>()))
            .unwrap();
        f.set_edge_role(EdgeRole::SourceLayer, "ls").unwrap();
        f.set_edge_role(EdgeRole::TargetLayer, "lt").unwrap();
        // The expected layer order: the table's, else first appearance
        // (source before target within a row).
        let with_table = rng.flag();
        let order: Vec<usize> = if with_table {
            let mut o: Vec<usize> = (0..nl).collect();
            rng.shuffle(&mut o);
            f.add_column(
                Which::Layers,
                "layer",
                text(&o.iter().map(|&k| name(k)).collect::<Vec<_>>()),
            )
            .unwrap();
            f.set_layer_role(LayerRole::Layer, "layer").unwrap();
            o
        } else {
            let mut o = Vec::new();
            for e in 0..m {
                for k in [ls[e], lt[e]] {
                    if !o.contains(&k) {
                        o.push(k);
                    }
                }
            }
            o
        };
        let v = validate(&f);
        assert!(v.report.is_valid(), "case {case}: {:?}", v.report);
        let built = build(&f, 2).unwrap();
        let s = summarize(&f, &built);
        let labels: Vec<String> = order.iter().map(|&k| name(k)).collect();
        assert_eq!(s.layer_labels, labels);
        assert_eq!(s.layer_table, with_table);
        let classes: Vec<usize> = (0..m)
            .map(|e| class(src[e], tgt[e], ls[e], lt[e]))
            .collect();
        for (c, got) in [(0, s.intralayer), (1, s.coupling), (2, s.interlayer)] {
            assert_eq!(got, classes.iter().filter(|&&x| x == c).count());
        }
        // Loops are one node in one layer; a coupling is not a loop.
        let loops = (0..m)
            .filter(|&e| src[e] == tgt[e] && ls[e] == lt[e])
            .count();
        assert_eq!(s.loops, loops);
        // Parallel edges share endpoints together with their layers.
        let key = |e: usize| {
            let (a, b) = ((src[e], ls[e]), (tgt[e], lt[e]));
            if directed || a <= b {
                (a, b)
            } else {
                (b, a)
            }
        };
        let brute: Vec<u32> = (0..m)
            .map(|e| (0..m).filter(|&g| key(g) == key(e)).count() as u32)
            .collect();
        assert_eq!(multiplicity(built.resolved.groups.as_ref().unwrap()), brute);
        // Per-layer counts and subgraphs hold the intralayer edges.
        for (l, &k) in order.iter().enumerate() {
            let rows: Vec<usize> = (0..m).filter(|&e| classes[e] == 0 && ls[e] == k).collect();
            assert_eq!(s.layer_counts[l], rows.len());
            let sub = layer_subgraph(&built, l).unwrap();
            assert_eq!(sub.edges, rows);
            assert_eq!(sub.graph.node_count(), n);
            let got: Vec<u32> = sub.graph.edge_sources().to_vec();
            let want: Vec<u32> = rows.iter().map(|&e| src[e] as u32).collect();
            assert_eq!(got, want);
        }
        assert!(layer_subgraph(&built, order.len()).is_err());
        // Declaring couplings only rejects exactly the interlayer edges.
        let mut g = f.clone();
        g.declarations_mut().couplings_only = Some(true);
        let v = validate(&g);
        let rows: Vec<usize> = (0..m).filter(|&e| classes[e] == 2).collect();
        if rows.is_empty() {
            assert!(v.report.is_valid(), "{:?}", v.report);
        } else {
            only_error(&v, Code::E_INTERLAYER_DECLARED, &rows);
        }
        // A layer absent from the table is reported on its row alone.
        if with_table && m > 0 {
            let e = rng.below(m as u64) as usize;
            let mut bad: Vec<String> = lt.iter().map(|&k| name(k)).collect();
            bad[e] = "nowhere".to_string();
            let mut g = f.clone();
            g.clear_edge_role(EdgeRole::TargetLayer);
            g.add_edge_column("lt2", text(&bad)).unwrap();
            g.set_edge_role(EdgeRole::TargetLayer, "lt2").unwrap();
            only_error(&validate(&g), Code::E_LAYER_UNKNOWN, &[e]);
        }
    }
}

#[test]
fn unobserved_pairs_are_kept_and_counted_but_not_built() {
    let mut rng = Rng(22);
    for case in 0..300 {
        let directed = rng.flag();
        let n = 1 + rng.below(15) as usize;
        let m = rng.below(40) as usize;
        let (src, tgt) = endpoints(&mut rng, n, m);
        let flags: Vec<bool> = (0..m).map(|_| rng.below(4) == 0).collect();
        let non_respondent: Vec<bool> = (0..n).map(|_| rng.below(5) == 0).collect();
        let mut f = base(directed, n, &src, &tgt);
        f.add_edge_column("na", Column::Bool(flags.clone(), None))
            .unwrap();
        f.set_edge_role(EdgeRole::EdgeMissing, "na").unwrap();
        f.add_node_column("nr", Column::Bool(non_respondent.clone(), None))
            .unwrap();
        f.set_node_role(NodeRole::NodeMissing, "nr").unwrap();
        let key = |e: usize| {
            if directed || src[e] <= tgt[e] {
                (src[e], tgt[e])
            } else {
                (tgt[e], src[e])
            }
        };
        let mut seen: HashMap<(usize, usize), (bool, bool)> = HashMap::new();
        for (e, &flag) in flags.iter().enumerate() {
            let s = seen.entry(key(e)).or_default();
            if flag {
                s.1 = true;
            } else {
                s.0 = true;
            }
        }
        let conflicts: Vec<usize> = (0..m).filter(|&e| flags[e] && seen[&key(e)].0).collect();
        let v = validate(&f);
        if !conflicts.is_empty() {
            only_error(&v, Code::E_MISSING_CONFLICT, &conflicts);
            continue;
        }
        assert!(v.report.is_valid(), "case {case}: {:?}", v.report);
        let built = build(&f, 2).unwrap();
        let observed: Vec<usize> = (0..m).filter(|&e| !flags[e]).collect();
        assert_eq!(built.edges.count, observed.len());
        assert_eq!(built.graph.edge_count(), observed.len());
        let rows: Vec<usize> = (0..built.edges.count).map(|e| built.edges.row(e)).collect();
        assert_eq!(rows, observed);
        let want: Vec<u32> = observed.iter().map(|&e| src[e] as u32).collect();
        assert_eq!(built.graph.edge_sources(), &want[..]);
        let pairs = missing_pairs(&built);
        let want: Vec<(u32, u32, usize)> = (0..m)
            .filter(|&e| flags[e])
            .map(|e| (src[e] as u32, tgt[e] as u32, e))
            .collect();
        assert_eq!(pairs, want);
        let s = summarize(&f, &built);
        assert_eq!(s.edges, m);
        assert_eq!(s.graph_edges, observed.len());
        assert_eq!(s.missing_edges, m - observed.len());
        assert_eq!(
            s.missing_nodes,
            non_respondent.iter().filter(|&&x| x).count()
        );
    }
}

#[test]
fn ego_networks_are_counted_and_split_by_ego() {
    let mut rng = Rng(23);
    for case in 0..300 {
        let directed = rng.flag();
        let n = 1 + rng.below(15) as usize;
        let ne = 1 + rng.below(n.min(4) as u64) as usize;
        // Ego k is node k (and belongs to itself) or a respondent who is
        // not a node.
        let is_node: Vec<bool> = (0..ne).map(|_| rng.flag()).collect();
        let label = |k: usize| {
            if is_node[k] {
                format!("{}", 10 + k)
            } else {
                format!("x{k}")
            }
        };
        let ego_of: Vec<usize> = (0..n)
            .map(|i| {
                if i < ne && is_node[i] {
                    i
                } else {
                    rng.below(ne as u64) as usize
                }
            })
            .collect();
        let members: Vec<Vec<usize>> = (0..ne)
            .map(|k| (0..n).filter(|&i| ego_of[i] == k).collect())
            .collect();
        let mut src = Vec::new();
        let mut tgt = Vec::new();
        for _ in 0..rng.below(30) {
            let k = rng.below(ne as u64) as usize;
            if members[k].is_empty() {
                continue;
            }
            let pick = |r: &mut Rng| members[k][r.below(members[k].len() as u64) as usize];
            src.push(pick(&mut rng));
            tgt.push(pick(&mut rng));
        }
        let m = src.len();
        // Text identifiers so ego labels and node identifiers compare.
        let mut f = Frame::new(directed);
        let ids: Vec<String> = (0..n).map(|i| format!("{}", 10 + i)).collect();
        f.add_node_column("id", text(&ids)).unwrap();
        f.set_node_role(NodeRole::Id, "id").unwrap();
        f.add_node_column(
            "ego",
            text(&ego_of.iter().map(|&k| label(k)).collect::<Vec<_>>()),
        )
        .unwrap();
        f.set_node_role(NodeRole::Ego, "ego").unwrap();
        f.add_edge_column(
            "from",
            text(&src.iter().map(|&i| ids[i].clone()).collect::<Vec<_>>()),
        )
        .unwrap();
        f.add_edge_column(
            "to",
            text(&tgt.iter().map(|&i| ids[i].clone()).collect::<Vec<_>>()),
        )
        .unwrap();
        f.set_edge_role(EdgeRole::Source, "from").unwrap();
        f.set_edge_role(EdgeRole::Target, "to").unwrap();
        let v = validate(&f);
        assert!(v.report.is_valid(), "case {case}: {:?}", v.report);
        let built = build(&f, 2).unwrap();
        let s = summarize(&f, &built);
        // Labels in order of first appearance in the node table.
        let mut order = Vec::new();
        for &k in &ego_of {
            if !order.contains(&k) {
                order.push(k);
            }
        }
        let labels: Vec<String> = order.iter().map(|&k| label(k)).collect();
        assert_eq!(s.ego_labels, labels);
        assert_eq!(
            s.egos_as_nodes,
            order.iter().filter(|&&k| is_node[k]).count()
        );
        for (c, &k) in order.iter().enumerate() {
            let rows: Vec<usize> = (0..m).filter(|&e| ego_of[src[e]] == k).collect();
            assert_eq!(s.ego_node_counts[c], members[k].len());
            assert_eq!(s.ego_edge_counts[c], rows.len());
            let sub = ego_subgraph(&built, c).unwrap();
            let nodes: Vec<u32> = members[k].iter().map(|&i| i as u32).collect();
            assert_eq!(sub.nodes, nodes);
            assert_eq!(sub.edges, rows);
            assert_eq!(sub.graph.node_count(), members[k].len());
        }
        assert!(ego_subgraph(&built, order.len()).is_err());
        // Declaring that egos are nodes rejects respondents who are not.
        let mut g = f.clone();
        g.declarations_mut().egos_are_nodes = Some(true);
        let v = validate(&g);
        if order.iter().all(|&k| is_node[k]) {
            assert!(v.report.is_valid(), "{:?}", v.report);
        } else {
            assert_eq!(error_codes(&v), vec![Code::E_EGO_UNKNOWN]);
        }
        // An edge between two egos' networks is reported on its row.
        let firsts: Vec<usize> = order.iter().map(|&k| members[k][0]).collect();
        if firsts.len() > 1 {
            let mut s2 = src.clone();
            let mut t2 = tgt.clone();
            s2.push(firsts[0]);
            t2.push(firsts[1]);
            let mut h = Frame::new(directed);
            h.add_node_column("id", text(&ids)).unwrap();
            h.set_node_role(NodeRole::Id, "id").unwrap();
            h.add_node_column(
                "ego",
                text(&ego_of.iter().map(|&k| label(k)).collect::<Vec<_>>()),
            )
            .unwrap();
            h.set_node_role(NodeRole::Ego, "ego").unwrap();
            h.add_edge_column(
                "from",
                text(&s2.iter().map(|&i| ids[i].clone()).collect::<Vec<_>>()),
            )
            .unwrap();
            h.add_edge_column(
                "to",
                text(&t2.iter().map(|&i| ids[i].clone()).collect::<Vec<_>>()),
            )
            .unwrap();
            h.set_edge_role(EdgeRole::Source, "from").unwrap();
            h.set_edge_role(EdgeRole::Target, "to").unwrap();
            only_error(&validate(&h), Code::E_EGO_MISMATCH, &[m]);
        }
    }
}

#[test]
fn node_spell_and_change_tables_are_validated() {
    let mut rng = Rng(24);
    for case in 0..300 {
        let n = 1 + rng.below(10) as usize;
        let m = rng.below(25) as usize;
        let (src, tgt) = endpoints(&mut rng, n, m);
        // Disjoint integer spells [s, e) with s < e per node (touching
        // allowed); some nodes have none and are present throughout.
        let mut spells: Vec<(usize, f64, f64)> = Vec::new();
        for i in 0..n {
            let mut t = rng.below(3) as f64;
            for _ in 0..rng.below(4) {
                let len = 1 + rng.below(4);
                spells.push((i, t, t + len as f64));
                t += (len + rng.below(3)) as f64;
            }
        }
        rng.shuffle(&mut spells);
        let present = |i: usize, t: f64| {
            let own: Vec<&(usize, f64, f64)> = spells.iter().filter(|s| s.0 == i).collect();
            own.is_empty() || own.iter().any(|s| s.1 <= t && t < s.2)
        };
        let times: Vec<f64> = (0..m).map(|_| rng.below(16) as f64).collect();
        let spell_frame = |spells: &[(usize, f64, f64)], timed: bool| {
            let mut f = base(true, n, &src, &tgt);
            if timed {
                f.add_edge_column("t", Column::Float(times.clone()))
                    .unwrap();
                f.set_edge_role(EdgeRole::Start, "t").unwrap();
            }
            f.add_column(
                Which::NodeSpells,
                "id",
                Column::Int(spells.iter().map(|s| 10 + s.0 as i64).collect(), None),
            )
            .unwrap();
            f.add_column(
                Which::NodeSpells,
                "s",
                Column::Float(spells.iter().map(|s| s.1).collect()),
            )
            .unwrap();
            f.add_column(
                Which::NodeSpells,
                "e",
                Column::Float(spells.iter().map(|s| s.2).collect()),
            )
            .unwrap();
            f.set_spell_role(SpellRole::Id, "id").unwrap();
            f.set_spell_role(SpellRole::Start, "s").unwrap();
            f.set_spell_role(SpellRole::End, "e").unwrap();
            f
        };
        // A contact must fall inside a spell of each endpoint.
        let f = spell_frame(&spells, true);
        let v = validate(&f);
        let outside: Vec<usize> = (0..m)
            .filter(|&e| !present(src[e], times[e]) || !present(tgt[e], times[e]))
            .collect();
        if !outside.is_empty() {
            only_error(&v, Code::E_SPELL_OUTSIDE_NODE, &outside);
        } else {
            assert!(v.report.is_valid(), "case {case}: {:?}", v.report);
            let built = build(&f, 2).unwrap();
            let t = built.resolved.spell_table.as_ref().unwrap();
            assert_eq!(t.offsets.len(), n + 1);
            assert_eq!(t.rows.len(), spells.len());
            for i in 0..n {
                let k = t.offsets[i]..t.offsets[i + 1];
                for j in k.clone() {
                    let row = t.rows[j];
                    assert_eq!(spells[row], (i, t.start[j], t.end[j]));
                    if j > k.start {
                        assert!(t.start[j - 1] <= t.start[j]);
                    }
                }
                assert_eq!(k.len(), spells.iter().filter(|s| s.0 == i).count());
            }
            let s = summarize(&f, &built);
            assert!(s.node_spell_table);
            assert_eq!(s.node_spell_rows, spells.len());
        }
        // An overlapping spell is an error unless declared allowed.
        if let Some(&(i, a, b)) = spells.first() {
            let mut more = spells.clone();
            more.push((i, (a + b) / 2.0, b + 1.0));
            let f = spell_frame(&more, false);
            let v = validate(&f);
            assert_eq!(error_codes(&v), vec![Code::E_NODE_SPELL_OVERLAP]);
            let rows = &v.report.find(Code::E_NODE_SPELL_OVERLAP).unwrap().rows;
            assert!(rows.contains(&spells.len()));
            let mut g = f.clone();
            g.declarations_mut().overlapping_spells = Some(true);
            assert!(validate(&g).report.is_valid());
        }
        // A spell of an unknown node is reported on its row.
        let mut more = spells.clone();
        more.push((n + 5, 0.0, 1.0));
        let v = validate(&spell_frame(&more, false));
        only_error(&v, Code::E_UNKNOWN_NODE, &[spells.len()]);
        assert_eq!(
            v.report.find(Code::E_UNKNOWN_NODE).unwrap().table,
            Some(Which::NodeSpells)
        );
        // Node changes: known identifiers and finite times.
        let k = rng.below(8) as usize;
        let mut ids: Vec<i64> = (0..k).map(|_| 10 + rng.below(n as u64) as i64).collect();
        let mut ts: Vec<f64> = (0..k).map(|_| rng.below(10) as f64).collect();
        let bad = k > 0 && rng.flag();
        let fault = if bad { rng.below(k as u64) as usize } else { 0 };
        let infinite = rng.flag();
        if bad {
            if infinite {
                ts[fault] = f64::INFINITY;
            } else {
                ids[fault] = 9;
            }
        }
        let mut f = base(true, n, &src, &tgt);
        f.add_column(Which::NodeChanges, "id", Column::Int(ids, None))
            .unwrap();
        f.add_column(Which::NodeChanges, "t", Column::Float(ts))
            .unwrap();
        f.set_change_role(ChangeRole::Id, "id").unwrap();
        f.set_change_role(ChangeRole::Time, "t").unwrap();
        let v = validate(&f);
        if !bad {
            assert!(v.report.is_valid(), "{:?}", v.report);
            let built = build(&f, 2).unwrap();
            assert_eq!(summarize(&f, &built).node_change_rows, k);
        } else if infinite {
            only_error(&v, Code::E_TIME_NOT_FINITE, &[fault]);
        } else {
            only_error(&v, Code::E_UNKNOWN_NODE, &[fault]);
        }
    }
}
