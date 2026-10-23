//! Randomized property tests for the frame: valid random frames pass
//! validation and build the same graph as the reference construction
//! (node order from the node table, then first appearance; endpoints
//! matched by identifier); a single injected fault is reported with
//! exactly its code and rows and nothing else; diagnostics do not depend
//! on the thread count; parallel-edge groups and simplified weights agree
//! with a quadratic brute force.

extern crate knappnet as engine;

use std::collections::HashMap;

use engine::frame::diagnostics::ROW_CAP;
use engine::frame::multiplicity::{is_simple, multiplicity};
use engine::frame::{
    build, simplify, summarize, Code, Column, Combine, CoordinateKind, EdgeRole, Frame, NodeRole,
    Severity, Strings, WeightType,
};
use engine::Graph;

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
}

/// A random network in identifier form, with every optional column.
struct Case {
    directed: bool,
    text_ids: bool,
    has_nodes: bool,
    /// Node identifiers as integers (rendered as text when `text_ids`).
    ids: Vec<i64>,
    src: Vec<i64>,
    tgt: Vec<i64>,
    weights: Vec<f64>,
    relations: Vec<i64>,
    modes: Vec<i64>,
    start: Vec<f64>,
    end: Vec<f64>,
    node_start: Vec<f64>,
    node_end: Vec<f64>,
}

fn random_case(rng: &mut Rng) -> Case {
    let n = 1 + rng.below(25) as usize;
    let m = rng.below((n * n / 2 + 2) as u64) as usize;
    // Identifiers are distinct, not contiguous, and not 0-based.
    let base = 1 + rng.below(100) as i64;
    let ids: Vec<i64> = (0..n).map(|i| base + 3 * i as i64).collect();
    let mut src = Vec::new();
    let mut tgt = Vec::new();
    let mut weights = Vec::new();
    let mut relations = Vec::new();
    let mut start = Vec::new();
    let mut end = Vec::new();
    let mut node_start = Vec::new();
    let mut node_end = Vec::new();
    for _ in 0..n {
        let a = rng.below(20) as f64;
        node_start.push(a);
        node_end.push(a + 20.0 + rng.below(20) as f64);
    }
    for _ in 0..m {
        let a = rng.below(n as u64) as usize;
        let b = rng.below(n as u64) as usize;
        src.push(ids[a]);
        tgt.push(ids[b]);
        weights.push(1.0 + rng.below(5) as f64);
        relations.push(rng.below(3) as i64);
        // A spell inside both endpoints' presence, [lo, hi) with hi > lo.
        let lo = node_start[a].max(node_start[b]);
        let hi = node_end[a].min(node_end[b]);
        let s = lo + rng.below((hi - lo) as u64) as f64;
        let e = (s + rng.below(3) as f64).min(hi);
        start.push(s);
        end.push(e);
    }
    Case {
        directed: rng.flag(),
        text_ids: rng.flag(),
        has_nodes: rng.flag(),
        ids,
        src,
        tgt,
        weights,
        relations,
        modes: (0..n).map(|_| rng.below(2) as i64).collect(),
        start,
        end,
        node_start,
        node_end,
    }
}

fn id_column(values: &[i64], text: bool) -> Column {
    if text {
        let v: Vec<Option<String>> = values.iter().map(|x| Some(format!("n{x}"))).collect();
        Column::Str(Strings::from_strs(&v))
    } else {
        Column::Int(values.to_vec(), None)
    }
}

/// The frame of a case, with the given optional roles.
struct Opts {
    weighted: bool,
    relations: bool,
    modes: bool,
    temporal: bool,
    node_spells: bool,
}

fn frame(c: &Case, o: &Opts) -> Frame {
    let mut f = Frame::new(c.directed);
    if c.has_nodes {
        f.add_node_column("id", id_column(&c.ids, c.text_ids))
            .unwrap();
        f.set_node_role(NodeRole::Id, "id").unwrap();
        if o.modes {
            f.add_node_column("mode", Column::Int(c.modes.clone(), None))
                .unwrap();
            f.set_node_role(NodeRole::Mode, "mode").unwrap();
        }
        if o.node_spells {
            f.add_node_column("arrive", Column::Float(c.node_start.clone()))
                .unwrap();
            f.add_node_column("leave", Column::Float(c.node_end.clone()))
                .unwrap();
            f.set_node_role(NodeRole::NodeStart, "arrive").unwrap();
            f.set_node_role(NodeRole::NodeEnd, "leave").unwrap();
        }
    }
    f.add_edge_column("from", id_column(&c.src, c.text_ids))
        .unwrap();
    f.add_edge_column("to", id_column(&c.tgt, c.text_ids))
        .unwrap();
    f.set_edge_role(EdgeRole::Source, "from").unwrap();
    f.set_edge_role(EdgeRole::Target, "to").unwrap();
    if o.weighted {
        f.add_edge_column("w", Column::Float(c.weights.clone()))
            .unwrap();
        f.set_edge_role(EdgeRole::Weight, "w").unwrap();
        f.set_weight_type(Some(WeightType::Distance));
    }
    if o.relations {
        f.add_edge_column("rel", Column::Int(c.relations.clone(), None))
            .unwrap();
        f.set_edge_role(EdgeRole::Relation, "rel").unwrap();
    }
    if o.temporal {
        f.add_edge_column("s", Column::Float(c.start.clone()))
            .unwrap();
        f.add_edge_column("e", Column::Float(c.end.clone()))
            .unwrap();
        f.set_edge_role(EdgeRole::Start, "s").unwrap();
        f.set_edge_role(EdgeRole::End, "e").unwrap();
    }
    f
}

fn random_opts(rng: &mut Rng, c: &Case) -> Opts {
    Opts {
        weighted: rng.flag(),
        relations: rng.flag(),
        modes: c.has_nodes && rng.flag(),
        temporal: rng.flag(),
        node_spells: c.has_nodes && rng.flag(),
    }
}

/// The reference construction: nodes in table order (or first appearance
/// without a table), endpoints matched by identifier.
fn reference(c: &Case, weighted: bool) -> (Graph, Vec<i64>) {
    let mut order: Vec<i64> = if c.has_nodes {
        c.ids.clone()
    } else {
        Vec::new()
    };
    let mut pos: HashMap<i64, i64> = order
        .iter()
        .enumerate()
        .map(|(i, &x)| (x, i as i64))
        .collect();
    let mut src = Vec::new();
    let mut tgt = Vec::new();
    for e in 0..c.src.len() {
        for (which, x) in [(&mut src, c.src[e]), (&mut tgt, c.tgt[e])] {
            let i = *pos.entry(x).or_insert_with(|| {
                order.push(x);
                order.len() as i64 - 1
            });
            which.push(i);
        }
    }
    let w = weighted.then_some(&c.weights[..]);
    (
        Graph::from_edges(order.len(), &src, &tgt, w, c.directed).unwrap(),
        order,
    )
}

/// Brute-force parallel-edge key comparison.
fn brute_multiplicity(c: &Case, o: &Opts) -> Vec<u32> {
    let m = c.src.len();
    let key = |e: usize| {
        let (a, b) = if c.directed || c.src[e] <= c.tgt[e] {
            (c.src[e], c.tgt[e])
        } else {
            (c.tgt[e], c.src[e])
        };
        let r = if o.relations { c.relations[e] } else { 0 };
        let (s, t) = if o.temporal {
            (c.start[e].to_bits(), c.end[e].to_bits())
        } else {
            (0, 0)
        };
        (a, b, r, s, t)
    };
    (0..m)
        .map(|e| (0..m).filter(|&f| key(f) == key(e)).count() as u32)
        .collect()
}

fn errors(f: &Frame) -> Vec<Code> {
    f.validate(1)
        .report
        .diagnostics
        .iter()
        .filter(|d| d.severity() == Severity::Error)
        .map(|d| d.code)
        .collect()
}

#[test]
fn valid_frames_build_the_reference_graph() {
    let mut rng = Rng(11);
    for _ in 0..300 {
        let c = random_case(&mut rng);
        let o = random_opts(&mut rng, &c);
        let f = frame(&c, &o);
        let v1 = f.validate(1);
        let v4 = f.validate(4);
        assert_eq!(v1.report, v4.report, "diagnostics depend on threads");
        assert!(v1.report.is_valid(), "{:?}", v1.report);
        let built = build(&f, 2).unwrap();
        let (expected, order) = reference(&c, o.weighted);
        assert_eq!(built.graph.node_count(), expected.node_count());
        assert_eq!(built.graph.edge_sources(), expected.edge_sources());
        assert_eq!(built.graph.edge_targets(), expected.edge_targets());
        assert_eq!(built.graph.edge_weights(), expected.edge_weights());
        assert_eq!(built.graph.is_directed(), expected.is_directed());
        // The added identifiers come back in node order.
        let added: Vec<i64> = order[built.map.table_rows..].to_vec();
        for (k, id) in added.iter().enumerate() {
            let rendered = built.map.added_ids.render(k);
            assert_eq!(
                rendered,
                if c.text_ids {
                    format!("n{id}")
                } else {
                    id.to_string()
                }
            );
        }
        // Notes agree with the data.
        let r = &built.resolved;
        let loops = (0..c.src.len()).filter(|&e| c.src[e] == c.tgt[e]).count();
        assert_eq!(r.loops, loops);
        assert_eq!(built.report.has(Code::N_LOOPS), loops > 0);
        let mult = multiplicity(r.groups.as_ref().unwrap());
        assert_eq!(mult, brute_multiplicity(&c, &o));
        let parallel = mult.iter().filter(|&&k| k > 1).count() > 0;
        assert_eq!(built.report.has(Code::N_PARALLEL), parallel);
        assert_eq!(is_simple(r), loops == 0 && !parallel);
        let s = summarize(&f, &built);
        assert_eq!(s.loops, loops);
        assert_eq!(s.edges, c.src.len());
        assert_eq!(s.nodes, order.len());
        let isolates = (0..order.len())
            .filter(|&i| !c.src.contains(&order[i]) && !c.tgt.contains(&order[i]))
            .count();
        assert_eq!(s.isolates, isolates);
        assert_eq!(built.report.has(Code::N_ISOLATES), isolates > 0);
        if c.has_nodes {
            assert!(!built.report.has(Code::N_NODES_ADDED) || f.declarations().add_nodes);
        }
    }
}

#[test]
fn simplify_matches_brute_force() {
    let mut rng = Rng(23);
    for _ in 0..200 {
        let c = random_case(&mut rng);
        let o = Opts {
            weighted: true,
            relations: rng.flag(),
            modes: false,
            temporal: rng.flag(),
            node_spells: false,
        };
        let f = frame(&c, &o);
        let r = f.validate(1).resolved.unwrap();
        let mult = brute_multiplicity(&c, &o);
        for (combine, drop_loops) in [
            (Combine::Sum, false),
            (Combine::Mean, true),
            (Combine::Min, false),
            (Combine::Max, true),
            (Combine::First, false),
            (Combine::Count, false),
        ] {
            let s = simplify(&r, combine, drop_loops);
            let m = c.src.len();
            // Every edge maps to the group of its first row; loops drop
            // when asked.
            let mut seen = vec![false; s.keep.len()];
            for (e, &mult_e) in mult.iter().enumerate() {
                let is_loop = c.src[e] == c.tgt[e];
                match s.group[e] {
                    None => assert!(drop_loops && is_loop),
                    Some(g) => {
                        assert!(!(drop_loops && is_loop));
                        let first = s.keep[g as usize];
                        assert!(first <= e);
                        seen[g as usize] = true;
                        // Same key as the first row of its group.
                        let members: Vec<usize> =
                            (0..m).filter(|&x| s.group[x] == Some(g)).collect();
                        assert_eq!(members.len() as u32, mult_e);
                        let w = s.weight.as_ref().unwrap()[g as usize];
                        let ws: Vec<f64> = members.iter().map(|&x| c.weights[x]).collect();
                        let expect = match combine {
                            Combine::Sum => ws.iter().sum(),
                            Combine::Mean => ws.iter().sum::<f64>() / ws.len() as f64,
                            Combine::Min => ws.iter().cloned().fold(f64::INFINITY, f64::min),
                            Combine::Max => ws.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
                            Combine::First => ws[0],
                            Combine::Count => ws.len() as f64,
                        };
                        assert!((w - expect).abs() < 1e-9);
                    }
                }
            }
            assert!(seen.iter().all(|&x| x));
        }
    }
}

/// One fault per frame: the expected code, with exactly the injected
/// rows, and no other error.
#[test]
fn injected_faults_are_reported_exactly() {
    let mut rng = Rng(37);
    let mut injected = 0;
    for _ in 0..400 {
        let mut c = random_case(&mut rng);
        c.has_nodes = true;
        let m = c.src.len();
        let n = c.ids.len();
        if m == 0 || n < 2 {
            continue;
        }
        let o = Opts {
            weighted: true,
            relations: true,
            modes: true,
            temporal: true,
            node_spells: true,
        };
        let row = rng.below(m as u64) as usize;
        let node = rng.below(n as u64) as usize;
        let fault = rng.below(14);
        let (code, rows): (Code, Vec<usize>) = match fault {
            0 => {
                // An endpoint absent from the node table.
                c.src[row] = c.ids[n - 1] + 1;
                (Code::E_UNKNOWN_NODE, vec![row])
            }
            1 => {
                c.ids[node] = c.ids[(node + 1) % n];
                (Code::E_DUPLICATE_ID, vec![node.max((node + 1) % n)])
            }
            2 => {
                c.weights[row] = f64::NAN;
                (Code::E_WEIGHT_NOT_FINITE, vec![row])
            }
            3 => {
                c.weights[row] = 0.0;
                (Code::E_DISTANCE_NOT_POSITIVE, vec![row])
            }
            4 => {
                c.end[row] = c.start[row] - 1.0;
                (Code::E_SPELL_ORDER, vec![row])
            }
            5 => {
                // Leave before arriving: also breaks every edge of that
                // node, so those rows are reported too.
                c.node_end[node] = c.node_start[node] - 1.0;
                (Code::E_NODE_SPELL_ORDER, vec![node])
            }
            6 => {
                c.start[row] = f64::NAN;
                (Code::E_TIME_MISSING, vec![row])
            }
            7 => {
                c.start[row] = f64::INFINITY;
                (Code::E_TIME_NOT_FINITE, vec![row])
            }
            8 => {
                // An edge starting before its source arrives.
                let a = c.ids.iter().position(|&x| x == c.src[row]).unwrap();
                c.start[row] = c.node_start[a] - 1.0;
                c.end[row] = c.start[row];
                (Code::E_SPELL_OUTSIDE_NODE, vec![row])
            }
            _ => {
                // Faults that change the frame rather than the data are
                // injected below.
                (Code::E_DECLARATION_CONFLICT, vec![])
            }
        };
        let mut f = frame(&c, &o);
        let (code, rows) = match fault {
            9 => {
                // Mode missing for one node.
                let mask: Vec<bool> = (0..n).map(|i| i == node).collect();
                f.add_node_column("mode2", Column::Int(c.modes.clone(), Some(mask)))
                    .unwrap();
                f.set_node_role(NodeRole::Mode, "mode2").unwrap();
                (Code::E_MODE_MISSING, vec![node])
            }
            10 => {
                let mask: Vec<bool> = (0..m).map(|i| i == row).collect();
                f.add_edge_column("rel2", Column::Int(c.relations.clone(), Some(mask)))
                    .unwrap();
                f.set_edge_role(EdgeRole::Relation, "rel2").unwrap();
                (Code::E_RELATION_MISSING, vec![row])
            }
            11 => {
                let mut sign = vec![1i64; m];
                sign[row] = 2;
                f.add_edge_column("sign", Column::Int(sign, None)).unwrap();
                f.set_edge_role(EdgeRole::Sign, "sign").unwrap();
                (Code::E_SIGN_VALUE, vec![row])
            }
            12 => {
                let mut y = vec![0.0; n];
                y[node] = 95.0;
                f.add_node_column("x", Column::Float(vec![0.0; n])).unwrap();
                f.add_node_column("y", Column::Float(y)).unwrap();
                f.set_node_role(NodeRole::X, "x").unwrap();
                f.set_node_role(NodeRole::Y, "y").unwrap();
                f.geometry_mut().coordinates = CoordinateKind::LonLat;
                (Code::E_LATITUDE_RANGE, vec![node])
            }
            13 => {
                f.set_weight_type(None);
                (Code::E_WEIGHT_TYPE_MISSING, vec![])
            }
            _ => (code, rows),
        };
        let v = f.validate(2);
        let found: Vec<&Code> = v.report.errors().map(|d| &d.code).collect();
        assert!(
            found.contains(&&code),
            "fault {fault}: expected {code}, got {found:?}"
        );
        let d = v.report.find(code).unwrap();
        if !rows.is_empty() {
            assert_eq!(d.rows, rows, "fault {fault}");
            assert_eq!(d.count, rows.len());
        }
        // The only other errors allowed are knock-on effects: a duplicated
        // identifier replaces one the edges still use, a broken node
        // spell breaks that node's edges, and a missing start time fails
        // the spell order check.
        for other in &found {
            let allowed = **other == code
                || (fault == 1 && **other == Code::E_UNKNOWN_NODE)
                || (fault == 5 && **other == Code::E_SPELL_OUTSIDE_NODE)
                || (fault == 6 && **other == Code::E_SPELL_ORDER);
            assert!(
                allowed,
                "fault {fault}: unexpected {other} alongside {code}"
            );
        }
        assert_eq!(
            errors(&frame(&c, &o)).is_empty(),
            fault >= 9,
            "fault {fault}"
        );
        injected += 1;
    }
    assert!(injected > 300);
}

#[test]
fn declared_properties_are_enforced() {
    let mut rng = Rng(41);
    for _ in 0..200 {
        let c = random_case(&mut rng);
        let o = Opts {
            weighted: false,
            relations: rng.flag(),
            modes: c.has_nodes,
            temporal: false,
            node_spells: false,
        };
        let mut f = frame(&c, &o);
        let m = c.src.len();
        let loops: Vec<usize> = (0..m).filter(|&e| c.src[e] == c.tgt[e]).collect();
        let mult = brute_multiplicity(&c, &o);
        let r = f.validate(1).resolved.unwrap();
        let first = &r.groups.as_ref().unwrap().first;
        let parallel: Vec<usize> = (0..m)
            .filter(|&e| mult[e] > 1 && !first.contains(&e))
            .collect();
        f.declarations_mut().loops = Some(false);
        f.declarations_mut().multiple = Some(false);
        if c.has_nodes {
            f.declarations_mut().between_modes_only = Some(true);
        }
        let v = f.validate(1);
        let check = |code: Code, rows: &[usize]| match v.report.find(code) {
            Some(d) => {
                assert_eq!(d.rows, &rows[..rows.len().min(ROW_CAP)]);
                assert_eq!(d.count, rows.len());
            }
            None => assert!(rows.is_empty(), "{code} missing for {rows:?}"),
        };
        check(Code::E_LOOP_DECLARED, &loops);
        check(Code::E_PARALLEL_DECLARED, &parallel);
        if c.has_nodes {
            let within: Vec<usize> = (0..m)
                .filter(|&e| {
                    let a = c.ids.iter().position(|&x| x == c.src[e]).unwrap();
                    let b = c.ids.iter().position(|&x| x == c.tgt[e]).unwrap();
                    c.modes[a] == c.modes[b]
                })
                .collect();
            check(Code::E_MODE_WITHIN, &within);
        }
        assert!(!v.report.has(Code::N_LOOPS) && !v.report.has(Code::N_PARALLEL));
        // Allowing them explicitly removes both the errors and the notes.
        f.declarations_mut().loops = Some(true);
        f.declarations_mut().multiple = Some(true);
        f.declarations_mut().between_modes_only = None;
        let v = f.validate(1);
        assert!(v.report.is_valid());
        assert!(!v.report.has(Code::N_LOOPS) && !v.report.has(Code::N_PARALLEL));
    }
}
