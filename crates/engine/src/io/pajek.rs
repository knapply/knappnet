//! Pajek `.net` files (Batagelj and Mrvar 1998, `batagelj1998pajek`;
//! de Nooy, Mrvar, and Batagelj 2018, `denooy2018exploratory`).
//!
//! Supported sections: `*Network`, `*Vertices n` (and two-mode
//! `*Vertices n n1`), `*Arcs`, `*Edges`, `*Arcslist`, `*Edgeslist`, and
//! `*Matrix`. Vertex lines may carry a label and x, y, z coordinates; other
//! drawing parameters are ignored. Lines starting with `%` are comments.
//!
//! When every vertex has a unique label, labels become node identifiers;
//! otherwise vertex numbers are used. A file with both `*Arcs` and
//! `*Edges` is read as directed, with each undirected edge stored as two
//! opposite arcs. Two-mode networks get a `mode` node column (1 or 2).

use super::table::{format_float, parse_float, Network, NetworkBuilder, Value};
use crate::error::{Error, Result};

const FORMAT: &str = "Pajek";

/// Splits a line into tokens, keeping quoted strings together.
fn tokens(line: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut chars = line.chars().peekable();
    while let Some(&c) = chars.peek() {
        if c.is_whitespace() {
            chars.next();
        } else if c == '"' {
            chars.next();
            let mut s = String::new();
            // Pajek has no escape syntax, but some writers (NetworkX among
            // them) emit \" for a quote inside a label; accept it.
            while let Some(c) = chars.next() {
                if c == '\\' && chars.peek() == Some(&'"') {
                    s.push('"');
                    chars.next();
                } else if c == '"' {
                    break;
                } else {
                    s.push(c);
                }
            }
            out.push(s);
        } else {
            let mut s = String::new();
            while let Some(&c) = chars.peek() {
                if c.is_whitespace() {
                    break;
                }
                s.push(c);
                chars.next();
            }
            out.push(s);
        }
    }
    out
}

#[derive(PartialEq, Clone, Copy)]
enum Section {
    None,
    Vertices,
    Arcs,
    Edges,
    ArcsList,
    EdgesList,
    Matrix,
}

/// Reads a Pajek network.
pub fn read(text: &str) -> Result<Network> {
    let err = |line: usize, message: String| Error::Parse {
        format: FORMAT,
        line: Some(line),
        message,
    };
    let mut n = 0usize;
    let mut n1: Option<usize> = None;
    let mut labels: Vec<Option<String>> = Vec::new();
    let mut coords: Vec<[Option<f64>; 3]> = Vec::new();
    // (from, to, weight, is_edge)
    let mut links: Vec<(usize, usize, Option<f64>, bool)> = Vec::new();
    let mut has_arcs = false;
    let mut has_edges = false;
    let mut section = Section::None;
    let mut matrix_row = 0usize;
    let mut name: Option<String> = None;
    let mut seen_vertices = false;

    for (i, raw) in text.lines().enumerate() {
        let line_no = i + 1;
        let line = raw.trim();
        if line.is_empty() || line.starts_with('%') {
            continue;
        }
        if let Some(rest) = line.strip_prefix('*') {
            let toks = tokens(rest);
            let head = toks
                .first()
                .map(|t| t.to_ascii_lowercase())
                .unwrap_or_default();
            section = match head.as_str() {
                "network" => {
                    name = Some(toks[1..].join(" "));
                    Section::None
                }
                "vertices" => {
                    if seen_vertices {
                        // Pajek project files hold several networks, each
                        // with its own *Vertices; links already read would
                        // refer to the wrong one.
                        return Err(err(
                            line_no,
                            "a second *Vertices section; only one network per file is supported"
                                .into(),
                        ));
                    }
                    seen_vertices = true;
                    n = toks
                        .get(1)
                        .and_then(|t| t.parse().ok())
                        .ok_or_else(|| err(line_no, "*Vertices needs a count".into()))?;
                    if n > super::MAX_DECLARED_NODES {
                        return Err(err(
                            line_no,
                            format!(
                                "*Vertices {n} declares more than the {} vertices supported",
                                super::MAX_DECLARED_NODES
                            ),
                        ));
                    }
                    n1 = toks.get(2).and_then(|t| t.parse().ok());
                    if n1.is_some_and(|k| k > n) {
                        return Err(err(
                            line_no,
                            format!("*Vertices {n} cannot have more than {n} vertices in its first mode"),
                        ));
                    }
                    labels = vec![None; n];
                    coords = vec![[None; 3]; n];
                    Section::Vertices
                }
                "arcs" => Section::Arcs,
                "edges" => Section::Edges,
                "arcslist" => Section::ArcsList,
                "edgeslist" => Section::EdgesList,
                "matrix" => {
                    matrix_row = 0;
                    Section::Matrix
                }
                other => return Err(err(line_no, format!("unsupported section *{other}"))),
            };
            match section {
                Section::Arcs | Section::ArcsList | Section::Matrix => has_arcs = true,
                Section::Edges | Section::EdgesList => has_edges = true,
                _ => {}
            }
            continue;
        }
        let toks = tokens(line);
        let vertex = |t: &str| -> Result<usize> {
            match t.parse::<usize>() {
                Ok(v) if v >= 1 && v <= n => Ok(v - 1),
                _ => Err(err(
                    line_no,
                    format!("{t:?} is not a vertex number between 1 and {n}"),
                )),
            }
        };
        match section {
            Section::None => return Err(err(line_no, "data before any section".into())),
            Section::Vertices => {
                let v = vertex(&toks[0])?;
                if let Some(l) = toks.get(1) {
                    labels[v] = Some(l.clone());
                }
                for (k, t) in toks.iter().skip(2).take(3).enumerate() {
                    match parse_float(t) {
                        Some(x) => coords[v][k] = Some(x),
                        None => break,
                    }
                }
            }
            Section::Arcs | Section::Edges => {
                if toks.len() < 2 {
                    return Err(err(line_no, "a link needs two vertices".into()));
                }
                let (a, b) = (vertex(&toks[0])?, vertex(&toks[1])?);
                let w = match toks.get(2) {
                    Some(t) => Some(
                        parse_float(t)
                            .ok_or_else(|| err(line_no, format!("invalid weight {t:?}")))?,
                    ),
                    None => None,
                };
                links.push((a, b, w, section == Section::Edges));
            }
            Section::ArcsList | Section::EdgesList => {
                let a = vertex(&toks[0])?;
                for t in &toks[1..] {
                    links.push((a, vertex(t)?, None, section == Section::EdgesList));
                }
            }
            Section::Matrix => {
                // Two-mode matrices are n1 x (n - n1); one-mode are n x n.
                let (rows, cols, col_offset) = match n1 {
                    Some(k) => (k, n - k, k),
                    None => (n, n, 0),
                };
                if matrix_row >= rows {
                    return Err(err(line_no, "too many matrix rows".into()));
                }
                if toks.len() > cols {
                    return Err(err(
                        line_no,
                        format!(
                            "a matrix row has {} values, but there are {cols} columns",
                            toks.len()
                        ),
                    ));
                }
                let row = matrix_row;
                for (j, t) in toks.iter().enumerate() {
                    let x = parse_float(t)
                        .ok_or_else(|| err(line_no, format!("invalid value {t:?}")))?;
                    if x != 0.0 {
                        links.push((row, col_offset + j, Some(x), false));
                    }
                }
                matrix_row += 1;
            }
        }
    }

    let directed = has_arcs || !has_edges;
    let use_labels = {
        let set: std::collections::HashSet<_> = labels.iter().flatten().collect();
        labels.iter().all(Option::is_some) && set.len() == n && n > 0
    };
    let ids: Vec<String> = (0..n)
        .map(|v| {
            if use_labels {
                labels[v].clone().unwrap()
            } else {
                (v + 1).to_string()
            }
        })
        .collect();

    let mut b = NetworkBuilder::new(directed);
    if let Some(name) = name.filter(|s| !s.is_empty()) {
        b.graph.push(("name".into(), Value::Str(name)));
    }
    let has_coords = [0, 1, 2].map(|k| coords.iter().any(|c| c[k].is_some()));
    let label_column = !use_labels && labels.iter().any(Option::is_some);
    for v in 0..n {
        let mut attrs = Vec::new();
        if label_column {
            attrs.push((
                "label".into(),
                labels[v].clone().map(Value::Str).unwrap_or(Value::Null),
            ));
        }
        for (k, axis) in ["x", "y", "z"].iter().enumerate() {
            if has_coords[k] {
                attrs.push((
                    axis.to_string(),
                    coords[v][k].map(Value::Float).unwrap_or(Value::Null),
                ));
            }
        }
        if let Some(k) = n1 {
            attrs.push(("mode".into(), Value::Int(if v < k { 1 } else { 2 })));
        }
        b.add_node(&ids[v], attrs);
    }
    let weighted = links.iter().any(|l| l.2.is_some());
    for (a, c, w, is_edge) in links {
        let attrs = if weighted {
            vec![("weight".to_string(), Value::Float(w.unwrap_or(1.0)))]
        } else {
            Vec::new()
        };
        b.add_edge(&ids[a], &ids[c], attrs.clone());
        if is_edge && directed && a != c {
            b.add_edge(&ids[c], &ids[a], attrs);
        }
    }
    b.finish()
}

/// Pajek labels cannot contain double quotes; they are written as single
/// quotes.
fn quote(s: &str) -> String {
    format!("\"{}\"", s.replace('"', "'"))
}

/// Writes a Pajek network. Vertex labels are the node identifiers, `x`,
/// `y`, and `z` columns become coordinates, and the `weight` edge column
/// (if numeric) becomes the link value.
pub fn write(net: &Network, weight: Option<&str>) -> Result<String> {
    net.validate()?;
    let ids = net.nodes.strings("id")?;
    let index = net.node_index()?;
    let mut out = format!("*Vertices {}\n", ids.len());
    let coord = |axis: &str, row: usize| match net.nodes.column(axis).map(|c| c.get(row)) {
        Some(Value::Float(f)) => Some(f),
        Some(Value::Int(i)) => Some(i as f64),
        _ => None,
    };
    for (row, id) in ids.iter().enumerate() {
        out.push_str(&format!("{} {}", row + 1, quote(id)));
        if let (Some(x), Some(y)) = (coord("x", row), coord("y", row)) {
            out.push_str(&format!(" {} {}", format_float(x), format_float(y)));
            if let Some(z) = coord("z", row) {
                out.push_str(&format!(" {}", format_float(z)));
            }
        }
        out.push('\n');
    }
    out.push_str(if net.directed { "*Arcs\n" } else { "*Edges\n" });
    let from = net.edges.strings("from")?;
    let to = net.edges.strings("to")?;
    let w = weight.or(Some("weight")).and_then(|n| net.edges.column(n));
    for row in 0..from.len() {
        out.push_str(&format!(
            "{} {}",
            index[&from[row]] + 1,
            index[&to[row]] + 1
        ));
        match w.map(|c| c.get(row)) {
            Some(Value::Float(f)) => out.push_str(&format!(" {}", format_float(f))),
            Some(Value::Int(i)) => out.push_str(&format!(" {i}")),
            _ => {}
        }
        out.push('\n');
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::table::Column;

    #[test]
    fn reads_labels_coordinates_and_weights() {
        let doc = "% comment\n*Network demo\n*Vertices 3\n1 \"Ann Lee\" 0.1 0.2\n2 \"Bob\" 0.3 0.4\n3 Cy\n*Arcs\n1 2 2.5\n*Edges\n2 3\n";
        let n = read(doc).unwrap();
        assert!(n.directed);
        assert_eq!(n.nodes.strings("id").unwrap(), vec!["Ann Lee", "Bob", "Cy"]);
        assert_eq!(
            n.nodes.column("x").unwrap(),
            &Column::Float(vec![Some(0.1), Some(0.3), None])
        );
        // The undirected edge becomes two arcs in a mixed file.
        assert_eq!(n.edges.nrow(), 3);
        assert_eq!(
            n.edges.column("weight").unwrap(),
            &Column::Float(vec![Some(2.5), Some(1.0), Some(1.0)])
        );
    }

    #[test]
    fn falls_back_to_numbers_and_reads_lists_and_matrices() {
        let q = read("*Vertices 1\n1 \"say \\\"hi\\\"\"\n").unwrap();
        assert_eq!(q.nodes.strings("id").unwrap(), vec!["say \"hi\""]);
        let n = read("*Vertices 3\n1 \"a\"\n2 \"a\"\n*Edgeslist\n1 2 3\n").unwrap();
        assert!(!n.directed);
        assert_eq!(n.nodes.strings("id").unwrap(), vec!["1", "2", "3"]);
        assert_eq!(n.edges.nrow(), 2);
        let m = read("*Vertices 2\n*Matrix\n0 3\n1 0\n").unwrap();
        assert_eq!(m.edges.strings("from").unwrap(), vec!["1", "2"]);
        assert!(read("*Vertices 2\n*Arcs\n1 5\n")
            .unwrap_err()
            .to_string()
            .contains("line 3"));
    }

    #[test]
    fn round_trips() {
        let doc = "*Vertices 3\n1 \"a\" 0.0 1.0\n2 \"b\" 1.0 0.0\n3 \"c\" 2.0 2.0\n*Edges\n1 2 2.0\n2 3 0.5\n";
        let n = read(doc).unwrap();
        let again = read(&write(&n, None).unwrap()).unwrap();
        assert_eq!(again, n);
    }
}
