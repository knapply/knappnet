//! UCINET DL files (Borgatti, Everett, and Freeman 2002,
//! `borgatti2002ucinet`).
//!
//! Supported: one-mode data (`n=`) in `fullmatrix` (the default),
//! `edgelist1`, and `nodelist1` formats, and two-mode data (`nr=` and
//! `nc=`) in `fullmatrix` (an nr x nc matrix), `edgelist2`, and
//! `nodelist2` formats; node labels given with `labels:`, `row labels:`,
//! and `column labels:` (or `col labels:`), or embedded in the data
//! (`labels embedded`). Keywords are case-insensitive and may be separated
//! by spaces, commas, or newlines. Multi-matrix and half-matrix files are
//! not yet supported and produce a clear error.
//!
//! # Two-mode data
//!
//! Rows and columns are different sets of nodes (people and events, say).
//! The node table holds the row nodes first, then the column nodes, with a
//! logical `type` column that is false for rows and true for columns (the
//! convention of igraph's bipartite graphs), and every edge joins a row to
//! a column. A matrix cell, an `edgelist2` line `row column [value]`, or a
//! `nodelist2` line `row column column ...` is one edge from the row to
//! the column. Identifiers come from the labels of each mode as for
//! one-mode data; unlabelled rows are numbered 1 to nr and unlabelled
//! columns nr + 1 to nr + nc, as in Pajek's two-mode files. Only when the
//! same identifier would appear in both modes are all identifiers prefixed
//! with `row::` and `col::`, with the unprefixed ones kept in a `label`
//! column. Two-mode data are read as undirected unless the caller asks
//! for directed edges (rows to columns). Whether igraph reads two-mode DL
//! files has not been checked; two-mode reading is tested against the
//! layouts in UCINET's documentation of the DL format.
//!
//! Labels become node identifiers when they are unique. When they are not
//! (UCINET does not require it), nodes are numbered from 1 and the labels
//! are kept in a `label` column, as for Pajek files; reading repeated
//! labels as identifiers would merge distinct nodes.
//!
//! DL files do not record direction. One-mode networks are read as
//! directed unless the caller says otherwise; when read as undirected, a
//! full matrix contributes only its upper triangle (diagonal included).
//!
//! # Writing
//!
//! One-mode networks are written in `edgelist1` format with node labels,
//! so isolates and labels survive. A network whose node table has a
//! logical `type` column without missing values, and whose every edge
//! joins a node of each type, is written as two-mode data: `nr=` and `nc=`
//! with row and column labels and an `edgelist2` data section (an edge
//! list keeps parallel edges and stays small for sparse data), rows being
//! the nodes whose `type` is false.

use super::table::{format_float, parse_float, Network, NetworkBuilder, Value};
use crate::error::{Error, Result};

const FORMAT: &str = "UCINET DL";

fn err<T>(message: impl Into<String>) -> Result<T> {
    Err(Error::Parse {
        format: FORMAT,
        line: None,
        message: message.into(),
    })
}

/// Splits on whitespace and commas, keeping quoted strings together.
fn words(text: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut quoted = false;
    for c in text.chars() {
        if quoted {
            if c == '"' {
                quoted = false;
                out.push(std::mem::take(&mut cur));
            } else {
                cur.push(c);
            }
        } else if c == '"' {
            quoted = true;
        } else if c.is_whitespace() || c == ',' {
            if !cur.is_empty() {
                out.push(std::mem::take(&mut cur));
            }
        } else {
            cur.push(c);
        }
    }
    if !cur.is_empty() {
        out.push(cur);
    }
    out
}

fn count(key: &str, v: &str) -> Result<usize> {
    let count: usize = v.parse().map_err(|_| Error::Parse {
        format: FORMAT,
        line: None,
        message: format!("invalid {key} {v:?}"),
    })?;
    if count > super::MAX_DECLARED_NODES {
        return err(format!(
            "{key}={count} declares more than the {} nodes supported",
            super::MAX_DECLARED_NODES
        ));
    }
    Ok(count)
}

fn value(t: &str) -> Result<f64> {
    parse_float(t).ok_or_else(|| Error::Parse {
        format: FORMAT,
        line: None,
        message: format!("invalid value {t:?}"),
    })
}

/// Reads a DL file. `directed` defaults to true for one-mode data and to
/// false for two-mode data.
pub fn read(text: &str, directed_opt: Option<bool>) -> Result<Network> {
    let directed = directed_opt.unwrap_or(true);
    let lower = text.to_ascii_lowercase();
    let data_at = match lower.find("data:") {
        Some(i) => i,
        None => return err("missing \"data:\" section"),
    };
    // Normalize "n = 5" to "n=5" before tokenizing the header.
    let header: String = text[..data_at]
        .split('=')
        .map(str::trim)
        .collect::<Vec<_>>()
        .join("=");
    let toks = words(&header);
    if toks.first().map(|t| t.to_ascii_lowercase()) != Some("dl".into()) {
        return err("file must start with \"dl\"");
    }
    let mut n: Option<usize> = None;
    let (mut nr, mut nc): (Option<usize>, Option<usize>) = (None, None);
    let mut format = "fullmatrix".to_string();
    let mut labels: Vec<String> = Vec::new();
    let mut row_labels: Vec<String> = Vec::new();
    let mut col_labels: Vec<String> = Vec::new();
    let mut embedded = false;
    let is_mode = |t: &str| matches!(t.to_ascii_lowercase().as_str(), "row" | "col" | "column");
    let mut i = 1;
    while i < toks.len() {
        let t = toks[i].to_ascii_lowercase();
        if let Some(v) = t.strip_prefix("n=") {
            n = Some(count("n", v)?);
        } else if let Some(v) = t.strip_prefix("nr=") {
            nr = Some(count("nr", v)?);
        } else if let Some(v) = t.strip_prefix("nc=") {
            nc = Some(count("nc", v)?);
        } else if let Some(v) = t.strip_prefix("format=") {
            format = v.to_string();
        } else if t.starts_with("nm=") && t != "nm=1" {
            return err("DL files with several matrices are not supported yet");
        } else if t.starts_with("diagonal=") && t != "diagonal=present" {
            return err("DL files with an absent diagonal are not supported yet");
        } else if t == "labels:" {
            let previous = toks[i - 1].to_ascii_lowercase();
            let target = match previous.as_str() {
                "row" => &mut row_labels,
                "col" | "column" => &mut col_labels,
                _ => &mut labels,
            };
            i += 1;
            while i < toks.len()
                && !toks[i].ends_with(':')
                && !toks[i].contains('=')
                && !(is_mode(&toks[i])
                    && toks
                        .get(i + 1)
                        .is_some_and(|x| x.to_ascii_lowercase().starts_with("labels")))
            {
                target.push(toks[i].clone());
                i += 1;
            }
            continue;
        } else if t == "embedded" || t == "embedded:" {
            embedded = true;
        }
        i += 1;
    }
    let data = words(&text[data_at + 5..]);

    let two_mode_format = matches!(format.as_str(), "edgelist2" | "el2" | "nodelist2" | "nl2");
    if nr.is_some() || nc.is_some() || two_mode_format {
        if n.is_some() {
            return err("give either n= (one-mode data) or nr= and nc= (two-mode data)");
        }
        let modes = TwoMode {
            nr,
            nc,
            row_labels,
            col_labels,
            embedded,
        };
        return read_two_mode(
            text,
            data_at,
            &data,
            &format,
            modes,
            directed_opt.unwrap_or(false),
        );
    }
    // One-mode files may label rows and columns separately; the row labels
    // name the nodes.
    if labels.is_empty() {
        labels = if row_labels.is_empty() {
            col_labels
        } else {
            row_labels
        };
    }

    let mut b = NetworkBuilder::new(directed);
    let mut links: Vec<(String, String, f64)> = Vec::new();
    let add_link = |_: &mut NetworkBuilder, a: &str, c: &str, w: f64| {
        links.push((a.to_string(), c.to_string(), w));
    };
    let mut add_link = add_link;

    match format.as_str() {
        "fullmatrix" | "fm" => {
            let n = match n {
                Some(n) => n,
                None => return err("fullmatrix data needs n="),
            };
            // Check the size before allocating anything for it.
            let needed = n
                .checked_mul(n + embedded as usize)
                .and_then(|cells| cells.checked_add(if embedded { n } else { 0 }));
            if needed.map_or(true, |k| k > data.len()) {
                return err("the matrix has too few values");
            }
            let mut pos = 0;
            if embedded {
                labels = data.get(..n).map(<[String]>::to_vec).unwrap_or_default();
                pos = n;
            }
            let names = add_nodes(&mut b, &labels, n)?;
            for r in 0..n {
                if embedded {
                    pos += 1; // row label
                }
                for c in 0..n {
                    let t = match data.get(pos) {
                        Some(t) => t,
                        None => return err("the matrix has too few values"),
                    };
                    pos += 1;
                    let x = parse_float(t).ok_or_else(|| Error::Parse {
                        format: FORMAT,
                        line: None,
                        message: format!("invalid value {t:?}"),
                    })?;
                    if x != 0.0 && (directed || c >= r) {
                        add_link(&mut b, &names[r], &names[c], x);
                    }
                }
            }
        }
        "edgelist1" | "el1" => {
            let names = match n {
                Some(n) if !embedded => add_nodes(&mut b, &labels, n)?,
                _ => Vec::new(),
            };
            let resolve = |t: &str| -> Result<String> {
                if embedded {
                    return Ok(t.to_string());
                }
                match t.parse::<usize>() {
                    Ok(k) if k >= 1 && k <= names.len() => Ok(names[k - 1].clone()),
                    _ => err(format!("{t:?} is not a node number")),
                }
            };
            // Rows are "from to [value]"; a value is present when the third
            // token is not itself a node reference on a new row. Use the
            // line structure to tell them apart.
            for line in text[data_at + 5..].lines() {
                let w = words(line);
                if w.is_empty() {
                    continue;
                }
                if w.len() < 2 {
                    return err(format!("edge row {line:?} needs two nodes"));
                }
                let value = match w.get(2) {
                    Some(t) => parse_float(t).ok_or_else(|| Error::Parse {
                        format: FORMAT,
                        line: None,
                        message: format!("invalid value {t:?}"),
                    })?,
                    None => 1.0,
                };
                add_link(&mut b, &resolve(&w[0])?, &resolve(&w[1])?, value);
            }
        }
        "nodelist1" | "nl1" => {
            let names = match n {
                Some(n) if !embedded => add_nodes(&mut b, &labels, n)?,
                _ => Vec::new(),
            };
            for line in text[data_at + 5..].lines() {
                let w = words(line);
                if w.is_empty() {
                    continue;
                }
                let get = |t: &str| -> Result<String> {
                    if embedded {
                        return Ok(t.to_string());
                    }
                    match t.parse::<usize>() {
                        Ok(k) if k >= 1 && k <= names.len() => Ok(names[k - 1].clone()),
                        _ => err(format!("{t:?} is not a node number")),
                    }
                };
                let a = get(&w[0])?;
                if !b.has_node(&a) {
                    b.add_node(&a, Vec::new());
                }
                for t in &w[1..] {
                    add_link(&mut b, &a, &get(t)?, 1.0);
                }
            }
        }
        other => return err(format!("format={other} is not supported yet")),
    }
    // Weights are recorded only when some value differs from 1.
    let weighted = links.iter().any(|l| l.2 != 1.0);
    for (a, c, w) in links {
        let attrs = if weighted {
            vec![("weight".to_string(), Value::Float(w))]
        } else {
            Vec::new()
        };
        b.add_edge(&a, &c, attrs);
    }
    b.finish()
}

/// Adds the `n` declared nodes and returns their identifiers: the labels
/// when they are unique, otherwise the numbers 1 to n with the labels as a
/// `label` attribute.
fn add_nodes(b: &mut NetworkBuilder, labels: &[String], n: usize) -> Result<Vec<String>> {
    if !labels.is_empty() && labels.len() != n {
        return err(format!("{} labels given for {n} nodes", labels.len()));
    }
    let unique = {
        let set: std::collections::HashSet<&String> = labels.iter().collect();
        set.len() == labels.len()
    };
    let ids: Vec<String> = if !labels.is_empty() && unique {
        labels.to_vec()
    } else {
        (1..=n).map(|i| i.to_string()).collect()
    };
    for (i, id) in ids.iter().enumerate() {
        let attrs = if unique {
            Vec::new()
        } else {
            vec![("label".to_string(), Value::Str(labels[i].clone()))]
        };
        b.add_node(id, attrs);
    }
    Ok(ids)
}

struct TwoMode {
    nr: Option<usize>,
    nc: Option<usize>,
    row_labels: Vec<String>,
    col_labels: Vec<String>,
    embedded: bool,
}

/// Names of one mode as they are read: given labels, or names added as
/// embedded labels appear.
struct Names {
    names: Vec<String>,
    index: std::collections::HashMap<String, usize>,
    labelled: bool,
}

impl Names {
    fn new(labels: Vec<String>) -> Names {
        let mut index = std::collections::HashMap::new();
        for (i, l) in labels.iter().enumerate() {
            index.entry(l.clone()).or_insert(i);
        }
        Names {
            labelled: !labels.is_empty(),
            names: labels,
            index,
        }
    }

    /// The position of an embedded label, adding it if new.
    fn embedded(&mut self, label: &str) -> usize {
        if let Some(&i) = self.index.get(label) {
            return i;
        }
        self.labelled = true;
        self.names.push(label.to_string());
        self.index.insert(label.to_string(), self.names.len() - 1);
        self.names.len() - 1
    }
}

/// A 1-based node number of one mode.
fn number(t: &str, count: usize, mode: &str) -> Result<usize> {
    match t.parse::<usize>() {
        Ok(k) if k >= 1 && k <= count => Ok(k - 1),
        _ => err(format!(
            "{t:?} is not a {mode} number between 1 and {count}"
        )),
    }
}

fn read_two_mode(
    text: &str,
    data_at: usize,
    data: &[String],
    format: &str,
    modes: TwoMode,
    directed: bool,
) -> Result<Network> {
    let TwoMode {
        nr,
        nc,
        row_labels,
        col_labels,
        embedded,
    } = modes;
    for (labels, count, mode) in [(&row_labels, nr, "row"), (&col_labels, nc, "column")] {
        if let Some(k) = count {
            if !labels.is_empty() && labels.len() != k {
                return err(format!(
                    "{} {mode} labels given for {k} {mode}s",
                    labels.len()
                ));
            }
        }
    }
    if nr.unwrap_or(0) + nc.unwrap_or(0) > super::MAX_DECLARED_NODES {
        return err(format!(
            "nr= and nc= declare more than the {} nodes supported",
            super::MAX_DECLARED_NODES
        ));
    }
    let mut rows = Names::new(row_labels);
    let mut cols = Names::new(col_labels);
    let mut links: Vec<(usize, usize, f64)> = Vec::new();
    let sizes = |what: &str| -> Result<(usize, usize)> {
        match (nr, nc) {
            (Some(r), Some(c)) => Ok((r, c)),
            _ => err(format!("{what} data need both nr= and nc=")),
        }
    };
    let (nrow, ncol) = match format {
        "fullmatrix" | "fm" => {
            let (nr, nc) = sizes("two-mode fullmatrix")?;
            let extra = if embedded { nr + nc } else { 0 };
            let needed = nr.checked_mul(nc).and_then(|k| k.checked_add(extra));
            if needed.map_or(true, |k| k > data.len()) {
                return err("the matrix has too few values");
            }
            let mut pos = 0;
            if embedded {
                cols = Names::new(data[..nc].to_vec());
                pos = nc;
            }
            let mut embedded_rows = Vec::new();
            for r in 0..nr {
                if embedded {
                    embedded_rows.push(data[pos].clone());
                    pos += 1;
                }
                for c in 0..nc {
                    let x = value(&data[pos])?;
                    pos += 1;
                    if x != 0.0 {
                        links.push((r, c, x));
                    }
                }
            }
            if embedded {
                rows = Names::new(embedded_rows);
            }
            (nr, nc)
        }
        "edgelist2" | "el2" | "nodelist2" | "nl2" => {
            let edge_list = matches!(format, "edgelist2" | "el2");
            let declared = if embedded {
                (nr.unwrap_or(0), nc.unwrap_or(0))
            } else {
                sizes(&format!("format={format}"))?
            };
            for line in text[data_at + 5..].lines() {
                let w = words(line);
                if w.is_empty() {
                    continue;
                }
                if w.len() < 2 && edge_list {
                    return err(format!("edge row {line:?} needs a row and a column"));
                }
                let r = if embedded {
                    rows.embedded(&w[0])
                } else {
                    number(&w[0], declared.0, "row")?
                };
                let ends: &[String] = if edge_list { &w[1..2] } else { &w[1..] };
                let x = match (edge_list, w.get(2)) {
                    (true, Some(t)) => value(t)?,
                    _ => 1.0,
                };
                for t in ends {
                    let c = if embedded {
                        cols.embedded(t)
                    } else {
                        number(t, declared.1, "column")?
                    };
                    links.push((r, c, x));
                }
            }
            if embedded {
                (
                    declared.0.max(rows.names.len()),
                    declared.1.max(cols.names.len()),
                )
            } else {
                declared
            }
        }
        "edgelist1" | "el1" | "nodelist1" | "nl1" => {
            return err(format!(
                "format={format} holds one-mode data; two-mode data use edgelist2 or nodelist2"
            ))
        }
        other => return err(format!("format={other} is not supported yet")),
    };

    // Identifiers of each mode, then prefixes if the modes share one.
    let mode_ids =
        |names: &Names, count: usize, offset: usize| -> (Vec<String>, Option<Vec<String>>) {
            let numbers = || {
                (offset + 1..=offset + count)
                    .map(|i| i.to_string())
                    .collect()
            };
            if !names.labelled {
                return (numbers(), None);
            }
            // Embedded labels can name fewer nodes than declared.
            let mut labels = names.names.clone();
            for i in labels.len()..count {
                labels.push((offset + i + 1).to_string());
            }
            let unique = labels
                .iter()
                .collect::<std::collections::HashSet<_>>()
                .len()
                == labels.len();
            if unique {
                (labels, None)
            } else {
                (numbers(), Some(labels))
            }
        };
    let (mut row_ids, mut row_label) = mode_ids(&rows, nrow, 0);
    let (mut col_ids, mut col_label) = mode_ids(&cols, ncol, nrow);
    let row_set: std::collections::HashSet<&String> = row_ids.iter().collect();
    if col_ids.iter().any(|c| row_set.contains(c)) {
        row_label = Some(row_label.unwrap_or_else(|| row_ids.clone()));
        col_label = Some(col_label.unwrap_or_else(|| col_ids.clone()));
        row_ids = row_ids.iter().map(|id| format!("row::{id}")).collect();
        col_ids = col_ids.iter().map(|id| format!("col::{id}")).collect();
    }

    let mut b = NetworkBuilder::new(directed);
    for (ids, label, is_col) in [(&row_ids, &row_label, false), (&col_ids, &col_label, true)] {
        for (i, id) in ids.iter().enumerate() {
            let mut attrs = vec![("type".to_string(), Value::Bool(is_col))];
            if let Some(l) = label {
                attrs.push(("label".to_string(), Value::Str(l[i].clone())));
            }
            b.add_node(id, attrs);
        }
    }
    let weighted = links.iter().any(|l| l.2 != 1.0);
    for (r, c, x) in links {
        let attrs = if weighted {
            vec![("weight".to_string(), Value::Float(x))]
        } else {
            Vec::new()
        };
        b.add_edge(&row_ids[r], &col_ids[c], attrs);
    }
    b.finish()
}

fn quote(s: &str) -> String {
    format!("\"{}\"", s.replace('"', "'"))
}

/// The node types of a two-mode network: a logical `type` column without
/// missing values, with every edge joining the two types.
fn node_types(
    net: &Network,
    index: &std::collections::HashMap<String, usize>,
) -> Result<Option<Vec<bool>>> {
    let types: Vec<bool> = match net.nodes.column("type") {
        Some(super::table::Column::Bool(v)) => {
            match v.iter().copied().collect::<Option<Vec<bool>>>() {
                Some(t) => t,
                None => return Ok(None),
            }
        }
        _ => return Ok(None),
    };
    let from = net.edges.texts("from")?;
    let to = net.edges.texts("to")?;
    let crossing = from
        .iter()
        .zip(&to)
        .all(|(a, b)| types[index[a.as_ref()]] != types[index[b.as_ref()]]);
    Ok(if crossing { Some(types) } else { None })
}

/// Writes a DL file: one-mode networks in `edgelist1` format with node
/// labels, so isolates and labels are preserved, and two-mode networks
/// (see the module documentation) with `nr=`, `nc=`, row and column
/// labels, and `edgelist2` data. The `weight` edge column (if numeric)
/// becomes the edge value.
pub fn write(net: &Network, weight: Option<&str>) -> Result<String> {
    net.validate()?;
    let ids = net.nodes.strings("id")?;
    let index = net.node_index()?;
    if let Some(types) = node_types(net, &index)? {
        return write_two_mode(net, weight, &ids, &index, &types);
    }
    let mut out = format!("dl n={} format=edgelist1\nlabels:\n", ids.len());
    out.push_str(&ids.iter().map(|s| quote(s)).collect::<Vec<_>>().join(","));
    out.push_str("\ndata:\n");
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

fn write_two_mode(
    net: &Network,
    weight: Option<&str>,
    ids: &[String],
    index: &std::collections::HashMap<String, usize>,
    types: &[bool],
) -> Result<String> {
    // Position of each node within its mode.
    let mut position = vec![0usize; ids.len()];
    let (mut rows, mut cols) = (Vec::new(), Vec::new());
    for (i, id) in ids.iter().enumerate() {
        let mode = if types[i] { &mut cols } else { &mut rows };
        position[i] = mode.len();
        mode.push(quote(id));
    }
    let mut out = format!(
        "dl nr={} nc={} format=edgelist2\nrow labels:\n{}\ncolumn labels:\n{}\ndata:\n",
        rows.len(),
        cols.len(),
        rows.join(","),
        cols.join(",")
    );
    let from = net.edges.strings("from")?;
    let to = net.edges.strings("to")?;
    let w = weight.or(Some("weight")).and_then(|n| net.edges.column(n));
    for row in 0..from.len() {
        let (mut a, mut c) = (index[&from[row]], index[&to[row]]);
        if types[a] {
            std::mem::swap(&mut a, &mut c);
        }
        out.push_str(&format!("{} {}", position[a] + 1, position[c] + 1));
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
    fn full_matrix_with_labels() {
        let doc = "DL N = 3\nFORMAT = FULLMATRIX\nLABELS:\nann,bob,\"cy d\"\nDATA:\n0 1 0\n1 0 2\n0 0 0\n";
        let n = read(doc, None).unwrap();
        assert_eq!(n.nodes.strings("id").unwrap(), vec!["ann", "bob", "cy d"]);
        assert_eq!(n.edges.strings("from").unwrap(), vec!["ann", "bob", "bob"]);
        assert_eq!(
            n.edges.column("weight").unwrap(),
            &Column::Float(vec![Some(1.0), Some(1.0), Some(2.0)])
        );
        let u = read(doc, Some(false)).unwrap();
        assert_eq!(u.edges.nrow(), 2);
    }

    #[test]
    fn embedded_labels_and_edge_lists() {
        let doc =
            "dl n=3 format=fullmatrix labels embedded\ndata:\na b c\na 0 1 0\nb 0 0 1\nc 0 0 0\n";
        let n = read(doc, None).unwrap();
        assert_eq!(n.edges.strings("to").unwrap(), vec!["b", "c"]);
        let el = read("dl n=4 format=edgelist1\ndata:\n1 2\n2 3 5\n", None).unwrap();
        assert_eq!(el.nodes.nrow(), 4);
        let nl = read(
            "dl format=nodelist1 labels embedded\ndata:\nx y z\ny z\n",
            None,
        )
        .unwrap();
        assert_eq!(nl.edges.nrow(), 3);
        assert!(read("dl nr=2 nc=3\ndata:\n", None)
            .unwrap_err()
            .to_string()
            .contains("too few values"));
    }

    #[test]
    fn round_trips() {
        let doc = "dl n=4 format=edgelist1\nlabels:\na,b,c,d\ndata:\n1 2\n2 3 0.5\n";
        let n = read(doc, None).unwrap();
        let again = read(&write(&n, None).unwrap(), None).unwrap();
        assert_eq!(again, n);
    }

    fn bools(v: &[bool]) -> Column {
        Column::Bool(v.iter().map(|&b| Some(b)).collect())
    }

    #[test]
    fn two_mode_full_matrix_with_row_and_column_labels() {
        let doc = "DL NR=2, NC=3\nFORMAT = FULLMATRIX\nROW LABELS:\nann bob\nCOLUMN LABELS:\nx,y,\"z w\"\nDATA:\n1 0 1\n0 2 1\n";
        let n = read(doc, None).unwrap();
        assert!(!n.directed);
        assert_eq!(
            n.nodes.strings("id").unwrap(),
            vec!["ann", "bob", "x", "y", "z w"]
        );
        assert_eq!(
            n.nodes.column("type").unwrap(),
            &bools(&[false, false, true, true, true])
        );
        assert!(n.nodes.column("label").is_none());
        assert_eq!(
            n.edges.strings("from").unwrap(),
            vec!["ann", "ann", "bob", "bob"]
        );
        assert_eq!(n.edges.strings("to").unwrap(), vec!["x", "z w", "y", "z w"]);
        assert_eq!(
            n.edges.column("weight").unwrap(),
            &Column::Float(vec![Some(1.0), Some(1.0), Some(2.0), Some(1.0)])
        );
        assert!(read(doc, Some(true)).unwrap().directed);
        // "col labels:" is accepted too.
        let col = read(&doc.replace("COLUMN LABELS", "col labels"), None).unwrap();
        assert_eq!(col.nodes, n.nodes);
    }

    #[test]
    fn two_mode_embedded_and_unlabelled() {
        let doc = "dl nr=2 nc=2 labels embedded\ndata:\n   e1 e2\nann 1 0\nbob 1 1\n";
        let n = read(doc, None).unwrap();
        assert_eq!(
            n.nodes.strings("id").unwrap(),
            vec!["ann", "bob", "e1", "e2"]
        );
        assert_eq!(n.edges.strings("to").unwrap(), vec!["e1", "e1", "e2"]);
        // Unlabelled rows are 1..nr and columns nr+1..nr+nc.
        let u = read("dl nr=2 nc=2\ndata:\n0 1\n1 0\n", None).unwrap();
        assert_eq!(u.nodes.strings("id").unwrap(), vec!["1", "2", "3", "4"]);
        assert_eq!(u.edges.strings("from").unwrap(), vec!["1", "2"]);
        assert_eq!(u.edges.strings("to").unwrap(), vec!["4", "3"]);
        assert!(u.edges.column("weight").is_none());
    }

    #[test]
    fn two_mode_edge_and_node_lists() {
        let el = "dl nr=3 nc=2 format=edgelist2\nrow labels:\na b c\ncol labels:\np q\ndata:\n1 2\n3 1 4\n";
        let n = read(el, None).unwrap();
        assert_eq!(
            n.nodes.strings("id").unwrap(),
            vec!["a", "b", "c", "p", "q"]
        );
        assert_eq!(n.edges.strings("from").unwrap(), vec!["a", "c"]);
        assert_eq!(n.edges.strings("to").unwrap(), vec!["q", "p"]);
        assert_eq!(
            n.edges.column("weight").unwrap(),
            &Column::Float(vec![Some(1.0), Some(4.0)])
        );
        let nl = "dl format=nodelist2 labels embedded\ndata:\nann party picnic\nbob picnic\n";
        let m = read(nl, None).unwrap();
        assert_eq!(
            m.nodes.strings("id").unwrap(),
            vec!["ann", "bob", "party", "picnic"]
        );
        assert_eq!(
            m.nodes.column("type").unwrap(),
            &bools(&[false, false, true, true])
        );
        assert_eq!(m.edges.nrow(), 3);
        assert!(read("dl nr=2 nc=2 format=edgelist2\ndata:\n1 3\n", None)
            .unwrap_err()
            .to_string()
            .contains("column number"));
        assert!(read("dl nr=2 format=edgelist2\ndata:\n1 1\n", None)
            .unwrap_err()
            .to_string()
            .contains("nr= and nc="));
        assert!(read("dl n=2 nr=2 nc=2\ndata:\n", None).is_err());
        assert!(read("dl nr=2 nc=2 format=edgelist1\ndata:\n1 1\n", None).is_err());
    }

    #[test]
    fn two_mode_labels_shared_by_both_modes_are_prefixed() {
        let doc = "dl format=edgelist2 labels embedded\ndata:\na a\nb a\n";
        let n = read(doc, None).unwrap();
        assert_eq!(
            n.nodes.strings("id").unwrap(),
            vec!["row::a", "row::b", "col::a"]
        );
        assert_eq!(n.nodes.strings("label").unwrap(), vec!["a", "b", "a"]);
        assert_eq!(n.edges.strings("to").unwrap(), vec!["col::a", "col::a"]);
        // Repeated labels within a mode fall back to numbers.
        let r = read("dl nr=2 nc=1\nrow labels:\nx x\ndata:\n1\n1\n", None).unwrap();
        assert_eq!(r.nodes.strings("id").unwrap(), vec!["1", "2", "3"]);
    }

    #[test]
    fn two_mode_networks_round_trip() {
        let doc = "dl nr=3 nc=2\nrow labels:\na b c\ncolumn labels:\np q\ndata:\n1 0\n0 2.5\n1 1\n";
        let n = read(doc, None).unwrap();
        let text = write(&n, None).unwrap();
        assert!(text.starts_with("dl nr=3 nc=2 format=edgelist2"), "{text}");
        assert_eq!(read(&text, None).unwrap(), n);
        // An edge stored column to row is written row to column.
        let mut b = NetworkBuilder::new(false);
        b.add_node("p", vec![("type".into(), Value::Bool(true))]);
        b.add_node("a", vec![("type".into(), Value::Bool(false))]);
        b.add_edge("p", "a", vec![]);
        let text = write(&b.finish().unwrap(), None).unwrap();
        assert!(text.contains("data:\n1 1\n"), "{text}");
        // An edge within one mode makes the network one-mode.
        let mut b = NetworkBuilder::new(false);
        b.add_node("a", vec![("type".into(), Value::Bool(false))]);
        b.add_node("b", vec![("type".into(), Value::Bool(false))]);
        b.add_edge("a", "b", vec![]);
        let text = write(&b.finish().unwrap(), None).unwrap();
        assert!(text.starts_with("dl n=2"), "{text}");
    }
}
