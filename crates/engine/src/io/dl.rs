//! UCINET DL files (Borgatti, Everett, and Freeman 2002,
//! `borgatti2002ucinet`).
//!
//! Supported: one-mode data in `fullmatrix` (the default), `edgelist1`,
//! and `nodelist1` formats; `n=` headers; node labels given with `labels:`
//! or embedded in the data (`labels embedded`). Keywords are
//! case-insensitive and may be separated by spaces, commas, or newlines.
//! Two-mode (`nr=`/`nc=`), multi-matrix, and half-matrix files are not yet
//! supported and produce a clear error.
//!
//! Labels become node identifiers when they are unique. When they are not
//! (UCINET does not require it), nodes are numbered from 1 and the labels
//! are kept in a `label` column, as for Pajek files; reading repeated
//! labels as identifiers would merge distinct nodes.
//!
//! DL files do not record direction. Networks are read as directed unless
//! the caller says otherwise; when read as undirected, a full matrix
//! contributes only its upper triangle (diagonal included).

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

/// Reads a DL file. `directed` defaults to true.
pub fn read(text: &str, directed: Option<bool>) -> Result<Network> {
    let directed = directed.unwrap_or(true);
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
    let mut format = "fullmatrix".to_string();
    let mut labels: Vec<String> = Vec::new();
    let mut embedded = false;
    let mut i = 1;
    while i < toks.len() {
        let t = toks[i].to_ascii_lowercase();
        if let Some(v) = t.strip_prefix("n=") {
            let count: usize = v.parse().map_err(|_| Error::Parse {
                format: FORMAT,
                line: None,
                message: format!("invalid n {v:?}"),
            })?;
            if count > super::MAX_DECLARED_NODES {
                return err(format!(
                    "n={count} declares more than the {} nodes supported",
                    super::MAX_DECLARED_NODES
                ));
            }
            n = Some(count);
        } else if let Some(v) = t.strip_prefix("format=") {
            format = v.to_string();
        } else if t.starts_with("nr=") || t.starts_with("nc=") {
            return err("two-mode DL files are not supported yet");
        } else if t.starts_with("nm=") && t != "nm=1" {
            return err("DL files with several matrices are not supported yet");
        } else if t.starts_with("diagonal=") && t != "diagonal=present" {
            return err("DL files with an absent diagonal are not supported yet");
        } else if t == "labels:" || t == "row" || t == "col" || t == "column" {
            if t == "labels:" {
                i += 1;
                while i < toks.len() && !toks[i].ends_with(':') && !toks[i].contains('=') {
                    labels.push(toks[i].clone());
                    i += 1;
                }
                continue;
            }
        } else if t == "embedded" || t == "embedded:" {
            embedded = true;
        }
        i += 1;
    }
    let data = words(&text[data_at + 5..]);

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

fn quote(s: &str) -> String {
    format!("\"{}\"", s.replace('"', "'"))
}

/// Writes a DL file in `edgelist1` format with node labels, so isolates
/// and labels are preserved. The `weight` edge column (if numeric)
/// becomes the edge value.
pub fn write(net: &Network, weight: Option<&str>) -> Result<String> {
    net.validate()?;
    let ids = net.nodes.strings("id")?;
    let index = net.node_index()?;
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
            .contains("two-mode"));
    }

    #[test]
    fn round_trips() {
        let doc = "dl n=4 format=edgelist1\nlabels:\na,b,c,d\ndata:\n1 2\n2 3 0.5\n";
        let n = read(doc, None).unwrap();
        let again = read(&write(&n, None).unwrap(), None).unwrap();
        assert_eq!(again, n);
    }
}
