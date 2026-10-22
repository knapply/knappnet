//! JSON node-link data, the format of NetworkX's `node_link_data` and of
//! D3 force layouts (`bostock2011d3`).
//!
//! ```json
//! {"directed": true, "graph": {...},
//!  "nodes": [{"id": "a", ...}],
//!  "edges": [{"source": "a", "target": "b", ...}]}
//! ```
//!
//! Edges are written under `edges` (the NetworkX default) and read from
//! either `edges` or `links` (the D3 convention). Endpoints may be node
//! identifiers or, as in some D3 data, 0-based positions in `nodes` when
//! nodes have no `id`. Nested attribute values are kept as JSON text.

use super::json::{self, Json};
use super::table::{Column, Network, NetworkBuilder, Value};
use crate::error::{Error, Result};

const FORMAT: &str = "JSON node-link";

fn err<T>(message: impl Into<String>) -> Result<T> {
    Err(Error::Parse {
        format: FORMAT,
        line: None,
        message: message.into(),
    })
}

fn to_value(v: &Json) -> Value {
    match v {
        Json::Null => Value::Null,
        Json::Bool(b) => Value::Bool(*b),
        Json::Int(i) => Value::Int(*i),
        Json::Float(f) => Value::Float(*f),
        Json::Str(s) => Value::Str(s.clone()),
        other => {
            let mut s = String::new();
            other.write(&mut s);
            Value::Str(s)
        }
    }
}

fn id_text(v: &Json) -> Option<String> {
    match v {
        Json::Str(s) => Some(s.clone()),
        Json::Int(i) => Some(i.to_string()),
        Json::Float(f) => Some(super::table::format_float(*f)),
        Json::Bool(b) => Some(b.to_string()),
        _ => None,
    }
}

fn members(v: &Json) -> &[(String, Json)] {
    match v {
        Json::Object(m) => m,
        _ => &[],
    }
}

/// Reads node-link JSON.
pub fn read(text: &str) -> Result<Network> {
    let doc = json::parse(text, FORMAT)?;
    if !matches!(doc, Json::Object(_)) {
        return err("expected a JSON object");
    }
    let directed = matches!(doc.get("directed"), Some(Json::Bool(true)));
    let mut b = NetworkBuilder::new(directed);
    if let Some(g) = doc.get("graph") {
        for (k, v) in members(g) {
            b.graph.push((k.clone(), to_value(v)));
        }
    }
    let nodes = match doc.get("nodes") {
        Some(Json::Array(n)) => n.as_slice(),
        Some(_) => return err("\"nodes\" must be an array"),
        None => &[],
    };
    let mut ids = Vec::with_capacity(nodes.len());
    for (i, node) in nodes.iter().enumerate() {
        let id = match node.get("id") {
            Some(v) => id_text(v).ok_or_else(|| Error::Parse {
                format: FORMAT,
                line: None,
                message: format!("node {i} has an invalid id"),
            })?,
            None => i.to_string(),
        };
        let attrs = members(node)
            .iter()
            .filter(|(k, _)| k != "id")
            .map(|(k, v)| (k.clone(), to_value(v)))
            .collect();
        if !b.add_node(&id, attrs) {
            return err(format!("node id {id:?} appears twice"));
        }
        ids.push(id);
    }
    let positional = nodes.iter().all(|n| n.get("id").is_none());
    let links = match doc.get("links").or_else(|| doc.get("edges")) {
        Some(Json::Array(l)) => l.as_slice(),
        Some(_) => return err("\"links\" must be an array"),
        None => &[],
    };
    for (i, link) in links.iter().enumerate() {
        let end = |key: &str| -> Result<String> {
            let v = match link.get(key) {
                Some(v) => v,
                None => return err(format!("link {i} has no {key}")),
            };
            match v {
                Json::Int(k) if positional => {
                    ids.get(*k as usize).cloned().ok_or_else(|| Error::Parse {
                        format: FORMAT,
                        line: None,
                        message: format!("link {i}: {key} {k} is out of range"),
                    })
                }
                Json::Object(_) => match v.get("id").and_then(id_text) {
                    Some(id) => Ok(id),
                    None => err(format!("link {i}: {key} object has no id")),
                },
                _ => id_text(v).ok_or_else(|| Error::Parse {
                    format: FORMAT,
                    line: None,
                    message: format!("link {i}: invalid {key}"),
                }),
            }
        };
        let (s, t) = (end("source")?, end("target")?);
        let attrs = members(link)
            .iter()
            .filter(|(k, _)| k != "source" && k != "target" && k != "key")
            .map(|(k, v)| (k.clone(), to_value(v)))
            .collect();
        b.add_edge(&s, &t, attrs);
    }
    b.finish()
}

fn cell(col: &Column, row: usize) -> Option<Json> {
    match col.get(row) {
        Value::Null => None,
        Value::Bool(b) => Some(Json::Bool(b)),
        Value::Int(i) => Some(Json::Int(i)),
        Value::Float(f) => Some(Json::Float(f)),
        Value::Str(s) => Some(Json::Str(s)),
    }
}

/// Writes node-link JSON readable by NetworkX and D3.
pub fn write(net: &Network) -> Result<String> {
    net.validate()?;
    let ids = net.nodes.strings("id")?;
    let from = net.edges.strings("from")?;
    let to = net.edges.strings("to")?;
    let mut seen = std::collections::HashSet::new();
    let multigraph = from.iter().zip(&to).any(|(a, b)| {
        let key = if net.directed || a <= b {
            (a, b)
        } else {
            (b, a)
        };
        !seen.insert(key)
    });
    let node_attrs: Vec<_> = net.node_attributes().collect();
    let edge_attrs: Vec<_> = net.edge_attributes().collect();
    let nodes = ids
        .iter()
        .enumerate()
        .map(|(row, id)| {
            let mut m = vec![("id".to_string(), Json::Str(id.clone()))];
            m.extend(
                node_attrs
                    .iter()
                    .filter_map(|(n, c)| cell(c, row).map(|v| (n.to_string(), v))),
            );
            Json::Object(m)
        })
        .collect();
    let links = (0..from.len())
        .map(|row| {
            let mut m = vec![
                ("source".to_string(), Json::Str(from[row].clone())),
                ("target".to_string(), Json::Str(to[row].clone())),
            ];
            m.extend(
                edge_attrs
                    .iter()
                    .filter_map(|(n, c)| cell(c, row).map(|v| (n.to_string(), v))),
            );
            Json::Object(m)
        })
        .collect();
    let graph = net
        .graph
        .iter()
        .map(|(k, v)| {
            let j = match v {
                Value::Null => Json::Null,
                Value::Bool(b) => Json::Bool(*b),
                Value::Int(i) => Json::Int(*i),
                Value::Float(f) => Json::Float(*f),
                Value::Str(s) => Json::Str(s.clone()),
            };
            (k.clone(), j)
        })
        .collect();
    let doc = Json::Object(vec![
        ("directed".into(), Json::Bool(net.directed)),
        ("multigraph".into(), Json::Bool(multigraph)),
        ("graph".into(), Json::Object(graph)),
        ("nodes".into(), Json::Array(nodes)),
        ("edges".into(), Json::Array(links)),
    ]);
    let mut out = String::new();
    doc.write(&mut out);
    out.push('\n');
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_networkx_and_d3_styles() {
        let nx = r#"{"directed": false, "multigraph": false, "graph": {"name": "g"},
            "nodes": [{"id": 0, "club": "A"}, {"id": 1, "club": "B", "pos": [1, 2]}],
            "edges": [{"source": 0, "target": 1, "weight": 4}]}"#;
        let n = read(nx).unwrap();
        assert!(!n.directed);
        assert_eq!(n.nodes.strings("id").unwrap(), vec!["0", "1"]);
        assert_eq!(
            n.nodes.column("pos").unwrap(),
            &Column::Str(vec![None, Some("[1, 2]".into())])
        );
        assert_eq!(
            n.edges.column("weight").unwrap(),
            &Column::Int(vec![Some(4)])
        );

        let d3 =
            r#"{"nodes": [{"name": "a"}, {"name": "b"}], "links": [{"source": 0, "target": 1}]}"#;
        let d = read(d3).unwrap();
        assert_eq!(d.edges.strings("to").unwrap(), vec!["1"]);
    }

    #[test]
    fn round_trips() {
        let src = r#"{"directed": true, "graph": {"k": 1.5}, "nodes": [{"id": "a", "x": 1}, {"id": "b"}],
            "links": [{"source": "a", "target": "b", "w": true}, {"source": "a", "target": "b"}]}"#;
        let n = read(src).unwrap();
        let text = write(&n).unwrap();
        assert!(text.contains("\"multigraph\": true"));
        assert_eq!(read(&text).unwrap(), n);
    }
}
