//! GML, the Graph Modelling Language (Himsolt 1997, `himsolt1997gml`).
//!
//! A GML file is a tree of `key value` pairs where values are integers,
//! reals, quoted strings, or bracketed lists. Nodes are `node [ id ... ]`
//! and edges `edge [ source ... target ... ]`; other top-level keys inside
//! `graph [ ]` become graph attributes. Nested lists inside a node or edge
//! (such as `graphics [ x 1 y 2 ]`) are flattened into columns named with
//! dots (`graphics.x`). Repeated keys keep the last value.
//!
//! When every node has a unique `label`, labels become the node
//! identifiers (as in NetworkX); otherwise the GML integer ids are used.
//! When writing, node identifiers go in `label` and nodes are numbered from
//! 0, so files round-trip and read naturally in NetworkX and igraph.

use super::table::{format_float, Network, NetworkBuilder, Value};
use crate::error::{Error, Result};

const FORMAT: &str = "GML";

#[derive(Debug, Clone)]
enum Gml {
    Int(i64),
    Float(f64),
    Str(String),
    List(Vec<(String, Gml)>),
}

struct Lexer<'a> {
    s: &'a str,
    pos: usize,
    line: usize,
}

#[derive(Debug, PartialEq)]
enum Tok {
    Key(String),
    Int(i64),
    Float(f64),
    Str(String),
    Open,
    Close,
}

impl<'a> Lexer<'a> {
    fn err<T>(&self, message: impl Into<String>) -> Result<T> {
        Err(Error::Parse {
            format: FORMAT,
            line: Some(self.line),
            message: message.into(),
        })
    }

    fn next(&mut self) -> Result<Option<Tok>> {
        let b = self.s.as_bytes();
        loop {
            while self.pos < b.len() && b[self.pos].is_ascii_whitespace() {
                if b[self.pos] == b'\n' {
                    self.line += 1;
                }
                self.pos += 1;
            }
            if self.pos < b.len() && b[self.pos] == b'#' {
                while self.pos < b.len() && b[self.pos] != b'\n' {
                    self.pos += 1;
                }
                continue;
            }
            break;
        }
        if self.pos >= b.len() {
            return Ok(None);
        }
        let c = b[self.pos];
        if c == b'[' {
            self.pos += 1;
            return Ok(Some(Tok::Open));
        }
        if c == b']' {
            self.pos += 1;
            return Ok(Some(Tok::Close));
        }
        if c == b'"' {
            let start = self.pos + 1;
            let end = match self.s[start..].find('"') {
                Some(e) => start + e,
                None => return self.err("unterminated string"),
            };
            let raw = &self.s[start..end];
            self.line += raw.matches('\n').count();
            self.pos = end + 1;
            return Ok(Some(Tok::Str(unescape(raw))));
        }
        let start = self.pos;
        while self.pos < b.len()
            && !b[self.pos].is_ascii_whitespace()
            && !matches!(b[self.pos], b'[' | b']' | b'"')
        {
            self.pos += 1;
        }
        let word = &self.s[start..self.pos];
        let first = word.as_bytes()[0];
        if first.is_ascii_alphabetic() || first == b'_' {
            match word {
                "NAN" | "NaN" => return Ok(Some(Tok::Float(f64::NAN))),
                "INF" => return Ok(Some(Tok::Float(f64::INFINITY))),
                _ => return Ok(Some(Tok::Key(word.to_string()))),
            }
        }
        if let Ok(i) = word.parse::<i64>() {
            return Ok(Some(Tok::Int(i)));
        }
        match word {
            "-INF" => return Ok(Some(Tok::Float(f64::NEG_INFINITY))),
            "+INF" => return Ok(Some(Tok::Float(f64::INFINITY))),
            _ => {}
        }
        match word.parse::<f64>() {
            Ok(f) => Ok(Some(Tok::Float(f))),
            Err(_) => self.err(format!("unexpected token {word:?}")),
        }
    }
}

/// GML strings use HTML-style character entities for quotes and
/// non-ASCII text, as written by NetworkX.
fn unescape(s: &str) -> String {
    if !s.contains('&') {
        return s.to_string();
    }
    let mut out = String::with_capacity(s.len());
    let mut rest = s;
    while let Some(i) = rest.find('&') {
        out.push_str(&rest[..i]);
        rest = &rest[i..];
        let end = rest.find(';').filter(|&e| e <= 10);
        let decoded = end.and_then(|e| {
            let ent = &rest[1..e];
            let ch = match ent {
                "quot" => Some('"'),
                "amp" => Some('&'),
                "lt" => Some('<'),
                "gt" => Some('>'),
                "apos" => Some('\''),
                _ if ent.starts_with("#x") => u32::from_str_radix(&ent[2..], 16)
                    .ok()
                    .and_then(char::from_u32),
                _ if ent.starts_with('#') => ent[1..].parse().ok().and_then(char::from_u32),
                _ => None,
            };
            ch.map(|c| (c, e))
        });
        match decoded {
            Some((c, e)) => {
                out.push(c);
                rest = &rest[e + 1..];
            }
            None => {
                out.push('&');
                rest = &rest[1..];
            }
        }
    }
    out.push_str(rest);
    out
}

fn escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '"' => out.push_str("&quot;"),
            '&' => out.push_str("&amp;"),
            c if !c.is_ascii() || (c as u32) < 0x20 => out.push_str(&format!("&#{};", c as u32)),
            c => out.push(c),
        }
    }
    out
}

/// The deepest nesting of lists accepted; files nest two or three levels.
const MAX_DEPTH: usize = 512;

fn parse_list(lex: &mut Lexer, depth: usize) -> Result<Vec<(String, Gml)>> {
    let nested = depth > 0;
    if depth > MAX_DEPTH {
        return lex.err(format!(
            "lists are nested more than {MAX_DEPTH} levels deep"
        ));
    }
    let mut items = Vec::new();
    loop {
        let key = match lex.next()? {
            Some(Tok::Key(k)) => k,
            Some(Tok::Close) if nested => return Ok(items),
            None if !nested => return Ok(items),
            None => return lex.err("missing closing bracket"),
            Some(t) => return lex.err(format!("expected a key, found {t:?}")),
        };
        let value = match lex.next()? {
            Some(Tok::Int(i)) => Gml::Int(i),
            Some(Tok::Float(f)) => Gml::Float(f),
            Some(Tok::Str(s)) => Gml::Str(s),
            Some(Tok::Open) => Gml::List(parse_list(lex, depth + 1)?),
            _ => return lex.err(format!("key {key:?} has no value")),
        };
        items.push((key, value));
    }
}

fn flatten(prefix: &str, items: &[(String, Gml)], skip: &[&str], out: &mut Vec<(String, Value)>) {
    for (k, v) in items {
        if prefix.is_empty() && skip.contains(&k.as_str()) {
            continue;
        }
        let name = if prefix.is_empty() {
            k.clone()
        } else {
            format!("{prefix}.{k}")
        };
        let value = match v {
            Gml::Int(i) => Value::Int(*i),
            Gml::Float(f) => Value::Float(*f),
            Gml::Str(s) => Value::Str(s.clone()),
            Gml::List(inner) => {
                flatten(&name, inner, &[], out);
                continue;
            }
        };
        match out.iter_mut().find(|(n, _)| *n == name) {
            Some(slot) => slot.1 = value,
            None => out.push((name, value)),
        }
    }
}

fn scalar_text(v: &Gml) -> Option<String> {
    match v {
        Gml::Int(i) => Some(i.to_string()),
        Gml::Str(s) => Some(s.clone()),
        Gml::Float(f) => Some(format_float(*f)),
        Gml::List(_) => None,
    }
}

/// Reads a GML document.
pub fn read(text: &str) -> Result<Network> {
    let mut lex = Lexer {
        s: text,
        pos: 0,
        line: 1,
    };
    let top = parse_list(&mut lex, 0)?;
    let graph = top
        .iter()
        .find_map(|(k, v)| match (k.as_str(), v) {
            ("graph", Gml::List(items)) => Some(items),
            _ => None,
        })
        .ok_or(Error::Parse {
            format: FORMAT,
            line: None,
            message: "no graph [ ... ] block".into(),
        })?;
    let directed = graph
        .iter()
        .any(|(k, v)| k == "directed" && matches!(v, Gml::Int(1)));
    let mut b = NetworkBuilder::new(directed);
    for (k, v) in graph {
        match (k.as_str(), v) {
            ("node", Gml::List(items)) => {
                let id = items
                    .iter()
                    .find(|(k, _)| k == "id")
                    .and_then(|(_, v)| scalar_text(v))
                    .ok_or(Error::Parse {
                        format: FORMAT,
                        line: None,
                        message: "a node has no id".into(),
                    })?;
                let mut attrs = Vec::new();
                flatten("", items, &["id"], &mut attrs);
                if !b.add_node(&id, attrs) {
                    return Err(Error::Parse {
                        format: FORMAT,
                        line: None,
                        message: format!("node id {id} appears twice"),
                    });
                }
            }
            ("edge", Gml::List(items)) => {
                let get = |key: &str| {
                    items
                        .iter()
                        .find(|(k, _)| k == key)
                        .and_then(|(_, v)| scalar_text(v))
                };
                let (s, t) = match (get("source"), get("target")) {
                    (Some(s), Some(t)) => (s, t),
                    _ => {
                        return Err(Error::Parse {
                            format: FORMAT,
                            line: None,
                            message: "an edge needs source and target".into(),
                        })
                    }
                };
                for id in [&s, &t] {
                    if !b.has_node(id) {
                        return Err(Error::Parse {
                            format: FORMAT,
                            line: None,
                            message: format!("edge refers to undefined node {id}"),
                        });
                    }
                }
                let mut attrs = Vec::new();
                flatten("", items, &["source", "target"], &mut attrs);
                b.add_edge(&s, &t, attrs);
            }
            ("directed", _) | ("multigraph", _) => {}
            (_, v) => {
                if let Some(text) = scalar_text(v) {
                    let value = match v {
                        Gml::Int(i) => Value::Int(*i),
                        Gml::Float(f) => Value::Float(*f),
                        _ => Value::Str(text),
                    };
                    b.graph.push((k.clone(), value));
                }
            }
        }
    }
    let mut net = b.finish()?;
    net.use_as_ids("label")?;
    Ok(net)
}

/// GML keys must start with a letter and contain only letters, digits,
/// and underscores.
fn key(name: &str) -> String {
    let mut k: String = name
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '_' {
                c
            } else {
                '_'
            }
        })
        .collect();
    if !k.starts_with(|c: char| c.is_ascii_alphabetic()) {
        k.insert(0, 'x');
    }
    k
}

fn value_text(v: &Value) -> Option<String> {
    match v {
        Value::Null => None,
        Value::Bool(b) => Some((*b as i32).to_string()),
        Value::Int(i) => Some(i.to_string()),
        Value::Float(f) if f.is_nan() => Some("NAN".into()),
        Value::Float(f) if f.is_infinite() => Some(if *f > 0.0 { "INF" } else { "-INF" }.into()),
        Value::Float(f) => Some(format_float(*f)),
        Value::Str(s) => Some(format!("\"{}\"", escape(s))),
    }
}

fn has_parallel_edges(net: &Network) -> Result<bool> {
    let from = net.edges.strings("from")?;
    let to = net.edges.strings("to")?;
    let mut seen = std::collections::HashSet::new();
    Ok(from.iter().zip(&to).any(|(a, b)| {
        let key = if net.directed || a <= b {
            (a, b)
        } else {
            (b, a)
        };
        !seen.insert(key)
    }))
}

/// Writes a GML document.
pub fn write(net: &Network) -> Result<String> {
    net.validate()?;
    let mut out = String::from("graph [\n");
    out.push_str(&format!("  directed {}\n", net.directed as i32));
    // Readers such as NetworkX reject parallel edges unless the file says
    // it is a multigraph.
    if has_parallel_edges(net)? {
        out.push_str("  multigraph 1\n");
    }
    for (k, v) in &net.graph {
        if let Some(t) = value_text(v) {
            out.push_str(&format!("  {} {t}\n", key(k)));
        }
    }
    let index = net.node_index()?;
    let ids = net.nodes.strings("id")?;
    let node_attrs: Vec<_> = net.node_attributes().collect();
    for (row, id) in ids.iter().enumerate() {
        out.push_str(&format!(
            "  node [\n    id {row}\n    label \"{}\"\n",
            escape(id)
        ));
        for (name, col) in &node_attrs {
            let k = if name.as_str() == "label" {
                "label_attr".to_string()
            } else {
                key(name)
            };
            if let Some(t) = value_text(&col.get(row)) {
                out.push_str(&format!("    {k} {t}\n"));
            }
        }
        out.push_str("  ]\n");
    }
    let from = net.edges.strings("from")?;
    let to = net.edges.strings("to")?;
    let edge_attrs: Vec<_> = net.edge_attributes().collect();
    for row in 0..from.len() {
        out.push_str(&format!(
            "  edge [\n    source {}\n    target {}\n",
            index[&from[row]], index[&to[row]]
        ));
        for (name, col) in &edge_attrs {
            if let Some(t) = value_text(&col.get(row)) {
                out.push_str(&format!("    {} {t}\n", key(name)));
            }
        }
        out.push_str("  ]\n");
    }
    out.push_str("]\n");
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::table::Column;

    #[test]
    fn reads_nested_and_typed_values() {
        let doc = r#"Creator "test"
graph [
  directed 1
  name "demo &quot;x&quot;"
  node [ id 1 label "a" graphics [ x 1.5 y 2 ] ]
  node [ id 2 label "b&#233;" ]
  edge [ source 1 target 2 weight 0.5 ]
  # a comment
  edge [ source 2 target 2 ]
]"#;
        let n = read(doc).unwrap();
        assert!(n.directed);
        assert_eq!(
            n.graph,
            vec![("name".into(), Value::Str("demo \"x\"".into()))]
        );
        assert_eq!(n.nodes.strings("id").unwrap(), vec!["a", "b\u{e9}"]);
        assert!(n.nodes.column("label").is_none());
        assert_eq!(n.edges.strings("from").unwrap(), vec!["a", "b\u{e9}"]);
        assert_eq!(
            n.nodes.column("graphics.x").unwrap(),
            &Column::Float(vec![Some(1.5), None])
        );
        assert_eq!(
            n.edges.column("weight").unwrap(),
            &Column::Float(vec![Some(0.5), None])
        );
    }

    #[test]
    fn writes_identifiers_as_labels() {
        let n =
            read("graph [ node [ id 7 ] node [ id 9 ] edge [ source 7 target 9 w 2 ] ]").unwrap();
        let text = write(&n).unwrap();
        assert!(text.contains("label \"7\""));
        let again = read(&text).unwrap();
        assert_eq!(again.nodes.strings("id").unwrap(), vec!["7", "9"]);
        assert_eq!(again.edges.strings("to").unwrap(), vec!["9"]);
        assert_eq!(
            again.edges.column("w").unwrap(),
            &Column::Int(vec![Some(2)])
        );
        assert!(read("graph [ edge [ source 1 target 2 ] ]").is_err());
    }
}
