//! Graphviz DOT output (Gansner and North 2000, `gansner2000open`), for
//! drawing networks with Graphviz and the many tools that read DOT.
//! Reading DOT is not supported.

use super::table::{Network, Value};
use crate::error::Result;

fn quote(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

fn attr_list(items: &[(String, String)]) -> String {
    if items.is_empty() {
        return String::new();
    }
    let body: Vec<String> = items
        .iter()
        .map(|(k, v)| format!("{}={}", quote(k), quote(v)))
        .collect();
    format!(" [{}]", body.join(", "))
}

/// Writes a DOT document. Attributes become DOT attributes, so Graphviz
/// styling columns such as `color`, `label`, or `penwidth` take effect.
pub fn write(net: &Network) -> Result<String> {
    net.validate()?;
    let (kw, arrow) = if net.directed {
        ("digraph", "->")
    } else {
        ("graph", "--")
    };
    let mut out = format!("{kw} {{\n");
    let from = net.edges.strings("from")?;
    let to = net.edges.strings("to")?;
    for (k, v) in &net.graph {
        if let Some(t) = v.to_text() {
            out.push_str(&format!("  {}={};\n", quote(k), quote(&t)));
        }
    }
    let node_attrs: Vec<_> = net.node_attributes().collect();
    for (row, id) in net.nodes.strings("id")?.iter().enumerate() {
        let items: Vec<(String, String)> = node_attrs
            .iter()
            .filter_map(|(n, c)| c.get(row).to_text().map(|t| (n.to_string(), t)))
            .collect();
        out.push_str(&format!("  {}{};\n", quote(id), attr_list(&items)));
    }
    let edge_attrs: Vec<_> = net.edge_attributes().collect();
    for row in 0..from.len() {
        let items: Vec<(String, String)> = edge_attrs
            .iter()
            .filter_map(|(n, c)| match c.get(row) {
                Value::Null => None,
                v => v.to_text().map(|t| (n.to_string(), t)),
            })
            .collect();
        out.push_str(&format!(
            "  {} {arrow} {}{};\n",
            quote(&from[row]),
            quote(&to[row]),
            attr_list(&items)
        ));
    }
    out.push_str("}\n");
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::table::NetworkBuilder;

    #[test]
    fn writes_quoted_dot() {
        let mut b = NetworkBuilder::new(true);
        b.add_node("a \"x\"", vec![("color".into(), Value::Str("red".into()))]);
        b.add_edge("a \"x\"", "b", vec![("w".into(), Value::Float(2.0))]);
        let text = write(&b.finish().unwrap()).unwrap();
        assert!(text.starts_with("digraph {"));
        assert!(text.contains("\"a \\\"x\\\"\" [\"color\"=\"red\"];"));
        assert!(text.contains("-> \"b\" [\"w\"=\"2.0\"];"));
    }
}
