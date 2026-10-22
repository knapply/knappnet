//! GraphML (Brandes et al. 2002, `brandes2002graphml`): the XML exchange
//! format used by igraph, NetworkX, Gephi, yEd, Cytoscape, and others.
//!
//! Reading supports typed `<key>` declarations with defaults, node, edge,
//! and graph attributes, `edgedefault`, and edge identifiers. Only the
//! first `<graph>` is read; nested graphs, hyperedges, and ports are
//! rejected with an error. Attributes declared for `all` apply to nodes,
//! edges, and the graph.

use super::table::{parse_bool, parse_float, Kind, Network, NetworkBuilder, Value};
use super::xml::{self, escape, Element};
use crate::error::{Error, Result};

const FORMAT: &str = "GraphML";

struct Key {
    name: String,
    kind: Kind,
    default: Option<String>,
}

fn err<T>(el: &Element, message: impl Into<String>) -> Result<T> {
    Err(Error::Parse {
        format: FORMAT,
        line: Some(el.line),
        message: message.into(),
    })
}

fn convert(raw: &str, kind: Kind, el: &Element) -> Result<Value> {
    let v = match kind {
        Kind::Str => Some(Value::Str(raw.to_string())),
        Kind::Bool => parse_bool(raw).map(Value::Bool),
        Kind::Int => raw.trim().parse().ok().map(Value::Int),
        Kind::Float => parse_float(raw).map(Value::Float),
    };
    match v {
        Some(v) => Ok(v),
        None => err(el, format!("value {raw:?} is not a valid {kind:?}")),
    }
}

/// Reads a GraphML document.
pub fn read(text: &str) -> Result<Network> {
    let root = xml::parse(text, FORMAT)?;
    if root.name != "graphml" {
        return err(&root, format!("expected <graphml>, found <{}>", root.name));
    }

    // Keys in document order, so column order is deterministic.
    let mut keys: Vec<(String, String, Key)> = Vec::new();
    for k in root.children_named("key") {
        let id = match k.attr("id") {
            Some(id) => id.to_string(),
            None => return err(k, "<key> has no id"),
        };
        let domain = k.attr("for").unwrap_or("all").to_string();
        let kind = match k.attr("attr.type").unwrap_or("string") {
            "boolean" => Kind::Bool,
            "int" | "long" => Kind::Int,
            "float" | "double" => Kind::Float,
            _ => Kind::Str,
        };
        let name = k.attr("attr.name").unwrap_or(&id).to_string();
        let default = k.child("default").map(|d| d.text.clone());
        keys.push((
            domain,
            id,
            Key {
                name,
                kind,
                default,
            },
        ));
    }
    let lookup = |domain: &str, id: &str| -> Option<&Key> {
        let find = |d: &str| keys.iter().find(|(kd, kid, _)| kd == d && kid == id);
        find(domain).or_else(|| find("all")).map(|(_, _, k)| k)
    };

    let graph = match root.child("graph") {
        Some(g) => g,
        None => return err(&root, "no <graph> element"),
    };
    let directed = match graph.attr("edgedefault") {
        Some("directed") | None => true,
        Some("undirected") => false,
        Some(other) => return err(graph, format!("unknown edgedefault {other:?}")),
    };
    let mut b = NetworkBuilder::new(directed);

    // Declared kinds and defaults, so columns keep their types even when
    // every value is missing.
    let mut node_defaults = Vec::new();
    let mut edge_defaults = Vec::new();
    for (domain, _, key) in &keys {
        let (for_nodes, for_edges) = match domain.as_str() {
            "node" => (true, false),
            "edge" => (false, true),
            "all" => (true, true),
            _ => (false, false),
        };
        if for_nodes {
            b.node_kind(&key.name, key.kind);
            if let Some(d) = &key.default {
                node_defaults.push((key.name.clone(), convert(d, key.kind, &root)?));
            }
        }
        if for_edges {
            b.edge_kind(&key.name, key.kind);
            if let Some(d) = &key.default {
                edge_defaults.push((key.name.clone(), convert(d, key.kind, &root)?));
            }
        }
    }

    let read_data = |el: &Element,
                     domain: &str,
                     defaults: &[(String, Value)]|
     -> Result<Vec<(String, Value)>> {
        let mut attrs: Vec<(String, Value)> = defaults.to_vec();
        for d in el.children_named("data") {
            let id = match d.attr("key") {
                Some(id) => id,
                None => return err(d, "<data> has no key"),
            };
            let key = match lookup(domain, id) {
                Some(k) => k,
                None => return err(d, format!("undeclared key {id:?}")),
            };
            let value = convert(&d.text, key.kind, d)?;
            match attrs.iter_mut().find(|(n, _)| *n == key.name) {
                Some(slot) => slot.1 = value,
                None => attrs.push((key.name.clone(), value)),
            }
        }
        Ok(attrs)
    };

    for d in graph.children_named("data") {
        if let Some(id) = d.attr("key") {
            if let Some(key) = lookup("graph", id) {
                b.graph
                    .push((key.name.clone(), convert(&d.text, key.kind, d)?));
            }
        }
    }

    for child in &graph.children {
        match child.name.as_str() {
            "node" => {
                if child.child("graph").is_some() {
                    return err(child, "nested graphs are not supported");
                }
                if child.child("port").is_some() {
                    return err(child, "ports are not supported");
                }
                let id = match child.attr("id") {
                    Some(id) => id,
                    None => return err(child, "<node> has no id"),
                };
                let attrs = read_data(child, "node", &node_defaults)?;
                if !b.add_node(id, attrs) {
                    return err(child, format!("node {id:?} is declared twice"));
                }
            }
            "edge" => {
                let (s, t) = match (child.attr("source"), child.attr("target")) {
                    (Some(s), Some(t)) => (s, t),
                    _ => return err(child, "<edge> needs source and target"),
                };
                if let Some(d) = child.attr("directed") {
                    if parse_bool(d) != Some(directed) {
                        return err(
                            child,
                            "edges whose direction differs from edgedefault are not supported",
                        );
                    }
                }
                // GraphML expects every node to be declared; nodes that only
                // appear in edges still receive the declared defaults.
                for id in [s, t] {
                    if !b.has_node(id) {
                        b.add_node(id, node_defaults.clone());
                    }
                }
                let mut attrs = read_data(child, "edge", &edge_defaults)?;
                if let Some(id) = child.attr("id") {
                    attrs.insert(0, ("id".to_string(), Value::Str(id.to_string())));
                }
                b.add_edge(s, t, attrs);
            }
            "hyperedge" => return err(child, "hyperedges are not supported"),
            _ => {}
        }
    }
    b.finish()
}

fn type_name(kind: Kind) -> &'static str {
    match kind {
        Kind::Bool => "boolean",
        Kind::Int => "long",
        Kind::Float => "double",
        Kind::Str => "string",
    }
}

/// Writes a GraphML document.
pub fn write(net: &Network) -> Result<String> {
    net.validate()?;
    let mut out = String::from(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
         <graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\" \
         xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" \
         xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns \
         http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">\n",
    );
    let mut key_id = 0;
    let mut node_keys = Vec::new();
    let mut edge_keys = Vec::new();
    for (name, value) in &net.graph {
        let kind = match value {
            Value::Bool(_) => Kind::Bool,
            Value::Int(_) => Kind::Int,
            Value::Float(_) => Kind::Float,
            _ => Kind::Str,
        };
        out.push_str(&format!(
            "  <key id=\"g{key_id}\" for=\"graph\" attr.name=\"{}\" attr.type=\"{}\"/>\n",
            escape(name),
            type_name(kind)
        ));
        key_id += 1;
    }
    for (name, col) in net.node_attributes() {
        out.push_str(&format!(
            "  <key id=\"d{key_id}\" for=\"node\" attr.name=\"{}\" attr.type=\"{}\"/>\n",
            escape(name),
            type_name(col.kind())
        ));
        node_keys.push((format!("d{key_id}"), col));
        key_id += 1;
    }
    let mut edge_ids = None;
    for (name, col) in net.edge_attributes() {
        if name == "id" && col.kind() == Kind::Str {
            edge_ids = Some(col);
            continue;
        }
        out.push_str(&format!(
            "  <key id=\"d{key_id}\" for=\"edge\" attr.name=\"{}\" attr.type=\"{}\"/>\n",
            escape(name),
            type_name(col.kind())
        ));
        edge_keys.push((format!("d{key_id}"), col));
        key_id += 1;
    }
    out.push_str(&format!(
        "  <graph edgedefault=\"{}\">\n",
        if net.directed {
            "directed"
        } else {
            "undirected"
        }
    ));
    for (i, (_, value)) in net.graph.iter().enumerate() {
        if let Some(text) = value.to_text() {
            out.push_str(&format!(
                "    <data key=\"g{i}\">{}</data>\n",
                escape(&text)
            ));
        }
    }
    let ids = net.nodes.strings("id")?;
    for (row, id) in ids.iter().enumerate() {
        write_element(
            &mut out,
            "node",
            &format!("id=\"{}\"", escape(id)),
            &node_keys,
            row,
        );
    }
    let from = net.edges.strings("from")?;
    let to = net.edges.strings("to")?;
    for row in 0..from.len() {
        let mut attrs = String::new();
        if let Some(Some(id)) = edge_ids.map(|c| c.get(row).to_text()) {
            attrs.push_str(&format!("id=\"{}\" ", escape(&id)));
        }
        attrs.push_str(&format!(
            "source=\"{}\" target=\"{}\"",
            escape(&from[row]),
            escape(&to[row])
        ));
        write_element(&mut out, "edge", &attrs, &edge_keys, row);
    }
    out.push_str("  </graph>\n</graphml>\n");
    Ok(out)
}

fn write_element(
    out: &mut String,
    tag: &str,
    attrs: &str,
    keys: &[(String, &super::table::Column)],
    row: usize,
) {
    let data: Vec<String> = keys
        .iter()
        .filter_map(|(key, col)| {
            col.get(row)
                .to_text()
                .map(|t| format!("      <data key=\"{key}\">{}</data>\n", escape(&t)))
        })
        .collect();
    if data.is_empty() {
        out.push_str(&format!("    <{tag} {attrs}/>\n"));
    } else {
        out.push_str(&format!("    <{tag} {attrs}>\n"));
        for d in data {
            out.push_str(&d);
        }
        out.push_str(&format!("    </{tag}>\n"));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::table::Column;

    const DOC: &str = r#"<?xml version="1.0"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns">
  <key id="d0" for="node" attr.name="color" attr.type="string"><default>gray</default></key>
  <key id="d1" for="edge" attr.name="weight" attr.type="double"/>
  <key id="d2" for="node" attr.name="age" attr.type="int"/>
  <key id="g" for="graph" attr.name="title" attr.type="string"/>
  <graph edgedefault="undirected">
    <data key="g">demo</data>
    <node id="a"><data key="d0">red</data><data key="d2">30</data></node>
    <node id="b"/>
    <edge id="e1" source="a" target="b"><data key="d1">1.5</data></edge>
    <edge source="b" target="c"/>
  </graph>
</graphml>"#;

    #[test]
    fn reads_types_defaults_and_implicit_nodes() {
        let n = read(DOC).unwrap();
        assert!(!n.directed);
        assert_eq!(n.graph, vec![("title".into(), Value::Str("demo".into()))]);
        assert_eq!(n.nodes.strings("id").unwrap(), vec!["a", "b", "c"]);
        assert_eq!(
            n.nodes.column("color").unwrap(),
            &Column::Str(vec![
                Some("red".into()),
                Some("gray".into()),
                Some("gray".into())
            ])
        );
        assert_eq!(
            n.nodes.column("age").unwrap(),
            &Column::Int(vec![Some(30), None, None])
        );
        assert_eq!(
            n.edges.column("weight").unwrap(),
            &Column::Float(vec![Some(1.5), None])
        );
        assert_eq!(
            n.edges.column("id").unwrap(),
            &Column::Str(vec![Some("e1".into()), None])
        );
    }

    #[test]
    fn round_trips() {
        let n = read(DOC).unwrap();
        let again = read(&write(&n).unwrap()).unwrap();
        assert_eq!(again.nodes, n.nodes);
        assert_eq!(again.edges, n.edges);
        assert_eq!(again.graph, n.graph);
        assert_eq!(again.directed, n.directed);
    }

    #[test]
    fn rejects_unsupported_features() {
        let hyper = DOC.replace("<edge source=\"b\" target=\"c\"/>", "<hyperedge/>");
        assert!(read(&hyper).unwrap_err().to_string().contains("hyperedges"));
        let bad = DOC.replace(">30<", ">thirty<");
        assert!(read(&bad).unwrap_err().to_string().contains("line"));
    }
}
