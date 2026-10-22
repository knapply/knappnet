//! GEXF, the native exchange format of Gephi (`bastian2009gephi`), as
//! specified by the GEXF Working Group (`gexf2009`).
//!
//! Reading supports versions 1.1 through 1.3: node labels, edge weights,
//! edge labels and types, declared attributes (`<attributes>` with
//! `<attvalues>`), `viz` position, size, and color, and the `start`/`end`
//! (or `startopen`/`endopen`) bounds of dynamic graphs, which are kept as
//! text columns. Spells and dynamic attribute values are not read. Graphs
//! with `defaultedgetype="mutual"` are read as undirected.

use super::table::{parse_bool, parse_float, Kind, Network, NetworkBuilder, Value};
use super::xml::{self, escape, Element};
use crate::error::{Error, Result};

const FORMAT: &str = "GEXF";

fn err<T>(el: &Element, message: impl Into<String>) -> Result<T> {
    Err(Error::Parse {
        format: FORMAT,
        line: Some(el.line),
        message: message.into(),
    })
}

struct Attribute {
    id: String,
    title: String,
    kind: Kind,
    default: Option<String>,
}

fn kind_of(t: &str) -> Kind {
    match t {
        "integer" | "long" | "short" | "byte" => Kind::Int,
        "double" | "float" | "bigdecimal" => Kind::Float,
        "boolean" => Kind::Bool,
        _ => Kind::Str,
    }
}

fn convert(raw: &str, kind: Kind, el: &Element) -> Result<Value> {
    let v = match kind {
        Kind::Str => Some(Value::Str(raw.to_string())),
        Kind::Bool => parse_bool(raw).map(Value::Bool),
        Kind::Int => raw.trim().parse().ok().map(Value::Int),
        Kind::Float => parse_float(raw).map(Value::Float),
    };
    v.map_or_else(
        || err(el, format!("value {raw:?} is not a valid {kind:?}")),
        Ok,
    )
}

fn declared(graph: &Element, class: &str) -> Vec<Attribute> {
    graph
        .children_named("attributes")
        .filter(|a| a.attr("class") == Some(class))
        .flat_map(|a| a.children_named("attribute"))
        .map(|a| {
            let id = a.attr("id").unwrap_or_default().to_string();
            Attribute {
                title: a
                    .attr("title")
                    .map(str::to_string)
                    .unwrap_or_else(|| id.clone()),
                id,
                kind: kind_of(a.attr("type").unwrap_or("string")),
                default: a.child("default").map(|d| d.text.clone()),
            }
        })
        .collect()
}

fn element_attrs(el: &Element, decl: &[Attribute], extra: &[&str]) -> Result<Vec<(String, Value)>> {
    let mut out = Vec::new();
    for name in extra {
        if let Some(v) = el.attr(name) {
            let value = if *name == "weight" {
                match parse_float(v) {
                    Some(f) => Value::Float(f),
                    None => return err(el, format!("invalid weight {v:?}")),
                }
            } else {
                Value::Str(v.to_string())
            };
            out.push((name.to_string(), value));
        }
    }
    for a in decl {
        if let Some(d) = &a.default {
            out.push((a.title.clone(), convert(d, a.kind, el)?));
        }
    }
    if let Some(values) = el.child("attvalues") {
        for v in values.children_named("attvalue") {
            let key = v.attr("for").or_else(|| v.attr("id")).unwrap_or_default();
            let a = match decl.iter().find(|a| a.id == key) {
                Some(a) => a,
                None => return err(v, format!("undeclared attribute {key:?}")),
            };
            let raw = v.attr("value").unwrap_or_default();
            let value = convert(raw, a.kind, v)?;
            match out.iter_mut().find(|(n, _)| *n == a.title) {
                Some(slot) => slot.1 = value,
                None => out.push((a.title.clone(), value)),
            }
        }
    }
    for viz in &el.children {
        let prefix = match viz.name.as_str() {
            "position" => "viz_",
            "size" => "viz_size",
            "color" => "viz_",
            "thickness" => "viz_thickness",
            _ => continue,
        };
        match viz.name.as_str() {
            "position" => {
                for axis in ["x", "y", "z"] {
                    if let Some(v) = viz.attr(axis).and_then(parse_float) {
                        out.push((format!("{prefix}{axis}"), Value::Float(v)));
                    }
                }
            }
            "color" => {
                for ch in ["r", "g", "b", "a"] {
                    if let Some(v) = viz.attr(ch).and_then(parse_float) {
                        out.push((format!("{prefix}{ch}"), Value::Float(v)));
                    }
                }
            }
            _ => {
                if let Some(v) = viz.attr("value").and_then(parse_float) {
                    out.push((prefix.to_string(), Value::Float(v)));
                }
            }
        }
    }
    Ok(out)
}

/// Reads a GEXF document.
pub fn read(text: &str) -> Result<Network> {
    let root = xml::parse(text, FORMAT)?;
    if root.name != "gexf" {
        return err(&root, format!("expected <gexf>, found <{}>", root.name));
    }
    let graph = match root.child("graph") {
        Some(g) => g,
        None => return err(&root, "no <graph> element"),
    };
    let directed = match graph.attr("defaultedgetype") {
        Some("directed") | None => true,
        Some("undirected") | Some("mutual") => false,
        Some(other) => return err(graph, format!("unknown defaultedgetype {other:?}")),
    };
    let node_decl = declared(graph, "node");
    let edge_decl = declared(graph, "edge");
    let mut b = NetworkBuilder::new(directed);
    for a in &node_decl {
        b.node_kind(&a.title, a.kind);
    }
    for a in &edge_decl {
        b.edge_kind(&a.title, a.kind);
    }
    let bounds = ["start", "end", "startopen", "endopen"];

    if let Some(nodes) = graph.child("nodes") {
        for node in nodes.children_named("node") {
            if node.child("nodes").is_some() {
                return err(node, "hierarchical graphs are not supported");
            }
            let id = match node.attr("id") {
                Some(id) => id,
                None => return err(node, "<node> has no id"),
            };
            let mut extra = vec!["label"];
            extra.extend(bounds);
            let attrs = element_attrs(node, &node_decl, &extra)?;
            if !b.add_node(id, attrs) {
                return err(node, format!("node {id:?} is declared twice"));
            }
        }
    }
    if let Some(edges) = graph.child("edges") {
        for edge in edges.children_named("edge") {
            let (s, t) = match (edge.attr("source"), edge.attr("target")) {
                (Some(s), Some(t)) => (s, t),
                _ => return err(edge, "<edge> needs source and target"),
            };
            if let Some(kind) = edge.attr("type") {
                let edge_directed = kind == "directed";
                if edge_directed != directed {
                    return err(
                        edge,
                        "edges whose type differs from defaultedgetype are not supported",
                    );
                }
            }
            let mut extra = vec!["id", "label", "weight"];
            extra.extend(bounds);
            let attrs = element_attrs(edge, &edge_decl, &extra)?;
            b.add_edge(s, t, attrs);
        }
    }
    b.finish()
}

fn gexf_type(kind: Kind) -> &'static str {
    match kind {
        Kind::Bool => "boolean",
        Kind::Int => "long",
        Kind::Float => "double",
        Kind::Str => "string",
    }
}

/// Writes a GEXF 1.2 document, the version read by default by Gephi and
/// NetworkX. A text column `label` becomes the node or
/// edge label and a numeric edge column `weight` becomes the edge weight;
/// other columns are written as declared attributes.
pub fn write(net: &Network) -> Result<String> {
    net.validate()?;
    let mut out = String::from(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
         <gexf xmlns=\"http://www.gexf.net/1.2draft\" version=\"1.2\">\n",
    );
    out.push_str(&format!(
        "  <graph defaultedgetype=\"{}\" mode=\"static\">\n",
        if net.directed {
            "directed"
        } else {
            "undirected"
        }
    ));

    let special = |name: &str, kind: Kind, edge: bool| {
        (name == "label" && kind == Kind::Str)
            || (edge && name == "weight" && matches!(kind, Kind::Float | Kind::Int))
            || (edge && name == "id" && kind == Kind::Str)
    };
    let node_attrs: Vec<_> = net
        .node_attributes()
        .filter(|(n, c)| !special(n, c.kind(), false))
        .collect();
    let edge_attrs: Vec<_> = net
        .edge_attributes()
        .filter(|(n, c)| !special(n, c.kind(), true))
        .collect();
    for (class, attrs) in [("node", &node_attrs), ("edge", &edge_attrs)] {
        if attrs.is_empty() {
            continue;
        }
        out.push_str(&format!("    <attributes class=\"{class}\">\n"));
        for (i, (name, col)) in attrs.iter().enumerate() {
            out.push_str(&format!(
                "      <attribute id=\"{i}\" title=\"{}\" type=\"{}\"/>\n",
                escape(name),
                gexf_type(col.kind())
            ));
        }
        out.push_str("    </attributes>\n");
    }

    let label = |table: &super::table::Table, row: usize| -> String {
        match table.column("label") {
            Some(c) if c.kind() == Kind::Str => c
                .get(row)
                .to_text()
                .map(|l| format!(" label=\"{}\"", escape(&l)))
                .unwrap_or_default(),
            _ => String::new(),
        }
    };
    let attvalues = |attrs: &[(&String, &super::table::Column)], row: usize| -> String {
        let vals: Vec<String> = attrs
            .iter()
            .enumerate()
            .filter_map(|(i, (_, c))| {
                c.get(row).to_text().map(|v| {
                    format!(
                        "          <attvalue for=\"{i}\" value=\"{}\"/>\n",
                        escape(&v)
                    )
                })
            })
            .collect();
        if vals.is_empty() {
            String::new()
        } else {
            format!(
                "        <attvalues>\n{}        </attvalues>\n",
                vals.concat()
            )
        }
    };

    out.push_str("    <nodes>\n");
    for (row, id) in net.nodes.strings("id")?.iter().enumerate() {
        let inner = attvalues(&node_attrs, row);
        let open = format!(
            "      <node id=\"{}\"{}",
            escape(id),
            label(&net.nodes, row)
        );
        if inner.is_empty() {
            out.push_str(&format!("{open}/>\n"));
        } else {
            out.push_str(&format!("{open}>\n{inner}      </node>\n"));
        }
    }
    out.push_str("    </nodes>\n    <edges>\n");
    let from = net.edges.strings("from")?;
    let to = net.edges.strings("to")?;
    let ids = net.edges.column("id").filter(|c| c.kind() == Kind::Str);
    let weights = net
        .edges
        .column("weight")
        .filter(|c| matches!(c.kind(), Kind::Float | Kind::Int));
    for row in 0..from.len() {
        let id = ids
            .and_then(|c| c.get(row).to_text())
            .unwrap_or_else(|| row.to_string());
        let mut open = format!(
            "      <edge id=\"{}\" source=\"{}\" target=\"{}\"",
            escape(&id),
            escape(&from[row]),
            escape(&to[row])
        );
        if let Some(w) = weights.and_then(|c| c.get(row).to_text()) {
            open.push_str(&format!(" weight=\"{w}\""));
        }
        open.push_str(&label(&net.edges, row));
        let inner = attvalues(&edge_attrs, row);
        if inner.is_empty() {
            out.push_str(&format!("{open}/>\n"));
        } else {
            out.push_str(&format!("{open}>\n{inner}      </edge>\n"));
        }
    }
    out.push_str("    </edges>\n  </graph>\n</gexf>\n");
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::table::Column;

    const DOC: &str = r#"<?xml version="1.0" encoding="UTF-8"?>
<gexf xmlns="http://gexf.net/1.3" xmlns:viz="http://gexf.net/1.3/viz" version="1.3">
  <graph defaultedgetype="directed">
    <attributes class="node">
      <attribute id="0" title="group" type="integer"><default>1</default></attribute>
      <attribute id="1" title="active" type="boolean"/>
    </attributes>
    <nodes>
      <node id="n0" label="Alice">
        <attvalues><attvalue for="0" value="3"/><attvalue for="1" value="true"/></attvalues>
        <viz:position x="1.5" y="-2.0" z="0.0"/>
      </node>
      <node id="n1" label="Bob" start="2001"/>
    </nodes>
    <edges>
      <edge id="e0" source="n0" target="n1" weight="2.5" label="knows"/>
    </edges>
  </graph>
</gexf>"#;

    #[test]
    fn reads_attributes_viz_and_weights() {
        let n = read(DOC).unwrap();
        assert!(n.directed);
        assert_eq!(
            n.nodes.column("group").unwrap(),
            &Column::Int(vec![Some(3), Some(1)])
        );
        assert_eq!(
            n.nodes.column("active").unwrap(),
            &Column::Bool(vec![Some(true), None])
        );
        assert_eq!(
            n.nodes.column("label").unwrap(),
            &Column::Str(vec![Some("Alice".into()), Some("Bob".into())])
        );
        assert_eq!(
            n.nodes.column("viz_x").unwrap(),
            &Column::Float(vec![Some(1.5), None])
        );
        assert_eq!(
            n.nodes.column("start").unwrap(),
            &Column::Str(vec![None, Some("2001".into())])
        );
        assert_eq!(
            n.edges.column("weight").unwrap(),
            &Column::Float(vec![Some(2.5)])
        );
    }

    #[test]
    fn round_trips() {
        let n = read(DOC).unwrap();
        let again = read(&write(&n).unwrap()).unwrap();
        assert_eq!(
            again.nodes.strings("id").unwrap(),
            n.nodes.strings("id").unwrap()
        );
        for name in ["group", "active", "label", "viz_x", "start"] {
            assert_eq!(again.nodes.column(name), n.nodes.column(name), "{name}");
        }
        for name in ["from", "to", "id", "weight", "label"] {
            assert_eq!(again.edges.column(name), n.edges.column(name), "{name}");
        }
    }
}
