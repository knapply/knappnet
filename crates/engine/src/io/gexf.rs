//! GEXF, the native exchange format of Gephi (`bastian2009gephi`), as
//! specified by the GEXF Working Group (`gexf2009`).
//!
//! Reading supports versions 1.1 through 1.3: node labels, edge weights,
//! edge labels and types, declared attributes (`<attributes>` with
//! `<attvalues>`), `viz` position, size, and color, and the `start`/`end`
//! (or `startopen`/`endopen`) bounds of dynamic graphs, which are kept as
//! text columns. Spells and dynamic attribute values are not read. Graphs
//! with `defaultedgetype="mutual"` are read as undirected.
//!
//! Files are read in one streaming pass over the document
//! ([`xml::Reader`]), without building a tree. Files that do not follow the
//! layout of the GEXF schema (declarations after the nodes or edges, edges
//! before nodes) and files with errors are read again by the tree-based
//! reader, so every file gives exactly the network, or the error message
//! and line, that the tree-based reader gives.

use std::borrow::Cow;

use super::table::{parse_bool, parse_float, Kind, Network, NetworkBuilder, Value};
use super::xml::{self, escape, Element, Event};
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

fn parse_value(raw: &str, kind: Kind) -> Option<Value> {
    match kind {
        Kind::Str => Some(Value::Str(raw.to_string())),
        Kind::Bool => parse_bool(raw).map(Value::Bool),
        Kind::Int => raw.trim().parse().ok().map(Value::Int),
        Kind::Float => parse_float(raw).map(Value::Float),
    }
}

fn convert(raw: &str, kind: Kind, el: &Element) -> Result<Value> {
    parse_value(raw, kind).map_or_else(
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
    match read_stream(text) {
        Ok(Some(net)) => Ok(net),
        // An unusual layout or an error: the tree-based reader gives the
        // result, or the error with its line.
        _ => read_tree(text),
    }
}

const BOUNDS: [&str; 4] = ["start", "end", "startopen", "endopen"];
const NODE_EXTRA: [&str; 5] = ["label", "start", "end", "startopen", "endopen"];
const EDGE_EXTRA: [&str; 7] = [
    "id",
    "label",
    "weight",
    "start",
    "end",
    "startopen",
    "endopen",
];

/// The values of `viz` elements (position, size, color, thickness), named
/// as columns.
fn viz_values(name: &str, r: &xml::Reader, out: &mut Vec<(String, Value)>) {
    match name {
        "position" => {
            for axis in ["x", "y", "z"] {
                if let Some(v) = r.attr(axis).and_then(parse_float) {
                    out.push((format!("viz_{axis}"), Value::Float(v)));
                }
            }
        }
        "color" => {
            for ch in ["r", "g", "b", "a"] {
                if let Some(v) = r.attr(ch).and_then(parse_float) {
                    out.push((format!("viz_{ch}"), Value::Float(v)));
                }
            }
        }
        "size" | "thickness" => {
            if let Some(v) = r.attr("value").and_then(parse_float) {
                out.push((format!("viz_{name}"), Value::Float(v)));
            }
        }
        _ => {}
    }
}

/// The node or edge being read.
struct Item<'a> {
    edge: bool,
    id: Option<Cow<'a, str>>,
    source: Option<Cow<'a, str>>,
    target: Option<Cow<'a, str>>,
    edge_type: Option<Cow<'a, str>>,
    extra: Vec<(String, Value)>,
    // (attribute id, raw value) from the first <attvalues>.
    values: Vec<(Cow<'a, str>, Cow<'a, str>)>,
    viz: Vec<(String, Value)>,
    hierarchical: bool,
    seen_attvalues: bool,
}

/// The streaming reader. `Ok(None)` (or an error) means the file must be
/// read by [`read_tree`].
fn read_stream(text: &str) -> Result<Option<Network>> {
    let mut r = xml::Reader::new(text, FORMAT);
    match r.next()? {
        Some(Event::Start { name: "gexf", .. }) => {}
        _ => return Ok(None),
    }
    loop {
        match r.next()? {
            Some(Event::Start { name: "graph", .. }) if r.depth() == 2 => break,
            Some(_) => {}
            None => return Ok(None),
        }
    }
    let directed = match r.attr("defaultedgetype") {
        Some("directed") | None => true,
        Some("undirected") | Some("mutual") => false,
        Some(_) => return Ok(None),
    };
    let mut b = NetworkBuilder::new(directed);
    let mut node_decl: Vec<Attribute> = Vec::new();
    let mut edge_decl: Vec<Attribute> = Vec::new();
    // Declared kinds are given to the builder when the first <nodes> or
    // <edges> starts; declarations after that need the tree reader.
    let mut declared = false;
    let (mut seen_nodes, mut seen_edges) = (false, false);

    #[derive(Clone, Copy, PartialEq, Eq)]
    enum Section {
        Other,
        NodeAttributes,
        EdgeAttributes,
        Nodes,
        Edges,
    }
    let mut section = Section::Other;
    let mut attribute: Option<Attribute> = None;
    let mut in_default = false;
    let mut item: Option<Item> = None;
    let mut in_attvalues = false;

    loop {
        let event = match r.next()? {
            Some(e) => e,
            None => return Ok(None),
        };
        match event {
            Event::Start { name, .. } => match (r.depth(), section) {
                (3, _) => {
                    section = match name {
                        "attributes" => {
                            if declared {
                                return Ok(None);
                            }
                            match r.attr("class") {
                                Some("node") => Section::NodeAttributes,
                                Some("edge") => Section::EdgeAttributes,
                                _ => Section::Other,
                            }
                        }
                        "nodes" if !seen_nodes => {
                            if seen_edges {
                                return Ok(None);
                            }
                            seen_nodes = true;
                            Section::Nodes
                        }
                        "edges" if !seen_edges => {
                            seen_edges = true;
                            Section::Edges
                        }
                        _ => Section::Other,
                    };
                    if matches!(section, Section::Nodes | Section::Edges) && !declared {
                        declared = true;
                        for a in &node_decl {
                            b.node_kind(&a.title, a.kind);
                        }
                        for a in &edge_decl {
                            b.edge_kind(&a.title, a.kind);
                        }
                    }
                }
                (4, Section::NodeAttributes | Section::EdgeAttributes) => {
                    if name == "attribute" {
                        let id = r.attr("id").unwrap_or_default().to_string();
                        attribute = Some(Attribute {
                            title: r
                                .attr("title")
                                .map(str::to_string)
                                .unwrap_or_else(|| id.clone()),
                            id,
                            kind: kind_of(r.attr("type").unwrap_or("string")),
                            default: None,
                        });
                    }
                }
                (5, Section::NodeAttributes | Section::EdgeAttributes) => {
                    if let Some(a) = attribute.as_mut() {
                        if name == "default" && a.default.is_none() {
                            a.default = Some(String::new());
                            in_default = true;
                        }
                    }
                }
                (4, Section::Nodes | Section::Edges) => {
                    let edge = section == Section::Edges;
                    if name != if edge { "edge" } else { "node" } {
                        continue;
                    }
                    let mut extra = Vec::new();
                    for key in if edge {
                        &EDGE_EXTRA[..]
                    } else {
                        &NODE_EXTRA[..]
                    } {
                        if let Some(v) = r.attr(key) {
                            let value = if *key == "weight" {
                                match parse_float(v) {
                                    Some(f) => Value::Float(f),
                                    None => return Ok(None),
                                }
                            } else {
                                Value::Str(v.to_string())
                            };
                            extra.push((key.to_string(), value));
                        }
                    }
                    item = Some(Item {
                        edge,
                        id: r.attr_cow("id"),
                        source: r.attr_cow("source"),
                        target: r.attr_cow("target"),
                        edge_type: r.attr_cow("type"),
                        extra,
                        values: Vec::new(),
                        viz: Vec::new(),
                        hierarchical: false,
                        seen_attvalues: false,
                    });
                }
                (5, Section::Nodes | Section::Edges) => {
                    if let Some(it) = item.as_mut() {
                        match name {
                            "nodes" => it.hierarchical = true,
                            "attvalues" if !it.seen_attvalues => {
                                it.seen_attvalues = true;
                                in_attvalues = true;
                            }
                            _ => viz_values(name, &r, &mut it.viz),
                        }
                    }
                }
                (6, Section::Nodes | Section::Edges) => {
                    if let Some(it) = item.as_mut() {
                        if in_attvalues && name == "attvalue" {
                            let key = r
                                .attr_cow("for")
                                .or_else(|| r.attr_cow("id"))
                                .unwrap_or_default();
                            let raw = r.attr_cow("value").unwrap_or_default();
                            it.values.push((key, raw));
                        }
                    }
                }
                _ => {}
            },
            Event::Text(t) => {
                if in_default && r.depth() == 5 {
                    if let Some(Attribute {
                        default: Some(d), ..
                    }) = attribute.as_mut()
                    {
                        d.push_str(&t);
                    }
                }
            }
            Event::End { .. } => match r.depth() {
                1 => break,
                2 => section = Section::Other,
                3 => {
                    if let Some(a) = attribute.take() {
                        match section {
                            Section::NodeAttributes => node_decl.push(a),
                            Section::EdgeAttributes => edge_decl.push(a),
                            _ => {}
                        }
                    }
                    if let Some(it) = item.take() {
                        let decl = if it.edge { &edge_decl } else { &node_decl };
                        let attrs = match assemble(&it, decl) {
                            Some(a) => a,
                            None => return Ok(None),
                        };
                        if it.edge {
                            let (s, t) = match (&it.source, &it.target) {
                                (Some(s), Some(t)) => (s, t),
                                _ => return Ok(None),
                            };
                            if let Some(kind) = &it.edge_type {
                                if (kind == "directed") != directed {
                                    return Ok(None);
                                }
                            }
                            b.add_edge(s, t, attrs);
                        } else {
                            if it.hierarchical {
                                return Ok(None);
                            }
                            let id = match &it.id {
                                Some(id) => id,
                                None => return Ok(None),
                            };
                            if !b.add_node(id, attrs) {
                                return Ok(None);
                            }
                        }
                    }
                }
                4 => {
                    in_default = false;
                    in_attvalues = false;
                }
                _ => {}
            },
        }
    }
    if !declared {
        for a in &node_decl {
            b.node_kind(&a.title, a.kind);
        }
        for a in &edge_decl {
            b.edge_kind(&a.title, a.kind);
        }
    }
    // The rest of the document must still be well formed.
    while r.next()?.is_some() {}
    b.finish().map(Some)
}

/// The attributes of a node or edge, in the order the tree-based reader
/// gives them, or `None` if one is invalid.
fn assemble(it: &Item, decl: &[Attribute]) -> Option<Vec<(String, Value)>> {
    let mut out = it.extra.clone();
    for a in decl {
        if let Some(d) = &a.default {
            out.push((a.title.clone(), parse_value(d, a.kind)?));
        }
    }
    for (key, raw) in &it.values {
        let a = decl.iter().find(|a| a.id == *key)?;
        let value = parse_value(raw, a.kind)?;
        match out.iter_mut().find(|(n, _)| *n == a.title) {
            Some(slot) => slot.1 = value,
            None => out.push((a.title.clone(), value)),
        }
    }
    out.extend(it.viz.iter().cloned());
    Some(out)
}

/// Reads a GEXF document through a tree of its elements. Gives the same
/// result as the streaming reader on every file the streaming reader
/// accepts, and reports errors.
fn read_tree(text: &str) -> Result<Network> {
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
    let bounds = BOUNDS;

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

    #[test]
    fn unusual_layouts_use_the_tree_reader() {
        // Edges before nodes: the tree reader adds the nodes first.
        let swapped = r#"<gexf><graph><edges><edge source="b" target="a"/></edges>
            <nodes><node id="a"/><node id="b"/></nodes></graph></gexf>"#;
        assert!(read_stream(swapped).unwrap().is_none());
        assert_eq!(
            read(swapped).unwrap().nodes.strings("id").unwrap(),
            vec!["a", "b"]
        );
        // Declarations after the nodes.
        let late = r#"<gexf><graph><nodes><node id="a"/></nodes>
            <attributes class="node"><attribute id="0" title="g" type="integer"/></attributes>
            </graph></gexf>"#;
        assert!(read_stream(late).unwrap().is_none());
        assert_eq!(
            read(late).unwrap().nodes.column("g").unwrap(),
            &Column::Int(vec![None])
        );
    }

    /// Random edits of GEXF documents: on every document the streaming
    /// reader accepts, it gives exactly the tree-based reader's network,
    /// and `read` always gives the tree-based reader's answer.
    #[test]
    fn streaming_agrees_with_the_tree_reader() {
        use crate::rng::Stream;
        let docs = [
            DOC.to_string(),
            DOC.replace("defaultedgetype=\"directed\"", "defaultedgetype=\"mutual\""),
            write(&read(DOC).unwrap()).unwrap(),
        ];
        const PIECES: [&str; 26] = [
            "<node id=\"n9\"/>",
            "<node id=\"n0\"/>",
            "<node>",
            "</node>",
            "<nodes>",
            "</nodes>",
            "<edge source=\"n0\" target=\"n9\" weight=\"x\"/>",
            "<edge source=\"n1\" target=\"n0\" type=\"undirected\"/>",
            "<edges>",
            "</edges>",
            "<attvalues><attvalue for=\"0\" value=\"7\"/></attvalues>",
            "<attvalue id=\"1\" value=\"no\"/>",
            "<attvalue for=\"9\" value=\"1\"/>",
            "<attributes class=\"node\"><attribute id=\"2\" title=\"group\"><default>z</default></attribute></attributes>",
            "<attributes class=\"edge\"><attribute id=\"0\" type=\"float\"/></attributes>",
            "<default>4</default>",
            "<viz:size value=\"2\"/>",
            "<viz:color r=\"1\" b=\"x\"/>",
            "<viz:position x=\"3\"/>",
            "<thickness value=\"1\"/>",
            " label=\"L\"",
            " start=\"1\"",
            " weight=\"2\"",
            "<!-- c -->",
            "&amp;",
            "\n",
        ];
        let mut accepted = 0;
        for case in 0..20_000u64 {
            let mut rng = Stream::new(13, case);
            let mut s = docs[rng.below(docs.len() as u64) as usize].clone();
            for _ in 0..1 + rng.below(3) {
                let mut at = rng.below(s.len() as u64 + 1) as usize;
                while !s.is_char_boundary(at) {
                    at -= 1;
                }
                match rng.below(4) {
                    0 => {
                        let end = (at + 1 + rng.below(20) as usize).min(s.len());
                        if s.is_char_boundary(end) {
                            s.replace_range(at..end, "");
                        }
                    }
                    _ => s.insert_str(at, PIECES[rng.below(PIECES.len() as u64) as usize]),
                }
            }
            // Debug output compares NaN values as equal.
            let tree = format!("{:?}", read_tree(&s));
            if let Ok(Some(net)) = read_stream(&s) {
                accepted += 1;
                assert_eq!(
                    tree,
                    format!("{:?}", Ok::<_, Error>(net)),
                    "case {case}: {s}"
                );
            }
            assert_eq!(format!("{:?}", read(&s)), tree, "case {case}");
        }
        assert!(accepted > 1_000, "only {accepted} documents were accepted");
    }
}
