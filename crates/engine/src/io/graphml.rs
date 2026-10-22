//! GraphML (Brandes et al. 2002, `brandes2002graphml`): the XML exchange
//! format used by igraph, NetworkX, Gephi, yEd, Cytoscape, and others.
//!
//! Reading supports typed `<key>` declarations with defaults, node, edge,
//! and graph attributes, `edgedefault`, and edge identifiers. Only the
//! first `<graph>` is read; nested graphs, hyperedges, and ports are
//! rejected with an error. Attributes declared for `all` apply to nodes,
//! edges, and the graph.
//!
//! Files are read in one streaming pass over the document
//! ([`xml::Reader`]), without building a tree. Files that do not follow the
//! layout of the GraphML schema (a `<key>` after the `<graph>`) and files
//! with errors are read again by the tree-based reader, so every file gives
//! exactly the network, or the error message and line, that the tree-based
//! reader gives.

use std::borrow::Cow;
use std::collections::HashMap;

use super::table::{parse_bool, parse_float, Column, Kind, Network, NetworkBuilder, Value};
use super::xml::{self, escape, escape_into, Element, Event};
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

fn parse_value(raw: &str, kind: Kind) -> Option<Value> {
    match kind {
        Kind::Str => Some(Value::Str(raw.to_string())),
        Kind::Bool => parse_bool(raw).map(Value::Bool),
        Kind::Int => raw.trim().parse().ok().map(Value::Int),
        Kind::Float => parse_float(raw).map(Value::Float),
    }
}

fn convert(raw: &str, kind: Kind, el: &Element) -> Result<Value> {
    match parse_value(raw, kind) {
        Some(v) => Ok(v),
        None => err(el, format!("value {raw:?} is not a valid {kind:?}")),
    }
}

fn key_kind(t: &str) -> Kind {
    match t {
        "boolean" => Kind::Bool,
        "int" | "long" => Kind::Int,
        "float" | "double" => Kind::Float,
        _ => Kind::Str,
    }
}

/// Reads a GraphML document.
pub fn read(text: &str) -> Result<Network> {
    match read_stream(text) {
        Ok(Some(net)) => Ok(net),
        // An unusual layout or an error: the tree-based reader gives the
        // result, or the error with its line.
        _ => read_tree(text),
    }
}

/// Key lookup for one domain: the first key declared for the domain with a
/// given id, or else the first declared for `all`.
fn domain_index(keys: &[(String, String, Key)], domain: &str) -> HashMap<String, usize> {
    let mut map = HashMap::new();
    for want in [domain, "all"] {
        for (i, (d, id, _)) in keys.iter().enumerate() {
            if d == want {
                map.entry(id.clone()).or_insert(i);
            }
        }
    }
    map
}

/// What is being read inside the `<graph>`.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Item {
    Other,
    Node,
    Edge,
    GraphData,
}

/// Appends text to a value that may still borrow from the document.
fn append<'a>(buf: &mut Cow<'a, str>, text: Cow<'a, str>) {
    if buf.is_empty() {
        *buf = text;
    } else {
        buf.to_mut().push_str(&text);
    }
}

/// The streaming reader. `Ok(None)` (or an error) means the file must be
/// read by [`read_tree`].
fn read_stream(text: &str) -> Result<Option<Network>> {
    let mut r = xml::Reader::new(text, FORMAT);
    match r.next()? {
        Some(Event::Start {
            name: "graphml", ..
        }) => {}
        _ => return Ok(None),
    }

    // Keys, up to the first <graph>.
    let mut keys: Vec<(String, String, Key)> = Vec::new();
    loop {
        let event = match r.next()? {
            Some(e) => e,
            None => return Ok(None),
        };
        let name = match event {
            Event::Start { name, .. } if r.depth() == 2 => name,
            _ => continue,
        };
        if name == "graph" {
            break;
        }
        if name != "key" {
            continue;
        }
        let id = match r.attr("id") {
            Some(id) => id.to_string(),
            None => return Ok(None),
        };
        let domain = r.attr("for").unwrap_or("all").to_string();
        let kind = key_kind(r.attr("attr.type").unwrap_or("string"));
        let name = r.attr("attr.name").unwrap_or(&id).to_string();
        // The direct text of the first <default> child.
        let mut default: Option<String> = None;
        let mut in_default = false;
        while r.depth() >= 2 {
            match r.next()? {
                Some(Event::Start { name, .. }) => {
                    if name == "default" && r.depth() == 3 && default.is_none() {
                        default = Some(String::new());
                        in_default = true;
                    }
                }
                Some(Event::Text(t)) => {
                    if in_default && r.depth() == 3 {
                        if let Some(d) = default.as_mut() {
                            d.push_str(&t);
                        }
                    }
                }
                Some(Event::End { .. }) => {
                    if r.depth() == 2 {
                        in_default = false;
                    }
                }
                None => return Ok(None),
            }
        }
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

    let directed = match r.attr("edgedefault") {
        Some("directed") | None => true,
        Some("undirected") => false,
        Some(_) => return Ok(None),
    };
    let mut b = NetworkBuilder::new(directed);
    let mut node_defaults = Vec::new();
    let mut edge_defaults = Vec::new();
    for (domain, _, key) in &keys {
        let (for_nodes, for_edges) = match domain.as_str() {
            "node" => (true, false),
            "edge" => (false, true),
            "all" => (true, true),
            _ => (false, false),
        };
        let default = match &key.default {
            Some(d) if for_nodes || for_edges => match parse_value(d, key.kind) {
                Some(v) => Some(v),
                None => return Ok(None),
            },
            _ => None,
        };
        if for_nodes {
            b.node_kind(&key.name, key.kind);
            if let Some(v) = &default {
                node_defaults.push((key.name.clone(), v.clone()));
            }
        }
        if for_edges {
            b.edge_kind(&key.name, key.kind);
            if let Some(v) = default {
                edge_defaults.push((key.name.clone(), v));
            }
        }
    }
    let node_keys = domain_index(&keys, "node");
    let edge_keys = domain_index(&keys, "edge");
    let graph_keys = domain_index(&keys, "graph");

    let mut item = Item::Other;
    let mut id: Option<Cow<str>> = None;
    let mut source: Option<Cow<str>> = None;
    let mut target: Option<Cow<str>> = None;
    let mut edge_directed: Option<Cow<str>> = None;
    // <data> children of the current node or edge: key and text.
    let mut data: Vec<(Option<Cow<str>>, Cow<str>)> = Vec::new();
    // Depth whose text belongs to the open <data>, if any.
    let mut capture: Option<usize> = None;

    loop {
        let event = match r.next()? {
            Some(e) => e,
            None => return Ok(None),
        };
        match event {
            Event::Start { name, .. } => {
                let depth = r.depth();
                if depth == 3 {
                    data.clear();
                    item = match name {
                        "node" => {
                            id = r.attr_cow("id");
                            Item::Node
                        }
                        "edge" => {
                            id = r.attr_cow("id");
                            source = r.attr_cow("source");
                            target = r.attr_cow("target");
                            edge_directed = r.attr_cow("directed");
                            Item::Edge
                        }
                        "hyperedge" => return Ok(None),
                        "data" => {
                            data.push((r.attr_cow("key"), Cow::Borrowed("")));
                            capture = Some(3);
                            Item::GraphData
                        }
                        _ => Item::Other,
                    };
                } else if depth == 4 && matches!(item, Item::Node | Item::Edge) {
                    match name {
                        "data" => {
                            data.push((r.attr_cow("key"), Cow::Borrowed("")));
                            capture = Some(4);
                        }
                        "graph" | "port" if item == Item::Node => return Ok(None),
                        _ => {}
                    }
                }
            }
            Event::Text(t) => {
                if capture == Some(r.depth()) {
                    if let Some((_, buf)) = data.last_mut() {
                        append(buf, t);
                    }
                }
            }
            Event::End { .. } => {
                let depth = r.depth();
                if capture == Some(depth + 1) {
                    capture = None;
                }
                if depth == 1 {
                    break;
                }
                if depth != 2 {
                    continue;
                }
                match item {
                    Item::Node => {
                        let node = match id.take() {
                            Some(id) => id,
                            None => return Ok(None),
                        };
                        let attrs = match collect(&data, &node_keys, &keys, &node_defaults) {
                            Some(a) => a,
                            None => return Ok(None),
                        };
                        if !b.add_node(&node, attrs) {
                            return Ok(None);
                        }
                    }
                    Item::Edge => {
                        let (s, t) = match (source.take(), target.take()) {
                            (Some(s), Some(t)) => (s, t),
                            _ => return Ok(None),
                        };
                        if let Some(d) = edge_directed.take() {
                            if parse_bool(&d) != Some(directed) {
                                return Ok(None);
                            }
                        }
                        let mut attrs = match collect(&data, &edge_keys, &keys, &edge_defaults) {
                            Some(a) => a,
                            None => return Ok(None),
                        };
                        if let Some(edge_id) = id.take() {
                            attrs.insert(0, ("id".to_string(), Value::Str(edge_id.into_owned())));
                        }
                        // Nodes that only appear in edges get the declared
                        // defaults.
                        b.add_edge_with(&s, &t, attrs, &node_defaults);
                    }
                    Item::GraphData => {
                        if let Some((Some(k), text)) = data.first() {
                            if let Some(&i) = graph_keys.get(k.as_ref()) {
                                let key = &keys[i].2;
                                match parse_value(text, key.kind) {
                                    Some(v) => b.graph.push((key.name.clone(), v)),
                                    None => return Ok(None),
                                }
                            }
                        }
                    }
                    Item::Other => {}
                }
                item = Item::Other;
            }
        }
    }

    // The rest of the document: a key here applies to the whole file, which
    // only the tree-based reader handles; other graphs are ignored.
    while let Some(event) = r.next()? {
        if let Event::Start { name: "key", .. } = event {
            if r.depth() == 2 {
                return Ok(None);
            }
        }
    }
    b.finish().map(Some)
}

/// The attributes of a node or edge from its `<data>` children, or `None`
/// if one is invalid.
fn collect(
    data: &[(Option<Cow<str>>, Cow<str>)],
    index: &HashMap<String, usize>,
    keys: &[(String, String, Key)],
    defaults: &[(String, Value)],
) -> Option<Vec<(String, Value)>> {
    let mut attrs: Vec<(String, Value)> = defaults.to_vec();
    for (key, text) in data {
        let key = &keys[*index.get(key.as_deref()?)?].2;
        let value = parse_value(text, key.kind)?;
        match attrs.iter_mut().find(|(n, _)| *n == key.name) {
            Some(slot) => slot.1 = value,
            None => attrs.push((key.name.clone(), value)),
        }
    }
    Some(attrs)
}

/// Reads a GraphML document through a tree of its elements. Gives the same
/// result as the streaming reader on every file the streaming reader
/// accepts, and reports errors.
fn read_tree(text: &str) -> Result<Network> {
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
        let kind = key_kind(k.attr("attr.type").unwrap_or("string"));
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
    let mut out = String::with_capacity(64 * (net.nodes.nrow() + net.edges.nrow()) + 1024);
    out.push_str(
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
    let ids = text_column(&net.nodes, "id")?;
    for (row, id) in ids.iter().enumerate() {
        out.push_str("    <node id=\"");
        escape_into(&mut out, id.as_deref().unwrap_or_default());
        out.push('"');
        write_data(&mut out, "node", &node_keys, row);
    }
    let from = text_column(&net.edges, "from")?;
    let to = text_column(&net.edges, "to")?;
    for (row, (from, to)) in from.iter().zip(to).enumerate() {
        out.push_str("    <edge ");
        if let Some(Some(id)) = edge_ids.map(|c| c.get(row).to_text()) {
            out.push_str("id=\"");
            escape_into(&mut out, &id);
            out.push_str("\" ");
        }
        out.push_str("source=\"");
        escape_into(&mut out, from.as_deref().unwrap_or_default());
        out.push_str("\" target=\"");
        escape_into(&mut out, to.as_deref().unwrap_or_default());
        out.push('"');
        write_data(&mut out, "edge", &edge_keys, row);
    }
    out.push_str("  </graph>\n</graphml>\n");
    Ok(out)
}

/// The values of a text column, which [`Network::validate`] has checked
/// exists and has no missing values.
fn text_column<'a>(table: &'a super::table::Table, name: &str) -> Result<&'a [Option<String>]> {
    match table.column(name) {
        Some(Column::Str(v)) => Ok(v),
        _ => Err(Error::InvalidArgument(format!(
            "column {name:?} must hold text"
        ))),
    }
}

/// Finishes the start tag opened by the caller and writes the `<data>`
/// children of one row.
fn write_data(out: &mut String, tag: &str, keys: &[(String, &Column)], row: usize) {
    let mut open = false;
    for (key, col) in keys {
        if let Some(t) = col.get(row).to_text() {
            if !open {
                out.push_str(">\n");
                open = true;
            }
            out.push_str("      <data key=\"");
            out.push_str(key);
            out.push_str("\">");
            escape_into(out, &t);
            out.push_str("</data>\n");
        }
    }
    if open {
        out.push_str("    </");
        out.push_str(tag);
        out.push_str(">\n");
    } else {
        out.push_str("/>\n");
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::table::Column;
    use crate::rng::Stream;

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

    #[test]
    fn streaming_handles_text_split_by_markup() {
        // Text broken by comments, CDATA, entities, and child elements, as
        // yEd and other editors write it.
        let doc = DOC.replace(
            "<data key=\"d0\">red</data>",
            "<data key=\"d0\">r<!-- x -->e<![CDATA[d]]>&amp;<y:Shape>skip</y:Shape>!</data>",
        );
        let n = read_stream(&doc).unwrap().unwrap();
        assert_eq!(
            n.nodes.column("color").unwrap().get(0),
            Value::Str("red&!".into())
        );
        assert_eq!(read_tree(&doc).unwrap(), n);
    }

    #[test]
    fn keys_after_the_graph_use_the_tree_reader() {
        let doc = DOC.replace(
            "</graph>\n",
            "</graph>\n<key id=\"late\" for=\"node\" attr.name=\"late\"/>\n",
        );
        assert!(read_stream(&doc).unwrap().is_none());
        assert!(read(&doc).unwrap().nodes.column("late").is_some());
    }

    /// Random edits of GraphML documents: on every document the streaming
    /// reader accepts, it gives exactly the tree-based reader's network,
    /// and it never accepts a document the tree-based reader rejects.
    #[test]
    fn streaming_agrees_with_the_tree_reader() {
        let docs = [
            DOC.to_string(),
            DOC.replace("undirected", "directed"),
            r#"<graphml><key id="k" attr.name="x" attr.type="boolean"><default>1</default><default>0</default></key>
<key id="k" for="edge" attr.name="w" attr.type="float"/><graph>
<node id="1"><data key="k">0</data><data key="k">yes</data></node>
<edge source="1" target="2" directed="true"><data key="k">2.5e3</data><port/></edge>
<data key="k">true</data><x><node id="ignored"/></x></graph><graph><node id="z"/></graph></graphml>"#
                .to_string(),
        ];
        const PIECES: [&str; 24] = [
            "<data key=\"d0\">x</data>",
            "<data key=\"d2\">7</data>",
            "<data key=\"d1\">nan</data>",
            "<data key=\"g\">t</data>",
            "<data>",
            "</data>",
            "<node id=\"a\"/>",
            "<node id=\"q\">",
            "</node>",
            "<edge source=\"a\" target=\"q\"/>",
            "<edge source=\"a\">",
            "</edge>",
            "<key id=\"d9\" for=\"all\" attr.name=\"color\"/>",
            "<default>5</default>",
            "<graph>",
            "</graph>",
            "<port/>",
            "<hyperedge/>",
            " directed=\"false\"",
            " id=\"e\"",
            "&amp;",
            "<![CDATA[z]]>",
            "<!-- c -->",
            "\n",
        ];
        let mut accepted = 0;
        for case in 0..20_000u64 {
            let mut rng = Stream::new(11, case);
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
            // Whatever path is taken, `read` is the tree reader's answer.
            assert_eq!(format!("{:?}", read(&s)), tree, "case {case}");
        }
        assert!(accepted > 1_000, "only {accepted} documents were accepted");
    }
}
