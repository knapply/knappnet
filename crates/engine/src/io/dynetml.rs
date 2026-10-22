//! DyNetML, the XML format of ORA (CASOS, Carnegie Mellon University) for
//! multimode "meta-networks" (Tsvetovat, Reminga, and Carley 2004,
//! `tsvetovat2004dynetml`; ORA User's Guide, Altman, Carley, and Reminga
//! 2020, `altman2020ora`).
//!
//! # Structure
//!
//! A meta-network holds node classes (sets of agents, organizations,
//! locations, events, resources, tasks, knowledge, and so on) and several
//! networks among them. A dynamic meta-network holds one meta-network per
//! time period. Two generations of element names are read:
//!
//! | Role | ORA (current) | 2004 report |
//! |---|---|---|
//! | dynamic meta-network (root) | `DynamicMetaNetwork`, `DynamicNetwork` | `DynamicNetwork` |
//! | meta-network | `MetaNetwork` (`id`) | `MetaMatrix` (`timePeriod`) |
//! | node class | `nodeclass` (`id`, `type`) | `nodeset` (`id`, `type`) |
//! | network | `network` | `graph` |
//! | link | `link` (`source`, `target`, `value`) | `edge` (`source`, `target`, `type`, `value`) |
//! | attribute | `property` (`id`, `value`), typed by `propertyIdentity` (`id`, `type`, `singleValued`) | `property` / `measure` (`name`, `type`, `value`) |
//!
//! A network names its node classes with `source` and `target` (the node
//! class `id`) and `sourceType` and `targetType` (the node class `type`,
//! such as `Agent`); either is enough. The root may also be a single
//! `MetaNetwork` or `MetaMatrix`.
//!
//! # Mapping onto tables
//!
//! Nodes: `id`, then `class` (the node class `id`), `type` (the node class
//! `type`), `title` when any node has one, and one column per node
//! property or measure. Other XML attributes of `<node>` become text
//! columns. Property types come from `propertyIdentities`: `number` and
//! `numberCategory` are read as doubles; `text`, `categoryText`, `URI`, and
//! `date` as text (dates are kept verbatim). The 2004 `type` attribute is
//! honored the same way (`double`, `string`, `binary` as logical). A
//! property declared with different types in different node classes is
//! read as text, and an undeclared property is text. An empty value is
//! missing. Names that collide with the structural columns get an `_attr`
//! suffix.
//!
//! Node identifiers: ORA identifies a node by its class and its `id`, so
//! the same `id` may appear in two classes. When it does, every node
//! identifier becomes `class::id`; when the file holds several
//! meta-networks, every identifier is prefixed with the period,
//! `period::id` (or `period::class::id`). Only then, the original `id` is
//! kept in a `node_id` column. Files with unique identifiers keep them
//! unchanged.
//!
//! Edges: `from`, `to`, `network` (the network `id`), and `weight` (the
//! link `value` as a double) when any link has a value; a link without a
//! value then gets weight 1, as for binary links. A link value of 0 is kept
//! as an edge with weight 0. Other link attributes and link properties
//! become edge columns. Links whose value is text are rejected, since the
//! weight must be numeric.
//!
//! Direction: each network records `isDirected` (missing means directed).
//! If every network is undirected, the result is undirected. If every
//! network is directed, it is directed. A file that mixes the two is read
//! as directed, with each undirected link stored as two opposite edges
//! (self-loops once) and a logical `directed` edge column marking which
//! network each edge came from, as Pajek files mixing `*Arcs` and `*Edges`
//! are read.
//!
//! Several meta-networks: node and edge tables gain a `period` column (the
//! `MetaNetwork` `id`, the 2004 `timePeriod`, or the 1-based position when
//! neither is given), and every period's nodes are separate rows, so the
//! result is the disjoint union of the periods. Edges never cross periods.
//!
//! Graph attributes: `name` is the meta-network `id` (or, with several
//! periods, the root `id`). Meta-network properties become graph
//! attributes; with several periods they are named `period::property`.
//! Other XML attributes of the root and meta-network elements (such as
//! namespace declarations) are ignored.
//!
//! Not supported, rejected with an error naming the element and line:
//! ports (`port`, `sourcePort`, `targetPort`), AnthroPac survey sections,
//! properties or measures on networks or node classes, repeated values of
//! one property on a node, two networks or node classes with the same `id`
//! in one meta-network, and any element not listed above.
//!
//! # Writing
//!
//! Writes a `DynamicNetwork` root holding one `MetaNetwork` per distinct
//! value of the node `period` column (one meta-network without it), the
//! structure of the DyNetML files distributed with ORA. Node classes come
//! from the `class` and `type` columns: a missing class is `Agent`, and a
//! missing type is the class name when it is one of ORA's node class types
//! (`Agent`, `Organization`, `Knowledge`, `Resource`, `Belief`, `Event`,
//! `Task`, `Location`, `Role`, `Action`), otherwise `Unknown`. The `title`
//! column becomes the node title, `node_id` (when present) the XML node
//! `id`, and other columns node properties with `propertyIdentities`
//! (numbers as `number`, text as `text`, logical as `categoryText`, so
//! logical columns read back as text and integers as doubles). Links are
//! grouped into one network per `network` value; without that column, one
//! network per pair of node classes, named `source x target` as ORA does.
//! The link value is the weight column. Graph attributes become
//! meta-network properties and `name` the meta-network `id` (the root `id`
//! with several periods, where a graph attribute `period::x` becomes
//! property `x` of that period and other graph attributes are written in
//! every period). In a directed network, an edge
//! whose `directed` column is false belongs to an undirected network; the
//! opposite copies made when reading are written once. Only networks with
//! links are written, except that an undirected network without links is
//! written with one empty undirected network so its direction survives.
//! Other edge columns
//! are written as link properties with a `type` attribute; whether ORA
//! reads link properties has not been checked.

use std::collections::{HashMap, HashSet};

use super::table::{
    format_float, parse_bool, parse_float, Column, Kind, Network, NetworkBuilder, Value,
};
use super::xml::{self, escape, Element};
use crate::error::{Error, Result};

const FORMAT: &str = "DyNetML";

/// Root element names that identify a DyNetML document.
pub const ROOTS: [&str; 4] = [
    "DynamicMetaNetwork",
    "DynamicNetwork",
    "MetaNetwork",
    "MetaMatrix",
];

/// ORA's node class types (ORA User's Guide 2020, `altman2020ora`).
const ORA_TYPES: [&str; 11] = [
    "Agent",
    "Organization",
    "Knowledge",
    "Resource",
    "Belief",
    "Event",
    "Task",
    "Location",
    "Role",
    "Action",
    "Unknown",
];

fn err<T>(el: &Element, message: impl Into<String>) -> Result<T> {
    Err(Error::Parse {
        format: FORMAT,
        line: Some(el.line),
        message: message.into(),
    })
}

fn err_at<T>(line: usize, message: impl Into<String>) -> Result<T> {
    Err(Error::Parse {
        format: FORMAT,
        line: Some(line),
        message: message.into(),
    })
}

/// The column kind for a DyNetML property type name.
fn kind_of(t: &str) -> Kind {
    match t.to_ascii_lowercase().as_str() {
        "number" | "numbercategory" | "double" | "float" => Kind::Float,
        "int" | "integer" | "long" => Kind::Int,
        "binary" | "boolean" | "bool" => Kind::Bool,
        _ => Kind::Str,
    }
}

/// A value as written in the file, converted once column types are known.
struct Raw {
    name: String,
    text: String,
    kind: Option<Kind>,
    line: usize,
}

struct RawNode {
    period: usize,
    class: usize,
    id: String,
    attrs: Vec<Raw>,
}

struct RawLink {
    period: usize,
    network: String,
    from: usize,
    to: usize,
    directed: bool,
    value: Option<f64>,
    attrs: Vec<Raw>,
}

struct NodeClass {
    id: String,
    ty: String,
}

struct Period {
    id: String,
    props: Vec<Raw>,
    classes: Vec<NodeClass>,
    /// (class index, node id) -> node index.
    nodes: HashMap<(usize, String), usize>,
    /// Node id -> node indices in every class.
    by_id: HashMap<String, Vec<usize>>,
}

#[derive(Default)]
struct Doc {
    periods: Vec<Period>,
    nodes: Vec<RawNode>,
    links: Vec<RawLink>,
    /// Declared direction of every network, including empty ones.
    directions: Vec<bool>,
    /// (period, network id) of every network.
    network_ids: Vec<(usize, String)>,
}

fn identities(el: &Element) -> Result<HashMap<String, Kind>> {
    let mut out = HashMap::new();
    for p in &el.children {
        if p.name != "propertyIdentity" {
            return err(
                p,
                format!("unexpected <{}> in <propertyIdentities>", p.name),
            );
        }
        let id = match p.attr("id").or_else(|| p.attr("name")) {
            Some(id) => id,
            None => return err(p, "<propertyIdentity> has no id"),
        };
        out.insert(id.to_string(), kind_of(p.attr("type").unwrap_or("text")));
    }
    Ok(out)
}

/// Reads `<property>` or `<measure>` elements (both carry a name, an
/// optional type, and a value) into `out`, rejecting repeated names.
fn read_values(
    el: &Element,
    item: &str,
    declared: &HashMap<String, Kind>,
    out: &mut Vec<Raw>,
) -> Result<()> {
    if el.name != item {
        return err(el, format!("unexpected <{}>; expected <{item}>", el.name));
    }
    let name = match el.attr("id").or_else(|| el.attr("name")) {
        Some(n) => n,
        None => return err(el, format!("<{item}> has no id or name")),
    };
    let text = el
        .attr("value")
        .map(str::to_string)
        .unwrap_or_else(|| el.text.trim().to_string());
    if out.iter().any(|r| r.name == name) {
        return err(
            el,
            format!("{item} {name:?} appears more than once; repeated values are not supported"),
        );
    }
    let kind = el
        .attr("type")
        .map(kind_of)
        .or_else(|| declared.get(name).copied());
    out.push(Raw {
        name: name.to_string(),
        text,
        kind,
        line: el.line,
    });
    Ok(())
}

/// Reads the property and measure children of `el` (a node or link).
fn read_attached(el: &Element, declared: &HashMap<String, Kind>, out: &mut Vec<Raw>) -> Result<()> {
    for c in &el.children {
        match c.name.as_str() {
            "properties" => {
                for p in &c.children {
                    read_values(p, "property", declared, out)?;
                }
            }
            "measures" => {
                for p in &c.children {
                    read_values(p, "measure", declared, out)?;
                }
            }
            "property" => read_values(c, "property", declared, out)?,
            "port" => return err(c, "ports are not supported"),
            n if n.starts_with("anthropac") => {
                return err(c, format!("AnthroPac survey data (<{n}>) is not supported"))
            }
            n => return err(c, format!("unsupported element <{n}> in <{}>", el.name)),
        }
    }
    Ok(())
}

impl Doc {
    fn find_class(&self, p: usize, id: &str) -> Option<usize> {
        self.periods[p].classes.iter().position(|c| c.id == id)
    }

    fn add_class(&mut self, p: usize, id: &str, ty: &str) -> usize {
        let classes = &mut self.periods[p].classes;
        classes.push(NodeClass {
            id: id.to_string(),
            ty: ty.to_string(),
        });
        classes.len() - 1
    }

    fn add_node(&mut self, p: usize, class: usize, id: &str, attrs: Vec<Raw>) -> usize {
        let idx = self.nodes.len();
        self.nodes.push(RawNode {
            period: p,
            class,
            id: id.to_string(),
            attrs,
        });
        self.periods[p].nodes.insert((class, id.to_string()), idx);
        self.periods[p]
            .by_id
            .entry(id.to_string())
            .or_default()
            .push(idx);
        idx
    }

    /// The node class a network end refers to: by class id (`source`),
    /// then by `sourceType` as a class id, then by a unique class of that
    /// type. A class named but not declared is created.
    fn resolve_class(&mut self, p: usize, set: Option<&str>, ty: Option<&str>) -> Option<usize> {
        if let Some(s) = set {
            return Some(match self.find_class(p, s) {
                Some(c) => c,
                None => self.add_class(p, s, ty.unwrap_or(s)),
            });
        }
        let t = ty?;
        if let Some(c) = self.find_class(p, t) {
            return Some(c);
        }
        let matches: Vec<usize> = (0..self.periods[p].classes.len())
            .filter(|&i| self.periods[p].classes[i].ty.eq_ignore_ascii_case(t))
            .collect();
        match matches.len() {
            0 => Some(self.add_class(p, t, t)),
            1 => Some(matches[0]),
            // Several classes of this type: endpoints are looked up by id.
            _ => None,
        }
    }

    fn endpoint(
        &mut self,
        p: usize,
        class: Option<usize>,
        id: &str,
        link: &Element,
    ) -> Result<usize> {
        if let Some(c) = class {
            if let Some(&n) = self.periods[p].nodes.get(&(c, id.to_string())) {
                return Ok(n);
            }
            return Ok(self.add_node(p, c, id, Vec::new()));
        }
        let found: &[usize] = self.periods[p].by_id.get(id).map_or(&[], |v| v);
        match found.len() {
            1 => Ok(found[0]),
            0 => err(
                link,
                format!("node {id:?} is not declared and the network does not name its node class"),
            ),
            _ => err(
                link,
                format!(
                    "node {id:?} is in several node classes and the network does not say which"
                ),
            ),
        }
    }

    fn read_meta_network(&mut self, el: &Element, position: usize) -> Result<()> {
        let id = el
            .attr("id")
            .or_else(|| el.attr("timePeriod"))
            .map(str::to_string)
            .unwrap_or_else(|| (position + 1).to_string());
        if self.periods.iter().any(|p| p.id == id) {
            return err(el, format!("meta-network {id:?} appears more than once"));
        }
        let p = self.periods.len();
        self.periods.push(Period {
            id,
            props: Vec::new(),
            classes: Vec::new(),
            nodes: HashMap::new(),
            by_id: HashMap::new(),
        });
        let mut declared = HashMap::new();
        for c in el.children_named("propertyIdentities") {
            declared.extend(identities(c)?);
        }
        for c in &el.children {
            match c.name.as_str() {
                "propertyIdentities" | "nodes" | "networks" => {}
                "properties" | "measures" => {
                    let item = if c.name == "properties" {
                        "property"
                    } else {
                        "measure"
                    };
                    let mut props = std::mem::take(&mut self.periods[p].props);
                    for v in &c.children {
                        read_values(v, item, &declared, &mut props)?;
                    }
                    self.periods[p].props = props;
                }
                n if n.starts_with("anthropac") => {
                    return err(c, format!("AnthroPac survey data (<{n}>) is not supported"))
                }
                n => return err(c, format!("unsupported element <{n}> in <{}>", el.name)),
            }
        }
        for nodes in el.children_named("nodes") {
            for class in &nodes.children {
                if class.name != "nodeclass" && class.name != "nodeset" {
                    return err(
                        class,
                        format!("unsupported element <{}> in <nodes>", class.name),
                    );
                }
                self.read_class(p, class, &declared)?;
            }
        }
        for networks in el.children_named("networks") {
            for net in &networks.children {
                if net.name != "network" && net.name != "graph" {
                    return err(
                        net,
                        format!("unsupported element <{}> in <networks>", net.name),
                    );
                }
                self.read_network(p, net)?;
            }
        }
        Ok(())
    }

    fn read_class(&mut self, p: usize, el: &Element, outer: &HashMap<String, Kind>) -> Result<()> {
        let (id, ty) = match (el.attr("id"), el.attr("type")) {
            (Some(i), Some(t)) => (i, t),
            (Some(i), None) => (i, i),
            (None, Some(t)) => (t, t),
            (None, None) => return err(el, format!("<{}> has neither id nor type", el.name)),
        };
        if self.find_class(p, id).is_some() {
            return err(el, format!("node class {id:?} appears more than once"));
        }
        let c = self.add_class(p, id, ty);
        let mut declared = outer.clone();
        for child in el.children_named("propertyIdentities") {
            declared.extend(identities(child)?);
        }
        for node in &el.children {
            match node.name.as_str() {
                "propertyIdentities" => continue,
                "node" => {}
                n => return err(node, format!("unsupported element <{n}> in <{}>", el.name)),
            }
            let nid = match node.attr("id") {
                Some(i) => i,
                None => return err(node, "<node> has no id"),
            };
            if self.periods[p].nodes.contains_key(&(c, nid.to_string())) {
                return err(
                    node,
                    format!("node {nid:?} appears twice in node class {id:?}"),
                );
            }
            let mut attrs: Vec<Raw> = node
                .attrs
                .iter()
                .filter(|(k, _)| k != "id")
                .map(|(k, v)| Raw {
                    name: k.clone(),
                    text: v.clone(),
                    kind: Some(Kind::Str),
                    line: node.line,
                })
                .collect();
            read_attached(node, &declared, &mut attrs)?;
            self.add_node(p, c, nid, attrs);
        }
        Ok(())
    }

    fn read_network(&mut self, p: usize, el: &Element) -> Result<()> {
        let directed = match el.attr("isDirected") {
            None => true,
            Some(d) => match parse_bool(d) {
                Some(b) => b,
                None => return err(el, format!("isDirected must be true or false, not {d:?}")),
            },
        };
        let source = self.resolve_class(p, el.attr("source"), el.attr("sourceType"));
        let target = self.resolve_class(p, el.attr("target"), el.attr("targetType"));
        let class_name = |d: &Doc, c: Option<usize>| {
            c.map(|c| d.periods[p].classes[c].id.clone())
                .unwrap_or_default()
        };
        let id = match el.attr("id") {
            Some(i) => i.to_string(),
            None => format!(
                "{} x {}",
                class_name(self, source),
                class_name(self, target)
            ),
        };
        if self.network_ids.contains(&(p, id.clone())) {
            return err(el, format!("network {id:?} appears more than once"));
        }
        self.network_ids.push((p, id.clone()));
        self.directions.push(directed);
        for link in &el.children {
            match link.name.as_str() {
                "link" | "edge" => {}
                "properties" | "measures" => {
                    return err(
                        link,
                        format!("{} on <{}> are not supported", link.name, el.name),
                    )
                }
                n => return err(link, format!("unsupported element <{n}> in <{}>", el.name)),
            }
            for port in ["sourcePort", "targetPort"] {
                if link.attr(port).is_some() {
                    return err(link, "ports are not supported");
                }
            }
            let (s, t) = match (link.attr("source"), link.attr("target")) {
                (Some(s), Some(t)) => (s, t),
                _ => return err(link, format!("<{}> needs source and target", link.name)),
            };
            let from = self.endpoint(p, source, s, link)?;
            let to = self.endpoint(p, target, t, link)?;
            let value = match link.attr("value") {
                None => None,
                Some(v) => {
                    match parse_float(v).or_else(|| parse_bool(v).map(|b| b as i64 as f64)) {
                        Some(x) => Some(x),
                        None => {
                            return err(
                                link,
                                format!(
                            "link value {v:?} is not a number; text-valued links are not supported"
                        ),
                            )
                        }
                    }
                }
            };
            let mut attrs: Vec<Raw> = link
                .attrs
                .iter()
                .filter(|(k, _)| !matches!(k.as_str(), "source" | "target" | "value" | "type"))
                .map(|(k, v)| Raw {
                    name: k.clone(),
                    text: v.clone(),
                    kind: Some(Kind::Str),
                    line: link.line,
                })
                .collect();
            read_attached(link, &HashMap::new(), &mut attrs)?;
            self.links.push(RawLink {
                period: p,
                network: id.clone(),
                from,
                to,
                directed,
                value,
                attrs,
            });
        }
        Ok(())
    }
}

/// Column kinds: the declared kind when every occurrence agrees, text
/// otherwise.
fn column_kinds<'a>(raws: impl Iterator<Item = &'a Raw>) -> HashMap<String, Kind> {
    let mut kinds: HashMap<String, Option<Kind>> = HashMap::new();
    for r in raws {
        let k = r.kind.unwrap_or(Kind::Str);
        kinds
            .entry(r.name.clone())
            .and_modify(|e| {
                if *e != Some(k) {
                    *e = None;
                }
            })
            .or_insert(Some(k));
    }
    kinds
        .into_iter()
        .map(|(n, k)| (n, k.unwrap_or(Kind::Str)))
        .collect()
}

fn convert(raw: &Raw, kind: Kind) -> Result<Value> {
    if kind != Kind::Str && raw.text.trim().is_empty() {
        // ORA leaves numeric properties empty when unknown.
        return Ok(Value::Null);
    }
    let v = match kind {
        Kind::Str => Some(Value::Str(raw.text.clone())),
        Kind::Bool => parse_bool(&raw.text).map(Value::Bool),
        Kind::Int => raw.text.trim().parse().ok().map(Value::Int),
        Kind::Float => parse_float(&raw.text).map(Value::Float),
    };
    match v {
        Some(v) => Ok(v),
        None => err_at(
            raw.line,
            format!(
                "value {:?} of {:?} is not a valid {kind:?}",
                raw.text, raw.name
            ),
        ),
    }
}

fn rename(name: &str, reserved: &[&str]) -> String {
    if reserved.contains(&name) {
        format!("{name}_attr")
    } else {
        name.to_string()
    }
}

/// Reads a DyNetML document.
pub fn read(text: &str) -> Result<Network> {
    let root = xml::parse(text, FORMAT)?;
    let mut doc = Doc::default();
    let dynamic = match root.name.as_str() {
        "DynamicMetaNetwork" | "DynamicNetwork" => true,
        "MetaNetwork" | "MetaMatrix" => false,
        other => {
            return err(
                &root,
                format!("expected <DynamicMetaNetwork>, <DynamicNetwork>, or <MetaNetwork>, found <{other}>"),
            )
        }
    };
    if dynamic {
        for (i, mn) in root.children.iter().enumerate() {
            match mn.name.as_str() {
                "MetaNetwork" | "MetaMatrix" => doc.read_meta_network(mn, i)?,
                n => return err(mn, format!("unsupported element <{n}> in <{}>", root.name)),
            }
        }
        if doc.periods.is_empty() {
            return err(&root, format!("<{}> holds no meta-network", root.name));
        }
    } else {
        doc.read_meta_network(&root, 0)?;
    }
    build(doc, &root, dynamic)
}

fn build(doc: Doc, root: &Element, dynamic: bool) -> Result<Network> {
    let multi = doc.periods.len() > 1;
    let mut seen: HashSet<(usize, &str)> = HashSet::new();
    let mut shared = false;
    for n in &doc.nodes {
        if !seen.insert((n.period, n.id.as_str())) {
            shared = true;
            break;
        }
    }
    let qualified = multi || shared;
    let key = |n: &RawNode| -> String {
        let mut parts: Vec<&str> = Vec::new();
        if multi {
            parts.push(&doc.periods[n.period].id);
        }
        if shared {
            parts.push(&doc.periods[n.period].classes[n.class].id);
        }
        parts.push(&n.id);
        parts.join("::")
    };
    let keys: Vec<String> = doc.nodes.iter().map(key).collect();

    let any_directed = doc.directions.iter().any(|&d| d);
    let any_undirected = doc.directions.iter().any(|&d| !d);
    let directed = any_directed || !any_undirected;
    let mixed = any_directed && any_undirected;
    let mut b = NetworkBuilder::new(directed);

    // Graph attributes.
    let meta_id = if dynamic {
        root.children.first()
    } else {
        Some(root)
    }
    .and_then(|m| m.attr("id"));
    let name = if multi {
        root.attr("id")
    } else {
        meta_id.or_else(|| root.attr("id").filter(|_| dynamic))
    };
    if let Some(n) = name.map(str::to_string) {
        b.graph.push(("name".into(), Value::Str(n)));
    }
    for period in &doc.periods {
        let kinds = column_kinds(period.props.iter());
        for raw in &period.props {
            let value = convert(raw, kinds[&raw.name])?;
            let mut n = if multi {
                format!("{}::{}", period.id, raw.name)
            } else {
                raw.name.clone()
            };
            if b.graph.iter().any(|(g, _)| *g == n) {
                n.push_str("_attr");
            }
            b.graph.push((n, value));
        }
    }

    // Node columns.
    let reserved_nodes = ["id", "node_id", "class", "type", "period"];
    if qualified {
        b.node_kind("node_id", Kind::Str);
    }
    b.node_kind("class", Kind::Str);
    b.node_kind("type", Kind::Str);
    if multi {
        b.node_kind("period", Kind::Str);
    }
    let node_kinds = column_kinds(doc.nodes.iter().flat_map(|n| n.attrs.iter()));
    if doc
        .nodes
        .iter()
        .any(|n| n.attrs.iter().any(|a| a.name == "title"))
    {
        b.node_kind("title", node_kinds["title"]);
    }
    for n in &doc.nodes {
        for a in &n.attrs {
            b.node_kind(&rename(&a.name, &reserved_nodes), node_kinds[&a.name]);
        }
    }
    for (n, k) in doc.nodes.iter().zip(&keys) {
        let class = &doc.periods[n.period].classes[n.class];
        let mut row = Vec::new();
        if qualified {
            row.push(("node_id".to_string(), Value::Str(n.id.clone())));
        }
        row.push(("class".into(), Value::Str(class.id.clone())));
        row.push(("type".into(), Value::Str(class.ty.clone())));
        if multi {
            row.push((
                "period".into(),
                Value::Str(doc.periods[n.period].id.clone()),
            ));
        }
        for a in &n.attrs {
            row.push((
                rename(&a.name, &reserved_nodes),
                convert(a, node_kinds[&a.name])?,
            ));
        }
        if !b.add_node(k, row) {
            return err(
                root,
                format!("node identifier {k:?} is ambiguous; rename nodes so identifiers do not contain \"::\""),
            );
        }
    }

    // Edge columns.
    let reserved_edges = ["from", "to", "network", "period", "directed", "weight"];
    b.edge_kind("network", Kind::Str);
    if multi {
        b.edge_kind("period", Kind::Str);
    }
    if mixed {
        b.edge_kind("directed", Kind::Bool);
    }
    let weighted = doc.links.iter().any(|l| l.value.is_some());
    if weighted {
        b.edge_kind("weight", Kind::Float);
    }
    let edge_kinds = column_kinds(doc.links.iter().flat_map(|l| l.attrs.iter()));
    for l in &doc.links {
        for a in &l.attrs {
            b.edge_kind(&rename(&a.name, &reserved_edges), edge_kinds[&a.name]);
        }
    }
    for l in &doc.links {
        let mut row = vec![("network".to_string(), Value::Str(l.network.clone()))];
        if multi {
            row.push((
                "period".into(),
                Value::Str(doc.periods[l.period].id.clone()),
            ));
        }
        if mixed {
            row.push(("directed".into(), Value::Bool(l.directed)));
        }
        if weighted {
            row.push(("weight".into(), Value::Float(l.value.unwrap_or(1.0))));
        }
        for a in &l.attrs {
            row.push((
                rename(&a.name, &reserved_edges),
                convert(a, edge_kinds[&a.name])?,
            ));
        }
        let (f, t) = (&keys[l.from], &keys[l.to]);
        b.add_edge(f, t, row.clone());
        if mixed && !l.directed && l.from != l.to {
            b.add_edge(t, f, row);
        }
    }
    b.finish()
}

fn text_column(col: Option<&Column>, row: usize) -> Option<String> {
    col.and_then(|c| c.get(row).to_text())
}

fn type_name(kind: Kind) -> &'static str {
    match kind {
        Kind::Int | Kind::Float => "number",
        Kind::Bool => "categoryText",
        Kind::Str => "text",
    }
}

fn default_type(class: &str) -> String {
    ORA_TYPES
        .iter()
        .find(|t| t.eq_ignore_ascii_case(class))
        .map(|t| t.to_string())
        .unwrap_or_else(|| "Unknown".to_string())
}

struct OutNetwork {
    period: usize,
    id: String,
    source: usize,
    target: usize,
    directed: bool,
    rows: Vec<usize>,
}

/// Writes a DyNetML document. `weight` names the edge column written as
/// the link value (default `weight`).
pub fn write(net: &Network, weight: Option<&str>) -> Result<String> {
    net.validate()?;
    let n = net.nodes.nrow();
    let ids = net.nodes.strings("id")?;
    let col = |name: &str| net.nodes.column(name);
    let invalid = |m: String| Error::InvalidArgument(m);

    // Periods and node classes, in order of first appearance.
    let period_col = col("period");
    let mut periods: Vec<String> = Vec::new();
    let mut node_period = vec![0usize; n];
    for (row, slot) in node_period.iter_mut().enumerate() {
        if let Some(c) = period_col {
            let p = c
                .get(row)
                .to_text()
                .ok_or_else(|| invalid(format!("node {:?} has no period", ids[row])))?;
            *slot = match periods.iter().position(|q| *q == p) {
                Some(i) => i,
                None => {
                    periods.push(p);
                    periods.len() - 1
                }
            };
        }
    }
    if periods.is_empty() {
        periods.push(String::new());
    }
    // (period, class id, type)
    let mut classes: Vec<(usize, String, String)> = Vec::new();
    let mut node_class = vec![0usize; n];
    let mut xml_ids: HashSet<(usize, String)> = HashSet::new();
    let xml_id: Vec<String> = (0..n)
        .map(|row| text_column(col("node_id"), row).unwrap_or_else(|| ids[row].clone()))
        .collect();
    for row in 0..n {
        let class = text_column(col("class"), row).unwrap_or_else(|| "Agent".to_string());
        let ty = text_column(col("type"), row).unwrap_or_else(|| default_type(&class));
        let p = node_period[row];
        let c = match classes.iter().position(|(q, c, _)| *q == p && *c == class) {
            Some(i) => {
                if classes[i].2 != ty {
                    return Err(invalid(format!(
                        "node class {class:?} has two types, {:?} and {ty:?}",
                        classes[i].2
                    )));
                }
                i
            }
            None => {
                classes.push((p, class.clone(), ty));
                classes.len() - 1
            }
        };
        node_class[row] = c;
        if !xml_ids.insert((c, xml_id[row].clone())) {
            return Err(invalid(format!(
                "node id {:?} appears twice in node class {class:?}",
                xml_id[row]
            )));
        }
    }

    let special_nodes = ["id", "node_id", "class", "type", "period", "title"];
    let props: Vec<(&String, &Column)> = net
        .node_attributes()
        .filter(|(name, _)| !special_nodes.contains(&name.as_str()))
        .collect();

    // Networks.
    let from = net.edges.strings("from")?;
    let to = net.edges.strings("to")?;
    let index = net.node_index()?;
    let weight_name = weight.unwrap_or("weight");
    let weight_col = net
        .edges
        .column(weight_name)
        .filter(|c| matches!(c.kind(), Kind::Float | Kind::Int));
    let directed_col = net
        .edges
        .column("directed")
        .filter(|c| c.kind() == Kind::Bool);
    let network_col = net.edges.column("network");
    let edge_props: Vec<(&String, &Column)> = net
        .edge_attributes()
        .filter(|(name, _)| {
            let name = name.as_str();
            !(name == "network"
                || (name == "period" && period_col.is_some())
                || (name == "directed" && directed_col.is_some() && net.directed)
                || (name == weight_name && weight_col.is_some()))
        })
        .collect();
    let mut networks: Vec<OutNetwork> = Vec::new();
    for row in 0..from.len() {
        let (a, b) = (index[&from[row]], index[&to[row]]);
        let p = node_period[a];
        if node_period[b] != p {
            return Err(invalid(format!(
                "edge {:?} -> {:?} joins two periods",
                from[row], to[row]
            )));
        }
        let (ca, cb) = (node_class[a], node_class[b]);
        let id = text_column(network_col, row)
            .unwrap_or_else(|| format!("{} x {}", classes[ca].1, classes[cb].1));
        let directed =
            net.directed && !matches!(directed_col.map(|c| c.get(row)), Some(Value::Bool(false)));
        match networks.iter_mut().find(|w| w.period == p && w.id == id) {
            Some(w) => {
                // Links of an undirected network may run either way.
                let reversed = !w.directed && w.source == cb && w.target == ca;
                if (w.source != ca || w.target != cb) && !reversed {
                    return Err(invalid(format!(
                        "network {id:?} links more than one pair of node classes"
                    )));
                }
                if w.directed != directed {
                    return Err(invalid(format!(
                        "network {id:?} mixes directed and undirected edges"
                    )));
                }
                w.rows.push(row);
            }
            None => networks.push(OutNetwork {
                period: p,
                id,
                source: ca,
                target: cb,
                directed,
                rows: vec![row],
            }),
        }
    }

    // Direction is recorded on networks only, so an undirected network
    // without links is written with one empty undirected network.
    if networks.is_empty() && !net.directed {
        if classes.is_empty() {
            classes.push((0, "Agent".into(), "Agent".into()));
        }
        networks.push(OutNetwork {
            period: classes[0].0,
            id: format!("{} x {}", classes[0].1, classes[0].1),
            source: 0,
            target: 0,
            directed: false,
            rows: Vec::new(),
        });
    }

    let mut out = String::from("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n");
    let name = net
        .graph
        .iter()
        .find(|(k, _)| k == "name")
        .and_then(|(_, v)| v.to_text());
    let multi = periods.len() > 1;
    match (&name, multi) {
        (Some(nm), true) => out.push_str(&format!("<DynamicNetwork id=\"{}\">\n", escape(nm))),
        _ => out.push_str("<DynamicNetwork>\n"),
    }
    for (p, period) in periods.iter().enumerate() {
        let mn_id = if multi {
            Some(period.clone())
        } else {
            name.clone()
        };
        match mn_id {
            Some(i) => out.push_str(&format!("  <MetaNetwork id=\"{}\">\n", escape(&i))),
            None => out.push_str("  <MetaNetwork>\n"),
        }
        // Meta-network properties from graph attributes.
        let prefix = format!("{period}::");
        let graph_props: Vec<(String, &Value)> = net
            .graph
            .iter()
            .filter(|(k, _)| k != "name")
            .filter_map(|(k, v)| {
                if !multi {
                    return Some((k.clone(), v));
                }
                match k.strip_prefix(&prefix) {
                    Some(rest) => Some((rest.to_string(), v)),
                    None if !periods.iter().any(|q| k.starts_with(&format!("{q}::"))) => {
                        Some((k.clone(), v))
                    }
                    None => None,
                }
            })
            .filter(|(_, v)| v.to_text().is_some())
            .collect();
        if !graph_props.is_empty() {
            out.push_str("    <propertyIdentities>\n");
            for (k, v) in &graph_props {
                let kind = match v {
                    Value::Bool(_) => Kind::Bool,
                    Value::Int(_) => Kind::Int,
                    Value::Float(_) => Kind::Float,
                    _ => Kind::Str,
                };
                out.push_str(&format!(
                    "      <propertyIdentity id=\"{}\" type=\"{}\" singleValued=\"true\"/>\n",
                    escape(k),
                    type_name(kind)
                ));
            }
            out.push_str("    </propertyIdentities>\n    <properties>\n");
            for (k, v) in &graph_props {
                out.push_str(&format!(
                    "      <property id=\"{}\" value=\"{}\"/>\n",
                    escape(k),
                    escape(&v.to_text().unwrap_or_default())
                ));
            }
            out.push_str("    </properties>\n");
        }

        out.push_str("    <nodes>\n");
        for (c, (cp, class, ty)) in classes.iter().enumerate() {
            if *cp != p {
                continue;
            }
            let rows: Vec<usize> = (0..n).filter(|&r| node_class[r] == c).collect();
            // Declare the properties this class uses, and those no node
            // uses at all, so empty columns survive.
            let declared: Vec<&(&String, &Column)> = props
                .iter()
                .filter(|(_, column)| {
                    rows.iter().any(|&r| column.get(r).to_text().is_some())
                        || (0..n).all(|r| column.get(r).to_text().is_none())
                })
                .collect();
            out.push_str(&format!(
                "      <nodeclass type=\"{}\" id=\"{}\">\n",
                escape(ty),
                escape(class)
            ));
            if !declared.is_empty() {
                out.push_str("        <propertyIdentities>\n");
                for (name, column) in &declared {
                    out.push_str(&format!(
                        "          <propertyIdentity id=\"{}\" type=\"{}\" singleValued=\"true\"/>\n",
                        escape(name),
                        type_name(column.kind())
                    ));
                }
                out.push_str("        </propertyIdentities>\n");
            }
            for &r in &rows {
                let mut tag = format!("        <node id=\"{}\"", escape(&xml_id[r]));
                if let Some(t) = text_column(col("title"), r) {
                    tag.push_str(&format!(" title=\"{}\"", escape(&t)));
                }
                let values: Vec<String> = declared
                    .iter()
                    .filter_map(|(name, column)| {
                        column.get(r).to_text().map(|v| {
                            format!(
                                "            <property id=\"{}\" value=\"{}\"/>\n",
                                escape(name),
                                escape(&v)
                            )
                        })
                    })
                    .collect();
                if values.is_empty() {
                    out.push_str(&tag);
                    out.push_str("/>\n");
                } else {
                    out.push_str(&tag);
                    out.push_str(">\n          <properties>\n");
                    for v in values {
                        out.push_str(&v);
                    }
                    out.push_str("          </properties>\n        </node>\n");
                }
            }
            out.push_str("      </nodeclass>\n");
        }
        out.push_str("    </nodes>\n    <networks>\n");
        for w in networks.iter().filter(|w| w.period == p) {
            let (src, tgt) = (&classes[w.source], &classes[w.target]);
            let self_loops = w.rows.iter().any(|&r| from[r] == to[r]);
            out.push_str(&format!(
                "      <network sourceType=\"{}\" source=\"{}\" targetType=\"{}\" target=\"{}\" \
                 id=\"{}\" isBinary=\"{}\" isDirected=\"{}\" allowSelfLoops=\"{}\">\n",
                escape(&src.2),
                escape(&src.1),
                escape(&tgt.2),
                escape(&tgt.1),
                escape(&w.id),
                weight_col.is_none(),
                w.directed,
                self_loops
            ));
            // An undirected network stored in a directed table holds each
            // link twice (see the module notes); write each pair once.
            let mut pending: HashMap<(usize, usize), usize> = HashMap::new();
            let dedupe = net.directed && !w.directed;
            for &r in &w.rows {
                let (a, b) = (index[&from[r]], index[&to[r]]);
                if dedupe && a != b {
                    if let Some(k) = pending.get_mut(&(b, a)).filter(|k| **k > 0) {
                        *k -= 1;
                        continue;
                    }
                    *pending.entry((a, b)).or_insert(0) += 1;
                }
                // Each link runs from the source class to the target class.
                let (a, b) = if node_class[a] != w.source {
                    (b, a)
                } else {
                    (a, b)
                };
                let mut tag = format!(
                    "        <link source=\"{}\" target=\"{}\"",
                    escape(&xml_id[a]),
                    escape(&xml_id[b])
                );
                match weight_col.map(|c| c.get(r)) {
                    Some(Value::Float(f)) => {
                        tag.push_str(&format!(" type=\"double\" value=\"{}\"", format_float(f)))
                    }
                    Some(Value::Int(i)) => tag.push_str(&format!(
                        " type=\"double\" value=\"{}\"",
                        format_float(i as f64)
                    )),
                    _ => {}
                }
                let values: Vec<String> = edge_props
                    .iter()
                    .filter_map(|(name, column)| {
                        column.get(r).to_text().map(|v| {
                            format!(
                                "            <property id=\"{}\" type=\"{}\" value=\"{}\"/>\n",
                                escape(name),
                                type_name(column.kind()),
                                escape(&v)
                            )
                        })
                    })
                    .collect();
                if values.is_empty() {
                    out.push_str(&tag);
                    out.push_str("/>\n");
                } else {
                    out.push_str(&tag);
                    out.push_str(">\n          <properties>\n");
                    for v in values {
                        out.push_str(&v);
                    }
                    out.push_str("          </properties>\n        </link>\n");
                }
            }
            out.push_str("      </network>\n");
        }
        out.push_str("    </networks>\n  </MetaNetwork>\n");
    }
    out.push_str("</DynamicNetwork>\n");
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Modeled on the meta-networks ORA saves: a `DynamicNetwork` root
    /// with one `MetaNetwork`, typed node classes with property identities,
    /// and one- and two-mode networks with `link` values.
    const ORA: &str = r#"<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<DynamicNetwork xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <MetaNetwork id="embassy">
    <propertyIdentities>
      <propertyIdentity id="source" type="text" singleValued="true"/>
    </propertyIdentities>
    <properties>
      <property id="source" value="field notes"/>
    </properties>
    <nodes>
      <nodeclass type="Agent" id="Agent">
        <propertyIdentities>
          <propertyIdentity id="age" type="number" singleValued="true"/>
          <propertyIdentity id="gender" type="categoryText" singleValued="true"/>
        </propertyIdentities>
        <node id="A01" title="Ann">
          <properties>
            <property id="age" value="34"/>
            <property id="gender" value="female"/>
          </properties>
        </node>
        <node id="A02" title="Bo">
          <properties>
            <property id="age" value=""/>
          </properties>
        </node>
        <node id="Kabul"/>
      </nodeclass>
      <nodeclass type="Location" id="Location">
        <node id="Kabul"/>
        <node id="Herat"/>
      </nodeclass>
      <nodeclass type="Task" id="Task">
        <node id="T1" title="Logistics"/>
      </nodeclass>
    </nodes>
    <networks>
      <network sourceType="Agent" source="Agent" targetType="Agent" target="Agent" id="agent x agent" isBinary="false" isDirected="true" allowSelfLoops="false">
        <link source="A01" target="A02" type="double" value="2.5"/>
        <link source="A02" target="A01" value="1"/>
        <link source="A02" target="Kabul"/>
      </network>
      <network sourceType="Agent" source="Agent" targetType="Location" target="Location" id="agent x location" isDirected="false">
        <link source="A01" target="Kabul" value="1"/>
      </network>
      <network sourceType="Agent" source="Agent" targetType="Task" target="Task" id="agent x task" isDirected="true">
        <link source="A02" target="T1" value="0"/>
      </network>
    </networks>
  </MetaNetwork>
</DynamicNetwork>
"#;

    fn strs(net: &Network, table: &str, col: &str) -> Vec<Option<String>> {
        let t = if table == "nodes" {
            &net.nodes
        } else {
            &net.edges
        };
        let c = t.column(col).unwrap_or_else(|| panic!("no column {col}"));
        (0..c.len()).map(|i| c.get(i).to_text()).collect()
    }

    fn some(v: &[&str]) -> Vec<Option<String>> {
        v.iter().map(|s| Some(s.to_string())).collect()
    }

    #[test]
    fn reads_ora_meta_network() {
        let n = read(ORA).unwrap();
        // "Kabul" is both an agent and a location, so identifiers are
        // qualified by class and the original id is kept.
        assert_eq!(
            n.nodes.strings("id").unwrap(),
            vec![
                "Agent::A01",
                "Agent::A02",
                "Agent::Kabul",
                "Location::Kabul",
                "Location::Herat",
                "Task::T1"
            ]
        );
        assert_eq!(
            strs(&n, "nodes", "node_id"),
            some(&["A01", "A02", "Kabul", "Kabul", "Herat", "T1"])
        );
        assert_eq!(
            strs(&n, "nodes", "type"),
            some(&["Agent", "Agent", "Agent", "Location", "Location", "Task"])
        );
        assert_eq!(n.nodes.column("class"), n.nodes.column("type"));
        assert_eq!(
            n.nodes.column("age").unwrap(),
            &Column::Float(vec![Some(34.0), None, None, None, None, None])
        );
        assert_eq!(
            n.nodes.column("title").unwrap(),
            &Column::Str(vec![
                Some("Ann".into()),
                Some("Bo".into()),
                None,
                None,
                None,
                Some("Logistics".into())
            ])
        );
        assert_eq!(
            n.graph,
            vec![
                ("name".into(), Value::Str("embassy".into())),
                ("source".into(), Value::Str("field notes".into()))
            ]
        );

        // Mixed direction: read as directed, the undirected link twice.
        assert!(n.directed);
        assert_eq!(
            strs(&n, "edges", "from"),
            some(&[
                "Agent::A01",
                "Agent::A02",
                "Agent::A02",
                "Agent::A01",
                "Location::Kabul",
                "Agent::A02"
            ])
        );
        assert_eq!(
            strs(&n, "edges", "to"),
            some(&[
                "Agent::A02",
                "Agent::A01",
                "Agent::Kabul",
                "Location::Kabul",
                "Agent::A01",
                "Task::T1"
            ])
        );
        assert_eq!(
            n.edges.column("weight").unwrap(),
            &Column::Float(vec![
                Some(2.5),
                Some(1.0),
                Some(1.0),
                Some(1.0),
                Some(1.0),
                Some(0.0)
            ])
        );
        assert_eq!(
            n.edges.column("directed").unwrap(),
            &Column::Bool(vec![
                Some(true),
                Some(true),
                Some(true),
                Some(false),
                Some(false),
                Some(true)
            ])
        );
        assert_eq!(
            strs(&n, "edges", "network"),
            some(&[
                "agent x agent",
                "agent x agent",
                "agent x agent",
                "agent x location",
                "agent x location",
                "agent x task"
            ])
        );
        n.validate().unwrap();
    }

    /// The structure of the 2004 report: `MetaMatrix` time periods,
    /// `nodeset`, `graph`, and `edge` with a required `type`, and
    /// properties and measures carrying `name`, `type`, and `value`.
    const REPORT_2004: &str = r#"<?xml version="1.0" encoding="UTF-8"?>
<DynamicNetwork>
  <MetaMatrix timePeriod="1">
    <measures>
      <measure name="ora.density" type="double" value="0.5"/>
    </measures>
    <nodes>
      <nodeset id="people" type="agent">
        <node id="actor1" title="First">
          <properties>
            <property name="rank" type="string" value="major"/>
            <property name="cleared" type="binary" value="true"/>
          </properties>
          <measures>
            <measure name="ora.degree" type="double" value="1"/>
          </measures>
        </node>
        <node id="actor2"/>
      </nodeset>
      <nodeset id="skills" type="knowledge">
        <node id="k1"/>
      </nodeset>
    </nodes>
    <networks>
      <graph id="comm" sourceType="agent" source="people" targetType="agent" target="people" isDirected="false">
        <edge source="actor1" target="actor2" type="binary" value="true"/>
      </graph>
      <graph id="knows" sourceType="agent" source="people" targetType="knowledge" target="skills" isDirected="true">
        <edge source="actor2" target="k1" type="double" value="0.25"/>
      </graph>
    </networks>
  </MetaMatrix>
  <MetaMatrix timePeriod="2">
    <nodes>
      <nodeset id="people" type="agent">
        <node id="actor1"/>
        <node id="actor3"/>
      </nodeset>
    </nodes>
    <networks>
      <graph id="comm" source="people" target="people" isDirected="false">
        <edge source="actor3" target="actor1" type="double" value="3"/>
      </graph>
    </networks>
  </MetaMatrix>
</DynamicNetwork>
"#;

    #[test]
    fn reads_2004_time_periods() {
        let n = read(REPORT_2004).unwrap();
        // Several periods: identifiers carry the period, and every period's
        // nodes are separate rows.
        assert_eq!(
            n.nodes.strings("id").unwrap(),
            vec!["1::actor1", "1::actor2", "1::k1", "2::actor1", "2::actor3"]
        );
        assert_eq!(
            strs(&n, "nodes", "node_id"),
            some(&["actor1", "actor2", "k1", "actor1", "actor3"])
        );
        assert_eq!(
            strs(&n, "nodes", "period"),
            some(&["1", "1", "1", "2", "2"])
        );
        assert_eq!(
            strs(&n, "nodes", "class"),
            some(&["people", "people", "skills", "people", "people"])
        );
        assert_eq!(
            strs(&n, "nodes", "type"),
            some(&["agent", "agent", "knowledge", "agent", "agent"])
        );
        assert_eq!(
            n.nodes.column("cleared").unwrap(),
            &Column::Bool(vec![Some(true), None, None, None, None])
        );
        assert_eq!(
            n.nodes.column("ora.degree").unwrap(),
            &Column::Float(vec![Some(1.0), None, None, None, None])
        );
        assert_eq!(n.graph, vec![("1::ora.density".into(), Value::Float(0.5))]);
        // "comm" is undirected and "knows" directed: mixed.
        assert!(n.directed);
        assert_eq!(
            strs(&n, "edges", "from"),
            some(&[
                "1::actor1",
                "1::actor2",
                "1::actor2",
                "2::actor3",
                "2::actor1"
            ])
        );
        assert_eq!(
            strs(&n, "edges", "period"),
            some(&["1", "1", "1", "2", "2"])
        );
        assert_eq!(
            n.edges.column("weight").unwrap(),
            &Column::Float(vec![Some(1.0), Some(1.0), Some(0.25), Some(3.0), Some(3.0)])
        );
    }

    #[test]
    fn reads_minimal_files() {
        // As in the CASOS "company" sample: lower-case types, no
        // isDirected, networks naming only node class types.
        let doc = r#"<DynamicNetwork><MetaNetwork>
<nodes>
<nodeclass type="agent" id="agent"><node id="LDR" title="Project Manager"/><node id="Mgr1"/></nodeclass>
<nodeclass type="task" id="task"><node id="T1"/></nodeclass>
</nodes>
<networks>
<network sourceType="agent" targetType="agent" id="agent x agent">
<link source="LDR" target="Mgr1" type="double" value="1.0000"/>
</network>
<network sourceType="agent" targetType="task">
<link source="Mgr1" target="T1"/>
</network>
</networks>
</MetaNetwork></DynamicNetwork>"#;
        let n = read(doc).unwrap();
        assert!(n.directed);
        assert_eq!(n.nodes.strings("id").unwrap(), vec!["LDR", "Mgr1", "T1"]);
        assert!(n.nodes.column("node_id").is_none());
        assert_eq!(
            strs(&n, "edges", "network"),
            some(&["agent x agent", "agent x task"])
        );
        assert_eq!(
            n.edges.column("weight").unwrap(),
            &Column::Float(vec![Some(1.0), Some(1.0)])
        );

        // A bare MetaNetwork root, undirected, unweighted, with a node only
        // mentioned in a link.
        let doc = r#"<MetaNetwork id="m"><nodes><nodeclass type="Agent" id="Agent"><node id="a"/></nodeclass></nodes>
<networks><network source="Agent" target="Agent" id="g" isDirected="false"><link source="a" target="b"/></network></networks></MetaNetwork>"#;
        let n = read(doc).unwrap();
        assert!(!n.directed);
        assert_eq!(n.nodes.strings("id").unwrap(), vec!["a", "b"]);
        assert!(n.edges.column("weight").is_none());
        assert!(n.edges.column("directed").is_none());
        assert_eq!(n.graph, vec![("name".into(), Value::Str("m".into()))]);
    }

    #[test]
    fn rejects_unsupported_content_with_lines() {
        let cases = [
            (
                ORA.replace(
                    "<node id=\"Herat\"/>",
                    "<node id=\"Herat\"><port name=\"p\"/></node>",
                ),
                "ports",
                31,
            ),
            (
                ORA.replace("value=\"2.5\"", "value=\"high\""),
                "not a number",
                39,
            ),
            (
                ORA.replace(
                    "<node id=\"Herat\"/>",
                    "<node id=\"Herat\"/><node id=\"Herat\"/>",
                ),
                "twice",
                31,
            ),
            (ORA.replace("<nodes>", "<nodes><bogus/>"), "<bogus>", 10),
            (
                ORA.replace("value=\"34\"", "value=\"old\""),
                "not a valid Float",
                18,
            ),
            (
                ORA.replace(
                    "<property id=\"gender\" value=\"female\"/>",
                    "<property id=\"age\" value=\"35\"/>",
                ),
                "more than once",
                19,
            ),
            (
                ORA.replace("isDirected=\"false\"", "isDirected=\"maybe\""),
                "isDirected",
                43,
            ),
            (
                REPORT_2004.replace("<measures>", "<anthropac/><measures>"),
                "AnthroPac",
                4,
            ),
        ];
        for (doc, needle, line) in cases {
            let e = read(&doc).unwrap_err().to_string();
            assert!(e.contains(needle), "{e}");
            assert!(e.contains(&format!("line {line}:")), "{e}");
        }
        assert!(read("<graphml/>")
            .unwrap_err()
            .to_string()
            .contains("found <graphml>"));
        assert!(read("<DynamicNetwork/>").is_err());
        let twice = ORA.replace("id=\"agent x task\"", "id=\"agent x agent\"");
        assert!(read(&twice)
            .unwrap_err()
            .to_string()
            .contains("network \"agent x agent\" appears more than once"));
        // A node id in two classes of one type cannot be resolved when the
        // network names only the type.
        let doc = r#"<MetaNetwork><nodes><nodeclass id="x" type="Agent"><node id="a"/></nodeclass>
<nodeclass id="y" type="Agent"><node id="a"/></nodeclass></nodes>
<networks><network sourceType="Agent" targetType="Agent"><link source="a" target="a"/></network></networks></MetaNetwork>"#;
        assert!(read(doc)
            .unwrap_err()
            .to_string()
            .contains("several node classes"));
    }

    fn assert_same(a: &Network, b: &Network) {
        assert_eq!(a.directed, b.directed);
        assert_eq!(a.graph, b.graph);
        assert_eq!(a.nodes, b.nodes);
        assert_eq!(a.edges, b.edges);
    }

    #[test]
    fn round_trips_files() {
        let n = read(ORA).unwrap();
        assert_same(&read(&write(&n, None).unwrap()).unwrap(), &n);

        // Logical properties are written as ORA's categoryText and come
        // back as text; everything else is unchanged.
        let mut n = read(REPORT_2004).unwrap();
        let back = read(&write(&n, None).unwrap()).unwrap();
        let pos = n.nodes.names.iter().position(|c| c == "cleared").unwrap();
        n.nodes.columns[pos] = Column::Str(vec![Some("true".into()), None, None, None, None]);
        assert_same(&back, &n);
    }

    #[test]
    fn round_trips_tables() {
        // A network built elsewhere, with classes, properties of every
        // kind, two networks, and weights.
        let mut b = NetworkBuilder::new(true);
        b.graph
            .push(("name".into(), Value::Str("demo & co".into())));
        b.graph.push(("year".into(), Value::Float(2004.0)));
        let node = |class: &str, title: Option<&str>, score: Option<f64>, team: Option<&str>| {
            vec![
                ("class".to_string(), Value::Str(class.into())),
                (
                    "title".to_string(),
                    title.map_or(Value::Null, |t| Value::Str(t.into())),
                ),
                ("score".to_string(), score.map_or(Value::Null, Value::Float)),
                (
                    "team".to_string(),
                    team.map_or(Value::Null, |t| Value::Str(t.into())),
                ),
            ]
        };
        b.add_node(
            "ana",
            node("Agent", Some("Ana \"A\""), Some(0.5), Some("red <1>")),
        );
        b.add_node("bo", node("Agent", None, None, Some("")));
        b.add_node(
            "acme",
            node("Organization", Some("ACME"), Some(-1e-7), None),
        );
        b.add_node("q", node("Survey", None, None, None));
        let edge = |net: &str, w: f64| {
            vec![
                ("network".to_string(), Value::Str(net.into())),
                ("weight".to_string(), Value::Float(w)),
            ]
        };
        b.add_edge("ana", "bo", edge("talks", 2.0));
        b.add_edge("bo", "bo", edge("talks", 0.1 + 0.2));
        b.add_edge("ana", "acme", edge("member", 1.0));
        b.add_edge("bo", "acme", edge("member", 1.0));
        let n = b.finish().unwrap();

        let text = write(&n, None).unwrap();
        assert!(text.contains("<nodeclass type=\"Organization\" id=\"Organization\">"));
        assert!(text.contains("<nodeclass type=\"Unknown\" id=\"Survey\">"));
        assert!(text.contains(
            "sourceType=\"Agent\" source=\"Agent\" targetType=\"Organization\" \
             target=\"Organization\" id=\"member\" isBinary=\"false\" isDirected=\"true\" \
             allowSelfLoops=\"false\""
        ));
        assert!(text.contains("allowSelfLoops=\"true\""));
        let back = read(&text).unwrap();
        assert_eq!(back.directed, n.directed);
        assert_eq!(back.graph, n.graph);
        assert_eq!(back.nodes.column("id"), n.nodes.column("id"));
        assert_eq!(back.nodes.column("class"), n.nodes.column("class"));
        assert_eq!(
            strs(&back, "nodes", "type"),
            some(&["Agent", "Agent", "Organization", "Unknown"])
        );
        for col in ["title", "score", "team"] {
            assert_eq!(back.nodes.column(col), n.nodes.column(col), "{col}");
        }
        assert_eq!(back.edges, n.edges);

        // Without class, network, or weight columns: one Agent class, one
        // network, binary links; undirected networks stay undirected.
        let mut b = NetworkBuilder::new(false);
        b.add_edge("1", "2", vec![]);
        b.add_edge("2", "3", vec![("kind".into(), Value::Int(4))]);
        let n = b.finish().unwrap();
        let text = write(&n, None).unwrap();
        assert!(text.contains("id=\"Agent x Agent\" isBinary=\"true\" isDirected=\"false\""));
        let back = read(&text).unwrap();
        assert!(!back.directed);
        assert_eq!(back.edges.strings("to").unwrap(), vec!["2", "3"]);
        assert_eq!(
            back.edges.column("kind").unwrap(),
            &Column::Float(vec![None, Some(4.0)])
        );
    }

    #[test]
    fn writes_the_weight_column_asked_for() {
        let mut b = NetworkBuilder::new(true);
        b.add_edge("a", "b", vec![("w".into(), Value::Int(3))]);
        let n = b.finish().unwrap();
        let text = write(&n, Some("w")).unwrap();
        assert!(text.contains("type=\"double\" value=\"3.0\""), "{text}");
        let back = read(&text).unwrap();
        assert_eq!(
            back.edges.column("weight").unwrap(),
            &Column::Float(vec![Some(3.0)])
        );
        assert!(back.edges.column("w").is_none());
    }

    #[test]
    fn writer_rejects_what_dynetml_cannot_hold() {
        let mut b = NetworkBuilder::new(true);
        b.add_node("a", vec![("class".into(), Value::Str("Agent".into()))]);
        b.add_node("b", vec![("class".into(), Value::Str("Task".into()))]);
        b.add_edge("a", "b", vec![("network".into(), Value::Str("n".into()))]);
        b.add_edge("a", "a", vec![("network".into(), Value::Str("n".into()))]);
        let e = write(&b.finish().unwrap(), None).unwrap_err().to_string();
        assert!(e.contains("more than one pair of node classes"), "{e}");

        let mut b = NetworkBuilder::new(true);
        b.add_node("x::1", vec![("node_id".into(), Value::Str("1".into()))]);
        b.add_node("y::1", vec![("node_id".into(), Value::Str("1".into()))]);
        let e = write(&b.finish().unwrap(), None).unwrap_err().to_string();
        assert!(e.contains("appears twice"), "{e}");
    }
}
