//! Reading and writing network files.
//!
//! Every format reads into, and writes from, a [`Network`]: a node table
//! whose first column is `id`, an edge table whose first columns are `from`
//! and `to`, typed attribute columns, and graph-level attributes.
//!
//! | Format | Read | Write | Typical source |
//! |---|---|---|---|
//! | GraphML | yes | yes | igraph, NetworkX, Gephi, yEd, Cytoscape |
//! | GEXF | yes | yes | Gephi |
//! | GML | yes | yes | igraph, NetworkX, Cytoscape |
//! | Pajek `.net` | yes | yes | Pajek, igraph, NetworkX |
//! | UCINET DL | yes | yes | UCINET, igraph |
//! | JSON node-link | yes | yes | NetworkX, D3 |
//! | DOT | no | yes | Graphviz |
//!
//! Each format module documents what it supports and how it maps the
//! format onto tables.

pub mod dl;
pub mod dot;
pub mod gexf;
pub mod gml;
pub mod graphml;
pub mod json;
pub mod nodelink;
pub mod pajek;
pub mod table;
pub mod xml;

use std::path::Path;

pub use table::{Column, Kind, Network, NetworkBuilder, Table, Value};

use crate::error::{Error, Result};

/// A supported file format.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    /// GraphML.
    GraphMl,
    /// GEXF (Gephi).
    Gexf,
    /// Graph Modelling Language.
    Gml,
    /// Pajek `.net`.
    Pajek,
    /// UCINET DL.
    Dl,
    /// JSON node-link data.
    NodeLink,
    /// Graphviz DOT (write only).
    Dot,
}

impl Format {
    /// All formats.
    pub const ALL: [Format; 7] = [
        Format::GraphMl,
        Format::Gexf,
        Format::Gml,
        Format::Pajek,
        Format::Dl,
        Format::NodeLink,
        Format::Dot,
    ];

    /// The canonical short name, as accepted by [`Format::from_name`].
    pub fn name(self) -> &'static str {
        match self {
            Format::GraphMl => "graphml",
            Format::Gexf => "gexf",
            Format::Gml => "gml",
            Format::Pajek => "pajek",
            Format::Dl => "dl",
            Format::NodeLink => "json",
            Format::Dot => "dot",
        }
    }

    /// Looks up a format by name (case-insensitive; common aliases accepted).
    pub fn from_name(name: &str) -> Option<Format> {
        Some(match name.to_ascii_lowercase().as_str() {
            "graphml" | "xml" => Format::GraphMl,
            "gexf" => Format::Gexf,
            "gml" => Format::Gml,
            "pajek" | "net" | "paj" => Format::Pajek,
            "dl" | "ucinet" => Format::Dl,
            "json" | "nodelink" | "node-link" | "node_link" => Format::NodeLink,
            "dot" | "gv" | "graphviz" => Format::Dot,
            _ => return None,
        })
    }

    /// Guesses the format from a file extension.
    pub fn from_path(path: &Path) -> Option<Format> {
        let ext = path.extension()?.to_str()?;
        Format::from_name(ext)
    }

    /// Whether the format can be read.
    pub fn can_read(self) -> bool {
        self != Format::Dot
    }
}

/// Options for reading.
#[derive(Debug, Clone, Copy, Default)]
pub struct ReadOptions {
    /// Direction for formats that do not record it (UCINET DL). Ignored by
    /// formats that do.
    pub directed: Option<bool>,
}

/// Options for writing.
#[derive(Debug, Clone, Default)]
pub struct WriteOptions {
    /// Edge column written as the link value by formats with a single
    /// positional value (Pajek, UCINET DL). Defaults to `weight`.
    pub weight: Option<String>,
}

/// Parses network text in the given format.
pub fn read_str(format: Format, text: &str, opts: &ReadOptions) -> Result<Network> {
    let text = text.strip_prefix('\u{feff}').unwrap_or(text);
    let net = match format {
        Format::GraphMl => graphml::read(text)?,
        Format::Gexf => gexf::read(text)?,
        Format::Gml => gml::read(text)?,
        Format::Pajek => pajek::read(text)?,
        Format::Dl => dl::read(text, opts.directed)?,
        Format::NodeLink => nodelink::read(text)?,
        Format::Dot => {
            return Err(Error::InvalidArgument(
                "reading DOT files is not supported; DOT is write-only".into(),
            ))
        }
    };
    net.validate()?;
    Ok(net)
}

/// Serializes a network in the given format.
pub fn write_string(format: Format, net: &Network, opts: &WriteOptions) -> Result<String> {
    let weight = opts.weight.as_deref();
    match format {
        Format::GraphMl => graphml::write(net),
        Format::Gexf => gexf::write(net),
        Format::Gml => gml::write(net),
        Format::Pajek => pajek::write(net, weight),
        Format::Dl => dl::write(net, weight),
        Format::NodeLink => nodelink::write(net),
        Format::Dot => dot::write(net),
    }
}

fn resolve(path: &Path, format: Option<Format>) -> Result<Format> {
    format.or_else(|| Format::from_path(path)).ok_or_else(|| {
        Error::InvalidArgument(format!(
            "cannot tell the format of {:?} from its extension; supported formats are {}",
            path,
            Format::ALL.map(Format::name).join(", ")
        ))
    })
}

/// Reads a network file. The format is inferred from the extension when
/// `format` is `None`.
pub fn read_file(path: &Path, format: Option<Format>, opts: &ReadOptions) -> Result<Network> {
    let format = resolve(path, format)?;
    let bytes = std::fs::read(path).map_err(|e| Error::Io(format!("cannot read {path:?}: {e}")))?;
    let text = match String::from_utf8(bytes) {
        Ok(t) => t,
        // Older Pajek and UCINET files are often Latin-1.
        Err(e) => e.into_bytes().iter().map(|&b| b as char).collect(),
    };
    read_str(format, &text, opts)
}

/// Writes a network file. The format is inferred from the extension when
/// `format` is `None`.
pub fn write_file(
    path: &Path,
    net: &Network,
    format: Option<Format>,
    opts: &WriteOptions,
) -> Result<()> {
    let format = resolve(path, format)?;
    let text = write_string(format, net, opts)?;
    std::fs::write(path, text).map_err(|e| Error::Io(format!("cannot write {path:?}: {e}")))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> Network {
        let mut b = NetworkBuilder::new(false);
        b.graph.push(("name".into(), Value::Str("sample".into())));
        b.add_node(
            "ann",
            vec![
                ("age".into(), Value::Int(31)),
                ("score".into(), Value::Float(0.25)),
            ],
        );
        b.add_node("bob", vec![("age".into(), Value::Int(45))]);
        b.add_node("cy", vec![]);
        b.add_node("isolate", vec![]);
        b.add_edge("ann", "bob", vec![("weight".into(), Value::Float(2.0))]);
        b.add_edge("bob", "cy", vec![("weight".into(), Value::Float(0.5))]);
        b.add_edge("cy", "cy", vec![("weight".into(), Value::Float(1.0))]);
        b.finish().unwrap()
    }

    #[test]
    fn every_readable_format_preserves_structure() {
        let net = sample();
        for format in Format::ALL.into_iter().filter(|f| f.can_read()) {
            let text = write_string(format, &net, &WriteOptions::default()).unwrap();
            let back = read_str(
                format,
                &text,
                &ReadOptions {
                    directed: Some(false),
                },
            )
            .unwrap();
            assert_eq!(
                back.nodes.strings("id").unwrap(),
                net.nodes.strings("id").unwrap(),
                "{format:?}"
            );
            assert_eq!(
                back.edges.strings("from").unwrap(),
                net.edges.strings("from").unwrap(),
                "{format:?}"
            );
            assert_eq!(
                back.edges.strings("to").unwrap(),
                net.edges.strings("to").unwrap(),
                "{format:?}"
            );
            assert_eq!(back.directed, net.directed, "{format:?}");
            assert_eq!(
                back.edges.column("weight"),
                net.edges.column("weight"),
                "{format:?}"
            );
        }
    }

    #[test]
    fn attribute_rich_formats_preserve_attributes() {
        let net = sample();
        for format in [Format::GraphMl, Format::Gexf, Format::NodeLink] {
            let text = write_string(format, &net, &WriteOptions::default()).unwrap();
            let back = read_str(format, &text, &ReadOptions::default()).unwrap();
            for name in ["age", "score"] {
                assert_eq!(
                    back.nodes.column(name),
                    net.nodes.column(name),
                    "{format:?} {name}"
                );
            }
        }
    }

    #[test]
    fn formats_from_names_and_paths() {
        assert_eq!(
            Format::from_path(Path::new("x/net.GraphML")),
            Some(Format::GraphMl)
        );
        assert_eq!(Format::from_path(Path::new("a.net")), Some(Format::Pajek));
        assert_eq!(Format::from_name("UCINET"), Some(Format::Dl));
        assert!(read_str(Format::Dot, "", &ReadOptions::default()).is_err());
    }
}
