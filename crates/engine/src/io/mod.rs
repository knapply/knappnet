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
//! | DyNetML | yes | yes | ORA (CASOS) |
//! | DOT | no | yes | Graphviz |
//!
//! Each format module documents what it supports and how it maps the
//! format onto tables.
//!
//! Readers take time and memory proportional to the size of what they
//! build: the input, plus a node count declared without listing the nodes
//! (Pajek, UCINET DL; limited to [`MAX_DECLARED_NODES`]). Edge tables hold
//! each endpoint's identifier as text, so memory also grows with the number
//! of edges times the length of their identifiers. Every reader is fuzzed
//! (`tests/io_fuzz.rs`): malformed input gives an error, never a panic.
//!
//! GraphML and DyNetML files both commonly end in `.xml`. When a `.xml`
//! file is read without a format, its root element decides: `<graphml>` is
//! GraphML, and `<DynamicMetaNetwork>`, `<DynamicNetwork>`, `<MetaNetwork>`,
//! or `<MetaMatrix>` is DyNetML. Writing a `.xml` file without a format
//! produces GraphML, as before, since GraphML is the more widely read of
//! the two; pass the `dynetml` format to write DyNetML.

pub mod dl;
pub mod dot;
pub mod dynetml;
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

/// The largest node count a file may declare without listing the nodes
/// (Pajek `*Vertices n`, UCINET DL `n=`). Both formats let a count alone
/// create isolates, so a short file can ask for any number of nodes; the
/// limit keeps a corrupt or hostile count from exhausting memory. A
/// network this large would need tens of gigabytes as tables.
pub const MAX_DECLARED_NODES: usize = 100_000_000;

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
    /// DyNetML (ORA).
    DynetMl,
    /// Graphviz DOT (write only).
    Dot,
}

impl Format {
    /// All formats.
    pub const ALL: [Format; 8] = [
        Format::GraphMl,
        Format::Gexf,
        Format::Gml,
        Format::Pajek,
        Format::Dl,
        Format::NodeLink,
        Format::DynetMl,
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
            Format::DynetMl => "dynetml",
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
            "dynetml" | "ora" => Format::DynetMl,
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
        Format::DynetMl => dynetml::read(text)?,
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
        Format::DynetMl => dynetml::write(net, weight),
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

/// The name of the root element of an XML document, if one can be found:
/// skips the XML declaration, processing instructions, comments, and a
/// DOCTYPE, and drops any namespace prefix.
pub fn xml_root_name(text: &str) -> Option<&str> {
    let mut rest = text.strip_prefix('\u{feff}').unwrap_or(text);
    loop {
        rest = rest.trim_start();
        if rest.starts_with("<?") {
            rest = &rest[rest.find("?>")? + 2..];
        } else if rest.starts_with("<!--") {
            rest = &rest[rest.find("-->")? + 3..];
        } else if rest.starts_with("<!") {
            // DOCTYPE, possibly with an internal subset.
            let close = rest.find('>')?;
            let end = match rest[..close].find('[') {
                Some(_) => rest.find("]>").or_else(|| rest.find(']'))? + 1,
                None => close,
            };
            rest = &rest[end..];
            rest = &rest[rest.find('>')? + 1..];
        } else {
            let tag = rest.strip_prefix('<')?;
            let end = tag
                .find(|c: char| c.is_whitespace() || c == '>' || c == '/')
                .unwrap_or(tag.len());
            let name = &tag[..end];
            return Some(name.rsplit(':').next().unwrap_or(name));
        }
    }
}

/// Chooses between GraphML and DyNetML for an XML document from its root
/// element (GraphML when the root is not a DyNetML root).
pub fn sniff_xml(text: &str) -> Format {
    match xml_root_name(text) {
        Some(root) if dynetml::ROOTS.contains(&root) => Format::DynetMl,
        _ => Format::GraphMl,
    }
}

/// Reads a network file. The format is inferred from the extension when
/// `format` is `None`; a `.xml` file is GraphML or DyNetML according to its
/// root element.
pub fn read_file(path: &Path, format: Option<Format>, opts: &ReadOptions) -> Result<Network> {
    let sniff = format.is_none()
        && path
            .extension()
            .and_then(|e| e.to_str())
            .is_some_and(|e| e.eq_ignore_ascii_case("xml"));
    let format = resolve(path, format)?;
    let bytes = std::fs::read(path).map_err(|e| Error::Io(format!("cannot read {path:?}: {e}")))?;
    let text = decode(bytes);
    let format = if sniff { sniff_xml(&text) } else { format };
    read_str(format, &text, opts)
}

/// Decodes file contents as UTF-8, or as Latin-1 when they are not valid
/// UTF-8 (older Pajek and UCINET files are often Latin-1). Every byte
/// sequence decodes, so readers only ever see valid text.
pub fn decode(bytes: Vec<u8>) -> String {
    match String::from_utf8(bytes) {
        Ok(t) => t,
        Err(e) => e.into_bytes().iter().map(|&b| b as char).collect(),
    }
}

/// Parses file contents in the given format, decoding them as
/// [`read_file`] does.
pub fn read_bytes(format: Format, bytes: &[u8], opts: &ReadOptions) -> Result<Network> {
    read_str(format, &decode(bytes.to_vec()), opts)
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
        assert_eq!(Format::from_name("ORA"), Some(Format::DynetMl));
        assert_eq!(Format::from_name("dynetml"), Some(Format::DynetMl));
        // Writing .xml stays GraphML; reading sniffs (see below).
        assert_eq!(Format::from_path(Path::new("a.xml")), Some(Format::GraphMl));
        assert!(read_str(Format::Dot, "", &ReadOptions::default()).is_err());
    }

    #[test]
    fn xml_files_are_sniffed_by_root_element() {
        let roots = [
            (
                "<?xml version=\"1.0\"?>\n<!-- c -->\n<graphml xmlns=\"x\"/>",
                Format::GraphMl,
            ),
            ("\u{feff}<DynamicNetwork>", Format::DynetMl),
            (
                "<!DOCTYPE d [ <!ENTITY e \"x\"> ]><DynamicMetaNetwork id=\"a\">",
                Format::DynetMl,
            ),
            ("<?xml version=\"1.0\"?><MetaNetwork>", Format::DynetMl),
            ("<ora:MetaMatrix xmlns:ora=\"u\">", Format::DynetMl),
            ("<gexf>", Format::GraphMl),
            ("not xml", Format::GraphMl),
        ];
        for (text, expected) in roots {
            assert_eq!(sniff_xml(text), expected, "{text}");
        }

        let dir = std::env::temp_dir();
        let net = sample();
        for format in [Format::GraphMl, Format::DynetMl] {
            let path = dir.join(format!(
                "knappnet-sniff-{}-{}.xml",
                format.name(),
                std::process::id()
            ));
            write_file(&path, &net, Some(format), &WriteOptions::default()).unwrap();
            let back = read_file(&path, None, &ReadOptions::default()).unwrap();
            std::fs::remove_file(&path).unwrap();
            assert_eq!(
                back.nodes.strings("id").unwrap(),
                net.nodes.strings("id").unwrap()
            );
            // Only DyNetML adds node classes.
            assert_eq!(
                back.nodes.column("class").is_some(),
                format == Format::DynetMl
            );
        }
    }
}
