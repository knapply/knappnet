//! Mutation-based fuzz tests for every file reader.
//!
//! Each case starts from a valid file (written by this crate's writers, or
//! written by hand to cover features the writers never produce), applies
//! one to four random mutations drawn from [`engine::rng::Stream`], and
//! reads the result as raw bytes, as [`engine::io::read_file`] would. A
//! reader must:
//!
//! - never panic;
//! - never hang: every read runs on a worker thread and must finish within
//!   a time limit (inputs are capped at a few hundred kilobytes, so this
//!   bounds the work per byte);
//! - either return a network that passes [`Network::validate`], whose
//!   tables are rectangular, and that every writer can serialize and the
//!   GraphML reader can read back unchanged; or return an error with a
//!   message, whose line number (when it gives one) lies inside the input
//!   and whose format name is the reader's.
//!
//! Mutations: bit flips, interesting bytes, deletions, duplications,
//! random insertions, truncation, splices of two valid files (of any
//! format), out-of-range and malformed numbers, deeply nested markup,
//! invalid UTF-8, repeated lines, and format keywords.
//!
//! `cargo test` runs a quick pass. The long pass is ignored by default:
//!
//! ```text
//! FUZZ_ITERATIONS=20000 cargo test --release --test io_fuzz -- --ignored
//! ```
//!
//! `FUZZ_ITERATIONS` is the number of cases per format, `FUZZ_SEED` the
//! seed (default 1), `FUZZ_TIMEOUT_SECS` the time limit per read
//! (default 30), and `FUZZ_FORMAT` (such as `graphml`) limits the run to
//! one format. `FUZZ_CASE` runs a single case and saves its input. A failing input is saved in the system temporary
//! directory and its path is printed.

extern crate knappnet as engine;

use std::sync::mpsc;
use std::thread;
use std::time::{Duration, Instant};

use engine::io::{
    read_bytes, read_str, write_string, Column, Format, Network, NetworkBuilder, ReadOptions,
    Value, WriteOptions,
};
use engine::rng::Stream;
use engine::Error;

/// Formats that can be read (CSV is handled in R).
const READABLE: [Format; 7] = [
    Format::GraphMl,
    Format::Gexf,
    Format::Gml,
    Format::Pajek,
    Format::Dl,
    Format::NodeLink,
    Format::DynetMl,
];

/// The format name each reader puts in its parse errors.
fn error_format(format: Format) -> &'static str {
    match format {
        Format::GraphMl => "GraphML",
        Format::Gexf => "GEXF",
        Format::Gml => "GML",
        Format::Pajek => "Pajek",
        Format::Dl => "UCINET DL",
        Format::NodeLink => "JSON node-link",
        Format::DynetMl => "DyNetML",
        Format::Dot => "DOT",
    }
}

fn sample(directed: bool) -> Network {
    let mut b = NetworkBuilder::new(directed);
    b.graph
        .push(("name".into(), Value::Str("fuzz <sample> & co".into())));
    b.graph.push(("year".into(), Value::Int(2024)));
    b.add_node(
        "ann",
        vec![
            ("age".into(), Value::Int(31)),
            ("score".into(), Value::Float(0.25)),
            ("lead".into(), Value::Bool(true)),
            ("label".into(), Value::Str("Ann \"A\" Lee".into())),
        ],
    );
    b.add_node(
        "bob",
        vec![
            ("age".into(), Value::Int(-4)),
            ("score".into(), Value::Float(1e300)),
        ],
    );
    b.add_node("c\u{e9}line", vec![("lead".into(), Value::Bool(false))]);
    b.add_node("d e", vec![("score".into(), Value::Float(-0.0))]);
    b.add_node("isolate", vec![]);
    let w = |x: f64| vec![("weight".to_string(), Value::Float(x))];
    b.add_edge("ann", "bob", w(2.0));
    b.add_edge("bob", "c\u{e9}line", w(0.5));
    b.add_edge("c\u{e9}line", "c\u{e9}line", w(1.0));
    b.add_edge("ann", "bob", w(3.25));
    b.add_edge("d e", "ann", w(1e-7));
    b.finish().unwrap()
}

const GRAPHML: &str = r#"<?xml version="1.0" encoding="UTF-8"?>
<!-- comment -->
<!DOCTYPE graphml [ <!ENTITY e "ignored"> ]>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns" xmlns:y="urn:y">
  <key id="d0" for="node" attr.name="color" attr.type="string"><default>gray</default></key>
  <key id="d1" for="edge" attr.name="weight" attr.type="double"/>
  <key id="d2" for="node" attr.name="age" attr.type="int"/>
  <key id="d3" for="all" attr.name="flag" attr.type="boolean"><default>true</default></key>
  <key id="g" for="graph" attr.name="title" attr.type="string"/>
  <graph id="G" edgedefault="undirected">
    <data key="g">demo &amp; <![CDATA[raw <x>]]></data>
    <node id="a"><data key="d0">red</data><data key="d2"> 30 </data><y:ShapeNode/></node>
    <node id="b"/>
    <node id="&#233;&#x41;"><data key="d3">false</data></node>
    <edge id="e1" source="a" target="b" directed="false"><data key="d1">1.5</data></edge>
    <edge source="b" target="c"/>
    <edge source="&#233;&#x41;" target="a"><data key="d1">-INF</data></edge>
  </graph>
</graphml>
"#;

const GEXF: &str = r#"<?xml version="1.0" encoding="UTF-8"?>
<gexf xmlns="http://gexf.net/1.3" xmlns:viz="http://gexf.net/1.3/viz" version="1.3">
  <graph defaultedgetype="directed" mode="dynamic">
    <attributes class="node">
      <attribute id="0" title="group" type="integer"><default>1</default></attribute>
      <attribute id="1" title="active" type="boolean"/>
    </attributes>
    <attributes class="edge">
      <attribute id="0" title="kind" type="string"/>
    </attributes>
    <nodes>
      <node id="n0" label="Alice">
        <attvalues><attvalue for="0" value="3"/><attvalue for="1" value="true"/></attvalues>
        <viz:position x="1.5" y="-2.0" z="0.0"/>
        <viz:color r="255" g="0" b="10" a="0.5"/>
        <viz:size value="3"/>
      </node>
      <node id="n1" label="Bob" start="2001" end="2003"/>
    </nodes>
    <edges>
      <edge id="e0" source="n0" target="n1" weight="2.5" label="knows">
        <attvalues><attvalue for="0" value="work"/></attvalues>
        <viz:thickness value="2"/>
      </edge>
      <edge id="e1" source="n1" target="n1" type="directed"/>
    </edges>
  </graph>
</gexf>
"#;

const GML: &str = r#"Creator "test"
graph [
  directed 1
  name "demo &quot;x&quot; &#233;"
  node [ id 1 label "a" graphics [ x 1.5 y 2 w NAN ] ]
  node [ id 2 label "b" ]
  node [ id 3 ]
  edge [ source 1 target 2 weight 0.5 ]
  # a comment
  edge [ source 2 target 2 weight -INF ]
  edge [ source 3 target 1 ]
]
"#;

const PAJEK: &str = "% comment\n*Network demo net\n*Vertices 5 2\n1 \"Ann Lee\" 0.1 0.2 0.3\n2 \"Bob\" 0.3 0.4\n3 Cy\n*Arcs\n1 2 2.5\n*Edges\n2 3\n*Arcslist\n1 3 4 5\n*Edgeslist\n4 5\n";

const PAJEK_MATRIX: &str =
    "*Vertices 3\n1 \"a\"\n2 \"b\"\n3 \"c\"\n*Matrix\n0 1 0\n2.5 0 0\n0 0 1\n";

const DL_MATRIX: &str =
    "DL N = 3\nFORMAT = FULLMATRIX\nLABELS:\nann,bob,\"cy d\"\nDATA:\n0 1 0\n1 0 2\n0 0 0\n";

const DL_EMBEDDED: &str =
    "dl n=3 format=fullmatrix labels embedded\ndata:\na b c\na 0 1 0\nb 0 0 1\nc 0 0 0\n";

const DL_LISTS: &str = "dl n=4 format=edgelist1\nlabels:\na,b,c,d\ndata:\n1 2\n2 3 0.5\n4 4\n";

const DL_NODELIST: &str = "dl format=nodelist1 labels embedded\ndata:\nx y z\ny z\nz\n";

const JSON_NX: &str = r#"{"directed": false, "multigraph": false, "graph": {"name": "g", "k": [1, {"a": null}]},
    "nodes": [{"id": 0, "club": "Aé😀"}, {"id": 1, "club": "B", "pos": [1, 2.5e-3]}, {"id": "x"}],
    "edges": [{"source": 0, "target": 1, "weight": 4}, {"source": "x", "target": 1, "weight": -1.5E2, "key": 0}]}"#;

const JSON_D3: &str = r#"{"nodes": [{"name": "a"}, {"name": "b"}, {"name": "c"}],
    "links": [{"source": 0, "target": 1, "value": true}, {"source": {"id": 2}, "target": 0}]}"#;

const DYNETML: &str = r#"<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
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
        <node id="A02" title="Bo"/>
        <node id="Kabul"/>
      </nodeclass>
      <nodeclass type="Location" id="Location">
        <node id="Kabul"/>
        <node id="Herat"/>
      </nodeclass>
    </nodes>
    <networks>
      <network sourceType="Agent" source="Agent" targetType="Agent" target="Agent" id="agent x agent" isDirected="true">
        <link source="A01" target="A02" type="double" value="2.5"/>
        <link source="A02" target="Kabul"/>
      </network>
      <network sourceType="Agent" source="Agent" targetType="Location" target="Location" id="agent x location" isDirected="false">
        <link source="A01" target="Kabul" value="1"/>
      </network>
    </networks>
  </MetaNetwork>
  <MetaNetwork id="later">
    <nodes>
      <nodeclass type="Agent" id="Agent"><node id="A01"/></nodeclass>
    </nodes>
    <networks/>
  </MetaNetwork>
</DynamicNetwork>
"#;

const DYNETML_2004: &str = r#"<DynamicNetwork>
  <MetaMatrix timePeriod="1">
    <nodes>
      <nodeset id="agents" type="Agent">
        <node id="LDR"><properties><property name="cleared" type="binary" value="true"/></properties>
          <measures><measure name="ora.degree" type="double" value="0.5"/></measures></node>
        <node id="Mgr1"/>
      </nodeset>
      <nodeset id="tasks" type="Task"><node id="T1"/></nodeset>
    </nodes>
    <networks>
      <graph sourceType="agents" targetType="agents" isDirected="true">
        <edge source="LDR" target="Mgr1" type="double" value="3"/>
      </graph>
      <graph sourceType="agents" targetType="tasks">
        <edge source="Mgr1" target="T1"/>
      </graph>
    </networks>
  </MetaMatrix>
</DynamicNetwork>
"#;

/// Valid starting points: hand-written files plus the sample network
/// written by every writer.
fn corpus() -> Vec<(Format, Vec<u8>)> {
    let mut out: Vec<(Format, Vec<u8>)> = vec![
        (Format::GraphMl, GRAPHML.into()),
        (Format::Gexf, GEXF.into()),
        (Format::Gml, GML.into()),
        (Format::Pajek, PAJEK.into()),
        (Format::Pajek, PAJEK_MATRIX.into()),
        (Format::Dl, DL_MATRIX.into()),
        (Format::Dl, DL_EMBEDDED.into()),
        (Format::Dl, DL_LISTS.into()),
        (Format::Dl, DL_NODELIST.into()),
        (Format::NodeLink, JSON_NX.into()),
        (Format::NodeLink, JSON_D3.into()),
        (Format::DynetMl, DYNETML.into()),
        (Format::DynetMl, DYNETML_2004.into()),
    ];
    for directed in [true, false] {
        let net = sample(directed);
        for format in READABLE {
            let text = write_string(format, &net, &WriteOptions::default()).unwrap();
            out.push((format, text.into_bytes()));
        }
    }
    // Every starting point must itself be readable.
    for (format, bytes) in &out {
        if let Err(e) = read_bytes(*format, bytes, &ReadOptions::default()) {
            panic!(
                "corpus file for {format:?} does not read: {e}\n{}",
                String::from_utf8_lossy(bytes)
            );
        }
    }
    out
}

const INTERESTING: &[u8] = b"<>\"'&;=[]{}:,/\\\n\r\t #*%!?-.+eE019\0\x7f\x80\xc3\xff";

const NUMBERS: &[&str] = &[
    "99999999999999999999999",
    "18446744073709551616",
    "9223372036854775807",
    "-9223372036854775808",
    "4294967296",
    "-1",
    "0",
    "-0",
    "1e309",
    "-1e309",
    "1e-400",
    "NaN",
    "inf",
    "0x10",
    "1.5.5",
    "--3",
    "",
];

const KEYWORDS: &[&str] = &[
    "<graph>",
    "</graph>",
    "<node id=\"x\"/>",
    "<node>",
    "</node>",
    "<edge source=\"a\" target=\"zz\"/>",
    "<edge/>",
    "<data key=\"d1\">nan</data>",
    "<data>",
    "<key id=\"d1\" for=\"edge\" attr.name=\"weight\" attr.type=\"double\"/>",
    "<default>",
    "<hyperedge/>",
    "<port name=\"p\"/>",
    "<attvalue for=\"0\" value=\"x\"/>",
    "<nodes>",
    "<edges>",
    "<nodeclass id=\"Agent\" type=\"Agent\">",
    "<network source=\"Agent\" target=\"Agent\">",
    "<link source=\"A01\" target=\"A01\" value=\"x\"/>",
    "<MetaNetwork>",
    "<property id=\"age\" value=\"x\"/>",
    "<![CDATA[",
    "]]>",
    "<!--",
    "-->",
    "<?",
    "?>",
    "<!DOCTYPE x [",
    "&#x110000;",
    "&#0;",
    "&#xD800;",
    "&amp",
    "&",
    "edgedefault=\"mutual\"",
    "directed=\"true\"",
    "xmlns:a=\"b\" a:c=\"d\"",
    "graph [",
    "node [ id 1 ]",
    "edge [ source 1 target 1 ]",
    " [ ",
    " ] ",
    "\"",
    "label \"x\"",
    "*Vertices 3",
    "*Vertices 2 5",
    "*Arcs",
    "*Edges",
    "*Arcslist",
    "*Edgeslist",
    "*Matrix",
    "*Network",
    "*Partition",
    "data:",
    "labels:",
    "labels embedded",
    "format=fullmatrix",
    "format=edgelist1",
    "format=nodelist1",
    "n=3",
    "n=",
    "nr=2",
    "dl",
    "{\"id\": 1}",
    "\"nodes\": []",
    "\"links\": [{\"source\": 9, \"target\": -1}]",
    "\"directed\": true,",
    "\\u",
    "\\ud800",
    "[",
    "{",
    ":",
    ",",
];

/// A mutated copy of `base`. `others` supplies material for splices.
fn mutate(rng: &mut Stream, base: &[u8], others: &[(Format, Vec<u8>)]) -> Vec<u8> {
    let mut s = base.to_vec();
    let rounds = 1 + rng.below(4);
    for _ in 0..rounds {
        let pos = |rng: &mut Stream, s: &Vec<u8>| rng.below(s.len() as u64 + 1) as usize;
        match rng.below(14) {
            0 if !s.is_empty() => {
                let i = rng.below(s.len() as u64) as usize;
                s[i] ^= 1 << rng.below(8);
            }
            1 if !s.is_empty() => {
                let i = rng.below(s.len() as u64) as usize;
                s[i] = INTERESTING[rng.below(INTERESTING.len() as u64) as usize];
            }
            2 if !s.is_empty() => {
                let a = pos(rng, &s);
                let span = if rng.below(4) == 0 { 256 } else { 16 };
                let len = 1 + rng.below(span) as usize;
                let b = (a + len).min(s.len());
                s.drain(a..b);
            }
            3 if !s.is_empty() => {
                let a = pos(rng, &s);
                let len = 1 + rng.below(64) as usize;
                let b = (a + len).min(s.len());
                let piece = s[a..b].to_vec();
                let at = pos(rng, &s);
                s.splice(at..at, piece);
            }
            4 => {
                let at = pos(rng, &s);
                let len = 1 + rng.below(8) as usize;
                let piece: Vec<u8> = (0..len).map(|_| rng.below(256) as u8).collect();
                s.splice(at..at, piece);
            }
            5 => {
                let at = pos(rng, &s);
                s.truncate(at);
            }
            6 => {
                let other = &others[rng.below(others.len() as u64) as usize].1;
                let a = pos(rng, &s);
                let b = rng.below(other.len() as u64 + 1) as usize;
                s.truncate(a);
                s.extend_from_slice(&other[b..]);
            }
            7 => {
                // Replace a run of digits (or insert) with an odd number.
                let num = NUMBERS[rng.below(NUMBERS.len() as u64) as usize].as_bytes();
                let digits: Vec<usize> = (0..s.len()).filter(|&i| s[i].is_ascii_digit()).collect();
                if digits.is_empty() {
                    let at = pos(rng, &s);
                    s.splice(at..at, num.iter().copied());
                } else {
                    let a = digits[rng.below(digits.len() as u64) as usize];
                    let mut b = a;
                    while b < s.len() && (s[b].is_ascii_digit() || s[b] == b'.') {
                        b += 1;
                    }
                    s.splice(a..b, num.iter().copied());
                }
            }
            8 => {
                // Deep nesting.
                let unit: &[u8] = match rng.below(6) {
                    0 => b"<a>",
                    1 => b"<a x='1'>",
                    2 => b"k [ ",
                    3 => b"[",
                    4 => b"{\"a\":",
                    _ => b"<node id='n'><graph>",
                };
                let count = 1 + rng.below(200_000 / unit.len() as u64) as usize;
                let at = pos(rng, &s);
                s.splice(at..at, unit.repeat(count));
            }
            9 => {
                const BAD: &[&[u8]] = &[
                    b"\x80",
                    b"\xc3\x28",
                    b"\xff\xfe",
                    b"\xed\xa0\x80",
                    b"\xf4\x90\x80\x80",
                    b"\xe2\x82",
                    b"\xef\xbb\xbf",
                ];
                let at = pos(rng, &s);
                let bad = BAD[rng.below(BAD.len() as u64) as usize];
                s.splice(at..at, bad.iter().copied());
            }
            10 if !s.is_empty() => {
                // Repeat one line many times.
                let at = rng.below(s.len() as u64) as usize;
                let start = s[..at]
                    .iter()
                    .rposition(|&b| b == b'\n')
                    .map_or(0, |i| i + 1);
                let end = s[at..]
                    .iter()
                    .position(|&b| b == b'\n')
                    .map_or(s.len(), |i| at + i + 1);
                let line = s[start..end].to_vec();
                let times = 1 + rng.below(300) as usize;
                s.splice(end..end, line.repeat(times));
            }
            11 => {
                let word = KEYWORDS[rng.below(KEYWORDS.len() as u64) as usize].as_bytes();
                let at = pos(rng, &s);
                s.splice(at..at, word.iter().copied());
            }
            12 if !s.is_empty() => {
                // Swap the case of a run of letters.
                let a = pos(rng, &s);
                let b = (a + 1 + rng.below(12) as usize).min(s.len());
                for c in &mut s[a..b] {
                    if c.is_ascii_alphabetic() {
                        *c ^= 0x20;
                    }
                }
            }
            _ => {
                // Pure random bytes, occasionally.
                if rng.below(8) == 0 {
                    let len = rng.below(512) as usize;
                    s = (0..len).map(|_| rng.below(256) as u8).collect();
                }
            }
        }
        // Keep inputs bounded so the time limit bounds work per byte.
        s.truncate(400_000);
    }
    s
}

fn check_tables(net: &Network) -> Result<(), String> {
    for (what, table, first) in [
        ("nodes", &net.nodes, &["id"][..]),
        ("edges", &net.edges, &["from", "to"][..]),
    ] {
        if table.names.len() != table.columns.len() {
            return Err(format!("{what}: names and columns differ in number"));
        }
        let rows = table.nrow();
        for (name, col) in table.names.iter().zip(&table.columns) {
            if col.len() != rows {
                return Err(format!("{what}: column {name:?} has {} rows", col.len()));
            }
        }
        for (i, name) in first.iter().enumerate() {
            if table.names.get(i).map(String::as_str) != Some(*name) {
                return Err(format!("{what}: column {i} is not {name:?}"));
            }
            if !matches!(table.columns[i], Column::Str(_)) {
                return Err(format!("{what}: column {name:?} is not text"));
            }
        }
    }
    Ok(())
}

/// What a reader must guarantee about one input.
fn check(format: Format, input: &[u8], result: engine::Result<Network>) -> Result<(), String> {
    match result {
        Ok(net) => {
            net.validate()
                .map_err(|e| format!("invalid network: {e}"))?;
            check_tables(&net)?;
            for out in READABLE.into_iter().chain([Format::Dot]) {
                // Writers may refuse (DyNetML cannot hold every network),
                // but must not panic; GraphML holds everything.
                let written = write_string(out, &net, &WriteOptions::default());
                if out == Format::GraphMl {
                    let text = written.map_err(|e| format!("GraphML writer failed: {e}"))?;
                    let back = read_str(Format::GraphMl, &text, &ReadOptions::default())
                        .map_err(|e| format!("GraphML round trip failed: {e}"))?;
                    if back.directed != net.directed
                        || back.nodes.strings("id").ok() != net.nodes.strings("id").ok()
                        || back.edges.strings("from").ok() != net.edges.strings("from").ok()
                        || back.edges.strings("to").ok() != net.edges.strings("to").ok()
                    {
                        return Err("GraphML round trip changed the structure".into());
                    }
                }
            }
            Ok(())
        }
        Err(e) => {
            let message = e.to_string();
            if message.trim().is_empty() {
                return Err("error without a message".into());
            }
            if let Error::Parse {
                format: name, line, ..
            } = &e
            {
                if *name != error_format(format) {
                    return Err(format!("error names format {name:?}: {message}"));
                }
                if let Some(line) = line {
                    let text = engine::io::decode(input.to_vec());
                    let lines = text.matches('\n').count() + 1;
                    if *line == 0 || *line > lines {
                        return Err(format!(
                            "error at line {line} of a {lines}-line input: {message}"
                        ));
                    }
                }
            }
            Ok(())
        }
    }
}

/// Whether a Pajek or DL input declares more than a million nodes. Both
/// formats create isolates from a count alone, so reading such a file is
/// legitimately proportional to the count (up to
/// `MAX_DECLARED_NODES`); these inputs are skipped so a long run stays
/// within a CI machine's memory.
fn declares_many_nodes(input: &[u8]) -> bool {
    let text = engine::io::decode(input.to_vec()).to_ascii_lowercase();
    let big = |rest: &str| {
        let digits: String = rest
            .trim_start()
            .chars()
            .take_while(char::is_ascii_digit)
            .collect();
        digits.len() > 6
    };
    text.match_indices("*vertices")
        .any(|(i, k)| big(&text[i + k.len()..]))
        || text.match_indices('n').any(|(i, _)| {
            let rest = text[i + 1..].trim_start();
            rest.starts_with('=') && big(&rest[1..])
        })
}

fn env_u64(name: &str, default: u64) -> u64 {
    std::env::var(name)
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(default)
}

enum Outcome {
    Done(engine::Result<Network>, Duration),
    Panicked,
    TimedOut,
}

/// Runs one read on a worker thread with a time limit.
fn read_with_limit(format: Format, input: Vec<u8>, opts: ReadOptions, limit: Duration) -> Outcome {
    let (tx, rx) = mpsc::channel();
    let start = Instant::now();
    thread::Builder::new()
        .name("reader".into())
        .spawn(move || {
            let result = std::panic::catch_unwind(|| read_bytes(format, &input, &opts));
            let _ = tx.send(result);
        })
        .expect("spawn reader thread");
    match rx.recv_timeout(limit) {
        Ok(Ok(result)) => Outcome::Done(result, start.elapsed()),
        Ok(Err(_)) => Outcome::Panicked,
        Err(_) => Outcome::TimedOut,
    }
}

fn save(format: Format, case: u64, input: &[u8]) -> String {
    let path = std::env::temp_dir().join(format!(
        "fuzz-failure-{}-{case}.{}",
        format.name(),
        format.name()
    ));
    let _ = std::fs::write(&path, input);
    path.display().to_string()
}

fn fuzz(iterations: u64) {
    let seed = env_u64("FUZZ_SEED", 1);
    let limit = Duration::from_secs(env_u64("FUZZ_TIMEOUT_SECS", 30));
    let corpus = corpus();
    let only = std::env::var("FUZZ_FORMAT").ok().filter(|f| !f.is_empty());
    for (f, format) in READABLE.into_iter().enumerate() {
        if only.as_deref().is_some_and(|name| name != format.name()) {
            continue;
        }
        let bases: Vec<&Vec<u8>> = corpus
            .iter()
            .filter(|(g, _)| *g == format)
            .map(|(_, b)| b)
            .collect();
        let (mut ok, mut errors, mut skipped) = (0u64, 0u64, 0u64);
        let (mut slowest, mut slowest_case) = (Duration::ZERO, 0);
        let cases = match std::env::var("FUZZ_CASE").ok().and_then(|c| c.parse().ok()) {
            Some(case) => case..case + 1,
            None => 0..iterations,
        };
        for case in cases {
            let mut rng = Stream::new(seed, (f as u64) << 40 | case);
            let base = bases[rng.below(bases.len() as u64) as usize];
            let input = mutate(&mut rng, base, &corpus);
            if std::env::var("FUZZ_CASE").is_ok() {
                println!("input saved to {}", save(format, case, &input));
            }
            if matches!(format, Format::Pajek | Format::Dl) && declares_many_nodes(&input) {
                skipped += 1;
                continue;
            }
            let opts = ReadOptions {
                directed: match rng.below(3) {
                    0 => None,
                    1 => Some(true),
                    _ => Some(false),
                },
            };
            let (result, took) = match read_with_limit(format, input.clone(), opts, limit) {
                Outcome::Done(result, took) => (result, took),
                Outcome::Panicked => {
                    let path = save(format, case, &input);
                    panic!("{format:?} reader panicked on case {case} (seed {seed}); input saved to {path}");
                }
                Outcome::TimedOut => {
                    let path = save(format, case, &input);
                    panic!("{format:?} reader took more than {limit:?} on case {case} (seed {seed}); input saved to {path}");
                }
            };
            if took > slowest {
                slowest = took;
                slowest_case = case;
            }
            if result.is_ok() {
                ok += 1;
            } else {
                errors += 1;
            }
            if let Err(problem) = check(format, &input, result) {
                let path = save(format, case, &input);
                panic!("{format:?} case {case} (seed {seed}): {problem}; input saved to {path}");
            }
        }
        println!(
            "{format:?}: {iterations} cases, {ok} read, {errors} rejected, {skipped} skipped, \
             slowest {slowest:?} (case {slowest_case})"
        );
    }
}

#[test]
fn readers_survive_mutated_files() {
    fuzz(env_u64("FUZZ_QUICK_ITERATIONS", 150));
}

#[test]
#[ignore = "long fuzzing run; set FUZZ_ITERATIONS"]
fn readers_survive_mutated_files_long() {
    fuzz(env_u64("FUZZ_ITERATIONS", 20_000));
}
