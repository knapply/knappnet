//! Inputs that once crashed or overloaded a file reader (most found by
//! `io_fuzz.rs`). Each must now give a clear error, or read correctly.

extern crate knappnet as engine;

use engine::io::{read_bytes, read_str, Format, ReadOptions};

fn read(format: Format, text: &str) -> engine::Result<engine::io::Network> {
    read_str(format, text, &ReadOptions::default())
}

fn error(format: Format, text: &str) -> String {
    match read(format, text) {
        Ok(_) => panic!("{format:?} accepted {text:?}"),
        Err(e) => e.to_string(),
    }
}

#[test]
fn deeply_nested_xml_is_rejected_without_exhausting_the_stack() {
    // Each level used to be a recursive call: 100,000 levels overflowed
    // the stack and aborted the process.
    for (format, root) in [
        (Format::GraphMl, "<graphml>"),
        (Format::Gexf, "<gexf>"),
        (Format::DynetMl, "<DynamicNetwork>"),
    ] {
        let doc = format!("{root}\n{}", "<a>".repeat(100_000));
        let e = error(format, &doc);
        assert!(e.contains("nested more than 512 levels deep"), "{e}");
        assert!(e.contains("line 2"), "{e}");
    }
}

#[test]
fn deeply_nested_gml_is_rejected_without_exhausting_the_stack() {
    let doc = format!("graph [ {}", "x [ ".repeat(100_000));
    let e = error(Format::Gml, &doc);
    assert!(e.contains("nested more than 512 levels deep"), "{e}");
    // Ordinary nesting still reads.
    let n = read(Format::Gml, "graph [ node [ id 1 a [ b [ c [ d 1 ] ] ] ] ]").unwrap();
    assert!(n.nodes.column("a.b.c.d").is_some());
}

#[test]
fn pajek_matrix_rows_cannot_be_longer_than_the_network() {
    // A fourth value in a row of a 3 x 3 matrix indexed past the vertices.
    let e = error(Format::Pajek, "*Vertices 3\n*Matrix\n0 1 0\n2.5 0 0 1\n");
    assert!(e.contains("line 4"), "{e}");
    assert!(e.contains("3 columns"), "{e}");
}

#[test]
fn pajek_two_mode_matrices_have_the_declared_shape() {
    // Two-mode: 2 rows (first mode) by 1 column (second mode).
    let n = read(Format::Pajek, "*Vertices 3 2\n*Matrix\n1\n0\n").unwrap();
    assert_eq!(n.edges.strings("from").unwrap(), vec!["1"]);
    assert_eq!(n.edges.strings("to").unwrap(), vec!["3"]);
    let e = error(Format::Pajek, "*Vertices 3 2\n*Matrix\n1 1\n");
    assert!(e.contains("1 columns"), "{e}");
    let e = error(Format::Pajek, "*Vertices 3 2\n*Matrix\n1\n0\n1\n");
    assert!(e.contains("too many matrix rows"), "{e}");
    // A first mode larger than the network.
    let e = error(Format::Pajek, "*Vertices 2 5\n*Matrix\n0 1\n");
    assert!(e.contains("first mode"), "{e}");
}

#[test]
fn huge_declared_node_counts_are_rejected_before_allocating() {
    // These used to abort with "capacity overflow" or exhaust memory.
    let e = error(Format::Pajek, "*Vertices 9223372036854775807\n*Arcs\n1 2\n");
    assert!(e.contains("declares more than"), "{e}");
    let e = error(
        Format::Dl,
        "dl n=9223372036854775807 format=edgelist1\ndata:\n1 2\n",
    );
    assert!(e.contains("declares more than"), "{e}");
    // A full matrix is checked against the data before anything is built.
    let e = error(Format::Dl, "dl n=100000000\ndata:\n0 1\n1 0\n");
    assert!(e.contains("too few values"), "{e}");
    let e = error(
        Format::Dl,
        "dl n=3 format=fullmatrix labels embedded\ndata:\na b c\na 0 1 0\n",
    );
    assert!(e.contains("too few values"), "{e}");
}

#[test]
fn unlabelled_pajek_vertices_are_read_in_linear_time() {
    // Whether any vertex had a label was checked once per vertex: a file
    // declaring 200,000 unlabelled vertices took 4 x 10^10 steps.
    let n = read(Format::Pajek, "*Vertices 200000\n*Arcs\n1 2\n").unwrap();
    assert_eq!(n.nodes.nrow(), 200_000);
    assert!(n.nodes.column("label").is_none());
}

#[test]
fn a_second_pajek_vertices_section_is_rejected() {
    // A smaller second *Vertices left earlier links pointing past the end.
    let e = error(Format::Pajek, "*Vertices 5\n*Arcs\n1 5\n*Vertices 2\n");
    assert!(e.contains("line 4"), "{e}");
    assert!(e.contains("second *Vertices"), "{e}");
}

#[test]
fn invalid_utf8_is_read_as_latin1() {
    let bytes = b"*Vertices 2\n1 \"caf\xe9\"\n2 \"b\"\n*Edges\n1 2\n";
    let n = read_bytes(Format::Pajek, bytes, &ReadOptions::default()).unwrap();
    assert_eq!(n.nodes.strings("id").unwrap(), vec!["caf\u{e9}", "b"]);
}

#[test]
fn repeated_dl_labels_do_not_merge_nodes() {
    // Two nodes labelled "a" used to become one node, folding two matrix
    // rows into one.
    let n = read(
        Format::Dl,
        "dl n=3\nlabels:\na,a,b\ndata:\n0 1 0\n0 0 1\n1 0 0\n",
    )
    .unwrap();
    assert_eq!(n.nodes.strings("id").unwrap(), vec!["1", "2", "3"]);
    assert_eq!(n.nodes.strings("label").unwrap(), vec!["a", "a", "b"]);
    assert_eq!(n.edges.strings("from").unwrap(), vec!["1", "2", "3"]);
    assert_eq!(n.edges.strings("to").unwrap(), vec!["2", "3", "1"]);
    // Unique labels are still identifiers.
    let u = read(
        Format::Dl,
        "dl n=2 format=edgelist1\nlabels:\nx,y\ndata:\n1 2\n",
    )
    .unwrap();
    assert_eq!(u.nodes.strings("id").unwrap(), vec!["x", "y"]);
    assert!(u.nodes.column("label").is_none());
}
