//! A minimal XML reader and writer, sufficient for GraphML, GEXF, and
//! DyNetML.
//!
//! Supports elements, attributes, text, CDATA, comments, processing
//! instructions, a DOCTYPE declaration, and the predefined and numeric
//! character entities. Namespace prefixes are dropped from element and
//! attribute names. External entities are never resolved, so documents
//! cannot cause file or network access.
//!
//! [`Reader`] streams through a document as a sequence of [`Event`]s (start
//! tags with their attributes, end tags, and text) without building a tree,
//! in time linear in the length of the document. [`parse`] builds a tree of
//! [`Element`]s from the same events, for formats that are easier to read
//! as a tree. Both report the first problem in document order, with its
//! line. Elements may be nested at most [`MAX_DEPTH`] levels deep, which
//! bounds the memory and stack a hostile document can demand.

use std::borrow::Cow;

use crate::error::{Error, Result};

/// The deepest nesting of elements accepted. Network files nest a handful
/// of levels; the limit only stops pathological input.
pub const MAX_DEPTH: usize = 512;

/// An XML element.
#[derive(Debug, Clone, Default)]
pub struct Element {
    /// Local name (namespace prefix removed).
    pub name: String,
    /// Attributes in document order (namespace prefixes removed).
    pub attrs: Vec<(String, String)>,
    /// Child elements in document order.
    pub children: Vec<Element>,
    /// Concatenated text content directly inside this element.
    pub text: String,
    /// 1-based line where the element starts.
    pub line: usize,
}

impl Element {
    /// The value of attribute `name`.
    pub fn attr(&self, name: &str) -> Option<&str> {
        self.attrs
            .iter()
            .find(|(k, _)| k == name)
            .map(|(_, v)| v.as_str())
    }

    /// Child elements called `name`.
    pub fn children_named<'a>(&'a self, name: &'a str) -> impl Iterator<Item = &'a Element> {
        self.children.iter().filter(move |c| c.name == name)
    }

    /// The first child element called `name`.
    pub fn child(&self, name: &str) -> Option<&Element> {
        self.children.iter().find(|c| c.name == name)
    }
}

/// One step through a document.
#[derive(Debug, Clone, PartialEq)]
pub enum Event<'a> {
    /// A start tag, with its local name and the line where it starts. Its
    /// attributes are available from [`Reader::attrs`] until the next call
    /// to [`Reader::next`]. An empty-element tag (`<a/>`) is reported as a
    /// start tag with `empty` set, immediately followed by its end tag.
    Start {
        /// Local name.
        name: &'a str,
        /// 1-based line of the `<`.
        line: usize,
        /// Whether the tag was an empty-element tag.
        empty: bool,
    },
    /// An end tag (or the end of an empty-element tag).
    End {
        /// Local name.
        name: &'a str,
    },
    /// Text directly inside the current element, with entities decoded; a
    /// CDATA section is reported verbatim. Consecutive pieces of text may
    /// arrive as several events.
    Text(Cow<'a, str>),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum State {
    Prolog,
    Content,
    Epilog,
    Done,
}

/// A streaming XML reader. Call [`Reader::next`] until it returns
/// `Ok(None)`; only then has the whole document been checked.
#[derive(Debug)]
pub struct Reader<'a> {
    s: &'a [u8],
    src: &'a str,
    pos: usize,
    format: &'static str,
    // The 1-based line number at byte offset `line_pos`.
    line: usize,
    line_pos: usize,
    stack: Vec<&'a str>,
    attrs: Vec<(&'a str, Cow<'a, str>)>,
    state: State,
    pending_end: bool,
}

impl<'a> Reader<'a> {
    /// A reader over `input`; errors name `format`.
    pub fn new(input: &'a str, format: &'static str) -> Reader<'a> {
        Reader {
            s: input.as_bytes(),
            src: input,
            pos: 0,
            format,
            line: 1,
            line_pos: 0,
            stack: Vec::new(),
            attrs: Vec::new(),
            state: State::Prolog,
            pending_end: false,
        }
    }

    /// Attributes of the most recent start tag, in document order, with
    /// namespace prefixes removed from their names.
    pub fn attrs(&self) -> &[(&'a str, Cow<'a, str>)] {
        &self.attrs
    }

    /// The first attribute of the most recent start tag called `name`.
    pub fn attr(&self, name: &str) -> Option<&str> {
        self.attrs
            .iter()
            .find(|(k, _)| *k == name)
            .map(|(_, v)| v.as_ref())
    }

    /// A copy of the first attribute of the most recent start tag called
    /// `name`, borrowing from the document when it has no entities.
    pub fn attr_cow(&self, name: &str) -> Option<Cow<'a, str>> {
        self.attrs
            .iter()
            .find(|(k, _)| *k == name)
            .map(|(_, v)| v.clone())
    }

    /// Number of elements currently open (1 inside the root element,
    /// counting the element just started).
    pub fn depth(&self) -> usize {
        self.stack.len()
    }

    /// The line at byte offset `pos`. Offsets are requested in increasing
    /// order, so lines are counted once overall.
    fn line_at(&mut self, pos: usize) -> usize {
        if pos < self.line_pos {
            return 1 + self.s[..pos].iter().filter(|&&b| b == b'\n').count();
        }
        self.line += self.s[self.line_pos..pos]
            .iter()
            .filter(|&&b| b == b'\n')
            .count();
        self.line_pos = pos;
        self.line
    }

    fn error_at(&mut self, pos: usize, message: impl Into<String>) -> Error {
        Error::Parse {
            format: self.format,
            line: Some(self.line_at(pos)),
            message: message.into(),
        }
    }

    fn err<T>(&mut self, message: impl Into<String>) -> Result<T> {
        Err(self.error_at(self.pos, message))
    }

    fn starts_with(&self, pat: &str) -> bool {
        self.s[self.pos..].starts_with(pat.as_bytes())
    }

    fn skip_ws(&mut self) {
        while self.pos < self.s.len() && self.s[self.pos].is_ascii_whitespace() {
            self.pos += 1;
        }
    }

    fn skip_until(&mut self, pat: &str) -> Result<()> {
        match self.src[self.pos..].find(pat) {
            Some(i) => {
                self.pos += i + pat.len();
                Ok(())
            }
            None => self.err(format!("unterminated construct, expected {pat:?}")),
        }
    }

    /// Skips a processing instruction, comment, or DOCTYPE; returns whether
    /// one was skipped.
    fn skip_misc(&mut self) -> Result<bool> {
        if self.starts_with("<?") {
            self.skip_until("?>")?;
        } else if self.starts_with("<!--") {
            self.skip_until("-->")?;
        } else if self.starts_with("<!DOCTYPE") || self.starts_with("<!doctype") {
            // Skip an internal subset if present.
            let rest = &self.src[self.pos..];
            let close = rest.find('>').unwrap_or(rest.len());
            if let Some(open) = rest[..close].find('[') {
                self.pos += open;
                self.skip_until("]")?;
                self.skip_until(">")?;
            } else {
                self.skip_until(">")?;
            }
        } else {
            return Ok(false);
        }
        Ok(true)
    }

    fn parse_name(&mut self) -> Result<&'a str> {
        let start = self.pos;
        while self.pos < self.s.len() {
            let b = self.s[self.pos];
            if b.is_ascii_whitespace() || b == b'>' || b == b'/' || b == b'=' {
                break;
            }
            self.pos += 1;
        }
        if start == self.pos {
            return self.err("expected a name");
        }
        // Names start after an ASCII byte and stop at one (or at the end),
        // so both ends are character boundaries.
        let src: &'a str = self.src;
        Ok(local(&src[start..self.pos]))
    }

    fn start_tag(&mut self) -> Result<Event<'a>> {
        let line = self.line_at(self.pos);
        if self.stack.len() >= MAX_DEPTH {
            return self.err(format!(
                "elements are nested more than {MAX_DEPTH} levels deep"
            ));
        }
        self.pos += 1; // '<'
        let name = self.parse_name()?;
        self.attrs.clear();
        let src: &'a str = self.src;
        loop {
            self.skip_ws();
            if self.starts_with("/>") {
                self.pos += 2;
                self.stack.push(name);
                self.pending_end = true;
                return Ok(Event::Start {
                    name,
                    line,
                    empty: true,
                });
            }
            if self.starts_with(">") {
                self.pos += 1;
                self.stack.push(name);
                return Ok(Event::Start {
                    name,
                    line,
                    empty: false,
                });
            }
            if self.pos >= self.s.len() {
                return self.err(format!("unterminated <{name}> tag"));
            }
            let key = self.parse_name()?;
            self.skip_ws();
            if !self.starts_with("=") {
                return self.err(format!("attribute {key:?} has no value"));
            }
            self.pos += 1;
            self.skip_ws();
            let quote = match self.s.get(self.pos) {
                Some(&q @ (b'"' | b'\'')) => q,
                _ => return self.err(format!("attribute {key:?} value must be quoted")),
            };
            self.pos += 1;
            let end = match self.s[self.pos..].iter().position(|&b| b == quote) {
                Some(e) => e,
                None => return self.err(format!("unterminated value for attribute {key:?}")),
            };
            let raw = &src[self.pos..self.pos + end];
            let value = match decode_entities(raw) {
                Ok(v) => v,
                Err(m) => return self.err(m),
            };
            self.pos += end + 1;
            self.attrs.push((key, value));
        }
    }

    fn close(&mut self, name: &'a str) -> Event<'a> {
        self.stack.pop();
        if self.stack.is_empty() {
            self.state = State::Epilog;
        }
        Event::End { name }
    }

    /// The next event, `Ok(None)` at the end of a well-formed document, or
    /// the first problem found.
    #[allow(clippy::should_implement_trait)]
    pub fn next(&mut self) -> Result<Option<Event<'a>>> {
        if self.pending_end {
            self.pending_end = false;
            let name = self.stack.last().copied().unwrap_or_default();
            return Ok(Some(self.close(name)));
        }
        let src: &'a str = self.src;
        match self.state {
            State::Prolog => {
                if self.starts_with("\u{feff}") {
                    self.pos += 3;
                }
                loop {
                    self.skip_ws();
                    if !self.skip_misc()? {
                        break;
                    }
                }
                if !self.starts_with("<") {
                    return self.err("expected the root element");
                }
                self.state = State::Content;
                self.start_tag().map(Some)
            }
            State::Content => loop {
                let open = self.stack.last().copied().unwrap_or_default();
                if self.pos >= self.s.len() {
                    return self.err(format!("<{open}> is never closed"));
                }
                if self.s[self.pos] == b'<' {
                    if self.starts_with("</") {
                        self.pos += 2;
                        let close = self.parse_name()?;
                        if close != open {
                            return self.err(format!("expected </{open}>, found </{close}>"));
                        }
                        self.skip_ws();
                        if !self.starts_with(">") {
                            return self.err("malformed closing tag");
                        }
                        self.pos += 1;
                        return Ok(Some(self.close(open)));
                    }
                    if self.starts_with("<![CDATA[") {
                        self.pos += 9;
                        let end = match src[self.pos..].find("]]>") {
                            Some(e) => e,
                            None => return self.err("unterminated CDATA section"),
                        };
                        let text = &src[self.pos..self.pos + end];
                        self.pos += end + 3;
                        return Ok(Some(Event::Text(Cow::Borrowed(text))));
                    }
                    if self.skip_misc()? {
                        continue;
                    }
                    return self.start_tag().map(Some);
                }
                let end = self.s[self.pos..]
                    .iter()
                    .position(|&b| b == b'<')
                    .unwrap_or(self.s.len() - self.pos);
                let raw = &src[self.pos..self.pos + end];
                let text = match decode_entities(raw) {
                    Ok(t) => t,
                    Err(m) => return self.err(m),
                };
                self.pos += end;
                return Ok(Some(Event::Text(text)));
            },
            State::Epilog => {
                loop {
                    self.skip_ws();
                    if self.pos >= self.s.len() || !self.skip_misc()? {
                        break;
                    }
                }
                if self.pos < self.s.len() {
                    return self.err("unexpected content after the root element");
                }
                self.state = State::Done;
                Ok(None)
            }
            State::Done => Ok(None),
        }
    }
}

/// Parses a document and returns its root element.
pub fn parse(input: &str, format: &'static str) -> Result<Element> {
    let mut r = Reader::new(input, format);
    let mut open: Vec<Element> = Vec::new();
    let mut root = None;
    while let Some(event) = r.next()? {
        match event {
            Event::Start { name, line, .. } => open.push(Element {
                name: name.to_string(),
                attrs: r
                    .attrs()
                    .iter()
                    .map(|(k, v)| (k.to_string(), v.to_string()))
                    .collect(),
                children: Vec::new(),
                text: String::new(),
                line,
            }),
            Event::End { .. } => {
                if let Some(el) = open.pop() {
                    match open.last_mut() {
                        Some(parent) => parent.children.push(el),
                        None => root = Some(el),
                    }
                }
            }
            Event::Text(t) => {
                if let Some(el) = open.last_mut() {
                    el.text.push_str(&t);
                }
            }
        }
    }
    root.ok_or(Error::Parse {
        format,
        line: None,
        message: "expected the root element".into(),
    })
}

fn local(name: &str) -> &str {
    // Names are short; a plain scan beats a general substring search.
    match name.bytes().rposition(|b| b == b':') {
        Some(i) => &name[i + 1..],
        None => name,
    }
}

fn decode_entities(s: &str) -> std::result::Result<Cow<'_, str>, String> {
    if !s.contains('&') {
        return Ok(Cow::Borrowed(s));
    }
    let mut out = String::with_capacity(s.len());
    let mut rest = s;
    while let Some(i) = rest.find('&') {
        out.push_str(&rest[..i]);
        rest = &rest[i + 1..];
        let end = rest
            .find(';')
            .ok_or_else(|| "unterminated character entity".to_string())?;
        let ent = &rest[..end];
        let ch = match ent {
            "lt" => '<',
            "gt" => '>',
            "amp" => '&',
            "quot" => '"',
            "apos" => '\'',
            _ if ent.starts_with("#x") || ent.starts_with("#X") => {
                u32::from_str_radix(&ent[2..], 16)
                    .ok()
                    .and_then(char::from_u32)
                    .ok_or_else(|| format!("invalid character reference &{ent};"))?
            }
            _ if ent.starts_with('#') => ent[1..]
                .parse::<u32>()
                .ok()
                .and_then(char::from_u32)
                .ok_or_else(|| format!("invalid character reference &{ent};"))?,
            _ => return Err(format!("unknown entity &{ent};")),
        };
        out.push(ch);
        rest = &rest[end + 1..];
    }
    out.push_str(rest);
    Ok(Cow::Owned(out))
}

/// Escapes text for use in element content or a double-quoted attribute.
pub fn escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    escape_into(&mut out, s);
    out
}

/// Appends `s` to `out`, escaped as by [`escape`].
pub fn escape_into(out: &mut String, s: &str) {
    let mut start = 0;
    for (i, b) in s.bytes().enumerate() {
        let rep = match b {
            b'<' => "&lt;",
            b'>' => "&gt;",
            b'&' => "&amp;",
            b'"' => "&quot;",
            b'\n' => "&#10;",
            b'\r' => "&#13;",
            b'\t' => "&#9;",
            _ => continue,
        };
        out.push_str(&s[start..i]);
        out.push_str(rep);
        start = i + 1;
    }
    out.push_str(&s[start..]);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Stream;

    #[test]
    fn parses_typical_document() {
        let doc = r#"<?xml version="1.0" encoding="UTF-8"?>
<!-- comment -->
<!DOCTYPE x [ <!ENTITY e "ignored"> ]>
<g:root xmlns:g="urn:x" a='1' b="x &amp; y">
  <child id="c1">text &lt;here&gt; &#65;&#x42;</child>
  <empty/>
  <![CDATA[raw <stuff>]]>
</g:root>"#;
        let root = parse(doc, "XML").unwrap();
        assert_eq!(root.name, "root");
        assert_eq!(root.attr("b"), Some("x & y"));
        assert_eq!(root.children.len(), 2);
        assert_eq!(root.children[0].text, "text <here> AB");
        assert_eq!(root.children[1].name, "empty");
        assert!(root.text.contains("raw <stuff>"));
        assert_eq!(root.children[0].line, 5);
    }

    #[test]
    fn reports_errors_with_lines() {
        let err = parse("<a>\n<b></a>", "XML").unwrap_err().to_string();
        assert!(err.contains("line 2"), "{err}");
        assert!(parse("<a x=1/>", "XML").is_err());
        assert!(parse("<a>&bogus;</a>", "XML").is_err());
    }

    #[test]
    fn escape_round_trips() {
        let s = "a<b>&\"c\"\n\u{e9}\t";
        let doc = format!("<x v=\"{}\">{}</x>", escape(s), escape(s));
        let root = parse(&doc, "XML").unwrap();
        assert_eq!(root.attr("v"), Some(s));
        assert_eq!(root.text, s);
    }

    #[test]
    fn streams_events() {
        let doc = "<a k='v'>x<b/>y &amp; z<![CDATA[<c>]]></a>";
        let mut r = Reader::new(doc, "XML");
        let mut seen = Vec::new();
        while let Some(e) = r.next().unwrap() {
            let attrs: Vec<String> = r.attrs().iter().map(|(k, v)| format!("{k}={v}")).collect();
            seen.push(match e {
                Event::Start { name, empty, .. } => format!("<{name} {attrs:?} {empty}"),
                Event::End { name } => format!("</{name}"),
                Event::Text(t) => format!("text {t}"),
            });
        }
        assert_eq!(
            seen,
            vec![
                "<a [\"k=v\"] false",
                "text x",
                "<b [] true",
                "</b",
                "text y & z",
                "text <c>",
                "</a"
            ]
        );
    }

    #[test]
    fn deep_nesting_is_an_error_not_a_crash() {
        let deep = "<a>".repeat(100_000);
        let err = parse(&deep, "XML").unwrap_err().to_string();
        assert!(err.contains("nested more than 512"), "{err}");
        assert!(err.contains("line 1"), "{err}");
        let ok = "<a>".repeat(MAX_DEPTH) + &"</a>".repeat(MAX_DEPTH);
        assert!(parse(&ok, "XML").is_ok());
    }

    /// The tree-building parser this module replaced, kept to check that
    /// the streaming reader accepts and rejects exactly the same documents,
    /// with the same messages and lines.
    mod reference {
        use crate::error::{Error, Result};

        #[derive(Debug, Clone, Default, PartialEq)]
        pub struct Element {
            pub name: String,
            pub attrs: Vec<(String, String)>,
            pub children: Vec<Element>,
            pub text: String,
            pub line: usize,
        }

        pub fn parse(input: &str, format: &'static str) -> Result<Element> {
            let mut p = Parser {
                s: input.as_bytes(),
                src: input,
                pos: 0,
                line: 1,
                format,
            };
            p.parse_document()
        }

        struct Parser<'a> {
            s: &'a [u8],
            src: &'a str,
            pos: usize,
            line: usize,
            format: &'static str,
        }

        impl<'a> Parser<'a> {
            fn err<T>(&self, message: impl Into<String>) -> Result<T> {
                Err(Error::Parse {
                    format: self.format,
                    line: Some(self.line),
                    message: message.into(),
                })
            }

            fn starts_with(&self, pat: &str) -> bool {
                self.s[self.pos..].starts_with(pat.as_bytes())
            }

            fn advance(&mut self, n: usize) {
                let end = (self.pos + n).min(self.s.len());
                self.line += self.s[self.pos..end]
                    .iter()
                    .filter(|&&b| b == b'\n')
                    .count();
                self.pos = end;
            }

            fn skip_until(&mut self, pat: &str) -> Result<()> {
                match self.src[self.pos..].find(pat) {
                    Some(i) => {
                        self.advance(i + pat.len());
                        Ok(())
                    }
                    None => self.err(format!("unterminated construct, expected {pat:?}")),
                }
            }

            fn skip_ws(&mut self) {
                while self.pos < self.s.len() && self.s[self.pos].is_ascii_whitespace() {
                    self.advance(1);
                }
            }

            fn skip_misc(&mut self) -> Result<bool> {
                if self.starts_with("<?") {
                    self.skip_until("?>")?;
                } else if self.starts_with("<!--") {
                    self.skip_until("-->")?;
                } else if self.starts_with("<!DOCTYPE") || self.starts_with("<!doctype") {
                    let rest = &self.src[self.pos..];
                    let close = rest.find('>').unwrap_or(rest.len());
                    if let Some(open) = rest[..close].find('[') {
                        self.advance(open);
                        self.skip_until("]")?;
                        self.skip_until(">")?;
                    } else {
                        self.skip_until(">")?;
                    }
                } else {
                    return Ok(false);
                }
                Ok(true)
            }

            fn parse_document(&mut self) -> Result<Element> {
                if self.starts_with("\u{feff}") {
                    self.advance(3);
                }
                loop {
                    self.skip_ws();
                    if !self.skip_misc()? {
                        break;
                    }
                }
                if !self.starts_with("<") {
                    return self.err("expected the root element");
                }
                let root = self.parse_element()?;
                loop {
                    self.skip_ws();
                    if self.pos >= self.s.len() || !self.skip_misc()? {
                        break;
                    }
                }
                if self.pos < self.s.len() {
                    return self.err("unexpected content after the root element");
                }
                Ok(root)
            }

            fn parse_name(&mut self) -> Result<String> {
                let start = self.pos;
                while self.pos < self.s.len() {
                    let b = self.s[self.pos];
                    if b.is_ascii_whitespace() || b == b'>' || b == b'/' || b == b'=' {
                        break;
                    }
                    self.pos += 1;
                }
                if start == self.pos {
                    return self.err("expected a name");
                }
                let name = &self.src[start..self.pos];
                Ok(super::super::local(name).to_string())
            }

            fn parse_element(&mut self) -> Result<Element> {
                let line = self.line;
                self.advance(1);
                let name = self.parse_name()?;
                let mut el = Element {
                    name,
                    line,
                    ..Default::default()
                };
                loop {
                    self.skip_ws();
                    if self.starts_with("/>") {
                        self.advance(2);
                        return Ok(el);
                    }
                    if self.starts_with(">") {
                        self.advance(1);
                        break;
                    }
                    if self.pos >= self.s.len() {
                        return self.err(format!("unterminated <{}> tag", el.name));
                    }
                    let key = self.parse_name()?;
                    self.skip_ws();
                    if !self.starts_with("=") {
                        return self.err(format!("attribute {key:?} has no value"));
                    }
                    self.advance(1);
                    self.skip_ws();
                    let quote = match self.s.get(self.pos) {
                        Some(&q @ (b'"' | b'\'')) => q as char,
                        _ => return self.err(format!("attribute {key:?} value must be quoted")),
                    };
                    self.advance(1);
                    let rest = &self.src[self.pos..];
                    let end = match rest.find(quote) {
                        Some(e) => e,
                        None => {
                            return self.err(format!("unterminated value for attribute {key:?}"))
                        }
                    };
                    let raw = &rest[..end];
                    let value = super::super::decode_entities(raw)
                        .map(|v| v.into_owned())
                        .or_else(|m| self.err(m))?;
                    self.advance(end + 1);
                    el.attrs.push((key, value));
                }
                loop {
                    if self.pos >= self.s.len() {
                        return self.err(format!("<{}> is never closed", el.name));
                    }
                    if self.starts_with("</") {
                        self.advance(2);
                        let close = self.parse_name()?;
                        if close != el.name {
                            return self
                                .err(format!("expected </{}>, found </{}>", el.name, close));
                        }
                        self.skip_ws();
                        if !self.starts_with(">") {
                            return self.err("malformed closing tag");
                        }
                        self.advance(1);
                        return Ok(el);
                    }
                    if self.starts_with("<![CDATA[") {
                        self.advance(9);
                        let rest = &self.src[self.pos..];
                        let end = match rest.find("]]>") {
                            Some(e) => e,
                            None => return self.err("unterminated CDATA section"),
                        };
                        el.text.push_str(&rest[..end]);
                        self.advance(end + 3);
                        continue;
                    }
                    if self.skip_misc()? {
                        continue;
                    }
                    if self.starts_with("<") {
                        let child = self.parse_element()?;
                        el.children.push(child);
                        continue;
                    }
                    let rest = &self.src[self.pos..];
                    let end = rest.find('<').unwrap_or(rest.len());
                    let text = super::super::decode_entities(&rest[..end])
                        .map(|t| t.into_owned())
                        .or_else(|m| self.err(m))?;
                    el.text.push_str(&text);
                    self.advance(end);
                }
            }
        }
    }

    fn convert(el: Element) -> reference::Element {
        reference::Element {
            name: el.name,
            attrs: el.attrs,
            children: el.children.into_iter().map(convert).collect(),
            text: el.text,
            line: el.line,
        }
    }

    /// Random edits of small documents: the streaming parser and the
    /// reference parser agree on every result and every error.
    #[test]
    fn agrees_with_the_reference_parser() {
        const DOCS: [&str; 3] = [
            "<?xml version=\"1.0\"?>\n<!DOCTYPE g [ <!ENTITY e \"x\"> ]>\n<!-- c -->\n<g:r xmlns:g='u' a=\"1 &amp; 2\">\n  <n id='a'>t&#65;<![CDATA[<x>]]></n>\n  <?pi x?>\n  <e/>\n</g:r>\n<!-- end -->\n",
            "<a>\n<b c='d'\n e=\"f\">text</b>\n<!--x-->\n</a>",
            "\u{feff}<r><s>&lt;&gt;&quot;&apos;&#x263A;</s>\n</r>",
        ];
        const PIECES: [&str; 22] = [
            "<",
            ">",
            "/",
            "=",
            "'",
            "\"",
            "&",
            ";",
            "&#0;",
            "&#x110000;",
            "&bogus;",
            "<!--",
            "-->",
            "<?",
            "?>",
            "<![CDATA[",
            "]]>",
            "<!DOCTYPE",
            "[",
            "]",
            "\n",
            " ",
        ];
        for case in 0..20_000u64 {
            let mut rng = Stream::new(7, case);
            let mut s: Vec<char> = DOCS[rng.below(3) as usize].chars().collect();
            for _ in 0..1 + rng.below(3) {
                let at = rng.below(s.len() as u64 + 1) as usize;
                match rng.below(3) {
                    0 if at < s.len() => {
                        s.remove(at);
                    }
                    1 => {
                        let p = PIECES[rng.below(PIECES.len() as u64) as usize];
                        for (k, c) in p.chars().enumerate() {
                            s.insert(at + k, c);
                        }
                    }
                    _ => s.truncate(at),
                }
            }
            let doc: String = s.into_iter().collect();
            let new = parse(&doc, "XML").map(convert);
            let old = reference::parse(&doc, "XML");
            assert_eq!(new, old, "case {case}: {doc:?}");
        }
    }
}
