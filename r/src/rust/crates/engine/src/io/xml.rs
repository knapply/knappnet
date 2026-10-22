//! A minimal XML reader and writer, sufficient for GraphML and GEXF.
//!
//! Supports elements, attributes, text, CDATA, comments, processing
//! instructions, a DOCTYPE declaration, and the predefined and numeric
//! character entities. Namespace prefixes are dropped from element and
//! attribute names. External entities are never resolved, so documents
//! cannot cause file or network access.

use crate::error::{Error, Result};

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

/// Parses a document and returns its root element.
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

    /// Skips prolog-style constructs; returns true if one was skipped.
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
        Ok(local(name).to_string())
    }

    fn parse_element(&mut self) -> Result<Element> {
        let line = self.line;
        self.advance(1); // '<'
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
                None => return self.err(format!("unterminated value for attribute {key:?}")),
            };
            let raw = &rest[..end];
            let value = decode_entities(raw).or_else(|m| self.err(m))?;
            self.advance(end + 1);
            el.attrs.push((key, value));
        }

        // Content.
        loop {
            if self.pos >= self.s.len() {
                return self.err(format!("<{}> is never closed", el.name));
            }
            if self.starts_with("</") {
                self.advance(2);
                let close = self.parse_name()?;
                if close != el.name {
                    return self.err(format!("expected </{}>, found </{}>", el.name, close));
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
            let text = decode_entities(&rest[..end]).or_else(|m| self.err(m))?;
            el.text.push_str(&text);
            self.advance(end);
        }
    }
}

fn local(name: &str) -> &str {
    name.rsplit(':').next().unwrap_or(name)
}

fn decode_entities(s: &str) -> std::result::Result<String, String> {
    if !s.contains('&') {
        return Ok(s.to_string());
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
    Ok(out)
}

/// Escapes text for use in element content or a double-quoted attribute.
pub fn escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '&' => out.push_str("&amp;"),
            '"' => out.push_str("&quot;"),
            '\n' => out.push_str("&#10;"),
            '\r' => out.push_str("&#13;"),
            '\t' => out.push_str("&#9;"),
            _ => out.push(c),
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

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
        let s = "a<b>&\"c\"\n";
        let doc = format!("<x v=\"{}\">{}</x>", escape(s), escape(s));
        let root = parse(&doc, "XML").unwrap();
        assert_eq!(root.attr("v"), Some(s));
        assert_eq!(root.text, s);
    }
}
