#!/usr/bin/env python3
"""Check and tidy the bibliography (docs/references.bib) and its index.

Every algorithm and measure cites its sources by BibTeX key, written in
backticks in Rust module comments, R help pages, manifests, and docs
(`brandes2001faster`). This script checks that

  - every key is unique and the entries are sorted by key, so a new entry
    has exactly one place to go and merges never conflict;
  - every entry has the fields its type needs (author or editor, title,
    year; journal or booktitle), and DOIs are well-formed (10.xxxx/...),
    not URLs, and not used by two entries;
  - every key cited in the sources exists in the bibliography;
  - docs/references.md, the index of entries and where each is cited, is
    up to date (it is generated from the bibliography by --update).

Entries that are not cited by key anywhere are listed as notes (many are
cited by author and year in R help pages); --strict makes them errors.

    python3 tools/check_bib.py            # check; exit 1 on any error
    python3 tools/check_bib.py --strict   # uncited entries are errors too
    python3 tools/check_bib.py --update   # sort the .bib and rewrite the index
    python3 tools/check_bib.py --keys     # list every key with its title
"""

from __future__ import annotations

import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BIB = ROOT / "docs" / "references.bib"
INDEX = ROOT / "docs" / "references.md"

# Where citations are looked for. r/src/rust is a copy of crates/.
SOURCES = [
    "crates/engine/src", "crates/engine/tests", "crates/capi/src",
    "r/R", "r/man", "r/vignettes", "r/inst", "docs", "validation", "bench",
    "data-raw", "README.md", "CONTRIBUTING.md",
]
SUFFIXES = {".rs", ".R", ".Rd", ".Rmd", ".md", ".yaml", ".yml", ".py", ".c", ".h", ".txt"}
# A key cited in prose or code: two or more letters, a four-digit year,
# then letters, digits, or underscores, not inside a longer word, path, or
# DOI.
CITATION = re.compile(r"(?<![\w/.-])([a-z]{2,}[a-z-]*(?:1[5-9]|20)\d\d[a-z0-9_]*)(?![\w/-])")
DOI = re.compile(r"^10\.\d{4,9}/\S+$")
ENTRY = re.compile(r"^@(\w+)\s*\{\s*([^,\s]+)\s*,", re.M)
FIELD = re.compile(r"^\s*(\w+)\s*=\s*(\{(?:[^{}]|\{(?:[^{}]|\{[^{}]*\})*\})*\}|\"[^\"]*\"|[^,\n]+)",
                   re.M)
HEADER = """\
% References for the algorithms, conventions, datasets, file formats, and
% software behind this project. Every implementation cites its sources by
% key in a comment at the top of its module, and every R help page lists
% its key references with \\doi{}. Entries are sorted by key; a comment
% directly above an entry (for example on a missing DOI) belongs to it.
% Check and sort with `python3 tools/check_bib.py --update`, which also
% regenerates docs/references.md. Verify every new entry against Crossref
% (https://api.crossref.org/works/<DOI>) or the publisher.
"""


@dataclass
class Entry:
    type: str
    key: str
    fields: dict
    text: str            # the entry as written, without surrounding blank lines
    comment: str = ""    # comment lines directly above it
    line: int = 0


def parse(text: str) -> list[Entry]:
    """Entries in file order, each with the comment lines directly above it.

    Section banners (lines of dashes with a heading between them) and the
    file's leading comment block are not attached to any entry."""
    entries = []
    prev_end = 0
    for m in ENTRY.finditer(text):
        if m.start() < prev_end:
            continue  # an "@" inside a field value, not an entry
        # The entry ends at the closing brace of its outermost group.
        depth = 0
        end = len(text)
        for i in range(m.start(), len(text)):
            if text[i] == "{":
                depth += 1
            elif text[i] == "}":
                depth -= 1
                if depth == 0:
                    end = i + 1
                    break
        chunk = text[m.start():end]
        fields = {}
        for f in FIELD.finditer(chunk[m.end() - m.start():]):
            value = f.group(2).strip()
            if value[:1] in "{\"":
                value = value[1:-1]
            fields[f.group(1).lower()] = value.strip()
        # Comment lines directly above the entry (separated from whatever
        # precedes them by a blank line, a banner, or an entry) belong to
        # it; the comment block at the very top of the file does not.
        gap = text[prev_end:m.start()].rstrip("\n").split("\n")
        block: list[str] = []
        for ln in reversed(gap):
            if ln.startswith("%") and not ln.startswith("% ---"):
                block.insert(0, ln)
            else:
                break
        if prev_end == 0 and len(block) == len(gap):
            block = []
        entries.append(Entry(m.group(1).lower(), m.group(2), fields, chunk,
                             comment="\n".join(block),
                             line=text.count("\n", 0, m.start()) + 1))
        prev_end = end
    return entries


def source_files():
    for src in SOURCES:
        p = ROOT / src
        if p.is_file():
            yield p
        elif p.is_dir():
            for f in sorted(p.rglob("*")):
                if f.is_file() and f.suffix in SUFFIXES:
                    yield f


def citations() -> dict[str, set[str]]:
    """Key -> files (relative to the repository root) that mention it."""
    cited: dict[str, set[str]] = {}
    for f in source_files():
        if f in (BIB, INDEX):
            continue
        try:
            text = f.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        for key in set(CITATION.findall(text)):
            cited.setdefault(key, set()).add(f.relative_to(ROOT).as_posix())
    return cited


def check_entry(e: Entry) -> list[str]:
    f, key, where = e.fields, e.key, f"{BIB.name}:{e.line}"
    errors = []
    if not re.match(r"^[A-Za-z][A-Za-z0-9_:-]*$", key):
        errors.append(f"{where}: key {key!r} has characters other than letters, digits, _ and -")
    if "title" not in f:
        errors.append(f"{where}: {key} has no title")
    web = e.type in ("misc", "manual", "online", "techreport", "software")
    if "author" not in f and "editor" not in f and not web:
        errors.append(f"{where}: {key} has neither author nor editor")
    if "year" not in f:
        if not web:
            errors.append(f"{where}: {key} has no year")
    elif not re.match(r"^\d{4}$", f["year"]):
        errors.append(f"{where}: {key} has a year that is not four digits: {f['year']!r}")
    if e.type == "article" and "journal" not in f:
        errors.append(f"{where}: article {key} has no journal")
    if e.type in ("inproceedings", "incollection") and "booktitle" not in f:
        errors.append(f"{where}: {e.type} {key} has no booktitle")
    doi = f.get("doi")
    if doi is not None:
        if doi.lower().startswith(("http", "doi:")):
            errors.append(f"{where}: {key} writes the DOI as a URL or with a prefix: {doi!r}")
        elif not DOI.match(doi):
            errors.append(f"{where}: {key} has a malformed DOI: {doi!r}")
    return errors


ACCENTS = {"'": "\u0301", "`": "\u0300", '"': "\u0308", "^": "\u0302", "~": "\u0303",
           "v": "\u030c", "c": "\u0327", "u": "\u0306", "H": "\u030b", "=": "\u0304",
           ".": "\u0307", "k": "\u0328", "r": "\u030a"}


def detex(s: str) -> str:
    """Plain text for a BibTeX field: accents composed, braces and math
    delimiters removed."""
    import unicodedata
    s = re.sub(r"\\url\{([^}]*)\}", r"\1", s)
    s = re.sub(r"\\([\'`\"^~vcuH=.kr])\{?([A-Za-z])\}?",
               lambda m: unicodedata.normalize("NFC", m.group(2) + ACCENTS[m.group(1)]), s)
    s = s.replace("\\&", "&").replace("--", "\u2013").replace("\\log", "log")
    s = re.sub(r"\\[a-zA-Z]+\s*", "", s)
    return re.sub(r"[{}$]", "", s).strip()


def authors(e: Entry) -> str:
    names = e.fields.get("author") or e.fields.get("editor") or ""
    parts = [p.strip() for p in re.split(r"\s+and\s+", names) if p.strip()]
    last = []
    for p in parts:
        if p.startswith("{"):
            last.append(detex(p))          # a corporate author, kept whole
        elif "," in p:
            last.append(detex(p.split(",")[0]))
        else:
            last.append(detex(p.split()[-1]))
    if not last:
        return ""
    if len(last) == 1:
        return last[0]
    if len(last) == 2:
        return f"{last[0]} and {last[1]}"
    if len(last) == 3:
        return f"{last[0]}, {last[1]}, and {last[2]}"
    return f"{last[0]} et al."


def venue(e: Entry) -> str:
    f = e.fields
    for name in ("journal", "booktitle", "publisher", "institution", "organization", "howpublished"):
        if name in f:
            return detex(f[name])
    return ""


def index_text(entries: list[Entry], cited: dict[str, set[str]]) -> str:
    lines = [
        "# References",
        "",
        "Every entry of [references.bib](references.bib), with the files that cite",
        "it by key. Generated by `python3 tools/check_bib.py --update`; edit the",
        "`.bib` file, not this page. An entry cited by key nowhere is usually",
        "cited by author and year in an R help page or a vignette instead.",
        "",
        "| Key | Reference | Cited in |",
        "|---|---|---|",
    ]
    for e in entries:
        title = detex(e.fields.get("title", ""))
        ref = f"{authors(e)} ({e.fields.get('year', 'n.d.')}). {title}."
        v = venue(e)
        if v:
            ref += f" *{v}*."
        ref = ref.replace("|", "\\|")
        files = sorted(cited.get(e.key, ()))
        shown = ", ".join(f"`{f}`" for f in files[:4])
        if len(files) > 4:
            shown += f", and {len(files) - 4} more"
        lines.append(f"| `{e.key}` | {ref} | {shown} |")
    return "\n".join(lines) + "\n"


def bib_text(entries: list[Entry]) -> str:
    blocks = []
    for e in sorted(entries, key=lambda e: e.key.lower()):
        blocks.append((e.comment + "\n" if e.comment else "") + e.text.strip("\n"))
    return HEADER + "\n" + "\n\n".join(blocks) + "\n"


def check(strict: bool) -> tuple[list[str], list[str]]:
    errors, notes = [], []
    text = BIB.read_text(encoding="utf-8")
    entries = parse(text)
    seen: dict[str, int] = {}
    dois: dict[str, str] = {}
    for e in entries:
        if e.key in seen:
            errors.append(f"{BIB.name}:{e.line}: key {e.key} is already used at line {seen[e.key]}")
        seen[e.key] = e.line
        errors += check_entry(e)
        doi = e.fields.get("doi", "").lower()
        if doi:
            if doi in dois:
                (errors if strict else notes).append(
                    f"{BIB.name}:{e.line}: {e.key} has the same DOI as {dois[doi]}; merge them")
            dois[doi] = e.key
    keys = [e.key.lower() for e in entries]
    if keys != sorted(keys):
        first = next(i for i in range(1, len(keys)) if keys[i] < keys[i - 1])
        errors.append(f"{BIB.name}:{entries[first].line}: entries are not sorted by key "
                      f"({entries[first].key} after {entries[first - 1].key}); "
                      "run tools/check_bib.py --update")
    cited = citations()
    for key, files in sorted(cited.items()):
        if key not in seen:
            errors.append(f"cited key {key} is not in {BIB.name} ("
                          + ", ".join(sorted(files)[:3]) + ")")
    for e in entries:
        if e.key not in cited:
            (errors if strict else notes).append(f"{e.key} is in {BIB.name} but cited by key nowhere")
    if not INDEX.exists() or INDEX.read_text(encoding="utf-8") != index_text(entries, cited):
        errors.append(f"{INDEX.relative_to(ROOT)} is out of date; run tools/check_bib.py --update")
    return errors, notes


def update() -> None:
    entries = parse(BIB.read_text(encoding="utf-8"))
    new = bib_text(entries)
    if new != BIB.read_text(encoding="utf-8"):
        BIB.write_text(new, encoding="utf-8", newline="\n")
        print(f"sorted {BIB.relative_to(ROOT)}")
    entries = parse(BIB.read_text(encoding="utf-8"))
    idx = index_text(entries, citations())
    if not INDEX.exists() or INDEX.read_text(encoding="utf-8") != idx:
        INDEX.write_text(idx, encoding="utf-8", newline="\n")
        print(f"wrote {INDEX.relative_to(ROOT)}")


def main(argv: list[str]) -> int:
    if "--keys" in argv:
        for e in parse(BIB.read_text(encoding="utf-8")):
            print(f"{e.key:36} {re.sub(r'[{}]', '', e.fields.get('title', ''))}")
        return 0
    if "--update" in argv:
        update()
    errors, notes = check("--strict" in argv)
    for n in notes:
        print(f"note: {n}")
    for e in errors:
        print(f"error: {e}")
    n = len(parse(BIB.read_text(encoding="utf-8")))
    print(f"{BIB.relative_to(ROOT)}: {n} entries, {len(errors)} error(s), {len(notes)} note(s)")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
