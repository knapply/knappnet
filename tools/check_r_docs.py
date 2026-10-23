#!/usr/bin/env python3
"""Static checks of the R package that do not need R installed.

    python3 tools/check_r_docs.py

Checks, in order:

* brackets, parentheses, and braces balance in every file under r/R and
  r/tests/testthat (strings and comments are skipped);
* every function exported from r/NAMESPACE is defined in r/R and has an
  \\alias in some file under r/man;
* every S3 method registered in r/NAMESPACE is defined in r/R;
* every \\usage entry in r/man names a function defined in r/R whose
  formals (names, order, and defaults) match the usage exactly;
* every Rd topic appears in the reference index of r/_pkgdown.yml;
* the top-level sections of every Rd file come in the package's order.

Exits 1 and prints one line per problem when anything fails.
"""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
R_DIR = ROOT / "r" / "R"
MAN_DIR = ROOT / "r" / "man"
TEST_DIR = ROOT / "r" / "tests" / "testthat"

# Sections in the order every Rd file keeps them: the order R renders them
# in (tools:::prepare_Rd sorts sections this way, custom \section blocks
# after \value), so the file reads as the help page does.
RD_ORDER = ["name", "alias", "docType", "title", "description", "usage",
            "arguments", "format", "details", "value", "section", "references",
            "seealso", "examples", "keyword"]

problems = []


def problem(msg):
    problems.append(msg)


# ---------------------------------------------------------------- R source

def strip_r(src):
    """Blanks strings and comments in R source, keeping line breaks so that
    positions still map to lines."""
    out = []
    i = 0
    n = len(src)
    while i < n:
        c = src[i]
        if c == "#":
            while i < n and src[i] != "\n":
                i += 1
            continue
        if c in "\"'`":
            quote = c
            out.append(" ")
            i += 1
            while i < n and src[i] != quote:
                if src[i] == "\\":
                    i += 1
                out.append("\n" if i < n and src[i] == "\n" else " ")
                i += 1
            out.append(" ")
            i += 1
            continue
        out.append(c)
        i += 1
    return "".join(out)


def check_balance(path):
    src = strip_r(path.read_text(encoding="utf-8"))
    pairs = {")": "(", "]": "[", "}": "{"}
    stack = []
    for pos, c in enumerate(src):
        if c in "([{":
            stack.append((c, pos))
        elif c in ")]}":
            if not stack or stack[-1][0] != pairs[c]:
                line = src[:pos].count("\n") + 1
                problem(f"{path}:{line}: unexpected `{c}`")
                return
            stack.pop()
    if stack:
        c, pos = stack[-1]
        line = src[:pos].count("\n") + 1
        problem(f"{path}:{line}: `{c}` is never closed")


def split_top_level(text):
    """Splits `text` on commas that are not inside brackets or strings."""
    parts = []
    depth = 0
    cur = []
    i = 0
    while i < len(text):
        c = text[i]
        if c in "\"'`":
            quote = c
            cur.append(c)
            i += 1
            while i < len(text) and text[i] != quote:
                if text[i] == "\\":
                    cur.append(text[i])
                    i += 1
                cur.append(text[i])
                i += 1
            cur.append(quote)
            i += 1
            continue
        if c in "([{":
            depth += 1
        elif c in ")]}":
            depth -= 1
        if c == "," and depth == 0:
            parts.append("".join(cur))
            cur = []
        else:
            cur.append(c)
        i += 1
    if "".join(cur).strip():
        parts.append("".join(cur))
    return parts


def normalize(expr):
    """Whitespace-insensitive form of an R expression (outside strings)."""
    out = []
    i = 0
    while i < len(expr):
        c = expr[i]
        if c in "\"'":
            j = i + 1
            while j < len(expr) and expr[j] != c:
                if expr[j] == "\\":
                    j += 1
                j += 1
            out.append(expr[i:j + 1])
            i = j + 1
            continue
        if not c.isspace():
            out.append(c)
        i += 1
    return "".join(out)


def parse_formals(text):
    """`text` is the inside of `function(...)`; returns [(name, default)]."""
    out = []
    for part in split_top_level(text):
        part = part.strip()
        if not part:
            continue
        m = re.match(r"^(`[^`]*`|[A-Za-z.][\w.]*|\.\.\.)\s*(=\s*(.*))?$", part, re.S)
        if not m:
            raise ValueError(f"cannot parse formal {part!r}")
        name = m.group(1).strip("`")
        default = normalize(m.group(3)) if m.group(3) is not None else None
        out.append((name, default))
    return out


def matching_paren(src, start):
    """Index of the `)` matching the `(` at `start` (strings skipped)."""
    depth = 0
    i = start
    while i < len(src):
        c = src[i]
        if c in "\"'`":
            quote = c
            i += 1
            while i < len(src) and src[i] != quote:
                if src[i] == "\\":
                    i += 1
                i += 1
        elif c == "#":
            while i < len(src) and src[i] != "\n":
                i += 1
        elif c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
            if depth == 0:
                return i
        i += 1
    raise ValueError("unbalanced parentheses")


def r_functions():
    """Every top-level function in r/R: name -> (formals, file)."""
    funs = {}
    pat = re.compile(r"^`?([A-Za-z.][\w.]*)`?\s*(?:<-|=)\s*function\s*\(", re.M)
    for path in sorted(R_DIR.glob("*.R")):
        src = path.read_text(encoding="utf-8")
        for m in pat.finditer(src):
            open_at = m.end() - 1
            close_at = matching_paren(src, open_at)
            try:
                formals = parse_formals(src[open_at + 1:close_at])
            except ValueError as e:
                line = src[:m.start()].count("\n") + 1
                problem(f"{path}:{line}: {e}")
                continue
            funs[m.group(1)] = (formals, path)
    return funs


# ------------------------------------------------------------------- NAMESPACE

def namespace():
    text = (ROOT / "r" / "NAMESPACE").read_text(encoding="utf-8")
    exports = re.findall(r"^export\(([^)]*)\)", text, re.M)
    exports = [e.strip().strip("`\"") for line in exports for e in line.split(",")]
    methods = []
    for m in re.finditer(r"^S3method\(([^,]+),\s*([^)]+)\)", text, re.M):
        generic = m.group(1).strip().strip("`\"").split("::")[-1]
        cls = m.group(2).strip().strip("`\"")
        methods.append((generic, cls))
    return exports, methods


# --------------------------------------------------------------------- Rd

def rd_sections(text):
    """[(section, body, line)] for the top-level sections of an Rd file."""
    out = []
    i = 0
    while i < len(text):
        m = re.compile(r"\\([A-Za-z]+)\{").search(text, i)
        if not m:
            break
        depth = 0
        j = m.end() - 1
        while j < len(text):
            c = text[j]
            if c == "\\":
                j += 2
                continue
            if c == "{":
                depth += 1
            elif c == "}":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        line = text[:m.start()].count("\n") + 1
        body = text[m.end():j]
        if m.group(1) == "section":
            # \section{title}{body}: skip the body, which is a second group.
            k = j + 1
            while k < len(text) and text[k] != "{":
                k += 1
            depth = 0
            while k < len(text):
                c = text[k]
                if c == "\\":
                    k += 2
                    continue
                if c == "{":
                    depth += 1
                elif c == "}":
                    depth -= 1
                    if depth == 0:
                        break
                k += 1
            j = k
        out.append((m.group(1), body, line))
        i = j + 1
    return out


def rd_usage_calls(body):
    """[(function name, formals text)] for each call in a \\usage body."""
    body = re.sub(r"%.*", "", body)  # Rd comments
    calls = []
    i = 0
    while i < len(body):
        m = re.compile(r"(\\method\{([^}]*)\}\{([^}]*)\}|`?([A-Za-z.][\w.]*)`?)\s*\(").search(body, i)
        if not m:
            break
        if m.group(2):
            name = f"{m.group(2)}.{m.group(3)}"
        else:
            name = m.group(4)
        open_at = m.end() - 1
        close_at = matching_paren(body, open_at)
        calls.append((name, body[open_at + 1:close_at]))
        i = close_at + 1
    return calls


def check_rd(funs, exports):
    aliases = {}
    topics = []
    for path in sorted(MAN_DIR.glob("*.Rd")):
        text = path.read_text(encoding="utf-8")
        sections = rd_sections(text)
        names = [s for s, _, _ in sections]
        if "name" not in names:
            problem(f"{path}: no \\name")
            continue
        topics.append(path.stem)
        # Section order.
        ranks = []
        for s, _, line in sections:
            if s not in RD_ORDER:
                problem(f"{path}:{line}: unknown top-level section \\{s}")
                continue
            ranks.append((RD_ORDER.index(s), s, line))
        for (r1, s1, _), (r2, s2, line) in zip(ranks, ranks[1:]):
            if r2 < r1:
                problem(f"{path}:{line}: \\{s2} comes after \\{s1}; the order is "
                        + ", ".join("\\" + s for s in RD_ORDER))
                break
        for s, body, line in sections:
            if s == "alias":
                aliases.setdefault(body.strip(), path)
            elif s == "usage":
                try:
                    calls = rd_usage_calls(body)
                except ValueError as e:
                    problem(f"{path}:{line}: {e}")
                    continue
                for name, args in calls:
                    if name not in funs:
                        problem(f"{path}:{line}: \\usage names `{name}`, which is not "
                                "defined in r/R")
                        continue
                    try:
                        usage = parse_formals(args)
                    except ValueError as e:
                        problem(f"{path}:{line}: {name}: {e}")
                        continue
                    actual, _ = funs[name]
                    if usage != actual:
                        problem(f"{path}:{line}: \\usage of {name}() differs from its "
                                f"definition:\n    usage:    {fmt(usage)}\n"
                                f"    function: {fmt(actual)}")
            elif s == "title":
                title = " ".join(body.split())
                if title.endswith("."):
                    problem(f"{path}:{line}: \\title ends with a period")
            elif s == "examples":
                for m in re.finditer(r"threads\s*=\s*(\d+)", body):
                    if int(m.group(1)) > 2:
                        problem(f"{path}:{line}: example uses threads = {m.group(1)}; "
                                "at most 2")
    for name in exports:
        if name not in funs:
            problem(f"NAMESPACE exports `{name}`, which is not defined in r/R")
        if name not in aliases:
            problem(f"NAMESPACE exports `{name}`, which has no \\alias in r/man")
    return topics


def fmt(formals):
    return ", ".join(n if d is None else f"{n} = {d}" for n, d in formals)


def check_pkgdown(topics):
    path = ROOT / "r" / "_pkgdown.yml"
    if not path.exists():
        return
    text = path.read_text(encoding="utf-8")
    ref = text.split("\nreference:", 1)
    if len(ref) < 2:
        problem(f"{path}: no reference section")
        return
    listed = set(re.findall(r"^\s*-\s*([A-Za-z.][\w.]*)\s*$", ref[1], re.M))
    for topic in topics:
        if topic not in listed:
            problem(f"{path}: help topic `{topic}` is not listed under reference:")
    for item in sorted(listed):
        if item not in topics:
            problem(f"{path}: reference: lists `{item}`, which has no Rd file")


def main():
    for path in sorted(R_DIR.glob("*.R")) + sorted(TEST_DIR.glob("*.R")):
        check_balance(path)
    funs = r_functions()
    exports, methods = namespace()
    for generic, cls in methods:
        if f"{generic}.{cls}" not in funs:
            problem(f"NAMESPACE registers S3method({generic}, {cls}), but "
                    f"`{generic}.{cls}` is not defined in r/R")
    topics = check_rd(funs, exports)
    check_pkgdown(topics)
    for p in problems:
        print(p)
    print(f"{len(funs)} functions, {len(exports)} exports, {len(methods)} S3 methods, "
          f"{len(topics)} help topics; {len(problems)} problem(s)")
    sys.exit(1 if problems else 0)


if __name__ == "__main__":
    main()
