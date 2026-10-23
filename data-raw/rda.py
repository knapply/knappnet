"""A small reader for R data files (.rda/.RData, XDR format, versions 2 and
3) and for the igraph objects stored in them, so that data-raw/ scripts can
read copies of data sets distributed as R packages without R installed.

Only what those scripts need is implemented: vectors, lists, pairlists,
attributes, environments, symbols, references, ALTREP compact sequences and
wrappers, and the two layouts of igraph graph objects (the list layout of
igraph 0.x-1.x and the `igraph_graph` layout of igraph 2.x, in which the
data are stored under a "graph_version"-tagged list).

The serialization format is described in "R Internals", section 1.8
(https://cran.r-project.org/doc/manuals/r-release/R-ints.html#Serialization-Formats).

    import rda
    objects = rda.load("ht_advice.rda")     # {name: value}
    nodes, edges, directed = rda.igraph_tables(objects["ht_advice"])
"""

import bz2
import gzip
import lzma
import math
import struct

NA_INT = -2147483648


class RNull:
    def __repr__(self):
        return "NULL"


NULL = RNull()


class RObject:
    """A parsed R object: `value` plus `attrs` (dict) and `tag`."""

    __slots__ = ("type", "value", "attrs")

    def __init__(self, type_, value, attrs=None):
        self.type = type_
        self.value = value
        self.attrs = attrs or {}

    def attr(self, name, default=None):
        return self.attrs.get(name, default)

    def names(self):
        n = self.attrs.get("names")
        return list(n.value) if n is not None else None

    def __getitem__(self, key):
        if isinstance(key, str):
            return self.value[self.names().index(key)]
        return self.value[key]

    def get(self, key, default=None):
        names = self.names() or []
        return self.value[names.index(key)] if key in names else default

    def __len__(self):
        return len(self.value)

    def __repr__(self):
        if not isinstance(self.value, (list, dict, str, bytes)):
            return "R<%s>" % self.type
        return "R<%s len=%d names=%s>" % (self.type, len(self.value), self.names())


class Env:
    __slots__ = ("frame", "attrs")

    def __init__(self):
        self.frame = {}
        self.attrs = {}


def _open(path):
    with open(path, "rb") as fh:
        head = fh.read(6)
    if head[:2] == b"\x1f\x8b":
        return gzip.open(path, "rb").read()
    if head[:3] == b"BZh":
        return bz2.open(path, "rb").read()
    if head[:6] == b"\xfd7zXZ\x00":
        return lzma.open(path, "rb").read()
    return open(path, "rb").read()


class _Reader:
    def __init__(self, data):
        self.data = data
        self.pos = 0
        self.refs = []

    def int(self):
        v = struct.unpack_from(">i", self.data, self.pos)[0]
        self.pos += 4
        return v

    def double(self):
        v = struct.unpack_from(">d", self.data, self.pos)[0]
        self.pos += 8
        return v

    def length(self):
        n = self.int()
        if n == -1:
            hi, lo = self.int(), self.int()
            n = (hi << 32) + lo
        return n

    def string_vec(self):
        if self.int() != 0:
            raise ValueError("names in persistent strings are not supported")
        return [self.item().value for _ in range(self.int())]

    def item(self):
        flags = self.int()
        t = flags & 0xFF
        has_attr = bool(flags & (1 << 9))
        has_tag = bool(flags & (1 << 10))
        if t == 254:  # NILVALUE_SXP
            return RObject("NULL", NULL)
        if t == 255:  # REFSXP
            idx = flags >> 8
            if idx == 0:
                idx = self.int()
            return self.refs[idx - 1]
        if t in (253, 252, 251, 250, 242, 241):  # global/empty/base env, unbound, missing
            return RObject({253: "globalenv", 242: "emptyenv", 241: "baseenv"}.get(t, "special"), None)
        if t in (249, 248, 247):  # NAMESPACESXP, PACKAGESXP, PERSISTSXP
            obj = RObject({249: "namespace", 248: "package", 247: "persistent"}[t],
                          self.string_vec())
            self.refs.append(obj)
            return obj
        if t == 1:  # SYMSXP
            obj = RObject("symbol", self.item().value)
            self.refs.append(obj)
            return obj
        if t == 4:  # ENVSXP
            self.int()  # locked
            env = RObject("environment", Env())
            self.refs.append(env)
            enclos = self.item()  # noqa: F841
            frame = self.item()
            hashtab = self.item()
            attrib = self.item()
            for src in (frame, hashtab):
                if src.type == "pairlist":
                    env.value.frame.update(src.value)
                elif src.type == "list":
                    for bucket in src.value:
                        if bucket.type == "pairlist":
                            env.value.frame.update(bucket.value)
            if attrib.type == "pairlist":
                env.attrs = attrib.value
            return env
        if t in (2, 6, 3, 5, 17, 239, 240):  # pairlist-like
            items = {}
            order = []
            attrs = None
            while True:
                if has_attr:
                    attrs = self.item()
                tag = self.item().value if has_tag else None
                car = self.item()
                if tag is None:
                    tag = len(order)
                items[tag] = car
                order.append(tag)
                flags = self.int()
                t2 = flags & 0xFF
                if t2 == 254:
                    break
                if t2 not in (2, 6, 3, 5, 17, 239, 240):
                    # A tail that is not a pairlist (dotted pair); rewind and read it.
                    self.pos -= 4
                    items["..."] = self.item()
                    break
                has_attr = bool(flags & (1 << 9))
                has_tag = bool(flags & (1 << 10))
            obj = RObject("pairlist", items)
            if attrs is not None and attrs.type == "pairlist":
                obj.attrs = attrs.value
            return obj
        if t == 238:  # ALTREP_SXP
            info = self.item()
            state = self.item()
            attr = self.item()
            cls = info.value[0].value
            obj = self._altrep(cls, state)
            if attr.type == "pairlist":
                obj.attrs = attr.value
            return obj
        if t == 9:  # CHARSXP
            n = self.int()
            if n == -1:
                return RObject("char", None)
            s = self.data[self.pos:self.pos + n]
            self.pos += n
            enc = "latin-1" if flags & (1 << 14) else "utf-8"
            return RObject("char", s.decode(enc, errors="replace"))
        if t == 10:  # LGLSXP
            n = self.length()
            vals = [None if v == NA_INT else bool(v) for v in
                    struct.unpack_from(">%di" % n, self.data, self.pos)]
            self.pos += 4 * n
            obj = RObject("logical", vals)
        elif t == 13:  # INTSXP
            n = self.length()
            vals = [None if v == NA_INT else v for v in
                    struct.unpack_from(">%di" % n, self.data, self.pos)]
            self.pos += 4 * n
            obj = RObject("integer", vals)
        elif t == 14:  # REALSXP
            n = self.length()
            vals = list(struct.unpack_from(">%dd" % n, self.data, self.pos))
            self.pos += 8 * n
            obj = RObject("double", vals)
        elif t == 15:  # CPLXSXP
            n = self.length()
            raw = struct.unpack_from(">%dd" % (2 * n), self.data, self.pos)
            self.pos += 16 * n
            obj = RObject("complex", [complex(raw[2 * i], raw[2 * i + 1]) for i in range(n)])
        elif t == 16:  # STRSXP
            n = self.length()
            obj = RObject("character", [self.item().value for _ in range(n)])
        elif t in (19, 20):  # VECSXP, EXPRSXP
            n = self.length()
            obj = RObject("list", [self.item() for _ in range(n)])
        elif t == 24:  # RAWSXP
            n = self.length()
            obj = RObject("raw", self.data[self.pos:self.pos + n])
            self.pos += n
        elif t == 25:  # S4SXP
            obj = RObject("S4", None)
        elif t == 22:  # EXTPTRSXP
            obj = RObject("externalptr", None)
            self.refs.append(obj)
            self.item()  # prot
            self.item()  # tag
        elif t == 23:  # WEAKREFSXP
            obj = RObject("weakref", None)
            self.refs.append(obj)
        elif t in (7, 8):  # SPECIALSXP, BUILTINSXP
            n = self.int()
            self.pos += n
            obj = RObject("builtin", None)
        else:
            raise ValueError("unsupported SEXP type %d at offset %d" % (t, self.pos))
        if has_attr:
            attrs = self.item()
            if attrs.type == "pairlist":
                obj.attrs = attrs.value
        return obj

    def _altrep(self, cls, state):
        if cls == "compact_intseq":
            n, start, step = state.value
            return RObject("integer", [int(start + i * step) for i in range(int(n))])
        if cls == "compact_realseq":
            n, start, step = state.value
            return RObject("double", [start + i * step for i in range(int(n))])
        if cls in ("wrap_integer", "wrap_real", "wrap_logical", "wrap_string",
                   "wrap_complex", "wrap_raw", "wrap_list"):
            return state.value[0]
        if cls == "deferred_string":
            arg = state.value[0] if state.type == "pairlist" else state
            if arg.type == "pairlist":
                arg = list(arg.value.values())[0]
            vals = ["NA" if v is None else
                    (str(int(v)) if isinstance(v, float) and v == int(v) else str(v))
                    for v in arg.value]
            return RObject("character", vals)
        raise ValueError("unsupported ALTREP class %s" % cls)


def load(path):
    """Reads an .rda file and returns {object name: RObject}."""
    data = _open(path)
    if data[:5] not in (b"RDX2\n", b"RDX3\n"):
        raise ValueError("%s is not an XDR .rda file" % path)
    r = _Reader(data)
    r.pos = 5
    fmt = data[5:7]
    if fmt != b"X\n":
        raise ValueError("only XDR format is supported (got %r)" % fmt)
    r.pos = 7
    version = r.int()
    r.int()  # writer version
    r.int()  # minimal reader version
    if version == 3:
        n = r.int()
        r.pos += n  # native encoding
    top = r.item()
    if top.type != "pairlist":
        raise ValueError("unexpected top-level object %s" % top.type)
    return top.value


def _r_list_to_dict(obj):
    if obj is None or obj.type == "NULL":
        return {}
    if obj.type == "pairlist":
        return dict(obj.value)
    names = obj.names() or []
    return {n: v for n, v in zip(names, obj.value)}


def _scalar(v):
    if v is None:
        return ""
    if isinstance(v, float):
        if math.isnan(v):
            return ""
        if v == int(v) and abs(v) < 1e15:
            return int(v)
    if isinstance(v, bool):
        return "TRUE" if v else "FALSE"
    return v


def igraph_tables(g):
    """Returns (nodes, edges, directed) for an igraph object: `nodes` is a
    list of dicts (one per vertex, in vertex order, 1-based "id" plus the
    vertex attributes); `edges` a list of dicts with 1-based "from", "to"
    plus the edge attributes, in edge order."""
    if g.type == "list":
        names = g.names()
        if names and "graph_version" in names:
            layout = "v2"
        else:
            layout = "v1"
    elif g.type == "environment":
        raise ValueError("unsupported igraph object layout (environment)")
    else:
        raise ValueError("not an igraph object: %s" % g.type)
    if layout == "v1":
        n = int(g.value[0].value[0])
        directed = bool(g.value[1].value[0])
        frm = [int(v) for v in g.value[2].value]
        to = [int(v) for v in g.value[3].value]
        attr = g.value[8]
    else:
        n = int(g["n"].value[0])
        directed = bool(g["directed"].value[0])
        frm = [int(v) for v in g["from"].value]
        to = [int(v) for v in g["to"].value]
        attr = g["attr"]
    attr = _attr_lists(attr)
    vattr = _r_list_to_dict(attr[1])
    eattr = _r_list_to_dict(attr[2])
    nodes = []
    for i in range(n):
        row = {"id": i + 1}
        for k, col in vattr.items():
            row[k] = _scalar(col.value[i])
        nodes.append(row)
    edges = []
    for j in range(len(frm)):
        row = {"from": frm[j] + 1, "to": to[j] + 1}
        for k, col in eattr.items():
            row[k] = _scalar(col.value[j])
        edges.append(row)
    return nodes, edges, directed


def _attr_lists(attr):
    """The attribute list of an igraph object as [graph, vertex, edge]
    attribute lists. igraph 0.x stores exactly those three; igraph 1.x
    prepends the attribute-format version number."""
    items = attr.value if attr.type == "list" else []
    if len(items) == 4 and items[0].type == "double":
        items = items[1:]
    items = list(items) + [RObject("list", [])] * (3 - len(items))
    return items


def igraph_graph_attrs(g):
    attr = g.value[8] if g.names() is None or "graph_version" not in g.names() else g["attr"]
    out = {}
    for k, v in _r_list_to_dict(_attr_lists(attr)[0]).items():
        out[k] = v.value[0] if isinstance(v.value, list) and len(v.value) == 1 else v.value
    return out
