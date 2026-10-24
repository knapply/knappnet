/*
 * CPython extension over the engine's C ABI (crates/capi/include/netcapi.h).
 *
 * The module has two jobs:
 *
 * 1. Convert a NetworkX graph to an engine graph in one pass over its
 *    dict-of-dicts adjacency (`G._node`, `G._adj`), reading only the edge
 *    structure and, when asked, one numeric edge attribute.
 * 2. Call engine functions on that graph with the GIL released and return
 *    plain Python objects (lists of floats, tuples, lists of node labels)
 *    that the pure-Python layer arranges into NetworkX's return types.
 *
 * It is written against the limited API of CPython 3.10, so one build
 * serves every later CPython version. Engine errors raise EngineError, a
 * subclass of NotImplementedError: NetworkX's dispatcher then runs its own
 * implementation, which reports the problem in NetworkX's terms.
 */
#define Py_LIMITED_API 0x030A0000
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "netcapi.h"

#define ERR_LEN 1024

static PyObject *EngineError = NULL;

/* ---- Graph object ------------------------------------------------------ */

typedef struct {
    PyObject_HEAD
    GraphHandle *h;
    int64_t n;
    int64_t m;
    int directed;
    int weighted;
    int64_t *src; /* edge endpoints in input order (m each) */
    int64_t *tgt;
    double *w;    /* NULL when unweighted */
    int64_t loops; /* number of self-loops */
    double wmin, wmax; /* weight range (1 when unweighted) */
    int wfinite;  /* every weight finite */
    PyObject *nodes; /* list of node labels, index order */
    PyObject *index; /* dict: node label -> index */
} GraphObject;

static PyTypeObject *GraphType = NULL;

static void graph_dealloc(PyObject *self) {
    GraphObject *g = (GraphObject *)self;
    PyTypeObject *tp = Py_TYPE(self);
    if (g->h) net_graph_free(g->h);
    free(g->src);
    free(g->tgt);
    free(g->w);
    Py_XDECREF(g->nodes);
    Py_XDECREF(g->index);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    tp_free(self);
    Py_DECREF(tp);
}

static PyObject *graph_get_n(PyObject *self, void *closure) {
    (void)closure;
    return PyLong_FromLongLong(((GraphObject *)self)->n);
}

static PyObject *graph_get_m(PyObject *self, void *closure) {
    (void)closure;
    return PyLong_FromLongLong(((GraphObject *)self)->m);
}

static PyObject *graph_get_directed(PyObject *self, void *closure) {
    (void)closure;
    return PyBool_FromLong(((GraphObject *)self)->directed);
}

static PyObject *graph_get_weighted(PyObject *self, void *closure) {
    (void)closure;
    return PyBool_FromLong(((GraphObject *)self)->weighted);
}

static PyObject *graph_get_nodes(PyObject *self, void *closure) {
    (void)closure;
    GraphObject *g = (GraphObject *)self;
    Py_INCREF(g->nodes);
    return g->nodes;
}

static PyObject *graph_get_index(PyObject *self, void *closure) {
    (void)closure;
    GraphObject *g = (GraphObject *)self;
    Py_INCREF(g->index);
    return g->index;
}

static PyObject *graph_get_loops(PyObject *self, void *closure) {
    (void)closure;
    return PyLong_FromLongLong(((GraphObject *)self)->loops);
}

static PyObject *graph_get_wmin(PyObject *self, void *closure) {
    (void)closure;
    return PyFloat_FromDouble(((GraphObject *)self)->wmin);
}

static PyObject *graph_get_wmax(PyObject *self, void *closure) {
    (void)closure;
    return PyFloat_FromDouble(((GraphObject *)self)->wmax);
}

static PyObject *graph_get_wfinite(PyObject *self, void *closure) {
    (void)closure;
    return PyBool_FromLong(((GraphObject *)self)->wfinite);
}

static PyGetSetDef graph_getset[] = {
    {"index", graph_get_index, NULL, "node label -> index", NULL},
    {"loops", graph_get_loops, NULL, "number of self-loops", NULL},
    {"wmin", graph_get_wmin, NULL, "smallest weight (1 if unweighted)", NULL},
    {"wmax", graph_get_wmax, NULL, "largest weight (1 if unweighted)", NULL},
    {"wfinite", graph_get_wfinite, NULL, "whether every weight is finite", NULL},
    {"n", graph_get_n, NULL, "number of nodes", NULL},
    {"m", graph_get_m, NULL, "number of edges", NULL},
    {"directed", graph_get_directed, NULL, "whether arcs are directed", NULL},
    {"weighted", graph_get_weighted, NULL, "whether weights were read", NULL},
    {"nodes", graph_get_nodes, NULL, "node labels in index order", NULL},
    {NULL, NULL, NULL, NULL, NULL},
};

static PyType_Slot graph_slots[] = {
    {Py_tp_dealloc, (void *)graph_dealloc},
    {Py_tp_getset, (void *)graph_getset},
    {Py_tp_doc, (void *)"An engine graph converted from a NetworkX graph."},
    {0, NULL},
};

static PyType_Spec graph_spec = {
    "knappnet._core.Graph",
    sizeof(GraphObject),
    0,
    Py_TPFLAGS_DEFAULT,
    graph_slots,
};

static GraphObject *as_graph(PyObject *obj) {
    if (!PyObject_TypeCheck(obj, GraphType)) {
        PyErr_SetString(PyExc_TypeError, "expected a knappnet._core.Graph");
        return NULL;
    }
    return (GraphObject *)obj;
}

/* ---- Conversion --------------------------------------------------------- */

/* Growable edge buffers. */
typedef struct {
    int64_t *src, *tgt;
    double *w;
    int64_t len, cap;
    int weighted;
} EdgeBuf;

static int edgebuf_push(EdgeBuf *b, int64_t s, int64_t t, double w) {
    if (b->len == b->cap) {
        int64_t cap = b->cap ? b->cap * 2 : 1024;
        int64_t *ns = realloc(b->src, (size_t)cap * sizeof(int64_t));
        if (!ns) return -1;
        b->src = ns;
        int64_t *nt = realloc(b->tgt, (size_t)cap * sizeof(int64_t));
        if (!nt) return -1;
        b->tgt = nt;
        if (b->weighted) {
            double *nw = realloc(b->w, (size_t)cap * sizeof(double));
            if (!nw) return -1;
            b->w = nw;
        }
        b->cap = cap;
    }
    b->src[b->len] = s;
    b->tgt[b->len] = t;
    if (b->weighted) b->w[b->len] = w;
    b->len++;
    return 0;
}

/* Iterates over a mapping: exact dicts with PyDict_Next, anything else
 * (NetworkX's filtered views) through a list of items. */
typedef struct {
    PyObject *map;
    PyObject *items; /* NULL for exact dicts */
    Py_ssize_t pos;
} MapIter;

static int mapiter_init(MapIter *it, PyObject *map) {
    it->map = map;
    it->pos = 0;
    it->items = NULL;
    if (PyDict_CheckExact(map)) return 0;
    it->items = PyMapping_Items(map);
    return it->items ? 0 : -1;
}

/* Borrowed key and value; returns 1 while items remain. */
static int mapiter_next(MapIter *it, PyObject **key, PyObject **value) {
    if (!it->items) return PyDict_Next(it->map, &it->pos, key, value);
    if (it->pos >= PyList_Size(it->items)) return 0;
    PyObject *pair = PyList_GetItem(it->items, it->pos++);
    *key = PyTuple_GetItem(pair, 0);
    *value = PyTuple_GetItem(pair, 1);
    return 1;
}

static void mapiter_clear(MapIter *it) { Py_CLEAR(it->items); }

/* Reads one numeric edge attribute; missing attributes take dflt. */
static int read_weight(PyObject *data, PyObject *key, double dflt, double *out) {
    PyObject *v;
    if (PyDict_CheckExact(data)) {
        v = PyDict_GetItemWithError(data, key);
        if (!v) {
            if (PyErr_Occurred()) return -1;
            *out = dflt;
            return 0;
        }
        Py_INCREF(v);
    } else {
        v = PyObject_GetItem(data, key);
        if (!v) {
            if (!PyErr_ExceptionMatches(PyExc_KeyError)) return -1;
            PyErr_Clear();
            *out = dflt;
            return 0;
        }
    }
    double x = PyFloat_AsDouble(v);
    Py_DECREF(v);
    if (x == -1.0 && PyErr_Occurred()) {
        PyErr_Clear();
        PyErr_SetString(EngineError, "edge attribute is not a number");
        return -1;
    }
    *out = x;
    return 0;
}

/* from_networkx(node, adj, directed, multigraph, weight, default)
 *
 * node: G._node; adj: G._adj (G._succ for directed graphs). Node i is the
 * i-th key of node. An undirected edge is read once, from the endpoint
 * that comes first in node order, so edges are numbered exactly as
 * G.edges() lists them. weight None reads no attribute. */
static PyObject *from_networkx(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *node, *adj, *weight;
    int directed, multigraph;
    double dflt;
    if (!PyArg_ParseTuple(args, "OOppOd", &node, &adj, &directed, &multigraph,
                          &weight, &dflt))
        return NULL;
    int weighted = weight != Py_None;

    EdgeBuf b = {NULL, NULL, NULL, 0, 0, weighted};
    PyObject *nodes = PySequence_List(node);
    if (!nodes) return NULL;
    Py_ssize_t n = PyList_Size(nodes);
    PyObject *index = PyDict_New();
    if (!index) {
        Py_DECREF(nodes);
        return NULL;
    }
    for (Py_ssize_t i = 0; i < n; i++) {
        PyObject *k = PyLong_FromSsize_t(i);
        if (!k || PyDict_SetItem(index, PyList_GetItem(nodes, i), k) < 0) {
            Py_XDECREF(k);
            goto fail_index;
        }
        Py_DECREF(k);
    }

    MapIter outer;
    if (mapiter_init(&outer, adj) < 0) goto fail_index;
    PyObject *u, *nbrs;
    while (mapiter_next(&outer, &u, &nbrs)) {
        PyObject *ui_obj = PyDict_GetItemWithError(index, u);
        if (!ui_obj) {
            if (!PyErr_Occurred())
                PyErr_SetString(PyExc_KeyError, "adjacency names an unknown node");
            goto fail_outer;
        }
        int64_t ui = PyLong_AsLongLong(ui_obj);
        MapIter inner;
        if (mapiter_init(&inner, nbrs) < 0) goto fail_outer;
        PyObject *v, *data;
        while (mapiter_next(&inner, &v, &data)) {
            PyObject *vi_obj = PyDict_GetItemWithError(index, v);
            if (!vi_obj) {
                if (!PyErr_Occurred())
                    PyErr_SetString(PyExc_KeyError, "adjacency names an unknown node");
                mapiter_clear(&inner);
                goto fail_outer;
            }
            int64_t vi = PyLong_AsLongLong(vi_obj);
            if (!directed && vi < ui) continue;
            if (multigraph) {
                MapIter keys;
                if (mapiter_init(&keys, data) < 0) {
                    mapiter_clear(&inner);
                    goto fail_outer;
                }
                PyObject *key, *edata;
                while (mapiter_next(&keys, &key, &edata)) {
                    double w = 1.0;
                    if ((weighted && read_weight(edata, weight, dflt, &w) < 0) ||
                        edgebuf_push(&b, ui, vi, w) < 0) {
                        if (!PyErr_Occurred()) PyErr_NoMemory();
                        mapiter_clear(&keys);
                        mapiter_clear(&inner);
                        goto fail_outer;
                    }
                }
                mapiter_clear(&keys);
            } else {
                double w = 1.0;
                if ((weighted && read_weight(data, weight, dflt, &w) < 0) ||
                    edgebuf_push(&b, ui, vi, w) < 0) {
                    if (!PyErr_Occurred()) PyErr_NoMemory();
                    mapiter_clear(&inner);
                    goto fail_outer;
                }
            }
        }
        mapiter_clear(&inner);
    }
    mapiter_clear(&outer);
    if (weighted && !b.w) {
        /* No edges: keep the graph weighted, so weighted calls accept it. */
        b.w = malloc(sizeof(double));
        if (!b.w) {
            PyErr_NoMemory();
            goto fail_index;
        }
    }

    GraphObject *g = PyObject_New(GraphObject, GraphType);
    if (!g) {
        free(b.src);
        free(b.tgt);
        free(b.w);
        Py_DECREF(nodes);
        Py_DECREF(index);
        return NULL;
    }
    g->h = NULL;
    g->loops = 0;
    g->wmin = 1.0;
    g->wmax = 1.0;
    g->wfinite = 1;
    for (int64_t e = 0; e < b.len; e++)
        if (b.src[e] == b.tgt[e]) g->loops++;
    if (weighted && b.len > 0) {
        g->wmin = INFINITY;
        g->wmax = -INFINITY;
        for (int64_t e = 0; e < b.len; e++) {
            double x = b.w[e];
            if (!isfinite(x)) g->wfinite = 0;
            if (x < g->wmin) g->wmin = x;
            if (x > g->wmax) g->wmax = x;
        }
    }
    g->n = n;
    g->m = b.len;
    g->directed = directed;
    g->weighted = weighted;
    g->src = b.src;
    g->tgt = b.tgt;
    g->w = b.w;
    g->nodes = nodes;
    g->index = index;
    char err[ERR_LEN];
    int32_t status;
    Py_BEGIN_ALLOW_THREADS
    status = net_graph_new(n, b.src, b.tgt, b.w, b.len, directed, &g->h, err,
                           ERR_LEN);
    Py_END_ALLOW_THREADS
    if (status) {
        Py_DECREF(g);
        PyErr_SetString(EngineError, err);
        return NULL;
    }
    return (PyObject *)g;

fail_outer:
    mapiter_clear(&outer);
fail_index:
    free(b.src);
    free(b.tgt);
    free(b.w);
    Py_DECREF(index);
    Py_DECREF(nodes);
    return NULL;
}

/* ---- Helpers ------------------------------------------------------------ */

static PyObject *engine_error(const char *err) {
    PyErr_SetString(EngineError, err);
    return NULL;
}

static double *alloc_doubles(int64_t len) {
    double *p = malloc((size_t)(len > 0 ? len : 1) * sizeof(double));
    if (!p) PyErr_NoMemory();
    return p;
}

static int64_t *alloc_int64(int64_t len) {
    int64_t *p = malloc((size_t)(len > 0 ? len : 1) * sizeof(int64_t));
    if (!p) PyErr_NoMemory();
    return p;
}

static PyObject *list_from_doubles(const double *x, int64_t len) {
    PyObject *out = PyList_New(len);
    if (!out) return NULL;
    for (int64_t i = 0; i < len; i++) {
        PyObject *v = PyFloat_FromDouble(x[i]);
        if (!v) {
            Py_DECREF(out);
            return NULL;
        }
        PyList_SetItem(out, i, v);
    }
    return out;
}

static PyObject *list_from_int64(const int64_t *x, int64_t len) {
    PyObject *out = PyList_New(len);
    if (!out) return NULL;
    for (int64_t i = 0; i < len; i++) {
        PyObject *v = PyLong_FromLongLong(x[i]);
        if (!v) {
            Py_DECREF(out);
            return NULL;
        }
        PyList_SetItem(out, i, v);
    }
    return out;
}

/* Consumes x (freed here) and returns it as a list of floats. */
static PyObject *finish_doubles(double *x, int64_t len) {
    PyObject *out = list_from_doubles(x, len);
    free(x);
    return out;
}

/* Reads a sequence of numbers into a new buffer of exactly len values. */
static double *doubles_from_seq(PyObject *seq, int64_t len) {
    PyObject *fast = PySequence_List(seq);
    if (!fast) return NULL;
    if (PyList_Size(fast) != len) {
        Py_DECREF(fast);
        PyErr_SetString(PyExc_ValueError, "sequence has the wrong length");
        return NULL;
    }
    double *x = alloc_doubles(len);
    if (!x) {
        Py_DECREF(fast);
        return NULL;
    }
    for (int64_t i = 0; i < len; i++) {
        x[i] = PyFloat_AsDouble(PyList_GetItem(fast, i));
        if (x[i] == -1.0 && PyErr_Occurred()) {
            free(x);
            Py_DECREF(fast);
            return NULL;
        }
    }
    Py_DECREF(fast);
    return x;
}

/* Reads a sequence of integers; *len receives its length. */
static int64_t *int64_from_seq(PyObject *seq, int64_t *len) {
    PyObject *fast = PySequence_List(seq);
    if (!fast) return NULL;
    *len = PyList_Size(fast);
    int64_t *x = alloc_int64(*len);
    if (!x) {
        Py_DECREF(fast);
        return NULL;
    }
    for (int64_t i = 0; i < *len; i++) {
        x[i] = PyLong_AsLongLong(PyList_GetItem(fast, i));
        if (x[i] == -1 && PyErr_Occurred()) {
            free(x);
            Py_DECREF(fast);
            return NULL;
        }
    }
    Py_DECREF(fast);
    return x;
}

static int32_t *int32_from_seq(PyObject *seq, int64_t len) {
    int64_t got = 0;
    int64_t *x = int64_from_seq(seq, &got);
    if (!x) return NULL;
    if (got != len) {
        free(x);
        PyErr_SetString(PyExc_ValueError, "sequence has the wrong length");
        return NULL;
    }
    int32_t *y = malloc((size_t)(len > 0 ? len : 1) * sizeof(int32_t));
    if (!y) {
        free(x);
        PyErr_NoMemory();
        return NULL;
    }
    for (int64_t i = 0; i < len; i++) y[i] = (int32_t)x[i];
    free(x);
    return y;
}

/* Groups (offsets/members of node indices) as a list of lists of labels. */
static PyObject *groups_to_labels(GraphObject *g, const int64_t *offsets,
                                  const int64_t *members, int64_t count) {
    PyObject *out = PyList_New(count);
    if (!out) return NULL;
    for (int64_t i = 0; i < count; i++) {
        int64_t lo = offsets[i], hi = offsets[i + 1];
        PyObject *grp = PyList_New(hi - lo);
        if (!grp) {
            Py_DECREF(out);
            return NULL;
        }
        for (int64_t j = lo; j < hi; j++) {
            PyObject *lab = PyList_GetItem(g->nodes, (Py_ssize_t)members[j]);
            Py_INCREF(lab);
            PyList_SetItem(grp, j - lo, lab);
        }
        PyList_SetItem(out, i, grp);
    }
    return out;
}

#define BEGIN_CALL                                                             \
    char err[ERR_LEN];                                                         \
    int32_t status;                                                            \
    Py_BEGIN_ALLOW_THREADS
#define END_CALL Py_END_ALLOW_THREADS

/* ---- Graph accessors --------------------------------------------------- */

/* edge_keys(g): list of (u, v) label pairs in edge order. */
static PyObject *edge_keys(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    if (!PyArg_ParseTuple(args, "O", &go)) return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    PyObject *out = PyList_New(g->m);
    if (!out) return NULL;
    for (int64_t e = 0; e < g->m; e++) {
        PyObject *u = PyList_GetItem(g->nodes, (Py_ssize_t)g->src[e]);
        PyObject *v = PyList_GetItem(g->nodes, (Py_ssize_t)g->tgt[e]);
        PyObject *t = PyTuple_Pack(2, u, v);
        if (!t) {
            Py_DECREF(out);
            return NULL;
        }
        PyList_SetItem(out, e, t);
    }
    return out;
}

/* ---- Degree and shortest paths ----------------------------------------- */

static PyObject *py_degree(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int mode, loops, weighted;
    if (!PyArg_ParseTuple(args, "Oiii", &go, &mode, &loops, &weighted)) return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    double *out = alloc_doubles(g->n);
    if (!out) return NULL;
    BEGIN_CALL
    status = net_degree(g->h, mode, loops, weighted, out, err, ERR_LEN);
    END_CALL
    if (status) {
        free(out);
        return engine_error(err);
    }
    return finish_doubles(out, g->n);
}

static PyObject *py_betweenness(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int weighted, normalized, threads;
    if (!PyArg_ParseTuple(args, "Oiii", &go, &weighted, &normalized, &threads))
        return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    double *out = alloc_doubles(g->n);
    if (!out) return NULL;
    BEGIN_CALL
    status = net_betweenness(g->h, weighted, 0, normalized, threads, out, err,
                             ERR_LEN);
    END_CALL
    if (status) {
        free(out);
        return engine_error(err);
    }
    return finish_doubles(out, g->n);
}

static PyObject *py_edge_betweenness(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int weighted, normalized, threads;
    if (!PyArg_ParseTuple(args, "Oiii", &go, &weighted, &normalized, &threads))
        return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    double *out = alloc_doubles(g->m);
    if (!out) return NULL;
    BEGIN_CALL
    status = net_edge_betweenness(g->h, weighted, 0, normalized, threads, out,
                                  err, ERR_LEN);
    END_CALL
    if (status) {
        free(out);
        return engine_error(err);
    }
    return finish_doubles(out, g->m);
}

static PyObject *py_closeness(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int harmonic, mode, weighted, normalized, unreachable, threads;
    if (!PyArg_ParseTuple(args, "Oiiiiii", &go, &harmonic, &mode, &weighted,
                          &normalized, &unreachable, &threads))
        return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    double *out = alloc_doubles(g->n);
    if (!out) return NULL;
    BEGIN_CALL
    status = net_closeness(g->h, harmonic, mode, weighted, 0, normalized,
                           unreachable, threads, out, err, ERR_LEN);
    END_CALL
    if (status) {
        free(out);
        return engine_error(err);
    }
    return finish_doubles(out, g->n);
}

/* distances(g, sources, mode, weighted, threads, integral, cutoff)
 *
 * One dict per source, keyed by node label, holding the reachable nodes
 * (distance at most cutoff when cutoff >= 0). integral != 0 returns ints.
 * sources is a sequence of node indices. */
static PyObject *py_distances(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go, *src_seq;
    int mode, weighted, threads, integral;
    double cutoff;
    if (!PyArg_ParseTuple(args, "OOiiiid", &go, &src_seq, &mode, &weighted,
                          &threads, &integral, &cutoff))
        return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    int64_t k = 0;
    int64_t *sources = int64_from_seq(src_seq, &k);
    if (!sources) return NULL;
    double *out = alloc_doubles(k * g->n);
    if (!out) {
        free(sources);
        return NULL;
    }
    BEGIN_CALL
    status = net_distances(g->h, sources, k, mode, weighted, 0, threads, out,
                           err, ERR_LEN);
    END_CALL
    free(sources);
    if (status) {
        free(out);
        return engine_error(err);
    }
    PyObject *res = PyList_New(k);
    if (!res) goto fail;
    for (int64_t s = 0; s < k; s++) {
        PyObject *d = PyDict_New();
        if (!d) goto fail_res;
        PyList_SetItem(res, s, d);
        for (int64_t v = 0; v < g->n; v++) {
            /* column-major: row s, column v */
            double x = out[s + v * k];
            if (isinf(x) || isnan(x)) continue;
            if (cutoff >= 0 && x > cutoff) continue;
            PyObject *val = integral ? PyLong_FromLongLong((long long)x)
                                     : PyFloat_FromDouble(x);
            if (!val) goto fail_res;
            int rc = PyDict_SetItem(d, PyList_GetItem(g->nodes, (Py_ssize_t)v), val);
            Py_DECREF(val);
            if (rc < 0) goto fail_res;
        }
    }
    free(out);
    return res;
fail_res:
    Py_DECREF(res);
fail:
    free(out);
    return NULL;
}

/* distance_matrix(g, sources, mode, weighted, threads): flat list of
 * len(sources) * n doubles, row-major by source; inf when unreachable. */
static PyObject *py_distance_matrix(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go, *src_seq;
    int mode, weighted, threads;
    if (!PyArg_ParseTuple(args, "OOiii", &go, &src_seq, &mode, &weighted, &threads))
        return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    int64_t k = 0;
    int64_t *sources = int64_from_seq(src_seq, &k);
    if (!sources) return NULL;
    double *out = alloc_doubles(k * g->n);
    if (!out) {
        free(sources);
        return NULL;
    }
    BEGIN_CALL
    status = net_distances(g->h, sources, k, mode, weighted, 0, threads, out,
                           err, ERR_LEN);
    END_CALL
    free(sources);
    if (status) {
        free(out);
        return engine_error(err);
    }
    PyObject *res = PyList_New(k * g->n);
    if (!res) {
        free(out);
        return NULL;
    }
    for (int64_t s = 0; s < k; s++)
        for (int64_t v = 0; v < g->n; v++) {
            PyObject *val = PyFloat_FromDouble(out[s + v * k]);
            if (!val) {
                Py_DECREF(res);
                free(out);
                return NULL;
            }
            PyList_SetItem(res, s * g->n + v, val);
        }
    free(out);
    return res;
}

/* path_summary(g, mode, weighted, unconnected, threads): the 7 totals of
 * net_path_summary over every source. */
static PyObject *py_path_summary(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int mode, weighted, unconnected, threads;
    if (!PyArg_ParseTuple(args, "Oiiii", &go, &mode, &weighted, &unconnected,
                          &threads))
        return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    double out[7];
    BEGIN_CALL
    status = net_path_summary(g->h, NULL, -1, mode, weighted, 0, unconnected,
                              threads, out, err, ERR_LEN);
    END_CALL
    if (status) return engine_error(err);
    return list_from_doubles(out, 7);
}

static PyObject *py_eccentricity(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int mode, weighted, threads;
    if (!PyArg_ParseTuple(args, "Oiii", &go, &mode, &weighted, &threads)) return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    double *out = alloc_doubles(g->n);
    if (!out) return NULL;
    BEGIN_CALL
    status = net_eccentricity(g->h, mode, weighted, 0, threads, out, err, ERR_LEN);
    END_CALL
    if (status) {
        free(out);
        return engine_error(err);
    }
    return finish_doubles(out, g->n);
}

static PyObject *py_reachability(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int mode, threads;
    if (!PyArg_ParseTuple(args, "Oii", &go, &mode, &threads)) return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    double *out = alloc_doubles(g->n);
    if (!out) return NULL;
    BEGIN_CALL
    status = net_reachability(g->h, mode, threads, out, err, ERR_LEN);
    END_CALL
    if (status) {
        free(out);
        return engine_error(err);
    }
    return finish_doubles(out, g->n);
}

/* ---- Spectral and walk-based centralities ------------------------------- */

static PyObject *py_eigenvector(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int mode, weighted, scale, loops, threads;
    double tol;
    long long max_iter;
    if (!PyArg_ParseTuple(args, "OiiiidLi", &go, &mode, &weighted, &scale, &loops,
                          &tol, &max_iter, &threads))
        return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    double *out = alloc_doubles(g->n);
    if (!out) return NULL;
    double value = NAN;
    int32_t connected = 0;
    BEGIN_CALL
    status = net_eigenvector(g->h, mode, weighted, scale, loops, tol, max_iter,
                             threads, out, &value, &connected, err, ERR_LEN);
    END_CALL
    if (status) {
        free(out);
        return engine_error(err);
    }
    PyObject *vals = finish_doubles(out, g->n);
    if (!vals) return NULL;
    return Py_BuildValue("(Ndi)", vals, value, (int)connected);
}

static PyObject *py_pagerank(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go, *pers_seq;
    int directed, weighted, dangling, loops, threads;
    double damping, tol;
    long long max_iter;
    if (!PyArg_ParseTuple(args, "OdiiOiidLi", &go, &damping, &directed, &weighted,
                          &pers_seq, &dangling, &loops, &tol, &max_iter, &threads))
        return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    double *pers = NULL;
    if (pers_seq != Py_None) {
        pers = doubles_from_seq(pers_seq, g->n);
        if (!pers) return NULL;
    }
    double *out = alloc_doubles(g->n);
    if (!out) {
        free(pers);
        return NULL;
    }
    BEGIN_CALL
    status = net_pagerank(g->h, damping, directed, weighted, pers, dangling, loops,
                          tol, max_iter, threads, out, err, ERR_LEN);
    END_CALL
    free(pers);
    if (status) {
        free(out);
        return engine_error(err);
    }
    return finish_doubles(out, g->n);
}

static PyObject *py_katz(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go, *beta_seq;
    int mode, weighted, loops, normalized, threads;
    double alpha, tol;
    long long max_iter;
    if (!PyArg_ParseTuple(args, "OdOiiiidLi", &go, &alpha, &beta_seq, &mode,
                          &weighted, &loops, &normalized, &tol, &max_iter, &threads))
        return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    Py_ssize_t blen = PySequence_Size(beta_seq);
    if (blen < 0) return NULL;
    double *beta = doubles_from_seq(beta_seq, blen);
    if (!beta) return NULL;
    double *out = alloc_doubles(g->n);
    if (!out) {
        free(beta);
        return NULL;
    }
    BEGIN_CALL
    status = net_katz(g->h, alpha, beta, blen, mode, weighted, loops, normalized,
                      tol, max_iter, threads, out, err, ERR_LEN);
    END_CALL
    free(beta);
    if (status) {
        free(out);
        return engine_error(err);
    }
    return finish_doubles(out, g->n);
}

static PyObject *py_hits(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int weighted, loops, scale, threads;
    double tol;
    long long max_iter;
    if (!PyArg_ParseTuple(args, "OiiidLi", &go, &weighted, &loops, &scale, &tol,
                          &max_iter, &threads))
        return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    double *hubs = alloc_doubles(g->n);
    double *auth = alloc_doubles(g->n);
    if (!hubs || !auth) {
        free(hubs);
        free(auth);
        return NULL;
    }
    double value = NAN;
    int32_t unique = 0;
    BEGIN_CALL
    status = net_hits(g->h, weighted, loops, scale, tol, max_iter, threads, hubs,
                      auth, &value, &unique, err, ERR_LEN);
    END_CALL
    if (status) {
        free(hubs);
        free(auth);
        return engine_error(err);
    }
    PyObject *h = finish_doubles(hubs, g->n);
    PyObject *a = finish_doubles(auth, g->n);
    if (!h || !a) {
        Py_XDECREF(h);
        Py_XDECREF(a);
        return NULL;
    }
    return Py_BuildValue("(NNdi)", h, a, value, (int)unique);
}

static PyObject *py_information_centrality(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int weighted;
    if (!PyArg_ParseTuple(args, "Oi", &go, &weighted)) return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    double *out = alloc_doubles(g->n);
    if (!out) return NULL;
    BEGIN_CALL
    status = net_information_centrality(g->h, weighted, out, err, ERR_LEN);
    END_CALL
    if (status) {
        free(out);
        return engine_error(err);
    }
    return finish_doubles(out, g->n);
}

/* ---- Structure and description ------------------------------------------ */

static PyObject *py_density(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int loops;
    if (!PyArg_ParseTuple(args, "Oi", &go, &loops)) return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    double out = NAN;
    BEGIN_CALL
    status = net_density(g->h, loops, &out, err, ERR_LEN);
    END_CALL
    if (status) return engine_error(err);
    return PyFloat_FromDouble(out);
}

static PyObject *py_components(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int kind;
    if (!PyArg_ParseTuple(args, "Oi", &go, &kind)) return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    int32_t *memb = malloc((size_t)(g->n > 0 ? g->n : 1) * sizeof(int32_t));
    if (!memb) return PyErr_NoMemory();
    int64_t count = 0;
    BEGIN_CALL
    status = net_components(g->h, kind, memb, &count, err, ERR_LEN);
    END_CALL
    if (status) {
        free(memb);
        return engine_error(err);
    }
    /* Components as lists of labels, in order of their lowest node. */
    int64_t *sizes = calloc((size_t)(count + 1), sizeof(int64_t));
    if (!sizes) {
        free(memb);
        return PyErr_NoMemory();
    }
    for (int64_t v = 0; v < g->n; v++) sizes[memb[v] + 1]++;
    for (int64_t c = 0; c < count; c++) sizes[c + 1] += sizes[c];
    int64_t *members = alloc_int64(g->n);
    int64_t *fill = alloc_int64(count + 1);
    if (!members || !fill) {
        free(memb);
        free(sizes);
        free(members);
        free(fill);
        return NULL;
    }
    memcpy(fill, sizes, (size_t)(count + 1) * sizeof(int64_t));
    for (int64_t v = 0; v < g->n; v++) members[fill[memb[v]]++] = v;
    PyObject *res = groups_to_labels(g, sizes, members, count);
    free(memb);
    free(sizes);
    free(members);
    free(fill);
    return res;
}

static PyObject *py_component_count(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int kind;
    if (!PyArg_ParseTuple(args, "Oi", &go, &kind)) return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    int32_t *memb = malloc((size_t)(g->n > 0 ? g->n : 1) * sizeof(int32_t));
    if (!memb) return PyErr_NoMemory();
    int64_t count = 0;
    BEGIN_CALL
    status = net_components(g->h, kind, memb, &count, err, ERR_LEN);
    END_CALL
    free(memb);
    if (status) return engine_error(err);
    return PyLong_FromLongLong(count);
}

static PyObject *py_compactness(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int mode, weighted, threads;
    if (!PyArg_ParseTuple(args, "Oiii", &go, &mode, &weighted, &threads)) return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    double out = NAN;
    BEGIN_CALL
    status = net_compactness(g->h, mode, weighted, 0, threads, &out, err, ERR_LEN);
    END_CALL
    if (status) return engine_error(err);
    return PyFloat_FromDouble(out);
}

static PyObject *py_reciprocity(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int mode;
    if (!PyArg_ParseTuple(args, "Oi", &go, &mode)) return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    double out = NAN;
    BEGIN_CALL
    status = net_reciprocity(g->h, mode, &out, err, ERR_LEN);
    END_CALL
    if (status) return engine_error(err);
    return PyFloat_FromDouble(out);
}

static PyObject *py_transitivity(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int undefined, threads;
    if (!PyArg_ParseTuple(args, "Oii", &go, &undefined, &threads)) return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    double out = NAN;
    BEGIN_CALL
    status = net_transitivity(g->h, undefined, threads, &out, err, ERR_LEN);
    END_CALL
    if (status) return engine_error(err);
    return PyFloat_FromDouble(out);
}

static PyObject *py_clustering(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int average, weighting, undefined, threads;
    if (!PyArg_ParseTuple(args, "Oiiii", &go, &average, &weighting, &undefined,
                          &threads))
        return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    int64_t len = average ? 1 : g->n;
    double *out = alloc_doubles(len);
    if (!out) return NULL;
    BEGIN_CALL
    status = net_clustering(g->h, average, weighting, undefined, threads, out, err,
                            ERR_LEN);
    END_CALL
    if (status) {
        free(out);
        return engine_error(err);
    }
    if (average) {
        double x = out[0];
        free(out);
        return PyFloat_FromDouble(x);
    }
    return finish_doubles(out, len);
}

static PyObject *py_triad_census(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int threads;
    if (!PyArg_ParseTuple(args, "Oi", &go, &threads)) return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    double out[16];
    BEGIN_CALL
    status = net_triad_census(g->h, threads, out, err, ERR_LEN);
    END_CALL
    if (status) return engine_error(err);
    return list_from_doubles(out, 16);
}

/* cut_structure(g): (cut-point flags per node, bridge flags per edge). */
static PyObject *py_cut_structure(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    if (!PyArg_ParseTuple(args, "O", &go)) return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    int32_t *cut = malloc((size_t)(g->n > 0 ? g->n : 1) * sizeof(int32_t));
    int32_t *bridge = malloc((size_t)(g->m > 0 ? g->m : 1) * sizeof(int32_t));
    if (!cut || !bridge) {
        free(cut);
        free(bridge);
        return PyErr_NoMemory();
    }
    BEGIN_CALL
    status = net_cut_structure(g->h, cut, bridge, err, ERR_LEN);
    END_CALL
    if (status) {
        free(cut);
        free(bridge);
        return engine_error(err);
    }
    PyObject *c = PyList_New(g->n), *b = PyList_New(g->m);
    if (!c || !b) goto fail;
    for (int64_t v = 0; v < g->n; v++) PyList_SetItem(c, v, PyBool_FromLong(cut[v]));
    for (int64_t e = 0; e < g->m; e++)
        PyList_SetItem(b, e, PyBool_FromLong(bridge[e]));
    free(cut);
    free(bridge);
    return Py_BuildValue("(NN)", c, b);
fail:
    Py_XDECREF(c);
    Py_XDECREF(b);
    free(cut);
    free(bridge);
    return NULL;
}

/* ---- Mixing patterns ---------------------------------------------------- */

static PyObject *py_assortativity_nominal(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go, *types_seq;
    long long k;
    int directed, loops, normalized;
    if (!PyArg_ParseTuple(args, "OOLiii", &go, &types_seq, &k, &directed, &loops,
                          &normalized))
        return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    int64_t len = 0;
    int64_t *types = int64_from_seq(types_seq, &len);
    if (!types) return NULL;
    if (len != g->n) {
        free(types);
        PyErr_SetString(PyExc_ValueError, "one type per node is needed");
        return NULL;
    }
    double out = NAN;
    BEGIN_CALL
    status = net_assortativity_nominal(g->h, types, k, directed, loops, normalized,
                                       &out, err, ERR_LEN);
    END_CALL
    free(types);
    if (status) return engine_error(err);
    return PyFloat_FromDouble(out);
}

static PyObject *py_assortativity_values(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go, *x_seq, *y_seq;
    int directed, loops, normalized;
    if (!PyArg_ParseTuple(args, "OOOiii", &go, &x_seq, &y_seq, &directed, &loops,
                          &normalized))
        return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    double *x = doubles_from_seq(x_seq, g->n);
    if (!x) return NULL;
    double *y = NULL;
    if (y_seq != Py_None) {
        y = doubles_from_seq(y_seq, g->n);
        if (!y) {
            free(x);
            return NULL;
        }
    }
    double out = NAN;
    BEGIN_CALL
    status = net_assortativity_values(g->h, x, y, directed, loops, normalized, &out,
                                      err, ERR_LEN);
    END_CALL
    free(x);
    free(y);
    if (status) return engine_error(err);
    return PyFloat_FromDouble(out);
}

static PyObject *py_assortativity_degree(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int weighted, directed, loops, normalized;
    if (!PyArg_ParseTuple(args, "Oiiii", &go, &weighted, &directed, &loops,
                          &normalized))
        return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    double out = NAN;
    BEGIN_CALL
    status = net_assortativity_degree(g->h, weighted, directed, loops, normalized,
                                      &out, err, ERR_LEN);
    END_CALL
    if (status) return engine_error(err);
    return PyFloat_FromDouble(out);
}

/* ---- Cohesion and groups ------------------------------------------------ */

static PyObject *py_coreness(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int mode;
    if (!PyArg_ParseTuple(args, "Oi", &go, &mode)) return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    double *out = alloc_doubles(g->n);
    if (!out) return NULL;
    BEGIN_CALL
    status = net_coreness(g->h, mode, out, err, ERR_LEN);
    END_CALL
    if (status) {
        free(out);
        return engine_error(err);
    }
    return finish_doubles(out, g->n);
}

static PyObject *py_max_cliques(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    long long min_size;
    int threads;
    if (!PyArg_ParseTuple(args, "OLi", &go, &min_size, &threads)) return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    GroupsHandle *gh = NULL;
    BEGIN_CALL
    status = net_max_cliques(g->h, min_size, threads, &gh, err, ERR_LEN);
    END_CALL
    if (status) return engine_error(err);
    int64_t count = net_groups_count(gh), total = net_groups_total(gh);
    int64_t *offsets = alloc_int64(count + 1);
    int64_t *members = alloc_int64(total);
    PyObject *res = NULL;
    if (offsets && members) {
        if (net_groups_copy(gh, offsets, members, err, ERR_LEN))
            engine_error(err);
        else
            res = groups_to_labels(g, offsets, members, count);
    }
    net_groups_free(gh);
    free(offsets);
    free(members);
    return res;
}

static PyObject *py_clique_number(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int threads;
    if (!PyArg_ParseTuple(args, "Oi", &go, &threads)) return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    int64_t out = 0;
    BEGIN_CALL
    status = net_clique_number(g->h, threads, &out, err, ERR_LEN);
    END_CALL
    if (status) return engine_error(err);
    return PyLong_FromLongLong(out);
}

static PyObject *py_bicomponents(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    if (!PyArg_ParseTuple(args, "O", &go)) return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    NodeSetsHandle *sh = NULL;
    BEGIN_CALL
    status = net_bicomponents(g->h, &sh, err, ERR_LEN);
    END_CALL
    if (status) return engine_error(err);
    int64_t count = net_node_sets_count(sh), total = net_node_sets_total(sh);
    int64_t *offsets = alloc_int64(count + 1);
    int64_t *members = alloc_int64(total);
    PyObject *res = NULL;
    if (offsets && members) {
        if (net_node_sets_copy(sh, offsets, members, err, ERR_LEN))
            engine_error(err);
        else
            res = groups_to_labels(g, offsets, members, count);
    }
    net_node_sets_free(sh);
    free(offsets);
    free(members);
    return res;
}

/* ---- Communities -------------------------------------------------------- */

static PyObject *py_modularity(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go, *memb_seq;
    double resolution;
    int weighted;
    if (!PyArg_ParseTuple(args, "OOdi", &go, &memb_seq, &resolution, &weighted))
        return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    int64_t len = 0;
    int64_t *memb = int64_from_seq(memb_seq, &len);
    if (!memb) return NULL;
    double out = NAN;
    BEGIN_CALL
    status = net_modularity(g->h, memb, len, resolution, weighted, &out, err,
                            ERR_LEN);
    END_CALL
    free(memb);
    if (status) return engine_error(err);
    return PyFloat_FromDouble(out);
}

/* Shared tail for the community functions: membership as a list. */
static PyObject *finish_membership(int64_t *memb, int64_t n) {
    PyObject *out = list_from_int64(memb, n);
    free(memb);
    return out;
}

static PyObject *py_louvain(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    double resolution;
    int weighted;
    unsigned long long seed;
    if (!PyArg_ParseTuple(args, "OdiK", &go, &resolution, &weighted, &seed))
        return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    int64_t *memb = alloc_int64(g->n);
    if (!memb) return NULL;
    double q = NAN;
    BEGIN_CALL
    status = net_louvain(g->h, resolution, weighted, (uint64_t)seed, memb, &q, err,
                         ERR_LEN);
    END_CALL
    if (status) {
        free(memb);
        return engine_error(err);
    }
    return finish_membership(memb, g->n);
}

static PyObject *py_fast_greedy(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    double resolution;
    int weighted;
    if (!PyArg_ParseTuple(args, "Odi", &go, &resolution, &weighted)) return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    int64_t *memb = alloc_int64(g->n);
    if (!memb) return NULL;
    double q = NAN;
    BEGIN_CALL
    status = net_fast_greedy(g->h, resolution, weighted, memb, &q, err, ERR_LEN);
    END_CALL
    if (status) {
        free(memb);
        return engine_error(err);
    }
    return finish_membership(memb, g->n);
}

static PyObject *py_label_propagation(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int variant, mode, weighted;
    unsigned long long seed;
    if (!PyArg_ParseTuple(args, "OiiiK", &go, &variant, &mode, &weighted, &seed))
        return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    int64_t *memb = alloc_int64(g->n);
    if (!memb) return NULL;
    BEGIN_CALL
    status = net_label_propagation(g->h, variant, mode, weighted, NULL, NULL,
                                   (uint64_t)seed, memb, err, ERR_LEN);
    END_CALL
    if (status) {
        free(memb);
        return engine_error(err);
    }
    return finish_membership(memb, g->n);
}

/* ---- Group centrality and structural holes ----------------------------- */

static PyObject *py_group_centrality(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go, *grp_seq;
    int measure, mode, weighted, normalized, threads;
    if (!PyArg_ParseTuple(args, "OOiiiii", &go, &grp_seq, &measure, &mode,
                          &weighted, &normalized, &threads))
        return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    int64_t k = 0;
    int64_t *grp = int64_from_seq(grp_seq, &k);
    if (!grp) return NULL;
    double out = NAN;
    BEGIN_CALL
    status = net_group_centrality(g->h, grp, k, measure, mode, weighted, 0,
                                  normalized, threads, &out, err, ERR_LEN);
    END_CALL
    free(grp);
    if (status) return engine_error(err);
    return PyFloat_FromDouble(out);
}

static PyObject *py_structural_holes(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int measure, weighted, threads;
    if (!PyArg_ParseTuple(args, "Oiii", &go, &measure, &weighted, &threads))
        return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    double *out = alloc_doubles(g->n);
    if (!out) return NULL;
    BEGIN_CALL
    status = net_structural_holes(g->h, measure, weighted, threads, out, err,
                                  ERR_LEN);
    END_CALL
    if (status) {
        free(out);
        return engine_error(err);
    }
    return finish_doubles(out, g->n);
}

/* ---- Two-mode networks -------------------------------------------------- */

static PyObject *py_twomode_density(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go, *types_seq;
    if (!PyArg_ParseTuple(args, "OO", &go, &types_seq)) return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    int32_t *types = int32_from_seq(types_seq, g->n);
    if (!types) return NULL;
    double out = NAN;
    BEGIN_CALL
    status = net_twomode_density(g->h, types, &out, err, ERR_LEN);
    END_CALL
    free(types);
    if (status) return engine_error(err);
    return PyFloat_FromDouble(out);
}

static PyObject *py_twomode_centrality(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go, *types_seq;
    int measure, normalized, unreachable, threads;
    if (!PyArg_ParseTuple(args, "OOiiii", &go, &types_seq, &measure, &normalized,
                          &unreachable, &threads))
        return NULL;
    GraphObject *g = as_graph(go);
    if (!g) return NULL;
    int32_t *types = int32_from_seq(types_seq, g->n);
    if (!types) return NULL;
    double *out = alloc_doubles(g->n);
    if (!out) {
        free(types);
        return NULL;
    }
    BEGIN_CALL
    status = net_twomode_centrality(g->h, types, measure, normalized, unreachable,
                                    threads, out, err, ERR_LEN);
    END_CALL
    free(types);
    if (status) {
        free(out);
        return engine_error(err);
    }
    return finish_doubles(out, g->n);
}

/* ---- Module ------------------------------------------------------------- */

static PyMethodDef methods[] = {
    {"from_networkx", from_networkx, METH_VARARGS, "Convert a NetworkX graph."},
    {"edge_keys", edge_keys, METH_VARARGS, "Edge endpoints as label pairs."},
    {"degree", py_degree, METH_VARARGS, NULL},
    {"betweenness", py_betweenness, METH_VARARGS, NULL},
    {"edge_betweenness", py_edge_betweenness, METH_VARARGS, NULL},
    {"closeness", py_closeness, METH_VARARGS, NULL},
    {"distances", py_distances, METH_VARARGS, NULL},
    {"distance_matrix", py_distance_matrix, METH_VARARGS, NULL},
    {"path_summary", py_path_summary, METH_VARARGS, NULL},
    {"eccentricity", py_eccentricity, METH_VARARGS, NULL},
    {"reachability", py_reachability, METH_VARARGS, NULL},
    {"eigenvector", py_eigenvector, METH_VARARGS, NULL},
    {"pagerank", py_pagerank, METH_VARARGS, NULL},
    {"katz", py_katz, METH_VARARGS, NULL},
    {"hits", py_hits, METH_VARARGS, NULL},
    {"information_centrality", py_information_centrality, METH_VARARGS, NULL},
    {"density", py_density, METH_VARARGS, NULL},
    {"components", py_components, METH_VARARGS, NULL},
    {"component_count", py_component_count, METH_VARARGS, NULL},
    {"compactness", py_compactness, METH_VARARGS, NULL},
    {"reciprocity", py_reciprocity, METH_VARARGS, NULL},
    {"transitivity", py_transitivity, METH_VARARGS, NULL},
    {"clustering", py_clustering, METH_VARARGS, NULL},
    {"triad_census", py_triad_census, METH_VARARGS, NULL},
    {"cut_structure", py_cut_structure, METH_VARARGS, NULL},
    {"assortativity_nominal", py_assortativity_nominal, METH_VARARGS, NULL},
    {"assortativity_values", py_assortativity_values, METH_VARARGS, NULL},
    {"assortativity_degree", py_assortativity_degree, METH_VARARGS, NULL},
    {"coreness", py_coreness, METH_VARARGS, NULL},
    {"max_cliques", py_max_cliques, METH_VARARGS, NULL},
    {"clique_number", py_clique_number, METH_VARARGS, NULL},
    {"bicomponents", py_bicomponents, METH_VARARGS, NULL},
    {"modularity", py_modularity, METH_VARARGS, NULL},
    {"louvain", py_louvain, METH_VARARGS, NULL},
    {"fast_greedy", py_fast_greedy, METH_VARARGS, NULL},
    {"label_propagation", py_label_propagation, METH_VARARGS, NULL},
    {"group_centrality", py_group_centrality, METH_VARARGS, NULL},
    {"structural_holes", py_structural_holes, METH_VARARGS, NULL},
    {"twomode_density", py_twomode_density, METH_VARARGS, NULL},
    {"twomode_centrality", py_twomode_centrality, METH_VARARGS, NULL},
    {NULL, NULL, 0, NULL},
};

static int module_exec(PyObject *m) {
    GraphType = (PyTypeObject *)PyType_FromSpec(&graph_spec);
    if (!GraphType) return -1;
    if (PyModule_AddObject(m, "Graph", (PyObject *)GraphType) < 0) {
        Py_DECREF(GraphType);
        return -1;
    }
    Py_INCREF(GraphType);
    EngineError = PyErr_NewExceptionWithDoc(
        "knappnet._core.EngineError",
        "The engine declined the call; NetworkX runs its own implementation.",
        PyExc_NotImplementedError, NULL);
    if (!EngineError) return -1;
    Py_INCREF(EngineError);
    if (PyModule_AddObject(m, "EngineError", EngineError) < 0) {
        Py_DECREF(EngineError);
        return -1;
    }
    return 0;
}

static PyModuleDef_Slot module_slots[] = {
    {Py_mod_exec, (void *)module_exec},
    {0, NULL},
};

static struct PyModuleDef module_def = {
    PyModuleDef_HEAD_INIT, "_core",
    "Engine graphs and engine calls for the NetworkX backend.",
    0, methods, module_slots, NULL, NULL, NULL,
};

PyMODINIT_FUNC PyInit__core(void) { return PyModuleDef_Init(&module_def); }
