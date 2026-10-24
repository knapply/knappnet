/* ---- Shortest paths (path graphs) --------------------------------------- */

typedef struct {
    PyObject_HEAD
    PathGraphHandle *h;
    int64_t n;
} PathGraphObject;

static PyTypeObject *PathGraphType = NULL;

static void path_graph_dealloc(PyObject *self) {
    PathGraphObject *p = (PathGraphObject *)self;
    PyTypeObject *tp = Py_TYPE(self);
    if (p->h) net_path_graph_free(p->h);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot path_graph_slots[] = {
    {Py_tp_dealloc, (void *)path_graph_dealloc},
    {Py_tp_doc, (void *)"Arcs with lengths for shortest-path searches."},
    {0, NULL},
};

static PyType_Spec path_graph_spec = {
    "knappnet._core.PathGraph",
    sizeof(PathGraphObject),
    0,
    Py_TPFLAGS_DEFAULT,
    path_graph_slots,
};

static PathGraphObject *as_pg(PyObject *obj) {
    if (!PyObject_TypeCheck(obj, PathGraphType)) {
        PyErr_SetString(PyExc_TypeError, "expected a knappnet._core.PathGraph");
        return NULL;
    }
    return (PathGraphObject *)obj;
}

#define PGRAPH(var, obj)                                                       \
    PathGraphObject *var = as_pg(obj);                                         \
    if (!var) return NULL

/* path_graph(g, mode, weighted, weights_as, lengths, collapse) */
static PyObject *py_path_graph(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go, *lens;
    int mode, weighted, weights_as, collapse;
    PARSE("OiiiOi", &go, &mode, &weighted, &weights_as, &lens, &collapse);
    GRAPH(g, go);
    double *lengths = NULL;
    if (opt_doubles(lens, g->m, &lengths) < 0) return NULL;
    PathGraphHandle *h = NULL;
    BEGIN_CALL
    status = net_path_graph_new(g->h, mode, weighted, weights_as, lengths, lengths ? g->m : 0,
                                collapse, &h, err, ERR_LEN);
    END_CALL
    free(lengths);
    FAIL_IF(status);
    PathGraphObject *p = PyObject_New(PathGraphObject, PathGraphType);
    if (!p) {
        net_path_graph_free(h);
        return NULL;
    }
    p->h = h;
    p->n = g->n;
    return (PyObject *)p;
}

/* sp_tree(pg, sources, cutoff, target, reverse) -> (dist, parent, parent_edge, origin) */
static PyObject *py_sp_tree(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *po, *so;
    double cutoff;
    long long target;
    int reverse;
    PARSE("OOdLi", &po, &so, &cutoff, &target, &reverse);
    PGRAPH(p, po);
    int64_t k = 0, *sources = NULL;
    if (req_int64(so, &sources, &k) < 0) return NULL;
    int64_t n = p->n;
    double *dist = alloc_doubles(n);
    int64_t *parent = alloc_int64(n), *pedge = alloc_int64(n), *origin = alloc_int64(n);
    PyObject *res = NULL;
    if (dist && parent && pedge && origin) {
        BEGIN_CALL
        status = net_shortest_paths(p->h, sources, k, cutoff, target, reverse, dist, parent,
                                    pedge, origin, err, ERR_LEN);
        END_CALL
        if (status) {
            engine_error(err);
        } else {
            res = Py_BuildValue("(NNNN)", PyBytes_FromStringAndSize((char *)dist, n * 8),
                                PyBytes_FromStringAndSize((char *)parent, n * 8),
                                PyBytes_FromStringAndSize((char *)pedge, n * 8),
                                PyBytes_FromStringAndSize((char *)origin, n * 8));
        }
    }
    free(sources);
    free(dist);
    free(parent);
    free(pedge);
    free(origin);
    return res;
}

/* sp_predecessors(pg, sources, cutoff, reverse, tie) -> (dist bytes, lists) */
static PyObject *py_sp_predecessors(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *po, *so;
    double cutoff;
    int reverse, tie;
    PARSE("OOdii", &po, &so, &cutoff, &reverse, &tie);
    PGRAPH(p, po);
    int64_t k = 0, *sources = NULL;
    if (req_int64(so, &sources, &k) < 0) return NULL;
    double *dist = alloc_doubles(p->n);
    PathListsHandle *lh = NULL;
    if (!dist) {
        free(sources);
        return NULL;
    }
    BEGIN_CALL
    status = net_shortest_path_predecessors(p->h, sources, k, cutoff, reverse, tie, dist, &lh,
                                            err, ERR_LEN);
    END_CALL
    free(sources);
    if (status) {
        free(dist);
        return engine_error(err);
    }
    PyObject *lists = path_lists_result(lh);
    if (!lists) {
        free(dist);
        return NULL;
    }
    PyObject *b = bytes_doubles(dist, p->n);
    PyObject *res = b ? PyTuple_Pack(2, b, lists) : NULL;
    Py_XDECREF(b);
    Py_DECREF(lists);
    return res;
}

/* sp_all_pairs(pg, sources|None, cutoff, threads, want_parent) ->
 * (dist bytes, parent bytes|None); column-major, one row per source. */
static PyObject *py_sp_all_pairs(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *po, *so;
    double cutoff;
    int threads, want_parent;
    PARSE("OOdii", &po, &so, &cutoff, &threads, &want_parent);
    PGRAPH(p, po);
    int64_t k = 0, *sources = NULL;
    if (opt_int64(so, &sources, &k) < 0) return NULL;
    int64_t rows = sources ? k : p->n;
    double *dist = alloc_doubles(rows * p->n);
    int64_t *parent = want_parent ? alloc_int64(rows * p->n) : NULL;
    PyObject *res = NULL;
    if (dist && (!want_parent || parent)) {
        BEGIN_CALL
        status = net_all_pairs_shortest_paths(p->h, sources, sources ? k : -1, cutoff, threads,
                                              dist, parent, err, ERR_LEN);
        END_CALL
        if (status) {
            engine_error(err);
        } else {
            PyObject *a = PyBytes_FromStringAndSize((char *)dist, rows * p->n * 8);
            PyObject *b = parent ? PyBytes_FromStringAndSize((char *)parent, rows * p->n * 8)
                                 : Py_NewRef(Py_None);
            res = (a && b) ? PyTuple_Pack(2, a, b) : NULL;
            Py_XDECREF(a);
            Py_XDECREF(b);
        }
    }
    free(sources);
    free(dist);
    free(parent);
    return res;
}

/* sp_bidirectional(pg, s, t) -> [(path, length)] */
static PyObject *py_sp_bidirectional(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *po;
    long long s, t;
    PARSE("OLL", &po, &s, &t);
    PGRAPH(p, po);
    PathListsHandle *lh = NULL;
    BEGIN_CALL
    status = net_bidirectional_shortest_path(p->h, s, t, &lh, err, ERR_LEN);
    END_CALL
    FAIL_IF(status);
    return path_lists_result(lh);
}

/* sp_negative(pg, sources, method, tie) -> (dist, parent, preds, cycle) */
static PyObject *py_sp_negative(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *po, *so;
    int method, tie;
    PARSE("OOii", &po, &so, &method, &tie);
    PGRAPH(p, po);
    int64_t k = 0, *sources = NULL;
    if (req_int64(so, &sources, &k) < 0) return NULL;
    int64_t n = p->n;
    double *dist = alloc_doubles(n);
    int64_t *parent = alloc_int64(n);
    PathListsHandle *preds = NULL, *cycle = NULL;
    PyObject *res = NULL;
    if (dist && parent) {
        BEGIN_CALL
        status = net_negative_shortest_paths(p->h, sources, k, method, tie, dist, parent, &preds,
                                             &cycle, err, ERR_LEN);
        END_CALL
        if (status) {
            engine_error(err);
        } else {
            PyObject *d = PyBytes_FromStringAndSize((char *)dist, n * 8);
            PyObject *par = PyBytes_FromStringAndSize((char *)parent, n * 8);
            PyObject *pr = preds ? path_lists_result(preds) : Py_NewRef(Py_None);
            PyObject *cy = path_lists_result(cycle);
            res = (d && par && pr && cy) ? PyTuple_Pack(4, d, par, pr, cy) : NULL;
            Py_XDECREF(d);
            Py_XDECREF(par);
            Py_XDECREF(pr);
            Py_XDECREF(cy);
        }
    }
    free(sources);
    free(dist);
    free(parent);
    return res;
}

/* sp_negative_cycle(pg, source) -> [(closed walk, nan)] */
static PyObject *py_sp_negative_cycle(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *po;
    long long source;
    PARSE("OL", &po, &source);
    PGRAPH(p, po);
    PathListsHandle *lh = NULL;
    BEGIN_CALL
    status = net_negative_cycle(p->h, source, &lh, err, ERR_LEN);
    END_CALL
    FAIL_IF(status);
    return path_lists_result(lh);
}

/* sp_johnson(pg, threads) -> (dist, parent, cycle) over every source */
static PyObject *py_sp_johnson(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *po;
    int threads;
    PARSE("Oi", &po, &threads);
    PGRAPH(p, po);
    int64_t n = p->n;
    double *dist = alloc_doubles(n * n);
    int64_t *parent = alloc_int64(n * n);
    PathListsHandle *cycle = NULL;
    PyObject *res = NULL;
    if (dist && parent) {
        BEGIN_CALL
        status = net_johnson(p->h, NULL, -1, threads, dist, parent, &cycle, err, ERR_LEN);
        END_CALL
        if (status) {
            engine_error(err);
        } else {
            PyObject *d = PyBytes_FromStringAndSize((char *)dist, n * n * 8);
            PyObject *par = PyBytes_FromStringAndSize((char *)parent, n * n * 8);
            PyObject *cy = path_lists_result(cycle);
            res = (d && par && cy) ? PyTuple_Pack(3, d, par, cy) : NULL;
            Py_XDECREF(d);
            Py_XDECREF(par);
            Py_XDECREF(cy);
        }
    }
    free(dist);
    free(parent);
    return res;
}

/* sp_floyd_warshall(pg, threads) -> (dist, pred, negative_cycle) */
static PyObject *py_sp_floyd_warshall(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *po;
    int threads;
    PARSE("Oi", &po, &threads);
    PGRAPH(p, po);
    int64_t n = p->n;
    double *dist = alloc_doubles(n * n);
    int64_t *pred = alloc_int64(n * n);
    int32_t neg = 0;
    PyObject *res = NULL;
    if (dist && pred) {
        BEGIN_CALL
        status = net_floyd_warshall(p->h, threads, dist, pred, &neg, err, ERR_LEN);
        END_CALL
        if (status) {
            engine_error(err);
        } else {
            PyObject *d = PyBytes_FromStringAndSize((char *)dist, n * n * 8);
            PyObject *pr = PyBytes_FromStringAndSize((char *)pred, n * n * 8);
            res = (d && pr) ? Py_BuildValue("(OOi)", d, pr, (int)neg) : NULL;
            Py_XDECREF(d);
            Py_XDECREF(pr);
        }
    }
    free(dist);
    free(pred);
    return res;
}

/* sp_count(pg, s, tie) -> bytes of n doubles */
static PyObject *py_sp_count(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *po;
    long long s;
    int tie;
    PARSE("OLi", &po, &s, &tie);
    PGRAPH(p, po);
    double *out = alloc_doubles(p->n);
    if (!out) return NULL;
    BEGIN_CALL
    status = net_count_shortest_paths(p->h, s, tie, out, err, ERR_LEN);
    END_CALL
    if (status) {
        free(out);
        return engine_error(err);
    }
    return bytes_doubles(out, p->n);
}

/* sp_all_paths(pg, s, t, tie, max_paths) -> [(path, length)] */
static PyObject *py_sp_all_paths(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *po;
    long long s, t, maxp;
    int tie;
    PARSE("OLLiL", &po, &s, &t, &tie, &maxp);
    PGRAPH(p, po);
    PathListsHandle *lh = NULL;
    BEGIN_CALL
    status = net_all_shortest_paths(p->h, s, t, tie, maxp, &lh, err, ERR_LEN);
    END_CALL
    FAIL_IF(status);
    return path_lists_result(lh);
}

#define PATHS_METHODS                                                          \
    {"path_graph", py_path_graph, METH_VARARGS, NULL},                         \
    {"sp_tree", py_sp_tree, METH_VARARGS, NULL},                               \
    {"sp_predecessors", py_sp_predecessors, METH_VARARGS, NULL},               \
    {"sp_all_pairs", py_sp_all_pairs, METH_VARARGS, NULL},                     \
    {"sp_bidirectional", py_sp_bidirectional, METH_VARARGS, NULL},             \
    {"sp_negative", py_sp_negative, METH_VARARGS, NULL},                       \
    {"sp_negative_cycle", py_sp_negative_cycle, METH_VARARGS, NULL},           \
    {"sp_johnson", py_sp_johnson, METH_VARARGS, NULL},                         \
    {"sp_floyd_warshall", py_sp_floyd_warshall, METH_VARARGS, NULL},           \
    {"sp_count", py_sp_count, METH_VARARGS, NULL},                             \
    {"sp_all_paths", py_sp_all_paths, METH_VARARGS, NULL},
