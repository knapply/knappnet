/*
 * Shared helpers for the wrappers in the _core_*.inc files (included into
 * _core.c after the graph type and the basic conversions).
 *
 * Conventions: wrappers take a Graph (or PathGraph) and plain Python
 * sequences and return lists, tuples, or bytes holding doubles; the Python
 * layer maps node indices to labels. Optional arrays may be None.
 */

#define PARSE(...)                                                             \
    if (!PyArg_ParseTuple(args, __VA_ARGS__)) return NULL
#define GRAPH(var, obj)                                                        \
    GraphObject *var = as_graph(obj);                                          \
    if (!var) return NULL

/* Reads an optional sequence of doubles (None gives NULL, len -1). A
 * required length is checked when expected >= 0. */
static int opt_doubles(PyObject *obj, int64_t expected, double **out) {
    *out = NULL;
    if (obj == Py_None) return 0;
    int64_t len = expected;
    if (expected < 0) {
        PyObject *l = PySequence_List(obj);
        if (!l) return -1;
        len = PyList_Size(l);
        Py_DECREF(l);
    }
    *out = doubles_from_seq(obj, len);
    return *out ? 0 : -1;
}

/* Reads an optional sequence of integers (None gives NULL, *len = -1). */
static int opt_int64(PyObject *obj, int64_t **out, int64_t *len) {
    *out = NULL;
    *len = -1;
    if (obj == Py_None) return 0;
    *out = int64_from_seq(obj, len);
    return *out ? 0 : -1;
}

/* Reads a required sequence of integers (an empty one gives a valid
 * pointer and length 0). */
static int req_int64(PyObject *obj, int64_t **out, int64_t *len) {
    *out = int64_from_seq(obj, len);
    return *out ? 0 : -1;
}

/* Consumes x and returns its doubles as a bytes object (native order). */
static PyObject *bytes_doubles(double *x, int64_t len) {
    PyObject *b = PyBytes_FromStringAndSize((const char *)x, (Py_ssize_t)(len * (int64_t)sizeof(double)));
    free(x);
    return b;
}

static PyObject *bytes_int64(int64_t *x, int64_t len) {
    PyObject *b = PyBytes_FromStringAndSize((const char *)x, (Py_ssize_t)(len * (int64_t)sizeof(int64_t)));
    free(x);
    return b;
}

/* Consumes x and returns it as a list of ints. */
static PyObject *finish_int64(int64_t *x, int64_t len) {
    PyObject *out = list_from_int64(x, len);
    free(x);
    return out;
}

static PyObject *list_from_int32(const int32_t *x, int64_t len) {
    PyObject *out = PyList_New(len);
    if (!out) return NULL;
    for (int64_t i = 0; i < len; i++) {
        PyObject *v = PyLong_FromLong(x[i]);
        if (!v) {
            Py_DECREF(out);
            return NULL;
        }
        PyList_SetItem(out, i, v);
    }
    return out;
}

/* A list of (items, value) pairs from offsets/items/values. */
static PyObject *pack_lists(int64_t count, const int64_t *off, const int64_t *items,
                            const double *vals) {
    PyObject *out = PyList_New(count);
    if (!out) return NULL;
    for (int64_t i = 0; i < count; i++) {
        PyObject *lst = list_from_int64(items + off[i], off[i + 1] - off[i]);
        PyObject *val = lst ? PyFloat_FromDouble(vals ? vals[i] : NAN) : NULL;
        PyObject *pair = val ? PyTuple_Pack(2, lst, val) : NULL;
        Py_XDECREF(lst);
        Py_XDECREF(val);
        if (!pair) {
            Py_DECREF(out);
            return NULL;
        }
        PyList_SetItem(out, i, pair);
    }
    return out;
}

/* Consumes a ListsHandle: [(items, value)]. */
static PyObject *lists_result(ListsHandle *lh) {
    int64_t count = net_lists_count(lh), total = net_lists_total(lh);
    int64_t *off = alloc_int64(count + 1);
    int64_t *items = alloc_int64(total);
    double *vals = alloc_doubles(count);
    PyObject *res = NULL;
    char err[ERR_LEN];
    if (off && items && vals) {
        if (net_lists_copy(lh, off, items, vals, err, ERR_LEN))
            engine_error(err);
        else
            res = pack_lists(count, off, items, vals);
    }
    net_lists_free(lh);
    free(off);
    free(items);
    free(vals);
    return res;
}

/* Consumes a PathListsHandle: [(nodes, value)]. */
static PyObject *path_lists_result(PathListsHandle *lh) {
    int64_t count = net_path_lists_count(lh), total = net_path_lists_total(lh);
    int64_t *off = alloc_int64(count + 1);
    int64_t *items = alloc_int64(total);
    double *vals = alloc_doubles(count);
    PyObject *res = NULL;
    char err[ERR_LEN];
    if (off && items && vals) {
        if (net_path_lists_copy(lh, off, items, vals, err, ERR_LEN))
            engine_error(err);
        else
            res = pack_lists(count, off, items, vals);
    }
    net_path_lists_free(lh);
    free(off);
    free(items);
    free(vals);
    return res;
}

/* Two output doubles/ints from an engine call that needs a value count. */
#define FAIL_IF(status)                                                        \
    if (status) return engine_error(err)

/* edge_weights(g) -> bytes of the m converted edge weights */
static PyObject *py_edge_weights(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    PARSE("O", &go);
    GRAPH(g, go);
    if (!g->w) return PyBytes_FromStringAndSize("", 0);
    return PyBytes_FromStringAndSize((const char *)g->w, (Py_ssize_t)(g->m * 8));
}

#define UTIL_METHODS {"edge_weights", py_edge_weights, METH_VARARGS, NULL},
