/* ---- Centrality variants and distance measures --------------------------- */

static PyObject *done_bytes(int32_t status, const char *err, double *out, int64_t len) {
    if (status) {
        free(out);
        return engine_error(err);
    }
    return bytes_doubles(out, len);
}

/* A pair (a, b) of bytes objects or None, consuming the buffers. */
static PyObject *pair_bytes(double *a, int64_t alen, double *b, int64_t blen) {
    PyObject *x = a ? bytes_doubles(a, alen) : Py_NewRef(Py_None);
    PyObject *y = b ? bytes_doubles(b, blen) : Py_NewRef(Py_None);
    if (!x || !y) {
        Py_XDECREF(x);
        Py_XDECREF(y);
        return NULL;
    }
    return Py_BuildValue("(NN)", x, y);
}

/* betweenness_subset(g, sources, targets, weighted, weights_as, endpoints,
 * normalized, rescale, parallel, edge_split, tie, threads, want_edges)
 * -> (nodes, edges|None); sources/targets None = every node */
static PyObject *py_betweenness_subset(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go, *so, *to;
    int weighted, weights_as, endpoints, normalized, rescale, parallel, split, tie, threads, want_edges;
    PARSE("OOOiiiiiiiiii", &go, &so, &to, &weighted, &weights_as, &endpoints, &normalized,
          &rescale, &parallel, &split, &tie, &threads, &want_edges);
    GRAPH(g, go);
    int64_t ns = 0, nt = 0, *src = NULL, *tgt = NULL;
    if (opt_int64(so, &src, &ns) < 0) return NULL;
    if (opt_int64(to, &tgt, &nt) < 0) {
        free(src);
        return NULL;
    }
    double *nodes = alloc_doubles(g->n);
    double *edges = want_edges ? alloc_doubles(g->m) : NULL;
    PyObject *res = NULL;
    if (nodes && (!want_edges || edges)) {
        BEGIN_CALL
        status = net_betweenness_subset(g->h, src, ns, tgt, nt, weighted, weights_as, endpoints,
                                        normalized, rescale, parallel, split, tie, threads,
                                        nodes, edges, err, ERR_LEN);
        END_CALL
        if (status) {
            engine_error(err);
            free(nodes);
            free(edges);
        } else {
            res = pair_bytes(nodes, g->n, edges, g->m);
        }
    } else {
        free(nodes);
        free(edges);
    }
    free(src);
    free(tgt);
    return res;
}

/* percolation(g, states, weighted, weights_as, parallel, tie, threads) */
static PyObject *py_percolation(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go, *st;
    int weighted, weights_as, parallel, tie, threads;
    PARSE("OOiiiii", &go, &st, &weighted, &weights_as, &parallel, &tie, &threads);
    GRAPH(g, go);
    double *states = doubles_from_seq(st, g->n);
    if (!states) return NULL;
    double *out = alloc_doubles(g->n);
    if (!out) {
        free(states);
        return NULL;
    }
    BEGIN_CALL
    status = net_percolation_centrality(g->h, states, weighted, weights_as, parallel, tie,
                                        threads, out, err, ERR_LEN);
    END_CALL
    free(states);
    return done_bytes(status, err, out, g->n);
}

/* load(g, weighted, weights_as, normalized, cutoff, networkx_ties, tie, threads) */
static PyObject *py_load(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int weighted, weights_as, normalized, nxties, tie, threads;
    double cutoff;
    PARSE("Oiiidiii", &go, &weighted, &weights_as, &normalized, &cutoff, &nxties, &tie, &threads);
    GRAPH(g, go);
    double *out = alloc_doubles(g->n);
    if (!out) return NULL;
    BEGIN_CALL
    status = net_load_centrality(g->h, weighted, weights_as, normalized, cutoff, nxties, tie,
                                 threads, out, err, ERR_LEN);
    END_CALL
    return done_bytes(status, err, out, g->n);
}

/* edge_load(g, cutoff, threads) -> 2m doubles */
static PyObject *py_edge_load(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    double cutoff;
    int threads;
    PARSE("Odi", &go, &cutoff, &threads);
    GRAPH(g, go);
    double *out = alloc_doubles(2 * g->m);
    if (!out) return NULL;
    BEGIN_CALL
    status = net_edge_load(g->h, cutoff, threads, out, err, ERR_LEN);
    END_CALL
    return done_bytes(status, err, out, 2 * g->m);
}

/* cf_betweenness(g, weighted, weights_as, normalized, threads, want_edges) */
static PyObject *py_cf_betweenness(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int weighted, weights_as, normalized, threads, want_edges;
    PARSE("Oiiiii", &go, &weighted, &weights_as, &normalized, &threads, &want_edges);
    GRAPH(g, go);
    double *nodes = alloc_doubles(g->n);
    double *edges = want_edges ? alloc_doubles(g->m) : NULL;
    if (!nodes || (want_edges && !edges)) {
        free(nodes);
        free(edges);
        return NULL;
    }
    BEGIN_CALL
    status = net_current_flow_betweenness(g->h, weighted, weights_as, normalized, threads, nodes,
                                          edges, err, ERR_LEN);
    END_CALL
    if (status) {
        free(nodes);
        free(edges);
        return engine_error(err);
    }
    return pair_bytes(nodes, g->n, edges, g->m);
}

/* cf_betweenness_subset(g, sources, targets, networkx, weighted, weights_as,
 * normalized, threads, want_edges) */
static PyObject *py_cf_betweenness_subset(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go, *so, *to;
    int nxflag, weighted, weights_as, normalized, threads, want_edges;
    PARSE("OOOiiiiii", &go, &so, &to, &nxflag, &weighted, &weights_as, &normalized, &threads,
          &want_edges);
    GRAPH(g, go);
    int64_t ns = 0, nt = 0, *src = NULL, *tgt = NULL;
    if (req_int64(so, &src, &ns) < 0) return NULL;
    if (req_int64(to, &tgt, &nt) < 0) {
        free(src);
        return NULL;
    }
    double *nodes = alloc_doubles(g->n);
    double *edges = want_edges ? alloc_doubles(g->m) : NULL;
    PyObject *res = NULL;
    if (nodes && (!want_edges || edges)) {
        BEGIN_CALL
        status = net_current_flow_betweenness_subset(g->h, src, ns, tgt, nt, nxflag, weighted,
                                                     weights_as, normalized, threads, nodes,
                                                     edges, err, ERR_LEN);
        END_CALL
        if (status) {
            engine_error(err);
            free(nodes);
            free(edges);
        } else {
            res = pair_bytes(nodes, g->n, edges, g->m);
        }
    } else {
        free(nodes);
        free(edges);
    }
    free(src);
    free(tgt);
    return res;
}

/* cf_closeness(g, weighted, weights_as, threads) */
static PyObject *py_cf_closeness(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int weighted, weights_as, threads;
    PARSE("Oiii", &go, &weighted, &weights_as, &threads);
    GRAPH(g, go);
    double *out = alloc_doubles(g->n);
    if (!out) return NULL;
    BEGIN_CALL
    status = net_current_flow_closeness(g->h, weighted, weights_as, threads, out, err, ERR_LEN);
    END_CALL
    return done_bytes(status, err, out, g->n);
}

/* resistance_distance(g, weighted, weights_as, threads) -> n*n doubles */
static PyObject *py_resistance_distance(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int weighted, weights_as, threads;
    PARSE("Oiii", &go, &weighted, &weights_as, &threads);
    GRAPH(g, go);
    double *out = alloc_doubles(g->n * g->n);
    if (!out) return NULL;
    BEGIN_CALL
    status = net_resistance_distance(g->h, weighted, weights_as, threads, out, err, ERR_LEN);
    END_CALL
    return done_bytes(status, err, out, g->n * g->n);
}

static PyObject *py_effective_resistance(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int weighted, weights_as, threads;
    PARSE("Oiii", &go, &weighted, &weights_as, &threads);
    GRAPH(g, go);
    double out = 0;
    BEGIN_CALL
    status = net_effective_graph_resistance(g->h, weighted, weights_as, threads, &out, err, ERR_LEN);
    END_CALL
    FAIL_IF(status);
    return PyFloat_FromDouble(out);
}

static PyObject *py_kemeny(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int weighted;
    PARSE("Oi", &go, &weighted);
    GRAPH(g, go);
    double out = 0;
    BEGIN_CALL
    status = net_kemeny_constant(g->h, weighted, &out, err, ERR_LEN);
    END_CALL
    FAIL_IF(status);
    return PyFloat_FromDouble(out);
}

static PyObject *py_second_order(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int weighted;
    PARSE("Oi", &go, &weighted);
    GRAPH(g, go);
    double *out = alloc_doubles(g->n);
    if (!out) return NULL;
    BEGIN_CALL
    status = net_second_order_centrality(g->h, weighted, out, err, ERR_LEN);
    END_CALL
    return done_bytes(status, err, out, g->n);
}

static PyObject *py_communicability(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    PARSE("O", &go);
    GRAPH(g, go);
    double *out = alloc_doubles(g->n * g->n);
    if (!out) return NULL;
    BEGIN_CALL
    status = net_communicability(g->h, out, err, ERR_LEN);
    END_CALL
    return done_bytes(status, err, out, g->n * g->n);
}

static PyObject *py_communicability_betweenness(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int threads;
    PARSE("Oi", &go, &threads);
    GRAPH(g, go);
    double *out = alloc_doubles(g->n);
    if (!out) return NULL;
    BEGIN_CALL
    status = net_communicability_betweenness(g->h, threads, out, err, ERR_LEN);
    END_CALL
    return done_bytes(status, err, out, g->n);
}

/* dispersion(g, pairs (flat u0, v0, u1, v1, ...), normalized, alpha, b, c) */
static PyObject *py_dispersion(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go, *po;
    int normalized;
    double alpha, b, c;
    PARSE("OOiddd", &go, &po, &normalized, &alpha, &b, &c);
    GRAPH(g, go);
    int64_t len = 0, *pairs = NULL;
    if (req_int64(po, &pairs, &len) < 0) return NULL;
    int64_t k = len / 2;
    double *out = alloc_doubles(k);
    if (!out) {
        free(pairs);
        return NULL;
    }
    BEGIN_CALL
    status = net_dispersion(g->h, pairs, k, normalized, alpha, b, c, out, err, ERR_LEN);
    END_CALL
    free(pairs);
    return done_bytes(status, err, out, k);
}

static PyObject *py_voterank(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    long long count;
    PARSE("OL", &go, &count);
    GRAPH(g, go);
    int64_t *out = alloc_int64(g->n);
    if (!out) return NULL;
    int64_t found = 0;
    BEGIN_CALL
    status = net_voterank(g->h, count, out, &found, err, ERR_LEN);
    END_CALL
    if (status) {
        free(out);
        return engine_error(err);
    }
    return finish_int64(out, found);
}

static PyObject *py_laplacian_centrality(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int weighted, normalized;
    PARSE("Oii", &go, &weighted, &normalized);
    GRAPH(g, go);
    double *out = alloc_doubles(g->n);
    if (!out) return NULL;
    BEGIN_CALL
    status = net_laplacian_centrality(g->h, weighted, normalized, out, err, ERR_LEN);
    END_CALL
    return done_bytes(status, err, out, g->n);
}

/* trophic(g, what, weighted, cannibalism): what 0 -> n, 1 -> m, 2 -> 1 values */
static PyObject *py_trophic(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int what, weighted, cannibalism;
    PARSE("Oiii", &go, &what, &weighted, &cannibalism);
    GRAPH(g, go);
    int64_t len = what == 0 ? g->n : what == 1 ? g->m : 1;
    double *out = alloc_doubles(len);
    if (!out) return NULL;
    BEGIN_CALL
    status = net_trophic(g->h, what, weighted, cannibalism, out, err, ERR_LEN);
    END_CALL
    return done_bytes(status, err, out, len);
}

/* reaching(g, weighted, normalized, threads) -> (local bytes, global float) */
static PyObject *py_reaching(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int weighted, normalized, threads;
    PARSE("Oiii", &go, &weighted, &normalized, &threads);
    GRAPH(g, go);
    double *local = alloc_doubles(g->n);
    double global = 0;
    if (!local) return NULL;
    BEGIN_CALL
    status = net_reaching_centrality(g->h, weighted, normalized, threads, local, &global, err,
                                     ERR_LEN);
    END_CALL
    if (status) {
        free(local);
        return engine_error(err);
    }
    return Py_BuildValue("(Nd)", bytes_doubles(local, g->n), global);
}

/* center_periphery(g, which, reachable, mode, weighted, weights_as, threads) */
static PyObject *py_center_periphery(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int which, reachable, mode, weighted, weights_as, threads;
    PARSE("Oiiiiii", &go, &which, &reachable, &mode, &weighted, &weights_as, &threads);
    GRAPH(g, go);
    int64_t *out = alloc_int64(g->n);
    int64_t count = 0;
    if (!out) return NULL;
    BEGIN_CALL
    status = net_center_periphery(g->h, which, reachable, mode, weighted, weights_as, threads,
                                  out, &count, err, ERR_LEN);
    END_CALL
    if (status) {
        free(out);
        return engine_error(err);
    }
    return finish_int64(out, count);
}

/* barycenter(g, mode, weighted, weights_as, tie, threads) -> (nodes, totals bytes) */
static PyObject *py_barycenter(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int mode, weighted, weights_as, tie, threads;
    PARSE("Oiiiii", &go, &mode, &weighted, &weights_as, &tie, &threads);
    GRAPH(g, go);
    int64_t *out = alloc_int64(g->n);
    double *totals = alloc_doubles(g->n);
    int64_t count = 0;
    if (!out || !totals) {
        free(out);
        free(totals);
        return NULL;
    }
    BEGIN_CALL
    status = net_barycenter(g->h, mode, weighted, weights_as, tie, threads, out, &count, totals,
                            err, ERR_LEN);
    END_CALL
    if (status) {
        free(out);
        free(totals);
        return engine_error(err);
    }
    return Py_BuildValue("(NN)", finish_int64(out, count), bytes_doubles(totals, g->n));
}

/* distance_index(g, index, mode, weighted, weights_as, threads) */
static PyObject *py_distance_index(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int index, mode, weighted, weights_as, threads;
    PARSE("Oiiiii", &go, &index, &mode, &weighted, &weights_as, &threads);
    GRAPH(g, go);
    double out = 0;
    BEGIN_CALL
    status = net_distance_index(g->h, index, mode, weighted, weights_as, threads, &out, err,
                                ERR_LEN);
    END_CALL
    FAIL_IF(status);
    return PyFloat_FromDouble(out);
}

static PyObject *py_local_efficiency(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int loop_is_neighbor, threads;
    PARSE("Oii", &go, &loop_is_neighbor, &threads);
    GRAPH(g, go);
    double out = 0;
    BEGIN_CALL
    status = net_local_efficiency(g->h, loop_is_neighbor, threads, &out, err, ERR_LEN);
    END_CALL
    FAIL_IF(status);
    return PyFloat_FromDouble(out);
}

static PyObject *py_closeness_vitality(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int mode, weighted, weights_as, threads;
    PARSE("Oiiii", &go, &mode, &weighted, &weights_as, &threads);
    GRAPH(g, go);
    double *out = alloc_doubles(g->n);
    if (!out) return NULL;
    BEGIN_CALL
    status = net_closeness_vitality(g->h, mode, weighted, weights_as, threads, out, err, ERR_LEN);
    END_CALL
    return done_bytes(status, err, out, g->n);
}

#define CX_METHODS                                                             \
    {"betweenness_subset", py_betweenness_subset, METH_VARARGS, NULL},         \
    {"percolation", py_percolation, METH_VARARGS, NULL},                       \
    {"load", py_load, METH_VARARGS, NULL},                                     \
    {"edge_load", py_edge_load, METH_VARARGS, NULL},                           \
    {"cf_betweenness", py_cf_betweenness, METH_VARARGS, NULL},                 \
    {"cf_betweenness_subset", py_cf_betweenness_subset, METH_VARARGS, NULL},   \
    {"cf_closeness", py_cf_closeness, METH_VARARGS, NULL},                     \
    {"resistance_distance", py_resistance_distance, METH_VARARGS, NULL},       \
    {"effective_resistance", py_effective_resistance, METH_VARARGS, NULL},     \
    {"kemeny", py_kemeny, METH_VARARGS, NULL},                                 \
    {"second_order", py_second_order, METH_VARARGS, NULL},                     \
    {"communicability", py_communicability, METH_VARARGS, NULL},               \
    {"communicability_betweenness", py_communicability_betweenness, METH_VARARGS, NULL}, \
    {"dispersion", py_dispersion, METH_VARARGS, NULL},                         \
    {"voterank", py_voterank, METH_VARARGS, NULL},                             \
    {"laplacian_centrality", py_laplacian_centrality, METH_VARARGS, NULL},     \
    {"trophic", py_trophic, METH_VARARGS, NULL},                               \
    {"reaching", py_reaching, METH_VARARGS, NULL},                             \
    {"center_periphery", py_center_periphery, METH_VARARGS, NULL},             \
    {"barycenter", py_barycenter, METH_VARARGS, NULL},                         \
    {"distance_index", py_distance_index, METH_VARARGS, NULL},                 \
    {"local_efficiency", py_local_efficiency, METH_VARARGS, NULL},             \
    {"closeness_vitality", py_closeness_vitality, METH_VARARGS, NULL},
