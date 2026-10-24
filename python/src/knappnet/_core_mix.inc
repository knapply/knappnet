/* ---- Cores, mixing, cuts, matching, colouring, similarity --------------- */

/* A list of node indices from an engine buffer. */
static PyObject *indices_result(int32_t status, const char *err, int64_t *buf, int64_t found) {
    if (status) {
        free(buf);
        return engine_error(err);
    }
    return finish_int64(buf, found);
}

/* core_subset(g, which, k, has_k) -> [nodes] */
static PyObject *py_core_subset(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int which, has_k;
    long long k;
    PARSE("OiLi", &go, &which, &k, &has_k);
    GRAPH(g, go);
    int64_t *nodes = alloc_int64(g->n), count = 0;
    if (!nodes) return NULL;
    BEGIN_CALL
    status = net_core_subset(g->h, which, k, has_k, nodes, &count, err, ERR_LEN);
    END_CALL
    return indices_result(status, err, nodes, count);
}

/* onion_layers(g) -> [layer] */
static PyObject *py_onion_layers(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    PARSE("O", &go);
    GRAPH(g, go);
    int64_t *out = alloc_int64(g->n);
    if (!out) return NULL;
    BEGIN_CALL
    status = net_onion_layers(g->h, out, err, ERR_LEN);
    END_CALL
    return indices_result(status, err, out, g->n);
}

/* k_truss(g, k) -> (edges, nodes) */
static PyObject *py_k_truss(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    long long k;
    PARSE("OL", &go, &k);
    GRAPH(g, go);
    int64_t *edges = alloc_int64(g->m), *nodes = alloc_int64(g->n), ne = 0, nn = 0;
    if (!edges || !nodes) {
        free(edges);
        free(nodes);
        return NULL;
    }
    BEGIN_CALL
    status = net_k_truss(g->h, k, edges, &ne, nodes, &nn, err, ERR_LEN);
    END_CALL
    if (status) {
        free(edges);
        free(nodes);
        return engine_error(err);
    }
    return Py_BuildValue("(NN)", finish_int64(edges, ne), finish_int64(nodes, nn));
}

/* mixing(g, kind, x, y, weighted, labels|None, nodes|None) -> (values, counts, k)
 * kind 0 = degree mixing, 1 = attribute mixing (labels) */
static PyObject *py_mixing(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go, *lo, *no;
    int kind, x, y, weighted;
    PARSE("OiiiiOO", &go, &kind, &x, &y, &weighted, &lo, &no);
    GRAPH(g, go);
    double *labels = NULL;
    int64_t *nodes = NULL, nn = -1, k = 0;
    if (opt_doubles(lo, g->n, &labels) < 0) return NULL;
    if (opt_int64(no, &nodes, &nn) < 0) {
        free(labels);
        return NULL;
    }
    char err[ERR_LEN];
    int32_t status;
    Py_BEGIN_ALLOW_THREADS
    status = kind == 0 ? net_degree_mixing(g->h, x, y, weighted, nodes, nn < 0 ? 0 : nn, 0, NULL,
                                           NULL, &k, err, ERR_LEN)
                       : net_attribute_mixing(g->h, labels, nodes, nn < 0 ? 0 : nn, 0, NULL, NULL,
                                              &k, err, ERR_LEN);
    Py_END_ALLOW_THREADS
    PyObject *res = NULL;
    if (status) {
        engine_error(err);
    } else {
        double *vals = alloc_doubles(k > 0 ? k : 1), *counts = alloc_doubles(k > 0 ? k * k : 1);
        if (vals && counts) {
            int64_t k2 = 0;
            Py_BEGIN_ALLOW_THREADS
            status = kind == 0 ? net_degree_mixing(g->h, x, y, weighted, nodes, nn < 0 ? 0 : nn, k,
                                                   vals, counts, &k2, err, ERR_LEN)
                               : net_attribute_mixing(g->h, labels, nodes, nn < 0 ? 0 : nn, k, vals,
                                                      counts, &k2, err, ERR_LEN);
            Py_END_ALLOW_THREADS
            if (status)
                engine_error(err);
            else
                res = Py_BuildValue("(NNL)", bytes_doubles(vals, k), bytes_doubles(counts, k * k),
                                    (long long)k);
            vals = counts = NULL;
        }
        free(vals);
        free(counts);
    }
    free(labels);
    free(nodes);
    return res;
}

/* assort_subset(g, coefficient, x, y, weighted, values|None, nodes|None) -> float */
static PyObject *py_assort_subset(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go, *vo, *no;
    int coef, x, y, weighted;
    PARSE("OiiiiOO", &go, &coef, &x, &y, &weighted, &vo, &no);
    GRAPH(g, go);
    double *values = NULL, out = NAN;
    int64_t *nodes = NULL, nn = -1;
    if (opt_doubles(vo, g->n, &values) < 0) return NULL;
    if (opt_int64(no, &nodes, &nn) < 0) {
        free(values);
        return NULL;
    }
    BEGIN_CALL
    status = net_assortativity_subset(g->h, coef, x, y, weighted, values, nodes, nn < 0 ? 0 : nn,
                                      &out, err, ERR_LEN);
    END_CALL
    free(values);
    free(nodes);
    if (status) return engine_error(err);
    return PyFloat_FromDouble(out);
}

/* neighbor_degree(g, source, target, weighted, nodes|None) -> bytes */
static PyObject *py_neighbor_degree(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go, *no;
    int source, target, weighted;
    PARSE("OiiiO", &go, &source, &target, &weighted, &no);
    GRAPH(g, go);
    int64_t *nodes = NULL, nn = -1;
    if (opt_int64(no, &nodes, &nn) < 0) return NULL;
    int64_t len = nn < 0 ? g->n : nn;
    double *out = alloc_doubles(len > 0 ? len : 1);
    if (!out) {
        free(nodes);
        return NULL;
    }
    BEGIN_CALL
    status = net_average_neighbor_degree(g->h, source, target, weighted, nodes, nn < 0 ? 0 : nn,
                                         out, err, ERR_LEN);
    END_CALL
    free(nodes);
    return done_bytes(status, err, out, len);
}

/* degree_connectivity(g, source, target, weighted, nodes|None) -> (degrees, values) */
static PyObject *py_degree_connectivity(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go, *no;
    int source, target, weighted;
    PARSE("OiiiO", &go, &source, &target, &weighted, &no);
    GRAPH(g, go);
    int64_t *nodes = NULL, nn = -1, count = 0;
    if (opt_int64(no, &nodes, &nn) < 0) return NULL;
    int64_t cap = g->n > 0 ? g->n : 1;
    double *deg = alloc_doubles(cap), *val = alloc_doubles(cap);
    if (!deg || !val) {
        free(deg);
        free(val);
        free(nodes);
        return NULL;
    }
    BEGIN_CALL
    status = net_average_degree_connectivity(g->h, source, target, weighted, nodes,
                                             nn < 0 ? 0 : nn, deg, val, &count, err, ERR_LEN);
    END_CALL
    free(nodes);
    if (status) {
        free(deg);
        free(val);
        return engine_error(err);
    }
    return Py_BuildValue("(NN)", bytes_doubles(deg, count), bytes_doubles(val, count));
}

/* cut_measures(g, s, t|None, weighted) -> bytes of 9 doubles */
static PyObject *py_cut_measures(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go, *so, *to;
    int weighted;
    PARSE("OOOi", &go, &so, &to, &weighted);
    GRAPH(g, go);
    int64_t *s = NULL, ns = 0, *t = NULL, nt = -1;
    if (req_int64(so, &s, &ns) < 0) return NULL;
    if (opt_int64(to, &t, &nt) < 0) {
        free(s);
        return NULL;
    }
    double *out = alloc_doubles(9);
    if (!out) {
        free(s);
        free(t);
        return NULL;
    }
    BEGIN_CALL
    status = net_cut_measures(g->h, s, ns, t, nt < 0 ? 0 : nt, weighted, out, err, ERR_LEN);
    END_CALL
    free(s);
    free(t);
    return done_bytes(status, err, out, 9);
}

/* square_clustering(g, threads) -> bytes */
static PyObject *py_square_clustering(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int threads;
    PARSE("Oi", &go, &threads);
    GRAPH(g, go);
    double *out = alloc_doubles(g->n > 0 ? g->n : 1);
    if (!out) return NULL;
    BEGIN_CALL
    status = net_square_clustering(g->h, threads, out, err, ERR_LEN);
    END_CALL
    return done_bytes(status, err, out, g->n);
}

/* generalized_degree(g, threads) -> (offsets, ks, counts) */
static PyObject *py_generalized_degree(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int threads;
    PARSE("Oi", &go, &threads);
    GRAPH(g, go);
    int64_t cap = 2 * g->m + 1, total = 0;
    int64_t *off = alloc_int64(g->n + 1), *ks = alloc_int64(cap), *cn = alloc_int64(cap);
    if (!off || !ks || !cn) {
        free(off);
        free(ks);
        free(cn);
        return NULL;
    }
    BEGIN_CALL
    status = net_generalized_degree(g->h, threads, off, ks, cn, cap, &total, err, ERR_LEN);
    END_CALL
    if (status) {
        free(off);
        free(ks);
        free(cn);
        return engine_error(err);
    }
    return Py_BuildValue("(NNN)", finish_int64(off, g->n + 1), finish_int64(ks, total),
                         finish_int64(cn, total));
}

/* matching(g, kind, weighted, max_cardinality) -> [edge indices]
 * kind 0 = max weight, 1 = min weight, 2 = maximal */
static PyObject *py_matching(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int kind, weighted, maxcard;
    PARSE("Oiii", &go, &kind, &weighted, &maxcard);
    GRAPH(g, go);
    int64_t *edges = alloc_int64(g->m > 0 ? g->m : 1), found = 0;
    if (!edges) return NULL;
    BEGIN_CALL
    status = kind == 0   ? net_max_weight_matching(g->h, weighted, maxcard, edges, &found, err, ERR_LEN)
             : kind == 1 ? net_min_weight_matching(g->h, weighted, edges, &found, err, ERR_LEN)
                         : net_maximal_matching(g->h, edges, &found, err, ERR_LEN);
    END_CALL
    return indices_result(status, err, edges, found);
}

/* check_matching(g, us, vs) -> (is_matching, is_maximal, is_perfect) */
static PyObject *py_check_matching(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go, *uo, *vo;
    PARSE("OOO", &go, &uo, &vo);
    GRAPH(g, go);
    int64_t *us = NULL, *vs = NULL, ku = 0, kv = 0;
    if (req_int64(uo, &us, &ku) < 0) return NULL;
    if (req_int64(vo, &vs, &kv) < 0 || ku != kv) {
        free(us);
        free(vs);
        if (!PyErr_Occurred()) PyErr_SetString(PyExc_ValueError, "length mismatch");
        return NULL;
    }
    int32_t out[3] = {0, 0, 0};
    BEGIN_CALL
    status = net_check_matching(g->h, us, vs, ku, out, err, ERR_LEN);
    END_CALL
    free(us);
    free(vs);
    if (status) return engine_error(err);
    return Py_BuildValue("(OOO)", out[0] ? Py_True : Py_False, out[1] ? Py_True : Py_False,
                         out[2] ? Py_True : Py_False);
}

/* greedy_color(g, strategy) -> [colours] */
static PyObject *py_greedy_color(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int strategy;
    PARSE("Oi", &go, &strategy);
    GRAPH(g, go);
    int64_t *out = alloc_int64(g->n > 0 ? g->n : 1);
    if (!out) return NULL;
    BEGIN_CALL
    status = net_greedy_color(g->h, strategy, 0, out, err, ERR_LEN);
    END_CALL
    return indices_result(status, err, out, g->n);
}

/* link_prediction(g, method, alpha, delta, communities|None, us|None, vs|None, threads)
 * -> bytes (pairs) or (us, vs, scores) for every non-edge */
static PyObject *py_link_prediction(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go, *co, *uo, *vo;
    int method, threads;
    double alpha, delta;
    PARSE("OiddOOOi", &go, &method, &alpha, &delta, &co, &uo, &vo, &threads);
    GRAPH(g, go);
    int64_t *comm = NULL, nc = -1, *us = NULL, *vs = NULL, np = 0, nv = 0;
    if (opt_int64(co, &comm, &nc) < 0) return NULL;
    if (uo != Py_None) {
        if (req_int64(uo, &us, &np) < 0 || req_int64(vo, &vs, &nv) < 0 || np != nv) {
            free(comm);
            free(us);
            free(vs);
            if (!PyErr_Occurred()) PyErr_SetString(PyExc_ValueError, "length mismatch");
            return NULL;
        }
        double *out = alloc_doubles(np > 0 ? np : 1);
        if (!out) {
            free(comm);
            free(us);
            free(vs);
            return NULL;
        }
        BEGIN_CALL
        status = net_link_prediction(g->h, method, alpha, delta, comm, us, vs, np, threads, out,
                                     err, ERR_LEN);
        END_CALL
        free(comm);
        free(us);
        free(vs);
        return done_bytes(status, err, out, np);
    }
    int64_t len = 0;
    char err[ERR_LEN];
    int32_t status = net_non_edge_count(g->h, &len, err, ERR_LEN);
    if (status) {
        free(comm);
        return engine_error(err);
    }
    int64_t *a = alloc_int64(len > 0 ? len : 1), *b = alloc_int64(len > 0 ? len : 1);
    double *sc = alloc_doubles(len > 0 ? len : 1);
    if (!a || !b || !sc) {
        free(a);
        free(b);
        free(sc);
        free(comm);
        return NULL;
    }
    Py_BEGIN_ALLOW_THREADS
    status = net_link_prediction_non_edges(g->h, method, alpha, delta, comm, threads, a, b, sc, len,
                                           err, ERR_LEN);
    Py_END_ALLOW_THREADS
    free(comm);
    if (status) {
        free(a);
        free(b);
        free(sc);
        return engine_error(err);
    }
    return Py_BuildValue("(NNN)", finish_int64(a, len), finish_int64(b, len), bytes_doubles(sc, len));
}

/* Borrowed UTF-8 strings from a list of str (or None), for the hash calls. */
static const char **strings_from_list(PyObject *seq, int64_t len) {
    if (seq == Py_None) return NULL;
    const char **out = malloc((size_t)(len > 0 ? len : 1) * sizeof(char *));
    if (!out) {
        PyErr_NoMemory();
        return NULL;
    }
    for (int64_t i = 0; i < len; i++) {
        PyObject *s = PyList_GetItem(seq, i);
        out[i] = s ? PyUnicode_AsUTF8AndSize(s, NULL) : NULL;
        if (!out[i]) {
            free(out);
            return NULL;
        }
    }
    return out;
}

/* wl_hash(g, node_labels|None, edge_labels|None, iterations, digest_size,
 * subgraph, include_initial, threads) -> str (graph hash) or list of str */
static PyObject *py_wl_hash(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go, *nl, *el;
    long long iterations, digest;
    int subgraph, initial, threads;
    PARSE("OOOLLiii", &go, &nl, &el, &iterations, &digest, &subgraph, &initial, &threads);
    GRAPH(g, go);
    if ((nl != Py_None && !PyList_Check(nl)) || (el != Py_None && !PyList_Check(el))) {
        PyErr_SetString(PyExc_TypeError, "labels must be lists of str");
        return NULL;
    }
    const char **nlab = strings_from_list(nl, g->n), **elab = strings_from_list(el, g->m);
    if (PyErr_Occurred()) {
        free(nlab);
        free(elab);
        return NULL;
    }
    int64_t w = 2 * digest, per = iterations + (initial ? 1 : 0);
    size_t len = subgraph ? (size_t)(g->n * per * w) : (size_t)w;
    char *buf = malloc(len + 1);
    if (!buf) {
        free(nlab);
        free(elab);
        return PyErr_NoMemory();
    }
    BEGIN_CALL
    status = subgraph ? net_wl_subgraph_hashes(g->h, nlab, elab, iterations, digest, initial, threads,
                                               buf, len, err, ERR_LEN)
                      : net_wl_graph_hash(g->h, nlab, elab, iterations, digest, threads, buf, len,
                                          err, ERR_LEN);
    END_CALL
    free(nlab);
    free(elab);
    PyObject *res = NULL;
    if (status) {
        engine_error(err);
    } else if (!subgraph) {
        res = PyUnicode_FromStringAndSize(buf, (Py_ssize_t)len);
    } else {
        res = PyList_New(g->n);
        for (int64_t v = 0; res && v < g->n; v++) {
            PyObject *row = PyList_New(per);
            if (!row) {
                Py_CLEAR(res);
                break;
            }
            for (int64_t i = 0; i < per; i++) {
                PyObject *s = PyUnicode_FromStringAndSize(buf + (v * per + i) * w, (Py_ssize_t)w);
                if (!s) {
                    Py_DECREF(row);
                    Py_CLEAR(res);
                    break;
                }
                PyList_SetItem(row, i, s);
            }
            if (res) PyList_SetItem(res, v, row);
        }
    }
    free(buf);
    return res;
}

/* simrank(g, importance, max_iterations, tolerance, weighted, stop, threads) -> bytes (n*n) */
static PyObject *py_simrank(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    double importance, tol;
    long long maxit;
    int weighted, stop, threads;
    PARSE("OdLdiii", &go, &importance, &maxit, &tol, &weighted, &stop, &threads);
    GRAPH(g, go);
    double *out = alloc_doubles(g->n * g->n > 0 ? g->n * g->n : 1);
    if (!out) return NULL;
    BEGIN_CALL
    status = net_simrank(g->h, importance, maxit, tol, weighted, stop, threads, out, err, ERR_LEN);
    END_CALL
    return done_bytes(status, err, out, g->n * g->n);
}

/* isomorphism_count(g1, g2, labels1|None, labels2|None, problem, empty_matches, limit) -> int */
static PyObject *py_vf2pp_count(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go1, *go2, *lo1, *lo2;
    int problem, empty;
    long long limit;
    PARSE("OOOOiiL", &go1, &go2, &lo1, &lo2, &problem, &empty, &limit);
    GRAPH(g1, go1);
    GRAPH(g2, go2);
    int64_t *l1 = NULL, *l2 = NULL, n1 = -1, n2 = -1, out = 0;
    if (opt_int64(lo1, &l1, &n1) < 0) return NULL;
    if (opt_int64(lo2, &l2, &n2) < 0) {
        free(l1);
        return NULL;
    }
    BEGIN_CALL
    status = net_vf2pp_count(g1->h, g2->h, l1, l2, problem, empty, limit, &out, err, ERR_LEN);
    END_CALL
    free(l1);
    free(l2);
    if (status) return engine_error(err);
    return PyLong_FromLongLong(out);
}

/* could_be_isomorphic(g1, g2, properties, threads) -> bool */
static PyObject *py_could_be_isomorphic(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go1, *go2;
    int props, threads;
    PARSE("OOii", &go1, &go2, &props, &threads);
    GRAPH(g1, go1);
    GRAPH(g2, go2);
    int32_t out = 0;
    BEGIN_CALL
    status = net_could_be_isomorphic(g1->h, g2->h, props, threads, &out, err, ERR_LEN);
    END_CALL
    if (status) return engine_error(err);
    return PyBool_FromLong(out);
}

/* edit_distance(g1, g2, root1, root2, upper_bound, max_steps) -> float or None */
static PyObject *py_edit_distance(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go1, *go2;
    long long root1, root2, steps;
    double upper;
    PARSE("OOLLdL", &go1, &go2, &root1, &root2, &upper, &steps);
    GRAPH(g1, go1);
    GRAPH(g2, go2);
    double out = 0;
    int32_t found = 0;
    BEGIN_CALL
    status = net_graph_edit_distance(g1->h, g2->h, NULL, NULL, NULL, NULL, NULL, root1, root2,
                                     upper, steps, &out, &found, err, ERR_LEN);
    END_CALL
    if (status) return engine_error(err);
    if (!found) Py_RETURN_NONE;
    return PyFloat_FromDouble(out);
}

#define MIX_METHODS                                                                              \
    {"core_subset", py_core_subset, METH_VARARGS, NULL},                                         \
    {"onion_layers", py_onion_layers, METH_VARARGS, NULL},                                       \
    {"k_truss", py_k_truss, METH_VARARGS, NULL},                                                 \
    {"mixing", py_mixing, METH_VARARGS, NULL},                                                   \
    {"assort_subset", py_assort_subset, METH_VARARGS, NULL},                                     \
    {"neighbor_degree", py_neighbor_degree, METH_VARARGS, NULL},                                 \
    {"degree_connectivity", py_degree_connectivity, METH_VARARGS, NULL},                         \
    {"cut_measures", py_cut_measures, METH_VARARGS, NULL},                                       \
    {"square_clustering", py_square_clustering, METH_VARARGS, NULL},                             \
    {"generalized_degree", py_generalized_degree, METH_VARARGS, NULL},                           \
    {"matching", py_matching, METH_VARARGS, NULL},                                               \
    {"check_matching", py_check_matching, METH_VARARGS, NULL},                                   \
    {"greedy_color", py_greedy_color, METH_VARARGS, NULL},                                       \
    {"link_prediction", py_link_prediction, METH_VARARGS, NULL},                                 \
    {"wl_hash", py_wl_hash, METH_VARARGS, NULL},                                                 \
    {"simrank", py_simrank, METH_VARARGS, NULL},                                                 \
    {"vf2pp_count", py_vf2pp_count, METH_VARARGS, NULL},                                         \
    {"could_be_isomorphic", py_could_be_isomorphic, METH_VARARGS, NULL},                         \
    {"edit_distance", py_edit_distance, METH_VARARGS, NULL},
