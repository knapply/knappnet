/* ---- Flow, connectivity, trees, DAGs and cycles -------------------------- */

/* Required lists of node indices, or None. */
#define NODE_LIST(var, count, obj)                                                               \
    int64_t *var = NULL, count = 0;                                                              \
    if (req_int64(obj, &var, &count) < 0) return NULL

/* A pair (list, float) result from a status, consuming the buffers. */
static PyObject *int_list_double(int32_t status, const char *err, int64_t *buf, int64_t found,
                                 double value) {
    if (status) {
        free(buf);
        return engine_error(err);
    }
    return Py_BuildValue("(Nd)", finish_int64(buf, found), value);
}

static PyObject *lists_from(int32_t status, const char *err, ListsHandle *lh) {
    if (status) return engine_error(err);
    return lists_result(lh);
}

/* max_flow(g, sources, sinks, capacity|None, want_flow) -> (value, flow bytes|None) */
static PyObject *py_max_flow(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go, *so, *to, *co;
    int want;
    PARSE("OOOOi", &go, &so, &to, &co, &want);
    GRAPH(g, go);
    NODE_LIST(src, ns, so);
    int64_t *snk = NULL, nt = 0;
    double *cap = NULL;
    if (req_int64(to, &snk, &nt) < 0 || opt_doubles(co, g->m, &cap) < 0) {
        free(src);
        free(snk);
        free(cap);
        return NULL;
    }
    double value = 0, *flow = want ? alloc_doubles(g->m > 0 ? g->m : 1) : NULL;
    BEGIN_CALL
    status = net_max_flow(g->h, src, ns, snk, nt, cap, &value, flow, err, ERR_LEN);
    END_CALL
    free(src);
    free(snk);
    free(cap);
    if (status) {
        free(flow);
        return engine_error(err);
    }
    PyObject *f = flow ? bytes_doubles(flow, g->m) : Py_NewRef(Py_None);
    return f ? Py_BuildValue("(dN)", value, f) : NULL;
}

/* min_cut(g, sources, sinks, capacity|None, side) -> (value, source_side list, cut edges list) */
static PyObject *py_min_cut(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go, *so, *to, *co;
    int side;
    PARSE("OOOOi", &go, &so, &to, &co, &side);
    GRAPH(g, go);
    NODE_LIST(src, ns, so);
    int64_t *snk = NULL, nt = 0;
    double *cap = NULL;
    if (req_int64(to, &snk, &nt) < 0 || opt_doubles(co, g->m, &cap) < 0) {
        free(src);
        free(snk);
        free(cap);
        return NULL;
    }
    double value = 0;
    int32_t *sides = malloc(sizeof(int32_t) * (size_t)(g->n > 0 ? g->n : 1));
    int32_t *cuts = malloc(sizeof(int32_t) * (size_t)(g->m > 0 ? g->m : 1));
    PyObject *res = NULL;
    if (sides && cuts) {
        BEGIN_CALL
        status = net_min_cut(g->h, src, ns, snk, nt, cap, side, &value, sides, cuts, err, ERR_LEN);
        END_CALL
        if (status) {
            engine_error(err);
        } else {
            res = Py_BuildValue("(dNN)", value, list_from_int32(sides, g->n),
                                list_from_int32(cuts, g->m));
        }
    } else {
        PyErr_NoMemory();
    }
    free(sides);
    free(cuts);
    free(src);
    free(snk);
    free(cap);
    return res;
}

/* min_cost_flow(g, demand, capacity|None, cost|None, want_flow) -> (cost, flow bytes|None) */
static PyObject *py_min_cost_flow(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go, *dobj, *co, *wo;
    int want;
    PARSE("OOOOi", &go, &dobj, &co, &wo, &want);
    GRAPH(g, go);
    double *demand = doubles_from_seq(dobj, g->n), *cap = NULL, *cost = NULL;
    if (!demand) return NULL;
    if (opt_doubles(co, g->m, &cap) < 0 || opt_doubles(wo, g->m, &cost) < 0) {
        free(demand);
        free(cap);
        free(cost);
        return NULL;
    }
    double total = 0, *flow = want ? alloc_doubles(g->m > 0 ? g->m : 1) : NULL;
    BEGIN_CALL
    status = net_min_cost_flow(g->h, demand, cap, cost, &total, flow, err, ERR_LEN);
    END_CALL
    free(demand);
    free(cap);
    free(cost);
    if (status) {
        free(flow);
        return engine_error(err);
    }
    PyObject *f = flow ? bytes_doubles(flow, g->m) : Py_NewRef(Py_None);
    return f ? Py_BuildValue("(dN)", total, f) : NULL;
}

/* gomory_hu(g, capacity|None) -> (parent list, weight bytes) */
static PyObject *py_gomory_hu(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go, *co;
    PARSE("OO", &go, &co);
    GRAPH(g, go);
    double *cap = NULL;
    if (opt_doubles(co, g->m, &cap) < 0) return NULL;
    int64_t *parent = alloc_int64(g->n > 0 ? g->n : 1);
    double *weight = alloc_doubles(g->n > 0 ? g->n : 1);
    if (!parent || !weight) {
        free(parent);
        free(weight);
        free(cap);
        return NULL;
    }
    BEGIN_CALL
    status = net_gomory_hu_tree(g->h, cap, parent, weight, err, ERR_LEN);
    END_CALL
    free(cap);
    if (status) {
        free(parent);
        free(weight);
        return engine_error(err);
    }
    return Py_BuildValue("(NN)", finish_int64(parent, g->n), bytes_doubles(weight, g->n));
}

/* local_connectivity(g, s, t, kind) -> int */
static PyObject *py_local_connectivity(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    long long s, t;
    int kind;
    PARSE("OLLi", &go, &s, &t, &kind);
    GRAPH(g, go);
    int64_t out = 0;
    BEGIN_CALL
    status = net_local_connectivity(g->h, s, t, kind, &out, err, ERR_LEN);
    END_CALL
    if (status) return engine_error(err);
    return PyLong_FromLongLong(out);
}

/* connectivity(g, kind, threads) -> int */
static PyObject *py_connectivity(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int kind, threads;
    PARSE("Oii", &go, &kind, &threads);
    GRAPH(g, go);
    int64_t out = 0;
    BEGIN_CALL
    status = net_connectivity(g->h, kind, threads, &out, err, ERR_LEN);
    END_CALL
    if (status) return engine_error(err);
    return PyLong_FromLongLong(out);
}

/* all_pairs_connectivity(g, threads) -> (matrix bytes, average) */
static PyObject *py_all_pairs_connectivity(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int threads;
    PARSE("Oi", &go, &threads);
    GRAPH(g, go);
    double *out = alloc_doubles(g->n * g->n > 0 ? g->n * g->n : 1), avg = 0;
    if (!out) return NULL;
    BEGIN_CALL
    status = net_all_pairs_node_connectivity(g->h, threads, out, &avg, err, ERR_LEN);
    END_CALL
    if (status) {
        free(out);
        return engine_error(err);
    }
    return Py_BuildValue("(Nd)", bytes_doubles(out, g->n * g->n), avg);
}

/* min_edge_cut(g, s, t, threads) -> (from list, to list) */
static PyObject *py_min_edge_cut(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    long long s, t;
    int threads;
    PARSE("OLLi", &go, &s, &t, &threads);
    GRAPH(g, go);
    int64_t *a = alloc_int64(g->m > 0 ? g->m : 1), *b = alloc_int64(g->m > 0 ? g->m : 1), size = 0;
    if (!a || !b) {
        free(a);
        free(b);
        return NULL;
    }
    BEGIN_CALL
    status = net_min_edge_cut(g->h, s, t, threads, a, b, &size, err, ERR_LEN);
    END_CALL
    if (status) {
        free(a);
        free(b);
        return engine_error(err);
    }
    return Py_BuildValue("(NN)", finish_int64(a, size), finish_int64(b, size));
}

/* min_node_cut(g, s, t, convention, threads) -> [nodes] */
static PyObject *py_min_node_cut(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    long long s, t;
    int convention, threads;
    PARSE("OLLii", &go, &s, &t, &convention, &threads);
    GRAPH(g, go);
    int64_t *nodes = alloc_int64(g->n > 0 ? g->n : 1), size = 0;
    if (!nodes) return NULL;
    BEGIN_CALL
    status = net_min_node_cut(g->h, s, t, convention, threads, nodes, &size, err, ERR_LEN);
    END_CALL
    return indices_result(status, err, nodes, size);
}

/* k_components(g, threads) -> [(nodes, k)] */
static PyObject *py_k_components(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int threads;
    PARSE("Oi", &go, &threads);
    GRAPH(g, go);
    ListsHandle *lh = NULL;
    BEGIN_CALL
    status = net_k_components(g->h, threads, &lh, err, ERR_LEN);
    END_CALL
    return lists_from(status, err, lh);
}

/* k_edge_components(g, k, threads) -> membership list */
static PyObject *py_k_edge_components(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    long long k;
    int threads;
    PARSE("OLi", &go, &k, &threads);
    GRAPH(g, go);
    int32_t *m = malloc(sizeof(int32_t) * (size_t)(g->n > 0 ? g->n : 1));
    if (!m) return PyErr_NoMemory();
    BEGIN_CALL
    status = net_k_edge_components(g->h, k, threads, m, err, ERR_LEN);
    END_CALL
    PyObject *res = status ? engine_error(err) : list_from_int32(m, g->n);
    free(m);
    return res;
}

/* local_bridges(g) -> (flags list, span bytes) */
static PyObject *py_local_bridges(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    PARSE("O", &go);
    GRAPH(g, go);
    int32_t *flags = malloc(sizeof(int32_t) * (size_t)(g->m > 0 ? g->m : 1));
    double *span = alloc_doubles(g->m > 0 ? g->m : 1);
    if (!flags || !span) {
        free(flags);
        free(span);
        return PyErr_NoMemory();
    }
    BEGIN_CALL
    status = net_local_bridges(g->h, flags, span, err, ERR_LEN);
    END_CALL
    PyObject *res = NULL;
    if (status)
        engine_error(err);
    else
        res = Py_BuildValue("(NN)", list_from_int32(flags, g->m), bytes_doubles(span, g->m));
    free(flags);
    if (status) free(span);
    return res;
}

/* condensation(g, order) -> (membership list, count) */
static PyObject *py_condensation(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int order;
    PARSE("Oi", &go, &order);
    GRAPH(g, go);
    int32_t *m = malloc(sizeof(int32_t) * (size_t)(g->n > 0 ? g->n : 1));
    int64_t count = 0;
    if (!m) return PyErr_NoMemory();
    BEGIN_CALL
    status = net_condensation(g->h, order, m, &count, err, ERR_LEN);
    END_CALL
    PyObject *res = status ? engine_error(err) : Py_BuildValue("(NL)", list_from_int32(m, g->n), (long long)count);
    free(m);
    return res;
}

/* attracting(g, order) -> membership list (-1 outside) */
static PyObject *py_attracting(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int order;
    PARSE("Oi", &go, &order);
    GRAPH(g, go);
    int32_t *m = malloc(sizeof(int32_t) * (size_t)(g->n > 0 ? g->n : 1));
    if (!m) return PyErr_NoMemory();
    BEGIN_CALL
    status = net_attracting_components(g->h, order, m, err, ERR_LEN);
    END_CALL
    PyObject *res = status ? engine_error(err) : list_from_int32(m, g->n);
    free(m);
    return res;
}

/* semiconnected(g) -> bool */
static PyObject *py_semiconnected(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    PARSE("O", &go);
    GRAPH(g, go);
    int32_t out = 0;
    BEGIN_CALL
    status = net_is_semiconnected(g->h, &out, err, ERR_LEN);
    END_CALL
    if (status) return engine_error(err);
    return PyBool_FromLong(out);
}

/* spanning_tree(g, weighted, maximum, algorithm) -> (edges, weight) */
static PyObject *py_spanning_tree(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int weighted, maximum, algorithm;
    PARSE("Oiii", &go, &weighted, &maximum, &algorithm);
    GRAPH(g, go);
    int64_t *edges = alloc_int64(g->n > 0 ? g->n : 1), size = 0;
    double weight = 0;
    if (!edges) return NULL;
    BEGIN_CALL
    status = net_spanning_tree(g->h, weighted, maximum, algorithm, edges, &size, &weight, err, ERR_LEN);
    END_CALL
    return int_list_double(status, err, edges, size, weight);
}

/* branching(g, weighted, kind, root) -> (edges, weight) */
static PyObject *py_branching(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int weighted, kind;
    long long root;
    PARSE("OiiL", &go, &weighted, &kind, &root);
    GRAPH(g, go);
    int64_t *edges = alloc_int64(g->n > 0 ? g->n : 1), size = 0;
    double weight = 0;
    if (!edges) return NULL;
    BEGIN_CALL
    status = net_branching(g->h, weighted, kind, root, edges, &size, &weight, err, ERR_LEN);
    END_CALL
    return int_list_double(status, err, edges, size, weight);
}

/* topological_sort(g, order) -> (nodes, generation) */
static PyObject *py_topological_sort(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int order;
    PARSE("Oi", &go, &order);
    GRAPH(g, go);
    int64_t *nodes = alloc_int64(g->n > 0 ? g->n : 1), *gen = alloc_int64(g->n > 0 ? g->n : 1);
    if (!nodes || !gen) {
        free(nodes);
        free(gen);
        return NULL;
    }
    BEGIN_CALL
    status = net_topological_sort(g->h, order, nodes, gen, err, ERR_LEN);
    END_CALL
    if (status) {
        free(nodes);
        free(gen);
        return engine_error(err);
    }
    return Py_BuildValue("(NN)", finish_int64(nodes, g->n), finish_int64(gen, g->n));
}

/* lex_topological_sort(g) -> nodes */
static PyObject *py_lex_topological_sort(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    PARSE("O", &go);
    GRAPH(g, go);
    int64_t *nodes = alloc_int64(g->n > 0 ? g->n : 1);
    if (!nodes) return NULL;
    BEGIN_CALL
    status = net_lexicographical_topological_sort(g->h, nodes, err, ERR_LEN);
    END_CALL
    return indices_result(status, err, nodes, g->n);
}

/* all_topological_sorts(g, order, limit) -> [(nodes, nan)] */
static PyObject *py_all_topological_sorts(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int order;
    long long limit;
    PARSE("OiL", &go, &order, &limit);
    GRAPH(g, go);
    ListsHandle *lh = NULL;
    BEGIN_CALL
    status = net_all_topological_sorts(g->h, order, limit, &lh, err, ERR_LEN);
    END_CALL
    return lists_from(status, err, lh);
}

/* is_dag(g) -> bool */
static PyObject *py_is_dag(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    PARSE("O", &go);
    GRAPH(g, go);
    int32_t out = 0;
    BEGIN_CALL
    status = net_is_dag(g->h, &out, err, ERR_LEN);
    END_CALL
    if (status) return engine_error(err);
    return PyBool_FromLong(out);
}

/* transitive_closure(g, reflexive) -> [(nodes, nan)] */
static PyObject *py_transitive_closure(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int reflexive;
    PARSE("Oi", &go, &reflexive);
    GRAPH(g, go);
    ListsHandle *lh = NULL;
    BEGIN_CALL
    status = net_transitive_closure(g->h, reflexive, &lh, err, ERR_LEN);
    END_CALL
    return lists_from(status, err, lh);
}

/* transitive_reduction(g) -> [(nodes, nan)] */
static PyObject *py_transitive_reduction(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    PARSE("O", &go);
    GRAPH(g, go);
    ListsHandle *lh = NULL;
    BEGIN_CALL
    status = net_transitive_reduction(g->h, &lh, err, ERR_LEN);
    END_CALL
    return lists_from(status, err, lh);
}

/* dag_longest_path(g, weighted, order) -> (path, length) */
static PyObject *py_dag_longest_path(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int weighted, order;
    PARSE("Oii", &go, &weighted, &order);
    GRAPH(g, go);
    int64_t *path = alloc_int64(g->n > 0 ? g->n : 1), size = 0;
    double length = 0;
    if (!path) return NULL;
    BEGIN_CALL
    status = net_dag_longest_path(g->h, weighted, order, path, &size, &length, err, ERR_LEN);
    END_CALL
    return int_list_double(status, err, path, size, length);
}

/* dominators(g, start) -> idom list (-1 unreachable) */
static PyObject *py_dominators(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    long long start;
    PARSE("OL", &go, &start);
    GRAPH(g, go);
    int64_t *out = alloc_int64(g->n > 0 ? g->n : 1);
    if (!out) return NULL;
    BEGIN_CALL
    status = net_immediate_dominators(g->h, start, out, err, ERR_LEN);
    END_CALL
    return indices_result(status, err, out, g->n);
}

/* dominance_frontiers(g, start) -> [(nodes, reachable)] */
static PyObject *py_dominance_frontiers(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    long long start;
    PARSE("OL", &go, &start);
    GRAPH(g, go);
    ListsHandle *lh = NULL;
    BEGIN_CALL
    status = net_dominance_frontiers(g->h, start, &lh, err, ERR_LEN);
    END_CALL
    return lists_from(status, err, lh);
}

/* lca(g, a, b, root, order) -> [lca index or -1] */
static PyObject *py_lca(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go, *ao, *bo;
    long long root;
    int order;
    PARSE("OOOLi", &go, &ao, &bo, &root, &order);
    GRAPH(g, go);
    int64_t *a = NULL, *b = NULL, ka = 0, kb = 0;
    if (req_int64(ao, &a, &ka) < 0) return NULL;
    if (req_int64(bo, &b, &kb) < 0 || ka != kb) {
        free(a);
        free(b);
        if (!PyErr_Occurred()) PyErr_SetString(PyExc_ValueError, "length mismatch");
        return NULL;
    }
    int64_t *out = alloc_int64(ka > 0 ? ka : 1);
    if (!out) {
        free(a);
        free(b);
        return NULL;
    }
    BEGIN_CALL
    status = net_lowest_common_ancestors(g->h, a, b, ka, root, order, out, err, ERR_LEN);
    END_CALL
    free(a);
    free(b);
    return indices_result(status, err, out, ka);
}

/* cycle_basis(g, order, root) -> [(nodes, nan)] */
static PyObject *py_cycle_basis(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int order;
    long long root;
    PARSE("OiL", &go, &order, &root);
    GRAPH(g, go);
    ListsHandle *lh = NULL;
    BEGIN_CALL
    status = net_cycle_basis(g->h, order, root, &lh, err, ERR_LEN);
    END_CALL
    return lists_from(status, err, lh);
}

/* minimum_cycle_basis(g, weighted, threads) -> [(nodes, weight)] */
static PyObject *py_minimum_cycle_basis(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int weighted, threads;
    PARSE("Oii", &go, &weighted, &threads);
    GRAPH(g, go);
    ListsHandle *lh = NULL;
    BEGIN_CALL
    status = net_minimum_cycle_basis(g->h, weighted, threads, &lh, err, ERR_LEN);
    END_CALL
    return lists_from(status, err, lh);
}

/* cycles(g, kind, length_bound, limit) -> [(nodes, nan)]; kind 0 simple, 1 chordless */
static PyObject *py_cycles(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int kind;
    long long bound, limit;
    PARSE("OiLL", &go, &kind, &bound, &limit);
    GRAPH(g, go);
    ListsHandle *lh = NULL;
    BEGIN_CALL
    status = kind == 0 ? net_simple_cycles(g->h, bound, limit, &lh, err, ERR_LEN)
                       : net_chordless_cycles(g->h, bound, limit, &lh, err, ERR_LEN);
    END_CALL
    return lists_from(status, err, lh);
}

/* girth(g) -> float */
static PyObject *py_girth(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    PARSE("O", &go);
    GRAPH(g, go);
    double out = 0;
    BEGIN_CALL
    status = net_girth(g->h, &out, err, ERR_LEN);
    END_CALL
    if (status) return engine_error(err);
    return PyFloat_FromDouble(out);
}

/* find_cycle(g, source, orientation, order) -> (edges, reversed) */
static PyObject *py_find_cycle(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    long long source;
    int orientation, order;
    PARSE("OLii", &go, &source, &orientation, &order);
    GRAPH(g, go);
    int64_t *edges = alloc_int64(g->m > 0 ? g->m : 1), size = 0;
    int32_t *rev = malloc(sizeof(int32_t) * (size_t)(g->m > 0 ? g->m : 1));
    if (!edges || !rev) {
        free(edges);
        free(rev);
        return PyErr_NoMemory();
    }
    BEGIN_CALL
    status = net_find_cycle(g->h, source, orientation, order, edges, rev, &size, err, ERR_LEN);
    END_CALL
    PyObject *res = NULL;
    if (status)
        engine_error(err);
    else
        res = Py_BuildValue("(NN)", finish_int64(edges, size), list_from_int32(rev, size));
    if (status) free(edges);
    free(rev);
    return res;
}

/* rich_club(g, threads) -> bytes of the unnormalized coefficients */
static PyObject *py_rich_club(PyObject *mod, PyObject *args) {
    (void)mod;
    PyObject *go;
    int threads;
    PARSE("Oi", &go, &threads);
    GRAPH(g, go);
    double *out = alloc_doubles(g->n > 0 ? g->n : 1);
    int64_t count = 0;
    if (!out) return NULL;
    BEGIN_CALL
    status = net_rich_club(g->h, 0, 0, 0, 0, threads, out, &count, err, ERR_LEN);
    END_CALL
    return done_bytes(status, err, out, count);
}

#define GRAPHOPS_METHODS                                                                         \
    {"max_flow", py_max_flow, METH_VARARGS, NULL},                                               \
    {"min_cut", py_min_cut, METH_VARARGS, NULL},                                                 \
    {"min_cost_flow", py_min_cost_flow, METH_VARARGS, NULL},                                     \
    {"gomory_hu", py_gomory_hu, METH_VARARGS, NULL},                                             \
    {"local_connectivity", py_local_connectivity, METH_VARARGS, NULL},                           \
    {"connectivity", py_connectivity, METH_VARARGS, NULL},                                       \
    {"all_pairs_connectivity", py_all_pairs_connectivity, METH_VARARGS, NULL},                   \
    {"min_edge_cut", py_min_edge_cut, METH_VARARGS, NULL},                                       \
    {"min_node_cut", py_min_node_cut, METH_VARARGS, NULL},                                       \
    {"k_edge_components", py_k_edge_components, METH_VARARGS, NULL},                             \
    {"k_components", py_k_components, METH_VARARGS, NULL},                                       \
    {"local_bridges", py_local_bridges, METH_VARARGS, NULL},                                     \
    {"condensation", py_condensation, METH_VARARGS, NULL},                                       \
    {"attracting", py_attracting, METH_VARARGS, NULL},                                           \
    {"semiconnected", py_semiconnected, METH_VARARGS, NULL},                                     \
    {"spanning_tree", py_spanning_tree, METH_VARARGS, NULL},                                     \
    {"branching", py_branching, METH_VARARGS, NULL},                                             \
    {"topological_sort", py_topological_sort, METH_VARARGS, NULL},                               \
    {"lex_topological_sort", py_lex_topological_sort, METH_VARARGS, NULL},                       \
    {"all_topological_sorts", py_all_topological_sorts, METH_VARARGS, NULL},                     \
    {"is_dag", py_is_dag, METH_VARARGS, NULL},                                                   \
    {"transitive_closure", py_transitive_closure, METH_VARARGS, NULL},                           \
    {"transitive_reduction", py_transitive_reduction, METH_VARARGS, NULL},                       \
    {"dag_longest_path", py_dag_longest_path, METH_VARARGS, NULL},                               \
    {"dominators", py_dominators, METH_VARARGS, NULL},                                           \
    {"dominance_frontiers", py_dominance_frontiers, METH_VARARGS, NULL},                         \
    {"lca", py_lca, METH_VARARGS, NULL},                                                         \
    {"cycle_basis", py_cycle_basis, METH_VARARGS, NULL},                                         \
    {"minimum_cycle_basis", py_minimum_cycle_basis, METH_VARARGS, NULL},                         \
    {"cycles", py_cycles, METH_VARARGS, NULL},                                                   \
    {"girth", py_girth, METH_VARARGS, NULL},                                                     \
    {"find_cycle", py_find_cycle, METH_VARARGS, NULL},                                           \
    {"rich_club", py_rich_club, METH_VARARGS, NULL},
