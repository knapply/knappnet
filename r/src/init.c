#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP C_graph_new(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_graph_alive(SEXP);
extern SEXP C_degree(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_betweenness(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_triad_census(SEXP, SEXP);
extern SEXP C_read_network(SEXP, SEXP, SEXP);
extern SEXP C_write_network(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_coreness(SEXP, SEXP);
extern SEXP C_max_cliques(SEXP, SEXP, SEXP);
extern SEXP C_clique_number(SEXP, SEXP);
extern SEXP C_modularity(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_louvain(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_structural_holes(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_brokerage(SEXP, SEXP, SEXP);
extern SEXP C_structural_equivalence(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_concor(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_qap_cor(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_mrqap(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_graph_simplify(SEXP, SEXP, SEXP);
extern SEXP C_density(SEXP, SEXP);
extern SEXP C_components(SEXP, SEXP);
extern SEXP C_fragmentation(SEXP, SEXP);
extern SEXP C_centralization(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_closeness(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_eigenvector(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef call_methods[] = {
    {"C_graph_new", (DL_FUNC)&C_graph_new, 5},
    {"C_graph_alive", (DL_FUNC)&C_graph_alive, 1},
    {"C_degree", (DL_FUNC)&C_degree, 4},
    {"C_betweenness", (DL_FUNC)&C_betweenness, 5},
    {"C_triad_census", (DL_FUNC)&C_triad_census, 2},
    {"C_read_network", (DL_FUNC)&C_read_network, 3},
    {"C_write_network", (DL_FUNC)&C_write_network, 7},
    {"C_coreness", (DL_FUNC)&C_coreness, 2},
    {"C_max_cliques", (DL_FUNC)&C_max_cliques, 3},
    {"C_clique_number", (DL_FUNC)&C_clique_number, 2},
    {"C_modularity", (DL_FUNC)&C_modularity, 4},
    {"C_louvain", (DL_FUNC)&C_louvain, 4},
    {"C_structural_holes", (DL_FUNC)&C_structural_holes, 4},
    {"C_brokerage", (DL_FUNC)&C_brokerage, 3},
    {"C_structural_equivalence", (DL_FUNC)&C_structural_equivalence, 4},
    {"C_concor", (DL_FUNC)&C_concor, 6},
    {"C_qap_cor", (DL_FUNC)&C_qap_cor, 10},
    {"C_mrqap", (DL_FUNC)&C_mrqap, 15},
    {"C_graph_simplify", (DL_FUNC)&C_graph_simplify, 3},
    {"C_density", (DL_FUNC)&C_density, 2},
    {"C_components", (DL_FUNC)&C_components, 2},
    {"C_fragmentation", (DL_FUNC)&C_fragmentation, 2},
    {"C_centralization", (DL_FUNC)&C_centralization, 7},
    {"C_closeness", (DL_FUNC)&C_closeness, 8},
    {"C_eigenvector", (DL_FUNC)&C_eigenvector, 8},
    {NULL, NULL, 0}};

/* The symbol name must match the package name; tools/rename.py updates it. */
void R_init_knappnet(DllInfo *dll) {
  R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
