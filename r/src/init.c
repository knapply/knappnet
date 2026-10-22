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

static const R_CallMethodDef call_methods[] = {
    {"C_graph_new", (DL_FUNC)&C_graph_new, 5},
    {"C_graph_alive", (DL_FUNC)&C_graph_alive, 1},
    {"C_degree", (DL_FUNC)&C_degree, 4},
    {"C_betweenness", (DL_FUNC)&C_betweenness, 5},
    {"C_triad_census", (DL_FUNC)&C_triad_census, 2},
    {"C_read_network", (DL_FUNC)&C_read_network, 3},
    {"C_write_network", (DL_FUNC)&C_write_network, 7},
    {NULL, NULL, 0}};

/* The symbol name must match the package name; tools/rename.py updates it. */
void R_init_knappnet(DllInfo *dll) {
  R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
