#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

extern void C_huber(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void C_quantile(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void C_squared(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void C_huber_l2(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void C_quantile_l2(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void C_squared_l2(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"C_huber",      (DL_FUNC) &C_huber,       22},
    {"C_quantile",   (DL_FUNC) &C_quantile,    22},
    {"C_squared",    (DL_FUNC) &C_squared,     21},
    {"C_huber_l2",   (DL_FUNC) &C_huber_l2,    17},
    {"C_quantile_l2",(DL_FUNC) &C_quantile_l2, 17},
    {"C_squared_l2", (DL_FUNC) &C_squared_l2,  16},
    {NULL, NULL, 0}
};

void R_init_PRIMsrc(DllInfo *dll) {
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
