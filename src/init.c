#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP dwiener(SEXP q, SEXP alpha, SEXP tau, SEXP beta, SEXP delta, SEXP give_log);
extern SEXP pwiener(SEXP q, SEXP alpha, SEXP tau, SEXP beta, SEXP delta);
extern SEXP qwiener(SEXP p, SEXP alpha, SEXP tau, SEXP beta, SEXP delta);
extern SEXP rwiener(SEXP alpha, SEXP tau, SEXP beta, SEXP delta);
extern SEXP pwiener_full(SEXP q, SEXP alpha, SEXP tau, SEXP beta, SEXP delta);
extern SEXP qwiener_full(SEXP p, SEXP alpha, SEXP tau, SEXP beta, SEXP delta);

static const R_CallMethodDef CallMethods[] = {
    {"dwiener", (DL_FUNC) &dwiener, 6},
    {"pwiener", (DL_FUNC) &pwiener, 5},
    {"qwiener", (DL_FUNC) &qwiener, 5},
    {"rwiener", (DL_FUNC) &rwiener, 4},
    {"pwiener_full", (DL_FUNC) &pwiener_full, 5},
    {"qwiener_full", (DL_FUNC) &qwiener_full, 5},
    {NULL, NULL, 0}
};

void R_init_RWiener(DllInfo *info)
{
    R_registerRoutines(info, NULL, CallMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}
