#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _mobsim_rThomas_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _mobsim_sSAC1_C(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
   {"_mobsim_rThomas_rcpp", (DL_FUNC) &_mobsim_rThomas_rcpp, 8},
   {"_mobsim_sSAC1_C",      (DL_FUNC) &_mobsim_sSAC1_C,      3},
   {NULL, NULL, 0}
};

void R_init_mobsim(DllInfo *dll)
{
   R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
   R_useDynamicSymbols(dll, FALSE);
}
