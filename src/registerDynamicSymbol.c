#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void TOY(double *, double *, int *, int *);
extern void SBDE(double *, double *, int *, double *, double *,
                 int *, double *, double *, double *, double *,
                 int *, int *, double *, int *, double *,
                 double *, double *, int *);
extern void DEV(double *, double *, int *, double *, double *,
                int *, double *, double *, double *, double *,
                double *, double *, int *);
extern void PRED(double *, double *, double *, int *, double *,
                 double *, double *, int *);
extern void QUANT(double *, double *, double *, int *, double *,
                  double *, double *, int *);

static const R_CMethodDef CEntries[] = {
    {"TOY", (DL_FUNC) &TOY, 4},
    {"SBDE", (DL_FUNC) &SBDE, 18},
    {"DEV",  (DL_FUNC) &DEV,  13},
    {"PRED", (DL_FUNC) &PRED,  8},
    {"QUANT", (DL_FUNC) &QUANT, 8},
    {NULL, NULL, 0}
};

void R_init_sbde(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
