/* RUMP's Poisson chi-square objective, called directly.
 *
 * `EvalChiPoisson` (genplot/curfit.c:557) is the function PERT minimises. It
 * takes an NLS_DATA struct, so this wrapper builds a minimal one -- only the
 * data/yfit/valid/npt/dof fields are read on that path.
 */

#include <stdlib.h>
#include <string.h>

#include "mytypes.h"
#include "curfit.h"

/* Evaluate the objective over `n` channels.
 *
 * `valid` may be NULL to use every channel. Returns the number of channels
 * with non-positive theory (RUMP's "Poisson statistics invalid!" count);
 * `out_residuals` receives the signed per-channel residuals and `out_reduced`
 * the reduced chi-square RUMP reports.
 */
int OracleChiPoisson(const double *data, const double *theory, const int *valid,
                     int n, int nvars, double *out_residuals, double *out_reduced) {
	NLS_DATA nls;
	REAL *d, *y, *chi;
	BOOL *v = NULL;
	int i, nbad, used = 0;

	d   = calloc(n, sizeof(REAL));
	y   = calloc(n, sizeof(REAL));
	chi = calloc(n, sizeof(REAL));
	for (i = 0; i < n; i++) { d[i] = (REAL) data[i]; y[i] = (REAL) theory[i]; }
	if (valid != NULL) {
		v = calloc(n, sizeof(BOOL));
		for (i = 0; i < n; i++) { v[i] = (BOOL) valid[i]; if (valid[i]) used++; }
	} else {
		used = n;
	}

	memset(&nls, 0, sizeof(nls));
	nls.data     = d;
	nls.yfit     = y;
	nls.valid    = v;
	nls.errorbar = NULL;
	nls.outchi   = chi;
	nls.npt      = n;
	nls.dof      = (used - nvars > 0) ? used - nvars : 1;

	nbad = EvalChiPoisson(&nls);

	for (i = 0; i < n; i++) out_residuals[i] = (double) chi[i];
	*out_reduced = (double) nls.chisqr;

	free(d); free(y); free(chi); if (v) free(v);
	return nbad;
}
