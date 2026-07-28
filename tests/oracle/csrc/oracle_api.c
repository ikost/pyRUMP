/* Flat C API over the RUMP physics routines, for cffi.
 *
 * Everything here is a thin, allocation-free wrapper: no RUMP struct crosses the
 * boundary, so the Python side needs no knowledge of REAL, STOPPING_TABLE, etc.
 * That keeps the binding stable if the C headers shift.
 */

#include <string.h>
#include <math.h>

/* NOTE: the RUMP headers have no include guards, so pulling in rumpdata.h or
 * rumpproto.h directly double-defines every enum. rump.h already includes both
 * (rump.h:17-21); stopping.h is the only one it leaves out. */
#include "rump.h"
#include "stopping.h"
#include "sigma.h"

extern void OracleSetDataDir(const char *dir);

/* stop_type (the selector between polynomials in E and in sqrt(E)) is already
 * declared as STOPPING_TYPE by rumpdata.h:132. */

/* ------------------------------------------------------------------ setup --- */

/* Load the three startup tables. Returns 0 on success, or a bitmask of failures.
 * Mirrors the load order in config.c:187-191, which matters: Kalbitzer entries
 * take priority over Ziegler at lookup time.
 */
int OracleInit(const char *data_dir) {
	int status = 0;
	OracleSetDataDir(data_dir);
	if (! RbsLoadAtomicData("atom4.dat"))     status |= 1;
	if (! RbsLoadZieglerData("pscoef.dat"))   status |= 2;
	if (! RbsLoadKalbitzerData("newstop.kal")) status |= 4;
	return status;
}

int  OracleGetStopType(void)     { return (int) stop_type; }
void OracleSetStopType(int type) { stop_type = (STOPPING_TYPE) type; }

/* ------------------------------------------------- raw Ziegler evaluation --- */

/* Direct ZBL85 stopping, bypassing RUMP's polynomial refit.
 *
 * units follows ziegler.c's convention:
 *   1 = eV/(1e15 atoms/cm^2)   2 = keV/micron
 *   3 = eV/Angstrom            4 = MeV/(mg/cm^2)
 *
 * Returns electronic and nuclear stopping separately, as the C does.
 */
void OracleZStop(int z1, double m1, int z2, double energy_keV,
                 int units, double *se, double *sn) {
	REAL r_se = 0.0f, r_sn = 0.0f;
	zstop(z1, (REAL) m1, z2, (REAL) energy_keV, &r_se, &r_sn, units);
	*se = (double) r_se;
	*sn = (double) r_sn;
}

/* ------------------------------------------- fitted polynomial evaluation --- */

/* Build (or reuse) the stopping table RUMP would use for this beam, then return
 * the degree-5 polynomial coefficients for one target element.
 *
 * This is the M3 crux: RUMP never evaluates Ziegler during a simulation, it
 * evaluates *these* coefficients (stopping.c:418-568).
 *
 * Returns 1 on success. `coef` must have room for NDEG (6) doubles.
 * `e_scale` receives the He-3/D energy-scaling factor applied by RbsStpfind.
 */
int OracleStoppingCoefficients(int z_beam, double m_beam, double e_beam_MeV,
                               int z_target, double *coef, double *e_scale) {
	REAL scale = 1.0f;
	STOPPING_TABLE *table;
	STOPPING_POWER *stop;
	int i;

	table = RbsStpfind(z_beam, (REAL) m_beam, &scale, (REAL) e_beam_MeV);
	if (table == NULL) return 0;

	stop = RbsLookupStop(table, z_target);
	if (stop == NULL) return 0;

	for (i = 0; i < NDEG; i++) coef[i] = (double) stop->p[i];
	*e_scale = (double) scale;
	return 1;
}

/* Evaluate the fitted polynomial the way the simulation kernel does, i.e. via
 * S_XFORM + S_POWER. Also returns the first and second derivatives, which the
 * 3rd-order energy-loss expansion needs (creatr.c:1546-1554).
 *
 * Energy is in keV; result is eV/(1e15 atoms/cm^2).
 */
int OracleStoppingEvaluate(int z_beam, double m_beam, double e_beam_MeV,
                           int z_target, const double *energies, int n,
                           double *out_s, double *out_ds, double *out_dds) {
	REAL scale = 1.0f;
	STOPPING_TABLE *table;
	STOPPING_POWER *stop;
	int i;

	table = RbsStpfind(z_beam, (REAL) m_beam, &scale, (REAL) e_beam_MeV);
	if (table == NULL) return 0;
	stop = RbsLookupStop(table, z_target);
	if (stop == NULL) return 0;

	for (i = 0; i < n; i++) {
		double x = (double) S_XFORM(energies[i] * scale);
		if (out_s)   out_s[i]   = (double) S_POWER(stop, x);
		if (out_ds)  out_ds[i]  = (double) (DS_POWER(stop, x) * scale);
		if (out_dds) out_dds[i] = (double) (DDS_POWER(stop, x) * scale * scale);
	}
	return 1;
}

/* Report the fit window RUMP derives from the beam energy (stopping.c:308-321).
 *
 * For STOP_SQRT: emin = 0.04*E, emax = 1.15*E, cutoff = 0.03*E, all in MeV.
 * The fitted polynomial is only meaningful inside [emin, emax]; outside it the
 * degree-5 fit diverges quickly. Any port must reproduce these bounds, since
 * they depend on the beam energy and therefore change the depth scale.
 */
int OracleStoppingRange(int z_beam, double m_beam, double e_beam_MeV,
                        double *emin, double *emax, double *cutoff, int *type) {
	REAL scale = 1.0f;
	STOPPING_TABLE *table = RbsStpfind(z_beam, (REAL) m_beam, &scale, (REAL) e_beam_MeV);
	if (table == NULL) return 0;
	*emin   = (double) table->emin;
	*emax   = (double) table->emax;
	*cutoff = (double) table->cutoff;
	*type   = (int)    table->type;
	return 1;
}

/* --------------------------------------------------------------- elements --- */

int    OracleNumElements(void)          { return NumElements; }
double OracleElementMass(int z)         { return (double) atom[z-1].mass; }
double OracleElementDensity(int z)      { return (double) atom[z-1].dense; }
double OracleRealMass(int z, int iso)   { return (double) RbsGetRealMass(z, iso); }

void OracleElementSymbol(int z, char *out) { strcpy(out, atom[z-1].name); }

/* Ziegler block, for cross-checking the pscoef.dat parser. */
void OracleZieglerParams(int z, double *out) {
	out[0] = (double) atom[z-1].zmm1;
	out[1] = (double) atom[z-1].zm1;
	out[2] = (double) atom[z-1].zm2;
	out[3] = (double) atom[z-1].zrho;
	out[4] = (double) atom[z-1].zatrho;
	out[5] = (double) atom[z-1].zvferm;
	out[6] = (double) atom[z-1].zlfctr;
}

/* ------------------------------------------------------ deterministic tests --- */

extern STOPPING_TABLE *stop_tables;

/* Flush the session cache of fitted stopping tables.
 *
 * RbsStpfind keeps every table it has built and reuses one whenever the new
 * beam energy merely fits inside its window (stopping.c:274-279). That makes
 * results order-dependent, which is fine for an interactive session but useless
 * for tests. Calling this before each comparison forces a fresh fit.
 */
void OracleResetStoppingTables(void) {
	STOPPING_TABLE *table = stop_tables, *next;
	while (table != NULL) {
		next = table->next;
		free(table);
		table = next;
	}
	stop_tables = NULL;
}

/* ---------------------------------------------------------- cross sections --- */

/* Evaluate a cross-section over an energy array.
 *
 * `recoil` selects SetupSigmaRecoil (ERD) over SetupSigmaScatter (RBS).
 * `phi_deg` is the TRUE scattering angle, not RUMP's supplement -- creatr.c
 * converts with `sp.phi = 180 - samm->phi` before calling in (creatr.c:1629).
 *
 * Result is barns/sr.
 */
int OracleCrossSection(int recoil, int z1, double m1, int z2, double m2,
                       double phi_deg, double kev_max,
                       const double *energies, int n, double *out) {
	SP sp;
	int i;

	memset(&sp, 0, sizeof(sp));
	sp.phi     = phi_deg;
	sp.sinph   = sin(phi_deg * 3.14159265358979323846 / 180.0);
	sp.cosph   = cos(phi_deg * 3.14159265358979323846 / 180.0);
	sp.z1      = z1;
	sp.z2      = z2;
	sp.m1      = m1;
	sp.m2      = m2;
	sp.kev_max = kev_max;

	if (! (recoil ? SetupSigmaRecoil(&sp) : SetupSigmaScatter(&sp))) return 0;
	for (i = 0; i < n; i++) out[i] = (*sp.calc)(energies[i], &sp);
	return 1;
}

/* Expose the setup constants themselves, so a port can be checked term by term
 * rather than only on the final number. */
int OracleSigmaConstants(int recoil, int z1, double m1, int z2, double m2,
                         double phi_deg, double kev_max, double *out) {
	SP sp;

	memset(&sp, 0, sizeof(sp));
	sp.phi     = phi_deg;
	sp.sinph   = sin(phi_deg * 3.14159265358979323846 / 180.0);
	sp.cosph   = cos(phi_deg * 3.14159265358979323846 / 180.0);
	sp.z1      = z1;
	sp.z2      = z2;
	sp.m1      = m1;
	sp.m2      = m2;
	sp.kev_max = kev_max;

	if (! (recoil ? SetupSigmaRecoil(&sp) : SetupSigmaScatter(&sp))) return 0;
	out[0] = sp.csigma;   /* 1/E^2 coefficient */
	out[1] = sp.csig_0;   /* constant term     */
	out[2] = sp.csig_f;   /* screening rolloff */
	return 1;
}
