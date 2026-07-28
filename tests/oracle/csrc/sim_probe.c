/* Brick-level access to RUMP's simulation engine.
 *
 * `creatr.c` is the engine. Rather than patch it (which would make the oracle a
 * modified RUMP), we exploit the fact that its output stage is a **function
 * pointer**: `SimFillSpectrum` (sample.h) is called once per brick. Pointing it
 * at a capture routine records the exact (efront, eback, hfront, hback, qqq,
 * sigf, sigb) tuples the engine produces, with the C untouched.
 *
 * This file supplies the small amount of RUMP host state creatr.c expects --
 * the buffer array, spectrum allocation, and a few helpers -- reimplemented
 * where they are trivial and copied faithfully where the numbers matter
 * (`RbsNormK`, `SimThickConvert`).
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "rump.h"
#include "sample.h"
#include "stopping.h"

/* ------------------------------------------------------------ host state --- */

SPECTRUM **RbsBuffers = NULL;
SPECTRUM  *RbsActiveBuf = NULL;
SAMPLE    *SimDefaultSample = NULL;

/* IMPROVED = inverse-density (hard-ball packing) averaging, the post-1/97
 * default; COMPATIBLE restores the older direct density average (creatr.c:607). */
DENSITYCALC RbsDensityCalc  = IMPROVED;
int         RbsRecoilZLimit = 1;              /* config.c:108 */

/* creatr.c reaches the simulation state through `Rmp`, which rumpdata.h:201
 * defines as an alias for RbsDataBlock. Only `autsim` is read on our path. */
static RMPTYPE probe_rmp;
RMPTYPE *RbsDataBlock = &probe_rmp;

void (*SimPileup)(SPECTRUM *buf) = NULL;
void (*SimInitFillSpectrum)(SPECTRUM *buf) = NULL;
void (*SimTermFillSpectrum)(SPECTRUM *buf) = NULL;
void (*SimFillSpectrum)(int z, REAL mass, REAL efront, REAL eback,
                        REAL hfront, REAL hback, REAL qqq,
                        REAL sigf, REAL sigb) = NULL;

/* ------------------------------------------------- spectrum housekeeping --- */

SPECTRUM *RbsAllocateSpectrum(SPECTRUM *proto, int channels) {
	SPECTRUM *buf = calloc(1, sizeof(SPECTRUM));
	if (proto != NULL) memcpy(buf, proto, sizeof(SPECTRUM));
	buf->counts  = calloc(channels > 0 ? channels : 1, sizeof(REAL));
	buf->nptmax  = channels;
	buf->npt     = 0;
	buf->nspectra = 1;
	return buf;
}

SPECTRUM *RbsResizeSpectrum(SPECTRUM *buf, int size) {
	if (size <= buf->nptmax) return buf;
	buf->counts = realloc(buf->counts, size * sizeof(REAL));
	memset(buf->counts + buf->nptmax, 0, (size - buf->nptmax) * sizeof(REAL));
	buf->nptmax = size;
	return buf;
}

int RbsFreeSpectrum(SPECTRUM *buf) {
	if (buf == NULL) return 0;
	free(buf->counts);
	free(buf);
	return 0;
}

/* Copies parameters but keeps the destination's own counts array, which is what
 * creatr.c relies on at :283 (`RbsCopySpectrum(ALTBUF, ibuf)`). */
SPECTRUM *RbsCopySpectrum(SPECTRUM *dest, SPECTRUM *source) {
	REAL *counts = dest->counts;
	int   nptmax = dest->nptmax;
	memcpy(dest, source, sizeof(SPECTRUM));
	dest->counts = counts;
	dest->nptmax = nptmax;
	return dest;
}

/* bmanip.c:880 -- verbatim; this sets the absolute yield scale. */
REAL RbsNormK(SPECTRUM *ibf) {
	REAL x = ibf->omega * ibf->q * ibf->kevch;
	if (ibf->cbeam != 0 && x != 0.0 && ibf->corr != 0.0) {
		x = ibf->cbeam * ibf->corr / x;
	} else {
		x = 1.0f;
	}
	return x;
}

/* sim2.c:2349 -- verbatim. */
REAL SimThickConvert(UNITS *units, REAL density, REAL sum, REAL *pdensity) {
	REAL thick_to_cm2;
	switch (units->type) {
		case ANGSTROMS: thick_to_cm2 = density * units->density; break;
		case ATOMIC:    thick_to_cm2 = 1;                        break;
		case MOLECULAR: thick_to_cm2 = sum;                      break;
		case ABSOLUTE:
			thick_to_cm2 = units->density;
			if (pdensity != NULL) *pdensity = units->density;
			break;
		default:        thick_to_cm2 = 1;
	}
	return thick_to_cm2;
}

void SimSetIdent(void) { /* cosmetic only: names the theory buffer */ }



/* --------------------------------------------------------- brick capture --- */

#define MAX_BRICKS 200000

typedef struct {
	int   z;
	double mass, efront, eback, hfront, hback, qqq, sigf, sigb;
} BRICK;

static BRICK  brick_store[MAX_BRICKS];
static int    brick_count = 0;
static int    brick_overflow = 0;

static void CaptureBrick(int z, REAL mass, REAL efront, REAL eback,
                         REAL hfront, REAL hback, REAL qqq,
                         REAL sigf, REAL sigb) {
	BRICK *b;
	if (brick_count >= MAX_BRICKS) { brick_overflow = 1; return; }
	b = &brick_store[brick_count++];
	b->z = z;  b->mass = mass;
	b->efront = efront;  b->eback = eback;
	b->hfront = hfront;  b->hback = hback;
	b->qqq = qqq;  b->sigf = sigf;  b->sigb = sigb;
}

int OracleBrickCount(void)    { return brick_count; }
int OracleBrickOverflow(void) { return brick_overflow; }

/* Copy captured bricks out as a flat (n, 9) double array. */
void OracleBricks(double *out, int n) {
	int i;
	if (n > brick_count) n = brick_count;
	for (i = 0; i < n; i++) {
		BRICK *b = &brick_store[i];
		out[i*9 + 0] = b->z;
		out[i*9 + 1] = b->mass;
		out[i*9 + 2] = b->efront;
		out[i*9 + 3] = b->eback;
		out[i*9 + 4] = b->hfront;
		out[i*9 + 5] = b->hback;
		out[i*9 + 6] = b->qqq;
		out[i*9 + 7] = b->sigf;
		out[i*9 + 8] = b->sigb;
	}
}

/* ------------------------------------------------------- sample assembly --- */

static UNITS  unit_atomic = {"/CM2", ATOMIC, 1.0f};
static SAMPLE probe_sample;
static LAYER  probe_layers[64];

/* Build a sample of uniform layers, thicknesses in 1e15 atoms/cm^2.
 *
 * `composition` is row-major (n_layer, n_element) and need not be normalised --
 * RUMP normalises internally.
 */
int OracleSetSample(int n_layer, const double *thickness,
                    int n_element, const int *element_z,
                    const double *composition,
                    const int *sublayers,
                    double maxpth, double straggle, double multiple) {
	int i, j;
	if (n_layer < 1 || n_layer > 64 || n_element < 1 || n_element > MAXEL) return 0;

	memset(&probe_sample, 0, sizeof(probe_sample));
	memset(probe_layers, 0, sizeof(probe_layers));

	for (i = 0; i < n_layer; i++) {
		LAYER *layer = &probe_layers[i];
		layer->previous = (i == 0)          ? NULL : &probe_layers[i-1];
		layer->next     = (i == n_layer-1)  ? NULL : &probe_layers[i+1];
		layer->thick.magn  = (REAL) thickness[i];
		layer->thick.units = &unit_atomic;
		layer->num_sublayers = sublayers ? sublayers[i] : 0;
		layer->thisub.magn   = 0;
		layer->thisub.units  = &unit_atomic;
		layer->eqn = NULL;
		for (j = 0; j < n_element; j++)
			layer->matrix[j] = (REAL) composition[i*n_element + j];
	}

	probe_sample.first   = &probe_layers[0];
	probe_sample.layer   = &probe_layers[0];
	probe_sample.g_first = NULL;
	probe_sample.g_layer = NULL;
	probe_sample.nel     = n_element;
	for (j = 0; j < n_element; j++) {
		probe_sample.z2[j]    = element_z[j];
		probe_sample.nukem[j] = 0;              /* natural abundance */
	}
	probe_sample.absorber_layers = 0;
	probe_sample.maxpth   = (REAL) maxpth;
	probe_sample.straggle = (REAL) straggle;
	probe_sample.multiple = (REAL) multiple;
	probe_sample.noise    = 0;
	probe_sample.background = NULL;

	SimDefaultSample = &probe_sample;
	return 1;
}

/* Configure the beam/detector, i.e. the active data buffer the engine copies. */
int OracleSetBeam(double e0_MeV, int zbeam, double mbeam, int cbeam,
                  double q_uC, double current_nA,
                  double kevch, double kev0, double first, int npt,
                  double fwhm, double tau, int geom,
                  double phi, double theta, double psi,
                  double omega, double corr) {
	if (RbsBuffers == NULL) {
		RbsBuffers = calloc(4, sizeof(SPECTRUM *));
		RbsBuffers[0] = RbsAllocateSpectrum(NULL, CMAX);   /* ALTBUF */
		RbsBuffers[1] = RbsAllocateSpectrum(NULL, CMAX);   /* MAINBUF */
		RbsActiveBuf  = RbsBuffers[1];
	}
	RbsActiveBuf->e0      = (REAL) e0_MeV;
	RbsActiveBuf->zbeam   = zbeam;
	RbsActiveBuf->mbeam   = (REAL) mbeam;
	RbsActiveBuf->cbeam   = cbeam;
	RbsActiveBuf->q       = (REAL) q_uC;
	RbsActiveBuf->current = (REAL) current_nA;
	RbsActiveBuf->kevch   = (REAL) kevch;
	RbsActiveBuf->kev0    = (REAL) kev0;
	RbsActiveBuf->first   = (REAL) first;
	RbsActiveBuf->npt     = npt;
	RbsActiveBuf->fwhm    = (REAL) fwhm;
	RbsActiveBuf->tau     = (REAL) tau;
	RbsActiveBuf->geom    = (GEOMETRY_TYPE) geom;
	RbsActiveBuf->phi     = (REAL) phi;
	RbsActiveBuf->theta   = (REAL) theta;
	RbsActiveBuf->psi     = (REAL) psi;
	RbsActiveBuf->omega   = (REAL) omega;
	RbsActiveBuf->corr    = (REAL) corr;
	RbsActiveBuf->type    = RBS;
	return 1;
}

/* Run the engine, capturing bricks instead of filling a spectrum.
 *
 * `capture_only` swaps the fill routine for our recorder. With it disabled the
 * real SimAnlyz runs and the resulting spectrum can be read with OracleSpectrum.
 */
int OracleSimulate(int capture_only) {
	brick_count = 0;
	brick_overflow = 0;
	SimFillSpectrum = capture_only ? CaptureBrick : SimAnlyz;
	SimPileup = NULL;
	if (SimDefaultSample == NULL || RbsBuffers == NULL) return 0;
	return SimCreateDetails(SimDefaultSample, -1, -1) != NULL;
}

int OracleSpectrum(double *out, int n) {
	SPECTRUM *alt = RbsBuffers[0];
	int i;
	if (n > alt->npt) n = alt->npt;
	for (i = 0; i < n; i++) out[i] = (double) alt->counts[i];
	return alt->npt;
}

/* ------------------------------------------------------- depth equations --- */

/* Mirror of sim2.c:101's eqlist, restricted to the fields creatr.c reads:
 * the type, the parameter count, the dose units, and -- importantly --
 * rcmd_sublayers, which overrides maxpath whenever an equation is present
 * (creatr.c:700). */
static UNITS probe_angstroms = {"A",    ANGSTROMS, 1.0};
static UNITS probe_atoms_cm2 = {"/cm2", ATOMIC,    1.0};

static SIMEQN probe_eqns[] = {
	{"None",           2, EQ_NONE,      0, {-1,-1,-1,-1,-1}, NULL,             -1, ""},
	{"Constant",       2, EQ_CONST,     1, { 0,-1,-1,-1,-1}, NULL,              5, ""},
	{"Linear",         1, EQ_LINEAR,    2, {12,13,-1,-1,-1}, NULL,             10, ""},
	{"ERFC",           2, EQ_ERFC,      3, { 1, 2, 3,-1,-1}, NULL,             10, ""},
	{"Exponential",    2, EQ_EXP,       3, { 1, 2, 4,-1,-1}, NULL,             20, ""},
	{"Semi-infinite",  1, EQ_SEMI_INF,  4, { 5, 2, 3, 6,-1}, NULL,             20, ""},
	{"Thinfilm",       4, EQ_THINFILM,  3, { 7, 2, 3,-1,-1}, &probe_angstroms, 30, ""},
	{"BuriedThinFilm", 3, EQ_BURIED,    4, { 7, 2, 3, 6,-1}, &probe_angstroms, 30, ""},
	{"ThickFilm",      4, EQ_THICFILM,  5, { 9,10, 2, 3,11}, NULL,             20, ""},
	{"Timedependent",  2, EQ_TIMEDEPE,  4, { 1, 2, 3, 8,-1}, NULL,             20, ""},
	{"Gaussian",       1, EQ_GAUSS,     3, {14,15,16,-1,-1}, &probe_atoms_cm2, 20, ""},
	{"Edgeworth",      4, EQ_EDGEWORTH, 5, {17,18,19,20,21}, &probe_atoms_cm2, 30, ""},
};
#define N_PROBE_EQNS ((int)(sizeof(probe_eqns)/sizeof(probe_eqns[0])))

/* Attach an equation to one layer. `eqn_index` indexes probe_eqns; -1 clears.
 * `species` is the n_element composition the equation blends toward. */
int OracleSetLayerEquation(int layer_index, int eqn_index,
                           const double *params, int n_param,
                           const double *species, int n_element) {
	LAYER *layer;
	int j;

	if (layer_index < 0 || layer_index >= 64) return 0;
	layer = &probe_layers[layer_index];

	if (eqn_index < 0) { layer->eqn = NULL; return 1; }
	if (eqn_index >= N_PROBE_EQNS) return 0;

	layer->eqn = &probe_eqns[eqn_index];
	layer->eqn_units = probe_eqns[eqn_index].dflt_units;
	for (j = 0; j < MAXPAR; j++)
		layer->par[j] = (j < n_param) ? (REAL) params[j] : 0.0f;
	for (j = 0; j < n_element; j++)
		layer->species[j] = (REAL) species[j];
	return 1;
}

int OracleEquationSublayers(int eqn_index) {
	if (eqn_index < 0 || eqn_index >= N_PROBE_EQNS) return -1;
	return probe_eqns[eqn_index].rcmd_sublayers;
}

/* --------------------------------------------------- absorber / fuzz --- */

/* Mark the first `n` layers as absorber (dead layer / window). */
int OracleSetAbsorber(int n_layers) {
	if (n_layers < 0 || n_layers > 64) return 0;
	probe_sample.absorber_layers = n_layers;
	probe_sample.fres_only_absorber = 0;
	return 1;
}

/* Attach lateral thickness roughness to a layer. */
int OracleSetLayerFuzz(int layer_index, double amount, int steps) {
	if (layer_index < 0 || layer_index >= 64) return 0;
	probe_layers[layer_index].fuzzd = (REAL) amount;
	probe_layers[layer_index].fuzzs = steps;
	return 1;
}
