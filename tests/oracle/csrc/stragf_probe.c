/* Access to anlyz.c's file-static numerical helpers.
 *
 * `SimStragf` and `SimAnlyz3` are `static`, so they cannot be linked against
 * directly. Including the translation unit gives access without editing it;
 * anlyz.c is therefore compiled *here* and excluded from the source list.
 *
 * SimStragf is the integral of a unit triangle convolved with a Gaussian,
 * evaluated through a hand-tuned rational approximation with seven energy
 * regimes (anlyz.c:371-487). pyRUMP replaces it with a closed erf form; this
 * probe is what lets the difference between the two be measured rather than
 * assumed.
 */

#include "anlyz.c"

double OracleStragf(double x, double sig) {
	return (double) SimStragf((REAL) x, (REAL) sig);
}
