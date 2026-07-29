/* RUMP's inverse normal CDF, extracted verbatim from lexp/gvcalc.c.
 *
 * NDTRI is reachable from the FUZZ feature (creatr.c:679), which uses it to
 * place Gaussian-weighted thickness replicas. It lives in gvcalc.c -- the
 * expression evaluator -- and that translation unit will not compile under
 * clang: it declares srand48/drand48 non-static then defines them static
 * (gvcalc.c:8481), and the SDK headers declare them too.
 *
 * Rather than patch or link the whole evaluator, the three functions actually
 * needed are copied here unchanged, so the oracle uses RUMP's own quantiles
 * rather than a substitute.
 */

#include <math.h>

#define PRIVATE  static
#define TMPREAL  double
#define SQRT(x)  sqrt(x)
#define LOG(x)   log(x)

PRIVATE TMPREAL poly_e(TMPREAL x, double *numer, int iorder) {
	
	TMPREAL tmp;
	numer += iorder;								/* Go to end of the array */
	tmp = *numer;
	while (iorder--) tmp = tmp*x + *(--numer);
	return(tmp);
}

PRIVATE TMPREAL rpoly_e(TMPREAL x, double *numer, int i1, double *denom, int i2) {

	TMPREAL top, bot;
	
	numer += i1;
	denom += i2;
	top = *numer;
	bot = *denom;
	while (i1--) top = top*x + *(--numer);
	while (i2--) bot = bot*x + *(--denom);
	return(top/bot);
}

TMPREAL NDTRI(TMPREAL p) {
	
#define LIMIT_ARG_VAL	38.47465;

	static double num[]   = {2.515517,		  0.802853,			0.010328};
	static double denom[] = {1.0,				  1.432788,			0.189269,		 0.001308};

	static double r2[] = {0.00029859715701, 0.00046690101221, -0.00041094492583, 1.4836042132e-005, 9.0445398334e-005,
								 -5.9417052852e-005, 2.4083658236e-005, -7.7312843132e-006, 1.9257277017e-006};
	static double r3[] = {0.00041854218865, -0.00014162624904, -0.00017261376736, 8.1839729896e-005, -1.2808004467e-005,
							   -2.3016012831e-006, 2.164855319e-006, -8.930115E-007, 3.281998E-007};
	static double r5[]  = {2.6521271801e-005, -0.00028301695906, 2.7742027695e-005, 1.5858032378e-005, -6.1300664207e-006,
								  1.2089310091e-006, -1.3779458613e-007,  -7.8864e-009, 7.9375e-009};
	static double r10[] = {-0.00033499183423, 0.0003643588597, 0.00016385965131, -0.00020531893775, 9.8624665016e-005, -2.7121319943e-005, 
								  -5.944520413e-007, 3.709582797e-006, -2.2303484285e-006, 4.2741e-006, -2.7178e-006};
	static double r25[] = {0.0004295823574, 9.4136771678e-005, -0.00019633197782, 7.2794900177e-005, -1.2725426308e-005,
								  -2.1193301912e-006, 3.0003340843e-006, -1.6890972433e-006, 9.0911607484e-007, -4.5248e-007, 1.1065e-007};
            
	/* Ten elements, exactly as gvcalc.c:4764 declares them -- NOT padded.
	 *
	 * poly_e reads taylor[10], one past the end. That is undefined behaviour,
	 * but at -O0 the adjacent static (ctay[0], ~6e-9) is folded to nothing by
	 * Horner and the result is correct. The shipped RUMP passes no -O flag, so
	 * this reproduces what RUMP actually computes rather than what its source
	 * intends. build_oracle.py pins -O0 for this file specifically. */
	static double taylor[] = {-4.5935005102e-010, 0.16666665975, 0.058333524176, 0.025197330986, 0.012041674594,
									   0.0060996947887, 0.0032662097955, 0.0014721469231, 0.0017529695707, -0.00033785703859};
	static double ctay[] = {6.4439238424e-009, -2.6928790649e-007, 3.7038273173e-006, -2.339886523e-005, 6.9000506134e-005, -2.0661086273e-005, 2.1380683449e-006};
	TMPREAL t;

	t = (p <= 0.5) ? p : 1.0-p;		/* Work with values on [0,0.5] */

	if (t <= 0) {							/* Lower/Upper bounds */
		t = LIMIT_ARG_VAL;
	} else if (t > 0.15) {										/* Near origin - get it right */
		t = (0.5-t)*2.50662827463100004;						/* Taylor expansion */
		t += t*poly_e(t*t,taylor,10);
		if (t > 0.40) {											/* Final corrections */
			t += poly_e(t*t*t*t, ctay, 6);					/* To rediculous precision */
		} else {
			t -= 1.15E-012*t;
		}
	} else {
		t  =  SQRT(-2.0*LOG(t));		/* LOG(1.0/(t*t)) */
		t  -= rpoly_e(t,num,2,denom,3);
		if (t < 2) {
			t -= poly_e(t-1.5, r2, 8);
		} else if (t < 3) {
			t -= poly_e(t-2.5, r3, 8);
		} else if (t < 5) {
			t -= poly_e(t-4.0, r5, 8);
		} else if (t < 14.75) {
			t -= poly_e((t-10.0)/5.0, r10, 10);
		} else {
			t -= poly_e((t-25.0)/10.0, r25, 10);
		}
	}

	return( p>0.5 ? t : -t);
}

/* Exposed so the port's quantiles can be checked against RUMP's own. */
double OracleNdtri(double p) { return NDTRI(p); }
