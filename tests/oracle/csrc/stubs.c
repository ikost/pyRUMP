/* Minimal host environment for the RUMP physics translation units.
 *
 * The physics code (ziegler.c, stopping.c, sigma.c, ...) is written against the
 * GENPLOT application shell: terminal I/O, a file-search path, and an expression
 * evaluator. None of that is needed to evaluate a stopping power, so this file
 * supplies just enough to link.
 *
 * Deliberate behaviours:
 *   - TTYprintf/ERRprintf are captured into a ring buffer rather than printed, so
 *     pytest output stays clean and tests can assert on warnings.
 *   - SysFindFile resolves against a caller-supplied data directory.
 *   - The GV* expression hooks are inert: except.stp overrides are simply not
 *     applied. Tests must therefore not rely on them.
 */

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>

#define LOG_SIZE 65536
static char log_buffer[LOG_SIZE];
static size_t log_used = 0;

static void log_append(const char *text) {
	size_t n = strlen(text);
	if (log_used + n + 1 >= LOG_SIZE) return;      /* silently drop when full */
	memcpy(log_buffer + log_used, text, n);
	log_used += n;
	log_buffer[log_used] = '\0';
}

const char *OracleLog(void)   { return log_buffer; }
void        OracleLogClear(void) { log_used = 0; log_buffer[0] = '\0'; }

/* ---------------------------------------------------------------- output --- */
int TTYprintf(const char *fmt, ...) {
	char tmp[4096];
	va_list ap;
	va_start(ap, fmt);
	vsnprintf(tmp, sizeof(tmp), fmt, ap);
	va_end(ap);
	log_append(tmp);
	return 0;
}

int ERRprintf(const char *fmt, ...) {
	char tmp[4096];
	va_list ap;
	va_start(ap, fmt);
	vsnprintf(tmp, sizeof(tmp), fmt, ap);
	va_end(ap);
	log_append(tmp);
	return 0;
}

int  TTYputs(const char *s)    { log_append(s); return 0; }
int  TTYputsnl(const char *s)  { log_append(s); log_append("\n"); return 0; }
int  ERRputs(const char *s)    { log_append(s); return 0; }
void gen_warn(const char *s)   { log_append("WARNING: "); log_append(s); log_append("\n"); }
void gen_err(const char *s)    { log_append("ERROR: ");   log_append(s); log_append("\n"); }

/* ------------------------------------------------------------ file search --- */
/* RUMP's data directory, set from Python before loading any table. */
static char oracle_data_dir[1024] = ".";

void OracleSetDataDir(const char *dir) {
	strncpy(oracle_data_dir, dir, sizeof(oracle_data_dir) - 1);
	oracle_data_dir[sizeof(oracle_data_dir) - 1] = '\0';
}

int SysFindFile(char *result, const char *name, const char *path,
                const char *ext, int mode) {
	char candidate[2048];
	(void) path; (void) ext;
	snprintf(candidate, sizeof(candidate), "%s/%s", oracle_data_dir, name);
	if (access(candidate, mode) == 0) { strcpy(result, candidate); return 1; }
	if (access(name, mode) == 0)      { strcpy(result, name);      return 1; }
	*result = '\0';
	return 0;
}

/* CAREFUL: mytypes.h:316 does `#define fopen(path,mode) E_FOPEN(path,mode)`, so
 * this is not an auxiliary helper -- it is *every* file open in the RUMP sources.
 * It must genuinely open the file. (Undefining the macro first, since this
 * translation unit includes no RUMP headers, would be equivalent.) */
FILE *E_FOPEN(const char *path, const char *mode) {
	return fopen(path, mode);
}

/* --------------------------------------------- expression evaluator (inert) --- */
/* except.stp overrides go through these; with the stubs they never apply. */
int  GVLinkArray(const char *name, void *data, int npt, int type) {
	(void) name; (void) data; (void) npt; (void) type; return 0;
}
int  GVEvalArrayExpr(const char *expr, void *result, int npt) {
	(void) expr; (void) result; (void) npt; return 0;    /* 0 == failure */
}
void GVDeallocate(void *p) { (void) p; }
int  GVGetInfo(const char *name, void *a, void *b) {
	(void) name; (void) a; (void) b; return 0;
}
int  GVLinkFnc(const char *name, void *fnc, int nargs) {
	(void) name; (void) fnc; (void) nargs; return 0;
}

/* stricmp/strnicmp/nint are NOT stubbed here: sys/sys_2.c provides the real
 * implementations and is linked in, so defining them again would clash. */
