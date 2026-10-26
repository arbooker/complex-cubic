/* Minimal stand-in for crlibm, so Platt's int_double14.2.h can be used without
 * it.  The header needs crlibm only for correctly-rounded transcendentals; the
 * covering kernel uses nothing but +, -, *, abs and comparisons, which
 * int_double implements directly in SSE with the rounding mode set to
 * round-down.  The stubs below abort: they exist so that the non-inline
 * _fpu_rndd() in the header links, and they guarantee that any *use* of a
 * transcendental is a loud failure rather than a silently wrong interval.
 */
#ifndef CRLIBM_SHIM_H
#define CRLIBM_SHIM_H
#include <stdlib.h>
#include <stdio.h>
#ifdef __cplusplus
extern "C" {
#endif
static inline void crlibm_init(void) { }
#define CRLIBM_STUB(f) \
  static inline double f(double x) { (void)x; \
    fprintf(stderr, "crlibm stub " #f " called: this build has no crlibm\n"); \
    abort(); return 0; }
CRLIBM_STUB(exp_rn) CRLIBM_STUB(exp_rd) CRLIBM_STUB(exp_ru)
CRLIBM_STUB(log_rn) CRLIBM_STUB(log_rd) CRLIBM_STUB(log_ru)
CRLIBM_STUB(log1p_rd) CRLIBM_STUB(log1p_ru)
CRLIBM_STUB(sin_rn) CRLIBM_STUB(sin_rd) CRLIBM_STUB(sin_ru)
CRLIBM_STUB(cos_rn) CRLIBM_STUB(cos_rd) CRLIBM_STUB(cos_ru)
CRLIBM_STUB(sinpi_rd) CRLIBM_STUB(sinpi_ru)
CRLIBM_STUB(cospi_rd) CRLIBM_STUB(cospi_ru)
CRLIBM_STUB(atan_rn) CRLIBM_STUB(atan_rd) CRLIBM_STUB(atan_ru)
CRLIBM_STUB(asin_rn) CRLIBM_STUB(acos_rn)
#undef CRLIBM_STUB
static inline double pow_rn(double x, double y) { (void)x; (void)y; abort(); return 0; }
#ifdef __cplusplus
}
#endif
#endif
