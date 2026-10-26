/* Exercise the actual Platt rounding guards and kernel under all MXCSR
 * rounding modes, including callers with flush-to-zero/denormals-are-zero.
 * Build with the same flags as euclid_platt; see the Makefile test target. */
#define main euclid_program_main
#include "euclid.cpp"
#undef main

#if IALIB != PLATT
#error "This test exercises the Platt backend"
#endif

static ival guarded_operation(const double *a, const double *b, bool product)
{
  IA_RND_SAVE(a, b);
  ival result = product ? ival(*a) * ival(*b) : ival(*a) + ival(*b);
  IA_RND_RESTORE(result);
  return result;
}

int main()
{
  platt_init();
  unsigned int original;
  read_mxcsr(original);
  const double one = 1.0, half_ulp = ldexp(1.0, -53);
  const double next_one = 1.0000000000000002;
  const double product_lo = 1.0000000000000004;
  const double product_hi = 1.0000000000000007;
  // Synthetic exact-double embeddings give a nontrivial bound without a field
  // file.  The field-level tests separately check all 58 genuine embeddings.
  for (int j = 0; j < D; j++) {
    W1I[j] = ival(j == 0 ? 1.1 : 0.0);
    W2RI[j] = ival(j == 1 ? 1.3 : 0.0);
    W2II[j] = ival(j == 2 ? 1.7 : 0.0);
  }
  Cube cube = {1, 2, 3, 3};
  long gamma[D] = {0, 0, 0};
  double reference = norm_bound(cube, gamma);
  int failures = 0;
  for (unsigned int mode = 0; mode < 4; mode++) {
    for (unsigned int flush : {0U, 0x8040U}) {
      unsigned int cw = (original & ~0xE07FU) | (mode << 13) | flush;
      asm volatile ("ldmxcsr %0" : : "m" (cw) : "memory");
      ival sum = guarded_operation(&one, &half_ulp, false);
      unsigned int after;
      read_mxcsr(after);
      if (after != cw || sum.left != one || -sum.right != next_one) failures++;
      ival product = guarded_operation(&next_one, &next_one, true);
      read_mxcsr(after);
      if (after != cw || product.left != product_lo || -product.right != product_hi)
        failures++;
      double bound = norm_bound(cube, gamma);
      read_mxcsr(after);
      if (after != cw || !(bound == reference) || !(bound > 0.0)) failures++;
    }
  }
  asm volatile ("ldmxcsr %0" : : "m" (original) : "memory");
  printf("rounding guards: %d failures across all 4 modes, with/without FTZ/DAZ\n",
         failures);
  return failures ? 1 : 0;
}
