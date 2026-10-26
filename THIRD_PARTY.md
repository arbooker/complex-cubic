# Third-party dependencies and data provenance

## Platt interval header (external, unmodified)

The optional Platt backend uses David Platt's `int_double14.2.h` from
[djplatt/code](https://github.com/djplatt/code/blob/42b21426718e542daa2b006dc05ea2d7f26426e6/includes/int_double14.2.h).

- Repository revision: `42b21426718e542daa2b006dc05ea2d7f26426e6`.
- Header SHA-256: `37cc444d50998799e4783a3fb8b7ad038d379c6ddebf5827d41f320278ea308e`.
- Retrieval: `python3 src/fetch_platt.py`.

The header is **not distributed in this package**, and the project MIT
license does not apply to it. The helper retrieves the unmodified file from
the author's public repository into the ignored `work/deps/platt/` directory.
It performs no text substitutions. No license notice was found in the
public header or repository; consult its author about redistribution terms
if you intend to bundle it elsewhere.

Our rounding-control safeguards are in `src/euclid.cpp`: a volatile
control-register read, explicit dependencies around setting and restoring
the rounding mode, and memory clobbers. Compiler flags are in our Makefile.
None of these require modifying the external header.

`src/shim/crlibm.h` is our fail-fast shim, not crlibm code. Transcendental
functions abort if called; the covering kernel uses only arithmetic and
comparisons and does not require the crlibm library.

## MPFI and other numerical libraries

FLINT (including Arb), fplll, MPFR, GMP, PARI/GP, and MPFI are external
dependencies. Their licenses are not replaced by the project license.
No compiled library is included here. In particular, install MPFI separately;
neither its header nor the working directory's static MPFI library is included.

The MPFI source is the [Inria Git repository](https://gitlab.inria.fr/mpfi/mpfi).
The release specified for this package is **1.5.5**, tag `1.5.5`, commit
`918adeef5e938baa677fa5fd8a7d72757e8335ea`, dated 27 September 2026.
MPFI requires GMP >= 4.1.0 and MPFR >= 4.0.1 and is licensed under
LGPL-2.1-or-later. The exported `euclid_mpfi` was built and tested against
this release. Two sequential smoke tests passed on 22 October 2026:

| Discriminant | Verdict | Cubes | Maximum depth |
|---|---|---:|---:|
| -23 | PROVED | 57 | 2 |
| -876 | PROVED | 56001 | 10 |

Both counts and depths match the stored Platt results. These are targeted
compatibility tests for MPFI 1.5.5, not a full 58-field run with this release.
The exact MPFI revision used for the paper's earlier full cross-check was
not recorded.
The export accepts installed headers and libraries through
`MPFI_CFLAGS` and `MPFI_LIBS`.
Magma is proprietary software and is not distributed here.

## Field data

`data/fields.csv` was normalized from an
[LMFDB number-field download](https://www.lmfdb.org/NumberField/).
It selects complex cubic fields
with class number one and absolute discriminant less than 170,520. Labels
and polynomial coefficients follow the LMFDB conventions. Acknowledge
The LMFDB Collaboration, *The L-functions and Modular Forms Database*,
https://www.lmfdb.org/.

The LMFDB database is licensed under
[CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/), as stated on
its [data licence page](https://www.lmfdb.org/license) (checked 26 October
2026 against the [official page source](https://github.com/LMFDB/lmfdb/blob/main/lmfdb/templates/license.html)).
The field list is distributed under that licence. Our changes select the
stated range and normalize the records to CSV, retaining LMFDB labels and
defining polynomials; no endorsement by the LMFDB is implied.

`data/nf/pari.dat` and `data/nf/magma.dat` were computed from this field list
by the included input-generation scripts. `results/certificates.txt` is the
final size-reduced set, not a concatenation of historical runs. The supplied
`results/euclidean.jsonl` is the complete Platt run used for the paper's
22 October 2026 wall-clock measurement. Its per-field `secs` values are
elapsed subprocess times, not CPU times.

The authors' CC0 dedication covers their own rights only and does not
replace the CC BY-SA 4.0 terms for LMFDB-derived material.
