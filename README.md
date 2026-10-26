# Norm-Euclidean complex cubic fields

Supplementary code and data for *Determining norm-Euclidean complex cubic
fields*, by Gustav Kjærbye Bagger, Andrew R. Booker, Kevin J. McGown, and
Tim Trudgian.

Repository: https://github.com/arbooker/complex-cubic.

Of the 17,283 complex cubic fields with class number one and
absolute discriminant at most 170,520, this package proves 17,225 not
norm-Euclidean and the remaining 58 norm-Euclidean. There is no complex cubic
field of discriminant −170,520, so the supplied list uses the strict bound.

The package contains **one final certificate set**, used throughout the paper.
Running the search programs is not necessary to verify the results.

Run the commands below from the package root. Historical working files are
not part of this package; new outputs should be written under `work/`.

## Contents

| Path | Contents |
|---|---|
| `data/fields.csv` | All 17,283 fields, with LMFDB labels and defining polynomials |
| `data/nf/pari.dat`, `data/nf/magma.dat` | Integral bases and units computed independently in PARI/GP and Magma |
| `results/certificates.txt` | The 17,225 final, self-contained certificates |
| `results/euclidean.jsonl` | Covering results for the 58 norm-Euclidean fields (Platt backend) |
| `src/certify.cpp`, `src/run_certify.py` | Common certificate verifier and parallel driver |
| `src/euclid.cpp`, `src/run_euclid.py` | Norm-Euclidean covering prover and parallel driver |
| `src/ne_sweep.c`, `src/run_sweep.py` | Method I: residue marking |
| `src/ne_chain.c`, `src/sample.py` | Method II: sampling |
| `src/nfdata.gp`, `src/nfdata.m`, `src/cmp_nfdata.py` | Algebraic input preparation and basis comparison |
| `src/parse_lmfdb.py` | Conversion of an LMFDB field download to CSV |
| `src/primitive.py`, `src/primitive.gp`, `src/shrink.py` | Optional certificate size reduction |
| `src/check_data.py` | Release checksums, record counts, and field-coverage checks |
| `src/test_verification.py`, `src/test_rounding.cpp`, `src/test_known.py` | Regression tests |
| `SHA256SUMS` | Checksums of published source, documentation, and data |

The paper is the mathematical reference. A public paper link will be added
when available; older drafts and implementation notes are not included here.

## Build

The core verification programs need a C++14 compiler, FLINT 3 (including Arb),
fplll, MPFR, and GMP, with development headers. Python 3 drivers use only the
standard library. The search programs additionally need the PARI C library;
the algebraic-input scripts use PARI/GP or Magma. Magma is optional and is not
needed to verify the supplied certificates.

The tested environment is x86-64 Linux, GCC 13.3, FLINT 3.1.0, fplll 5.4.5,
MPFR 4.2.1, GMP 6.3.0, and PARI/GP 2.18.1. These are tested
versions, not a claim about the oldest supported releases.

For MPFI, use release **1.5.5** from the
[Inria Git repository](https://gitlab.inria.fr/mpfi/mpfi), tag `1.5.5`, commit
`918adeef5e938baa677fa5fd8a7d72757e8335ea` (27 September 2026).
MPFI 1.5.5 was tested with the exported backend on discriminants -23 and -876;
both proved norm-Euclideanity and matched the stored Platt cube counts and
depths. This is a targeted compatibility check, not a full 58-field run with
this release.
The Git checkout requires `autoreconf -v -i` before `./configure` if the
configure script is absent, as described in MPFI's README. This requires
Autoconf, Automake, Libtool, and M4 in addition to the compiler and libraries.

```sh
make -C src                         # certificate verifier and Arb covering prover
make -C src search                  # optional: Method I and Method II search
make -C src euclid_mpfi             # optional: installed MPFI headers and library
python3 src/fetch_platt.py          # optional: fetch pinned, unmodified Platt header
make -C src euclid_platt            # optional: x86-64/SSE Platt backend
```

The export Makefile uses `pkg-config` for FLINT, fplll, and PARI when available,
with ordinary library flags as a fallback. Custom installations can override
`FLINT_CFLAGS`, `FLINT_LIBS`, `FPLLL_CFLAGS`, `FPLLL_LIBS`, `PARI_CFLAGS`,
`PARI_LIBS`, `MPFI_CFLAGS`, and `MPFI_LIBS`. For example:

```sh
make -C src euclid_mpfi MPFI_CFLAGS=-I/opt/mpfi/include \
  MPFI_LIBS='-L/opt/mpfi/lib -Wl,-rpath,/opt/mpfi/lib -lmpfi'
```

No compiled libraries or executables are distributed. MPFI is an external
dependency, not a bundled static library. Platt's header is also an external
dependency: `fetch_platt.py` downloads an unmodified, revision-pinned copy
and verifies its SHA-256. It is stored under the ignored `work/deps/platt/`,
not among the published files. To use an existing copy, set `PLATT_INCLUDE`
to its absolute directory when invoking Make. See `THIRD_PARTY.md` for the
source revision and checksum. Do not use `-ffast-math`: the Platt build
requires the directed-rounding flags in the Makefile.

## Verify the data and certificates

```sh
python3 src/check_data.py
python3 src/cmp_nfdata.py data/nf/pari.dat data/nf/magma.dat
python3 src/run_certify.py --jobs 12 --argmin --out work/certified.jsonl
```

`check_data.py` verifies package checksums, unique field identities
(discriminant **and** defining polynomial), both algebraic-input lists, and
that the certificates and 58 covering fields partition the input list.
This is a consistency check, not a substitute for running the proof programs.

The final certificate file has 17,225 records: 16,825 explicit-modulus
(Method I) records and 400 implicit-modulus (Method II) records. Its size is
2,443,898 bytes; its longest record has 852 characters, excluding the newline.
Its SHA-256 is:

```text
67e9f2cc4c9c83bff7bc66f97770fc11b5c2a0b878f70bb820a6c64a679a328d
```

The `--argmin` run checks both non-norm-Euclideanity and the exact coset
minimum. It must produce 17,225 `CERTIFIED` verdicts and 17,225 `MINIMAL`
representatives, with no refusals. Omitting `--argmin` checks the certificates
without requiring an exact minimum. `REFUTED` means the witness fails;
`REFUSED` means the verifier could not complete a rigorous check, not that the
field is norm-Euclidean. Either is a failed verification run.

The driver distributes records among single-threaded processes. It is not
needed for the mathematics; the complete serial computation is:

```sh
mkdir -p work
src/certify --argmin < results/certificates.txt > work/certify.log
```

### Certificate format

Each non-comment line is self-contained:

```text
disc | a0,a1,a2,a3 | den | b00,b01,b02,b10,...,b22 | e0,e1,e2 | certificate
```

The polynomial is `a0 + a1*x + a2*x^2 + a3*x^3`, with ascending coefficients.
If `t` is its root, the integral basis is
`omega_j = (sum_k b_jk*t^k)/den`, and the unit is `sum_j e_j*omega_j`.
The last column is one of:

```text
b B0 B1 B2 t0 t1 t2       beta = sum_j Bj*omega_j
n n t0 t1 t2             beta = epsilon^n - 1
```

In both cases `alpha = sum_j tj*omega_j`. The witness pair is `(alpha,beta)`;
the format writes the modulus first for parsing. It asserts

$$m(\alpha+\beta\mathcal O_K)\geq |N\beta|,$$

and hence $M(K,\alpha/\beta)\geq1$. Every coefficient is an exact integer.
There is no fixed-width bound on a record or its coefficients.

The algebraic-input files use the same first five columns, followed by the
field discriminant rather than a certificate. PARI and Magma may give different
bases and units; `cmp_nfdata.py` checks that their bases span the same lattice.

## Prove norm-Euclideanity

The default backend is Arb, so the core build is sufficient:

```sh
python3 src/run_euclid.py --jobs 12 --out work/euclidean-arb.jsonl
```

To cross-check with the optional backends:

```sh
python3 src/run_euclid.py --jobs 12 --binary src/euclid_platt \
  --out work/euclidean-platt.jsonl
python3 src/run_euclid.py --jobs 12 --binary src/euclid_mpfi \
  --out work/euclidean-mpfi.jsonl
```

Each run must prove all 58 fields. The input list is fixed in `test_known.py`,
not inferred from old output. The driver writes results atomically and uses
the selected backend for its scheduling pre-pass as well as the proof pass.
All three builds use the same overall covering algorithm, with different
implementations of interval bounds; they are not three independent algorithms.
Agreement of verdicts is required. Subdivision counts can differ between Arb
and the endpoint-based backends.

### Timings reported in the paper

On 22 October 2026, an Intel Core i7-1265U laptop with twelve parallel
single-threaded workers took:

| Computation | Wall time |
|---|---:|
| All 17,225 certificates, including exact coset minima (`--argmin`) | 272.8 seconds |
| All 58 norm-Euclidean proofs, Platt backend, including driver setup and scheduling pre-pass | 6.53 seconds |

The paper rounds these to 4.5 minutes and 7 seconds. These are elapsed times,
not summed CPU times. The final certificates were selected from several
searches, so there is no single associated certificate-search time.

## Tests

```sh
make -C src test                    # Python bookkeeping tests; no Platt dependency
make -C src test-verification       # additionally tests SSE rounding guards
python3 src/test_known.py           # optional search regression; build search first
```

The tests check missing inputs, interrupted output handling, the fixed covering
list, and (on x86-64) rounding guards under all four hardware rounding modes.
Use a separate `work/` directory for new results; do not replace the published
certificate file when experimenting.

## Optional searches and algebraic input

```sh
python3 src/run_sweep.py --jobs 12 --outdir work/search
python3 src/sample.py --all --jobs 12 --results work/search/results.jsonl \
  --outdir work/sample
```

These runs search for new witnesses; they need not reproduce the particular
representatives in the final file. Method I marks residues modulo a chosen
small modulus. Method II samples residues modulo `epsilon^n - 1`. Search
output is provisional until checked by `certify`; an unsuccessful search is
not a norm-Euclideanity verdict. The size-reduction scripts accept certificate
files; use their `--help` options and write to a new output path.

To recompute algebraic input without overwriting the supplied records:

```sh
mkdir -p work
printf 'read("src/nfdata.gp"); ne_nfdata("data/fields.csv")\n' | gp -q > work/pari.dat
magma -b csv:=data/fields.csv src/nfdata.m > work/magma.dat
python3 src/cmp_nfdata.py work/pari.dat work/magma.dat
```

The field list was normalized from an LMFDB number-field download. The supplied
CSV fixes the input used here; regenerating algebraic input does not regenerate
or independently establish completeness of the field list.

## Correctness and notation

Given a correct basis of the full ring of integers, the certificate verifier
checks the order arithmetic, the supplied unit, unit stability of the coset,
`alpha` not belonging to `beta*O_K`, and an exhaustive enumeration with exact
norms and certified interval bounds. The unit need not be fundamental. No
class number or regulator is assumed by this verifier. Completeness of the
classification separately uses the field list and unconditional class-number
computations described in the paper.

LLL chooses a useful change of basis. Its unimodularity is checked exactly,
and the enumeration bound is reconstructed with interval arithmetic; no
guarantee of reducedness is trusted. A poor basis can increase the work or
cause a refusal, but is not grounds for accepting a certificate.

| Implementation | Paper |
|---|---|
| `theta`, `th` | Witness numerator $\alpha$ |
| Embedding order `(real, Re complex, Im complex)` | $(z,x,y)$, permuted from the paper's $(x,y,z)$ |
| `s1`, `s2` | Slice scales $a_i,b_i$ |
| `W` | Reduced basis $A_i$, with rows permuted to real-first order |
| `TG`, `tgt` | Coset coordinate centre $c_i$ |
| `nsl` | Number of slices |

Integral-basis coordinates in the records are unaffected by this embedding
permutation. The paper uses $\xi$ for a field element and $\gamma$ for an
algebraic integer. In `ne_chain.c`, the member `A` is `W` inverse times the
scaling matrix, not the paper's $A_i$.

## License and attribution

Original project code and documentation are under the MIT license (`LICENSE`).
The authors' rights in the data are dedicated under CC0 1.0 (`LICENSE-DATA`).
The LMFDB-derived field list, `data/fields.csv`, is distributed under
[CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/), with
attribution to The LMFDB Collaboration. See `THIRD_PARTY.md` for provenance,
changes to the source data, and external dependencies. Cite the accompanying
paper when using these computations.
