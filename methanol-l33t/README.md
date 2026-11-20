# Methanol L33T

Reproducible construction of the **L33T v3** A/E LAMDA-format files for
CH3OH and 13CH3OH. This source package calculates torsion–rotation energies
and transition strengths, converts them to Einstein A coefficients, normalizes
published collision matrices, applies the documented extrapolations, and writes
the four molecular files with row-level provenance.

Public release **v1.0** distributes the frozen numerical model `l33t-v3`, adopted
for numerical experiments on 6 November 2026. The four ready-to-use species
archives carry the original molecular and provenance checksums. The collision
extrapolations, vibrational prescriptions and finite basis set limits on their
use. [Download v1.0](https://github.com/blankhaar/methanol-l33t/releases/tag/v1.0)

See the [companion document](docs/companion.pdf) ([editable LaTeX](docs/companion.tex))
for the methods, source hierarchy and limitations.

| File | Levels | Radiative transitions |
|---|---:|---:|
| `a-12ch3oh-l33t-E2000-J46.dat` | 7,600 | 146,635 |
| `e-12ch3oh-l33t-E2000-J46.dat` | 7,669 | 163,932 |
| `a-13ch3oh-l33t-E2000-J46.dat` | 7,714 | 149,820 |
| `e-13ch3oh-l33t-E2000-J46.dat` | 7,794 | 167,491 |

The four files contain 7,811,683 collision rows in total. Levels have energy
below 2,000 cm^-1 relative to the ground of their A/E species, with J <= 46
and a native |K| <= 22 output box. Included manifolds are vt=0–4, CO stretch,
CH3 rock, and nu4, nu5, nu6 and nu10. Each file has separate para-H2, ortho-H2
and He blocks on the 10, 20, 30, 50, 75, 100, 150, 200, 300, 400 K grid.

The literal `[EXPERIMENTAL]` tag in the four frozen DAT titles is historical.
It is retained only to preserve the v3 checksums; the maintainer's 6 November
2026 adoption supersedes that label. Software displaying the raw title may
still show it. Current delivery manifests identify the adopted-use status.

## Requirements

- Python 3.11 or newer; the generation path uses its standard library.
- GNU Fortran >= 10 and LAPACK/BLAS; on macOS the compiler links Apple Accelerate.
- The exact external inputs listed in [inputs.lock.json](provenance/inputs.lock.json).
- For the companion PDF only: a LaTeX installation with `latexmk`, `pdflatex`
  and BibTeX.

Allow 10 GiB of free disk space and about 8 GiB of available memory as a
conservative working provision. The existing four DAT files occupy 1.39 GB;
sidecars, native spectra and temporary files add to this. Resource use is not
a formal cross-platform minimum. Native calculations run with one BLAS/OpenMP
thread. No cluster account or radiative-transfer code is required.

## Collision coverage and mixture policy

Coverage is unequal. The table counts **levels appearing in any row**, including
zero coefficients; it does not imply complete pair coverage or positive coupling.

| Product | p-H2 rows / levels touched | o-H2 rows / levels touched | He rows / levels touched |
|---|---:|---:|---:|
| A-12 | 1,528,128 / 7,600 | 31,350 / 836 | 332,615 / 2,684 |
| E-12 | 1,615,634 / 7,669 | 31,350 / 836 | 335,378 / 2,696 |
| A-13 | 1,556,888 / 7,714 | 31,350 / 836 | 332,615 / 2,684 |
| E-13 | 1,649,647 / 7,794 | 31,350 / 836 | 335,378 / 2,696 |

Only para-H2 receives the J/K donor extrapolation. Ortho-H2 leaves 6,764--6,958
levels without any ortho row; He leaves 4,916--5,098 without any He row.
[Measured coverage](provenance/collider-coverage.json) includes every manifold.

The file-level policy keeps the three physical colliders separate. For an ortho
fraction `f = OPR/(1+OPR)`, a consumer forms
`C = n(H2)*[(1-f)*k_p + f*k_o] + n(He)*k_He` for each pair and temperature.
An absent row contributes no term in this representation; it is **unknown
coupling, not an established physical zero**. Increasing OPR can therefore
remove much of the represented high-level H2 excitation. No automatic
ortho-from-para or H2-from-He fallback is supplied. Such substitutions require
an explicit, separately documented consumer model; physical He must not be
counted twice. Declare densities, OPR and any fallback in each calculation.

## Inputs and redistribution

Python stages the coefficients and compiles native Fortran prediction routines.
BELGI supplies the Hamiltonian matrices. Local overlays perform the A-species
parity-block eigensolve/state assignment and all cached line-strength contractions.
The GLOBALFIT program is not compiled or run; native unpruned INTEN is used
only by the separate comparison check. This distinction preserves upstream
BELGI-Cs credit and identifies the local numerical work. The source-only
archive contains the interface and overlays; the local review bundle also
contains all locked inputs. Neither archive contains generated molecular files.

Users of this package or its molecular files should acknowledge **BELGI-Cs**
and cite **Hougen, Kleiner & Godefroid (1994), Journal of Molecular
Spectroscopy 163, 559–586**, [doi:10.1006/jmsp.1994.1047](https://doi.org/10.1006/jmsp.1994.1047),
alongside the fitted methanol Hamiltonian and collision sources in the companion.
The project-authored code is licensed under EUPL-1.2. The upstream material and the derived data retain their source-specific terms. [License and acknowledgements](LICENSE-NOTICE.md)

A separate, privately held reconstruction bundle includes all **97 exact inputs (153,588,511 bytes)**,
including 58 BELGI Fortran files and the frozen identity/observation tables.
A holder of that bundle can reproduce the model with the following commands.

```sh
cd methanol-l33t
python -B scripts/l33t.py check-inputs --inputs inputs
python -B scripts/l33t.py rebuild --inputs inputs --output work/rebuild
```

The source-only archive remains suitable for reviewing code separately. To
prepare either archive from an existing checkout:

```sh
python scripts/export_source.py --output dist/methanol-l33t-source-1.0.zip
python scripts/export_source.py --inputs work/inputs \
  --output dist/methanol-l33t-review-1.0.zip
```

The input lock gives relative filenames, exact SHA-256 hashes, byte sizes and
roles. A user who already has the original
authorized source trees can stage them with:

```sh
python scripts/l33t.py inputs \
  --resource-root /path/to/methanol-resource \
  --repair-root /path/to/parity-repair-source \
  --dest work/inputs
```

These paths are **source roots for an import**, not runtime dependencies.
The command copies only the explicitly listed files after checking their hashes.
Alternatively, arrange separately obtained inputs beneath `work/inputs`
according to each record's `path`, then run:

```sh
python scripts/l33t.py check-inputs --inputs work/inputs
```

[Input acquisition and provenance](docs/INPUTS.md) identifies the source
families. The exact original BELGI snapshot and other raw inputs are held separately from the public repository.

## Reproduce all four files

From this directory:

```sh
python scripts/l33t.py rebuild --inputs work/inputs --output work/rebuild
```

The output directory must be new. The command:

1. Reconstructs all 419,322 normalized collision rows from the 18 BASECOL
   matrices and checks the frozen master checksum.
2. Compiles the native Hamiltonian with the parity/full-vector overlay and
   evaluates the permanent, a-axis and direct c-axis dipole operators for
   both isotope Hamiltonians.
3. Builds the levels, source-bound frequencies, Einstein A values, vibrational
   proxies, E-species collision hybrid and signed J/K extrapolation.
4. Writes the four DAT files and compressed level/radiation/collision ledgers.
5. Checks exact DAT hashes and decompressed sidecar hashes, then runs the
   separate streaming structural/format audit.

The numerical products are in `work/rebuild/products/12` and `products/13`.
`REPRODUCTION.json` and the build manifests record success only after all checks
pass. Shared manifests use role-based paths such as `@inputs` and `@output`. Native stages
retain compiler commands, source hashes and state-assignment diagnostics.
The collision table derived inside `work/inputs` is reusable only if its
reference hash matches. Original imported inputs are unchanged.

Stages can also be invoked separately:

```sh
python scripts/l33t.py normalize --inputs work/inputs
python scripts/l33t.py spectroscopy --inputs work/inputs \
  --isotope 12 --output work/spectroscopy-12
python scripts/l33t.py spectroscopy --inputs work/inputs \
  --isotope 13 --output work/spectroscopy-13
python scripts/l33t.py build --inputs work/inputs \
  --spectroscopy12 work/spectroscopy-12 \
  --spectroscopy13 work/spectroscopy-13 --output work/products
```

The low-level spectroscopy script permits diagnostic J/K limits. Reduced
spectra do **not** define the published-size v3 product: the builder's empirical
calibration requires its complete reference transitions, and the reference
check will reject changed products. Cross-platform floating-point differences
may also fail exact DAT checks. Preserve the outputs and compare them by quantum
labels with the command below; reference hashes are never changed automatically.

## Extract files for a radiative-transfer calculation

For the main-isotope A/E pair, with deterministic gzip compression:

```sh
python scripts/l33t.py extract --products work/rebuild/products \
  --isotope 12 --species both --gzip --output work/main-isotope-delivery
```

Use `--isotope both` for all four, or `--species A`/`E` for one symmetry.
Omit `--gzip` for ordinary `.dat` files. Extraction verifies all four reference
products and sidecars first, copies the selected complete files and their
provenance ledgers, and writes `EXTRACTION.json`. It does not truncate or
reindex the level system. To read a compressed DAT, use `gunzip` or Python's
`gzip` module.

The file syntax follows the [LAMDA format](https://home.strw.leidenuniv.nl/~moldata/molformat.html):
energies in cm^-1, Einstein A in s^-1, frequencies in GHz, and downward
collision coefficients in cm^3 s^-1. Weights are 2J+1; A and E populations
are separate. Large level counts may exceed a consumer's compiled limits.
Test its parser, statistical weights, interpolation and density conventions.
The files do not specify an H2 ortho:para ratio or a He abundance. He is a
separate collider; treating it as an H2 fallback is an explicit consumer model.

## Verification modes and scientific limits

`verify` requires exact DAT bytes and exact decompressed CSV content. Gzip
headers/compression levels may differ. `--archive-bytes` additionally requires
the original compressed bytes. A separate comparison matches quantum labels,
transition endpoints, collider IDs and temperature grids, independently of row
order. Its tolerances are numerical acceptance criteria, not model error bars:
energies use 1e-8 cm^-1 absolute plus 1e-10 relative; frequencies (GHz) and upper
energies (K) use 1e-9 absolute plus 1e-9 relative; A and collision coefficients
use 1e-7 relative with no absolute allowance. Zero rates must remain zero.
Weights, masses, grids and identity sets must match exactly. Sidecar provenance
is outside the numerical comparator and remains covered by `verify`.

```sh
python scripts/l33t.py verify --products work/rebuild/products --archive-bytes
python scripts/l33t.py compare --reference /path/to/reference-products \
  --products work/rebuild/products --report numerical-comparison.json
python scripts/l33t.py coverage --products work/rebuild/products \
  --report collider-coverage.json
python scripts/check_native_cache.py --inputs inputs --output work/native-cache
```

A numerical pass is reported separately from release identity. The structural
audit checks serialization, labels, selection rules and connectivity; it does
not independently validate physical rate values. The native/cache comparison
uses the same repaired eigensystem, tests PERMANENT and VCO_A through J=12,
and does not replace the inherited independent c-operator audit.


Run the rule tests and verify the generated files with:

```sh
python -m unittest discover -s tests -v
python scripts/l33t.py verify --products work/rebuild/products
```

[Validation evidence](provenance/validation.json) distinguishes the fresh
packaging checks from the inherited independent scientific audit. The exact
reference products and sidecars are in [products.json](provenance/products.json).
The prior audit checked the direct c operator through J=12 independently in
both isotopes, source-donor reconstruction and actual consumer rate assembly.
Those checks are described as inherited evidence, not repeated by unit tests.

The E hybrid inherits 37 RF-zero pairs per H2 collider at all ten temperatures
(370 pair-temperature nodes). Including identity copies, these give 333 all-zero
p-H2 rows and 407 all-zero o-H2 rows in each E file. Across all sources and
extrapolations, E-12 and E-13 contain 55,935 and 57,133 all-zero p-H2 rows,
respectively; [the measured breakdown](provenance/zero-rate-coverage.json)
separates their origins. These are represented zero rates, distinct from absent rows. Exact zero quartet ratios are now
honoured, with synthetic fixtures; none occurs in the locked input set, so
this correction leaves the v3 rates unchanged. The DAT display masses retain
the frozen conventions (12C: 32.04186; 13C: 33.02957 u); isotope collision
scaling uses the separate monoisotopic masses recorded in the code.

The dominant limitations are physical: high-J and high-|K| collision proxies,
the empirical E collision hybrid, copied excited-manifold rates, sparse
torsion-changing collisions, endpoint sqrt(T) continuation, frozen-ground
vibrational operators and isotope transfers. No new quantum-scattering
calculation is performed. The companion document gives measured backtest
errors; neither successful serialization nor exact reproduction bounds the
errors of a maser prediction.

The supported `l33t.py` workflow clears `PYTHONOPTIMIZE` for child generators.
The remaining generation invariants also use explicit exceptions, and direct
builder/spectroscopy help works under `-OO`.

## Source layout and maintenance

`scripts/l33t.py` is the supported interface. `src/methanol_lamda` contains
the frozen molecular-data implementation; `fortran` contains prediction
drivers and the parity/dipole overlays. The original compatibility modules
remain in the snapshot, but their historical command interfaces are outside
this package's tested release workflow. `scripts/build_l33t.py` is the v3
assembly, with portable path handling. `highj_scheme_v2.py` and
`hybrid_e_lowj.py` expose the collision prescriptions. `provenance` pins the
inputs, products and validation; `docs` contains the editable scientific guide.

Build the document with:

```sh
cd docs
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build companion.tex
```

## Citation and verification history

Cite Boy Lankhaar (2026), *L33T methanol molecular data for excitation calculations*, version 1.0. Identify the isotope, symmetry species and exact DAT checksum. Credit BELGI-Cₛ and the spectroscopy and collision authors listed in the companion. See the [builder citation](CITATION.cff), [source notices](docs/SOURCE-TERMS.md) and [source bibliography](docs/references.bib).

The full native reconstruction ran on macOS/arm64 for release candidate rc2 and reproduced all four DAT files and twelve provenance tables. The rc4 review checked the test suite and the main-isotope PERMANENT tables through J = 46. The v1.0 downloads retain those molecular bytes and document the same numerical prescriptions. Native reproduction on a second platform remains untested. [Validation evidence](provenance/validation.json) · [Release-candidate review](docs/REVIEW.md)
