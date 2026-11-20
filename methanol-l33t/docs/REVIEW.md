# Historical rc4 response to the independent appraisal

This dated record documents the November 6 review of the release candidate. The v1.0 distribution and credits are described in the [source notices](SOURCE-TERMS.md).

Source candidate 0.1.0rc4, 6 November 2026. The numerical target remains the
four frozen L33T v3 files. On 6 November 2026 the maintainer adopted this dataset
for numerical experiments. The present source revision updates that status and
preserves the molecular bytes and model approximations. Public release remains
under review; the previous rc2 numerical verification is retained as evidence.

| Finding | Response and evidence | Remaining limit |
|---|---|---|
| F1: unequal collider coverage | `provenance/collider-coverage.json` counts rows and actual level endpoints for every collider and manifold. README and companion give the measured table, OPR mixture formula and missing-coupling consequences. | No new ortho-H2/He extrapolation or fallback is introduced. Endpoint presence is not pair completeness or positive coupling. |
| F2: inaccessible exact inputs | The separate review ZIP includes all 97 locked inputs, including 58 native BELGI files and frozen derived maps. Import, export and input checks use exact hashes. | Public hosting and artifact-specific terms remain owner release decisions. Human source adjudication is not regenerated. |
| F3: gzip-coupled verification | Default verification checks decompressed CSV hashes/lengths; optional archive-byte mode retains the original compressed identity test. A disk-backed label-keyed DAT comparator reports explicit tolerances, missing/extra identities and field errors. Tests exercise changed compression, row order, identity/grid/rate defects and preserved zeros. | Numerical agreement is not exact release identity or physical validation. No second-platform result is claimed. |
| F4: BELGI/local implementation | README, companion, source notice and CITATION.cff distinguish BELGI matrices from local parity/assignment and cached contractions. `check_native_cache.py` repeats PERMANENT/VCO_A comparisons through J=12 for both isotopes; see `provenance/native-cache-validation.json`. | The comparison shares the repaired eigensystem; c-operator validation remains inherited independent evidence. GLOBALFIT is not run. |
| F5: builder test coverage | Synthetic fixtures now import the builder rules: absolute A prefactor, copy reversal, signed K translation, A parity, isotope factors, tie handling, observed-frequency precedence and E-hybrid ratios/zeros/fallback. | These are bounded fixtures plus regression against frozen products, not a fresh primary-paper or complete Hamiltonian audit. |
| F6: audit naming | README/companion explicitly call the streaming check a structural/format audit. Consumer detailed balance is identified as an arithmetic check. | No independent physical accuracy certification is implied. |
| F7: E-hybrid zero handling | Exact zero quartet ratios are honoured using explicit None checks. Locked inputs contain no such zero ratio, and all hybrid factors remain identical. RF E zeros (370 pair-temperature nodes per collider) and their propagation are documented. | RF zero inheritance remains part of the frozen model. |
| F8: provenance/code hygiene | Successful build manifests record passed verification; historical builder key corrected; known local roots become role-based paths; patch guards survive optimized Python; unreachable b-exchange and standalone hybrid rewrite paths removed. Header masses and the non-independent raw ledger are explained. | The display-mass inconsistency is retained to preserve v3 bytes. Source rates use the separately pinned monoisotopic masses. |
| F9: relocation/platform evidence | A full rebuild is exercised from an extracted complete bundle outside the research workspace. Commands and code identities are retained in the maintainer validation receipt. | No Linux/container runtime was available. Second-platform native build and comparison remain untested; no remote job was submitted. |
| F10: release records | BELGI citation metadata and upstream/local contributor boundaries are explicit; original vendor notices accompany the local bundle. `LICENSE-NOTICE.md` records what is known. | Owner selection of the package code license, final contributor list and artifact-specific input/product terms is still required before public release. No license grant is invented from a citation. |

## Reproduction sequence

Unzip the complete review bundle into a new directory. From `methanol-l33t`:

```sh
python -B -m unittest discover -s tests -v
python -B scripts/l33t.py check-inputs --inputs inputs
python -B scripts/l33t.py rebuild --inputs inputs --output work/rebuild
python -B scripts/l33t.py verify --products work/rebuild/products --archive-bytes
python -B scripts/check_native_cache.py --inputs inputs --output work/native-cache
python -B scripts/l33t.py extract --products work/rebuild/products \
  --isotope both --species both --gzip --output work/delivery
```

On another platform, retain a failed exact-check build for diagnosis, and run
`compare` against the reference products as documented in the README. Preserve
both the numerical report and any sidecar/identity differences. Do not replace
the frozen hashes simply to obtain a pass.

The tested environment and completed checks are in `provenance/validation.json`.
A source ZIP alone intentionally omits inputs. The complete local review ZIP
has an exact member inventory and contains no generated DAT or native output.

## Independent follow-up and rc4 corrections

Claude independently confirmed the F1 coverage counts, archive/input identities,
31 tests in normal, -O and -OO modes, source lineage and native/cache agreement
for 143,176 strengths. No new numerical defect was found. F1, F3-F7 are addressed;
F2 is resolved for bundle holders, while public acquisition remains a release
decision. F8's remaining optimization guards are corrected below. F9 still lacks
a second-platform result. F10 still requires owner release decisions.

The follow-up's N1-N5 prompted rc4: extraction manifests now say adopted for
numerical experiments; the immutable DAT title is explained as historical;
validation separates current metadata from the dated rc2 rebuild; twelve remaining
invariant asserts become explicit exceptions; generator help tolerates stripped
docstrings; child generators clear PYTHONOPTIMIZE; all-zero E rows are quantified;
and native/cache rounding and a command-block introduction are corrected. All
scientific file hashes and valid-input numerical operations remain unchanged.
SOURCE-MANIFEST.json is the complete exported inventory; source-lineage.json
records the files inherited from the parent source. Public release remains pending.
