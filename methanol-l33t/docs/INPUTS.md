# Input sources and exact custody

The numerical workflow uses the exact files in `../provenance/inputs.lock.json`.
`source_root`/`source_path` describe a maintainer's existing source trees;
`path` is the portable destination beneath the chosen input directory.
Changing a file requires a separately reviewed input lock and product version.
The local import command performs no network requests.

| Input family | Role | Acquisition/custody |
|---|---|---|
| BELGI Fortran matrix routines | Native torsion–rotation Hamiltonian, compiled by the Python interface | Exact pre-existing research snapshot (58 Fortran files); external in the source-only archive and included in the local review bundle. Cite BELGI-Cs: Hougen, Kleiner & Godefroid (1994), doi:10.1006/jmsp.1994.1047. Exact public retrieval is unverified and the snapshot contains no upstream license text. Compare every hash. |
| Xu et al. 2008 constants | Main-isotope Hamiltonian | Archived final fitted constant listing; [paper](https://doi.org/10.1016/j.jms.2008.03.017). |
| Xu et al. 2014 Table 2 | Carbon-13 Hamiltonian | Retained constant transcription; [paper](https://doi.org/10.1016/j.jms.2014.06.005). |
| BASECOL2023 and 2025 | 18 collision matrices and their level tables | Archived [BASECOL](https://basecol.vamdc.eu/) releases; exact filenames and IDs in the lock. The v1.0 downloads include the derived rates and their provenance. |
| Frozen identity basis | Restricts collision normalization to the historical vt<=2, J<=22, absolute K<=14 basis | Compact derived level table. This preserves v3's source coverage, including exclusion of 6,108 source rows; expanding it would change the collision model. |
| Xu2008 observation crosswalk | Accepted frequency replacements | Retained derived identity map; unresolved assignments are not inferred during release preparation. It is a frozen input, not regenerated from the original full fit listing by this pipeline. |
| CDMS 033502 | Source-bound observed isotope frequencies | Archived negative-tag catalogue records; [CDMS](https://cdms.astro.uni-koeln.de/). |
| Xu et al. 2004 V7/V8 | Empirical CO/rock energies and strengths | Exact archived HITRAN-format mirror; [source paper](https://doi.org/10.1016/j.jms.2004.05.017). No claim of a current HITRAN release. |
| Mukhopadhyay et al. 1995; Predoi-Cross 1997 | Carbon-13 CO and rock models | Retained numerical transcriptions with their source notes; references in the companion document. |
| Source configuration | Thirteen-feature CO calibration, band moments and frequency overrides | Frozen reviewed numerical choices. These are not newly fitted in the release preparation. |

No source PDF is necessary at runtime. The minimal source-plus-input package
includes the Hamiltonian implementation dependency, fitted parameters, collision
tables, fixed identity/observation maps and model configuration. Normalized
collision matrices, native spectra, DAT files and expanded provenance are
generated. A public website alone is not evidence that the exact locked bytes
can still be recovered; public acquisition of the complete input set remains
unverified.

## Complete local review bundle

`export_source.py --inputs <locked-input-root> --output <new.zip>` includes exactly
the 97 locked files beneath `methanol-l33t/inputs`, together with the source.
Every archive member has a size and SHA-256 in `SOURCE-MANIFEST.json`. Unzip into
a new directory and run `python -B scripts/l33t.py check-inputs --inputs inputs`
from the package root before rebuilding. Normalized rates and generated spectra
are deliberately regenerated. These exact inputs total 153,588,511 bytes before
compression. This local bundle is the review reproduction route; no external
workspace, hidden credentials or network access is required at runtime.

The raw SHA256SUMS compatibility ledger is derived from the same input lock,
not a second independent integrity claim. The frozen observation and identity
maps are retained inputs; the pipeline does not redo their source adjudication.
GitHub hosts the derived molecular products and the source-only tree. The raw
inputs remain separate. See [source notices](SOURCE-TERMS.md) for credits and scope.
