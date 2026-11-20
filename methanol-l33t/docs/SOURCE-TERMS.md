# Sources, credits and distribution in v1.0

The four species archives contain the constructed DAT files and row-level provenance. The source tree provides the builder, documentation, input checksums and reference product checksums. The 97 original build inputs are held separately.

| Material | Credit | Treatment in v1.0 |
|---|---|---|
| Project-authored Python, drivers and overlays | Boy Lankhaar | EUPL-1.2 for project contributions. Upstream fragments retain their notices. |
| BELGI-Cₛ Hamiltonian and source-matching fragments | Hougen, Kleiner & Godefroid (1994) and the authors in the upstream notices | The builder supplies prediction overlays. The 58 original Fortran files are checksum-identified external inputs. Their archived snapshot has no separate upstream license text. |
| BASECOL collision calculations | Rabli & Flower (2010, 2011), Dagdigian (2024), BASECOL2023 and VAMDC | State-to-state rates and their transformations are identified in the collision tables. The original 18 rate matrices and their level tables are external inputs. |
| CDMS 033502 catalogue | CDMS and the original laboratory measurements cited by the catalogue | Source-bound observed frequencies are identified in the radiative tables. The raw catalogue is an external input. |
| Xu V7/V8 HITRAN-format line lists | Xu et al. (2004) | The frozen CO/rock calculation uses a historical mirror. The raw line lists are external inputs. |
| Fitted constants, observation map and vibrational models | Xu et al. (2008, 2014), Mukhopadhyay et al. (1995), Predoi-Cross (1997) and sources in the companion | The input lock identifies the retained transcriptions and the frozen model configuration. |
| Derived molecular files and provenance | Boy Lankhaar and the underlying sources above | Cite version 1.0 and the exact DAT hash. Underlying data retain their source-specific terms. |

The [code notice](../LICENSE-NOTICE.md) defines the scope of the project license. Consult the provider records when reusing upstream material. [BASECOL](https://basecol.vamdc.eu/) · [BASECOL2023 paper](https://doi.org/10.1051/0004-6361/202348233) · [VAMDC citation policy](https://vamdc.org/structure/policy-citation/) · [CDMS citation guidance](https://cdms.ph1.uni-koeln.de/cdms/portal/general) · [Scientific bibliography](references.bib)
