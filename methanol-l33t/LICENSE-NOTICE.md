# License and acknowledgements

Boy Lankhaar licenses the project-authored reconstruction code under the European Union Public Licence, version 1.2 (SPDX `EUPL-1.2`). The license text accompanies the source. [LICENSE](LICENSE)

## Spectroscopy software

Acknowledge BELGI-Cₛ and cite J. T. Hougen, I. Kleiner and M. Godefroid (1994), *Journal of Molecular Spectroscopy* **163**, 559–586. The archived BELGI source also credits Jonathan Fisher, Ortigoso and Li-Hong Xu. Preserve the original notices wherever their material is used. [BELGI-Cₛ paper](https://doi.org/10.1006/jmsp.1994.1047)

BELGI supplies the native Hamiltonian matrices. The project provides prediction drivers, A-species parity assignment and cached dipole contractions. The interface and overlays include BELGI-derived material and source-matching fragments, including those in `fortran/cached_inten.f` and `src/methanol_lamda/parity_repair.py`. The project license covers the project's contributions to these files. Original material remains subject to its upstream terms. The 58 original BELGI files are separately obtained build inputs. Their archived snapshot carries author notices and contains no separate upstream license text. [Source lineage](provenance/source-lineage.json)

## Molecular data

Cite Boy Lankhaar, *L33T methanol molecular data for excitation calculations*, version 1.0 (2026). Identify the isotope, A/E symmetry and DAT checksum used. Credit the original Hamiltonian fits, collision calculations, line lists and catalogues listed in the companion and the provenance tables. For BASECOL rates, cite the relevant scattering authors and BASECOL2023. [Scientific references](docs/references.bib) · [BASECOL2023](https://doi.org/10.1051/0004-6361/202348233)

The molecular files and their CSV provenance tables incorporate data from several sources. The code license applies to project-authored software only. The release makes no separate blanket license grant for the molecular data, raw inputs or upstream code. The individual sources and their treatment are recorded in the [source notices](docs/SOURCE-TERMS.md).
