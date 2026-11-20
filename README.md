# L33T methanol molecular data v1.0

Four LAMDA-format files describe the A and E species of ¹²CH₃OH and ¹³CH₃OH for non-LTE excitation calculations. Choose one file for each isotope and symmetry species in your model. Each contains level energies, radiative transitions, and collisional rate coefficients for para-H₂, ortho-H₂ and He.

## Choose a file

Download the archive for your species from the table. Each archive contains one ready-to-use `.dat` file and three compressed CSV tables documenting its levels, radiative transitions and collision rates.

| Molecule | Species | Release asset | Levels | Radiative transitions |
|---|---|---|---:|---:|
| ¹²CH₃OH | A | [l33t-v1.0-a-12ch3oh.zip](https://github.com/blankhaar/methanol-l33t/releases/download/v1.0/l33t-v1.0-a-12ch3oh.zip) | 7,600 | 146,635 |
| ¹²CH₃OH | E | [l33t-v1.0-e-12ch3oh.zip](https://github.com/blankhaar/methanol-l33t/releases/download/v1.0/l33t-v1.0-e-12ch3oh.zip) | 7,669 | 163,932 |
| ¹³CH₃OH | A | [l33t-v1.0-a-13ch3oh.zip](https://github.com/blankhaar/methanol-l33t/releases/download/v1.0/l33t-v1.0-a-13ch3oh.zip) | 7,714 | 149,820 |
| ¹³CH₃OH | E | [l33t-v1.0-e-13ch3oh.zip](https://github.com/blankhaar/methanol-l33t/releases/download/v1.0/l33t-v1.0-e-13ch3oh.zip) | 7,794 | 167,491 |

For example, unpack the main-isotope A file and give the resulting `products/a-12ch3oh-l33t-E2000-J46.dat` to your excitation code.

```sh
unzip l33t-v1.0-a-12ch3oh.zip -d products
```

To check a single download, compare its `shasum -a 256` digest to the corresponding entry in [SHA256SUMS](SHA256SUMS). If you have all four downloads, put them and `SHA256SUMS` in one directory and run the following there. Run the last command from the repository root with the path to the unpacked directory. The verifier requires Python 3.11 or newer.

```sh
shasum -a 256 -c SHA256SUMS
for archive in l33t-v1.0-*.zip; do unzip "$archive" -d products; done
```

```sh
python3 -B methanol-l33t/scripts/l33t.py verify --products /path/to/products --archive-bytes
```

The [product manifest](methanol-l33t/provenance/products.json) gives SHA-256 values for each DAT file and for both the compressed and expanded provenance tables.

## Physical coverage

The level set extends to J = 46 and |K| = 22. Each file retains levels below 2,000 cm⁻¹ relative to the ground state of that symmetry species. It covers `vt = 0–4`, the CO stretch, the CH₃ rock and the ν₄, ν₅, ν₆ and ν₁₀ bands. Collision coefficients are tabulated at 10, 20, 30, 50, 75, 100, 150, 200, 300 and 400 K.

Only para-H₂ has an extrapolation to high J and high |K|. Ortho-H₂ rows reach 836 levels in each file. He rows reach 2,684 levels in each A file and 2,696 in each E file. These counts include represented zero coefficients. Missing collider pairs have unknown rates. Report the H₂ ortho-to-para ratio, He abundance and any substitution you use for missing rates. The [coverage tables](methanol-l33t/provenance/collider-coverage.json) give the counts by manifold.

Some excited-manifold rates are copied, the carbon-13 collision model uses isotope scaling, and the vibrational treatment uses frozen-ground operators. For maser work, the [scientific companion](methanol-l33t/docs/companion.pdf) describes where these approximations enter and gives the collision backtests. The `[EXPERIMENTAL]` text in the DAT headers dates from the original build. The maintainer adopted these exact files for numerical experiments on 6 November 2026.

## Reproduce and cite

Public release v1.0 includes builder version 1.0 and the frozen numerical model `l33t-v3`. Using a prebuilt DAT file requires only the selected archive. Regenerating the data requires GNU Fortran, LAPACK/BLAS and 97 separately obtained, checksum-locked inputs. See the [builder instructions](methanol-l33t/README.md), [input inventory](methanol-l33t/docs/INPUTS.md) and [editable companion](methanol-l33t/docs/companion.tex).

Cite Boy Lankhaar (2026), *L33T methanol molecular data for excitation calculations*, version 1.0, and give the exact DAT file and checksum used. Credit BELGI-Cs to Hougen, Kleiner & Godefroid (1994), the fitted Hamiltonians to Xu et al. (2008, 2014), and the collision calculations to Rabli & Flower (2010, 2011) and Dagdigian (2024), as applicable. BASECOL requests citation of the original rate calculations and its database paper. See the [dataset citation](CITATION.cff), [scientific references](methanol-l33t/docs/references.bib), [BASECOL citation policy](https://basecol.vamdc.eu/) and [BASECOL2023 paper](https://doi.org/10.1051/0004-6361/202348233).

The project-authored code is licensed under EUPL-1.2. The molecular data retain their source-specific terms. [License](LICENSE) · [Source notices](methanol-l33t/LICENSE-NOTICE.md) · [Source-by-source record](methanol-l33t/docs/SOURCE-TERMS.md) · [v1.0 release](https://github.com/blankhaar/methanol-l33t/releases/tag/v1.0)
