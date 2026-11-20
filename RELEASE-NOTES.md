# Methanol L33T v1.0

The first public release provides four LAMDA-format molecular files for excitation calculations with A- and E-symmetry ¹²CH₃OH and ¹³CH₃OH. Each species ZIP includes the DAT file and three compressed CSVs tracing the level energies, radiative transitions and collision rates.

| Isotopologue | A symmetry | E symmetry |
|---|---|---|
| ¹²CH₃OH | [A-12 ZIP](https://github.com/blankhaar/methanol-l33t/releases/download/v1.0/l33t-v1.0-a-12ch3oh.zip) | [E-12 ZIP](https://github.com/blankhaar/methanol-l33t/releases/download/v1.0/l33t-v1.0-e-12ch3oh.zip) |
| ¹³CH₃OH | [A-13 ZIP](https://github.com/blankhaar/methanol-l33t/releases/download/v1.0/l33t-v1.0-a-13ch3oh.zip) | [E-13 ZIP](https://github.com/blankhaar/methanol-l33t/releases/download/v1.0/l33t-v1.0-e-13ch3oh.zip) |

Across the four files there are 30,777 levels, 627,878 radiative transitions and 7,811,683 collision rows. Only para-H₂ has the high-J and high-|K| donor extrapolation. The companion describes the sparser ortho-H₂ and He coverage, isotope scaling and vibrational approximations. [Companion PDF](https://github.com/blankhaar/methanol-l33t/blob/v1.0/methanol-l33t/docs/companion.pdf) · [Editable LaTeX](https://github.com/blankhaar/methanol-l33t/blob/v1.0/methanol-l33t/docs/companion.tex)

Version 1.0 distributes the adopted numerical model `l33t-v3`. Its frozen molecular and provenance files retain their exact checksums. The source tree includes the builder, tests, bibliography and input lock. A native rebuild requires the 97 original inputs separately. [Download checksums](https://github.com/blankhaar/methanol-l33t/releases/download/v1.0/SHA256SUMS) · [Product checksums](https://github.com/blankhaar/methanol-l33t/blob/v1.0/methanol-l33t/provenance/products.json)

Please cite Boy Lankhaar (2026), *L33T methanol molecular data for excitation calculations*, version 1.0, and identify the DAT file used. Credit BELGI-Cₛ and the original spectroscopy and collision sources as applicable. The project-authored code is licensed under EUPL-1.2. [Dataset citation](https://github.com/blankhaar/methanol-l33t/blob/v1.0/CITATION.cff) · [Source notices](https://github.com/blankhaar/methanol-l33t/blob/v1.0/methanol-l33t/LICENSE-NOTICE.md)
