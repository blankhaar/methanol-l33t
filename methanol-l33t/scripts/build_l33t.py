#!/usr/bin/env python3
"""Build l33t v3 with the owner-selected nu10 cross-band floor of 1e-7 D2.

Derived from the audited D49 builder (SHA256 d9ae09de32418c73e79b93d4ab8ed21a679550fdf41a376a0204d60c0cb4639e).
Owner choices, 2026-11-01: lower only the direct-c nu10 cross-band cutoff;
retain the audited collision model. All other radiation keeps its 1e-4-D2
floor. Levels, source constants, band moments and all collision rules are
unchanged. The complete construction and limitations are documented in the
D49 report and the v3 validation note. Reuses existing stage-1 spectroscopy.

Usage: build_l33t.py --inputs DIR --spectroscopy DIR
       --vco-c-spectroscopy DIR --output DIR [--isotope 13 --main-spectroscopy DIR]
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import bisect
import os
from collections import Counter, defaultdict
from dataclasses import replace
from pathlib import Path
import sys
import time

PACKAGE = Path(__file__).resolve().parents[1]
GRID = (10.0, 20.0, 30.0, 50.0, 75.0, 100.0, 150.0, 200.0, 300.0, 400.0)
ECUT = 2000.0
JEDGE = 12
DJMAX = 10
FLOOR = 1e-4
NU10_FLOOR = 1e-7
UPSTREAM_BUILDER_SHA256 = "d9ae09de32418c73e79b93d4ab8ed21a679550fdf41a376a0204d60c0cb4639e"
AU2D2 = 2.541746473 ** 2
REF2 = 0.27 ** 2
C2 = 1.4387768775039338
DAG = "dagdigian2024-h2"
# Sunaga, Bartok & Matyus (2026), arXiv:2602.17219v2 Table III: gas-phase experimental band
# origins (A; E relative to the A zero-point level) and computed <mu_fi>^2 in e^2 a0^2.
NEW_BANDS = {
    "vnu6": dict(mode="nu6 COH bend", symmetry="A'", origin_A=1320.6, origin_E=1297.5, mu2_A=6.14e-4, mu2_E=2.42e-4, axis="a"),
    "vnu5": dict(mode="nu5 CH3 sym. deformation", symmetry="A'", origin_A=1453.3, origin_E=1462.1, mu2_A=1.48e-4, mu2_E=1.41e-4, axis="a"),
    "vnu4": dict(mode="nu4 CH3 asym. deformation", symmetry="A'", origin_A=1486.1, origin_E=1483.3, mu2_A=1.74e-4, mu2_E=1.31e-4, axis="a"),
    "vnu10": dict(mode="nu10 CH3 asym. deformation", symmetry="A''", origin_A=1481.5, origin_E=1473.9, mu2_A=1.29e-4, mu2_E=9.67e-5, axis="c"),
}
EXCLUDED = {
    "nu11 CH3 rock (A'')": "<mu>^2 = 1.77e-6 / 1.61e-6 e^2 a0^2 (|mu| ~ 0.0034 D) is below the 0.01-D minimum effective transition dipole policy",
    "nu1, nu2, nu3, nu9 stretches": "band origins 2844-3685 cm-1 lie above the 2,000 cm-1 cutoff",
}
VIB_TAGS = ("vCO1", "vrock1", *NEW_BANDS)


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def split(label):
    """label -> (species, manifold, j, k, parity)."""
    p = label.split(":")
    if p[1].startswith("vt"):
        man, rest = p[1], p[2:]
    else:
        man, rest = p[1], p[3:]
    j = int(rest[0][1:])
    k = int(rest[1][1:])
    par = rest[2] if p[0] == "A" else ""
    return p[0], man, j, k, par


def make_label(sp, man, j, k, par):
    if man.startswith("vt"):
        return f"A:{man}:J{j}:K{k}:{par}" if sp == "A" else f"E:{man}:J{j}:K{k:+d}"
    return f"A:{man}:vt0:J{j}:K{k}:{par}" if sp == "A" else f"E:{man}:vt0:J{j}:K{k:+d}"


def phys_parity(label):
    sp, man, j, k, par = split(label)
    if sp != "A":
        return None
    vt = int(man[2:]) if man.startswith("vt") else 0
    p = (1 if par == "+" else -1) * (-1) ** (j + vt)
    return -p if man == "vnu10" else p


def read_c_strengths(path, level_key):
    """VCO_C native lines -> [(key, key, S mu^2)]; fails unless every A row links equal parity."""
    rows, seen = [], set()
    with open(path) as f:
        for r in csv.DictReader(f):
            sp = r["species"].strip()
            x, y = (level_key(sp, int(r[e + "_vt"]), int(r[e + "_j"]), int(r[e + "_k"]), r[e + "_parity"].strip())
                    for e in ("upper", "lower"))
            s = float(r["smu2_debye2"])
            if not (s > 0 and math.isfinite(s)) or abs(x.j - y.j) > 1 or x.j == y.j == 0 or frozenset((x, y)) in seen:
                raise ValueError(f"invalid VCO_C row {x} {y}")
            if 0 not in (x.vt, y.vt):
                raise ValueError(f"VCO_C row without a vt0 endpoint {x} {y}")
            if sp == "A" and phys_parity(x.encode()) != phys_parity(y.encode()):
                raise ValueError(f"c operator linking opposite-parity ground functions {x} {y}")
            seen.add(frozenset((x, y)))
            rows.append((x, y, s))
    return rows


def gz_csv(path, header, rows):
    with path.open("xb") as raw, gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as gz, \
            io.TextIOWrapper(gz, encoding="utf-8", newline="") as tx:
        w = csv.writer(tx, lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def missing_pair_donor(kh, kl, net):
    """Choose the frozen signed J/K donor for a pair absent from the target network."""
    import highj_scheme_v2 as scheme
    man = kh.manifold
    kmax, dj = max(abs(kh.k), abs(kl.k)), abs(kh.j - kl.j)
    choice, rule_prefix = None, ""
    if kmax <= JEDGE:  # the pinned v4 rule, unchanged
        if max(kh.j, kl.j) > JEDGE:
            choice = scheme.select_donor(kh, kl, net, edge=JEDGE, method="nearest_shift")
    s_k = kmax - (JEDGE - dj)  # translation that lets a same-dJ donor fit below the J edge
    if choice is None and s_k > 0:
        big = kh.k if abs(kh.k) >= abs(kl.k) else kl.k
        sigma = 1 if big >= 0 else -1
        mapped = []
        for kk in (kh, kl):
            nk, nj = kk.k - s_k * sigma, kk.j - s_k
            if (kk.species == "A" and nk < 0) or abs(nk) > nj:
                mapped = None; break
            mapped.append(scheme.Key(kk.species, man, nj, nk, ("+" if nk == 0 else kk.parity) if kk.species == "A" else ""))
        if mapped is not None:
            mh, ml = mapped
            rule_prefix = f"k_translation_s{s_k}+"
            if max(mh.j, ml.j) <= JEDGE:
                donor = net.rates.get(frozenset((mh, ml)))
                if donor is not None and max(donor.upper.j, donor.lower.j) <= JEDGE:
                    choice = scheme.DonorChoice(mh, ml, 1.0, "direct", donor.upper != mh, donor.upper, donor.lower)
            else:
                choice = scheme.select_donor(mh, ml, net, edge=JEDGE, method="nearest_shift")
    return choice, rule_prefix


def copy_orientation(u, l, tu, tl, rates, phys, weight, grid):
    """Copy a downward source row, reversing with source quantities when needed."""
    if phys[tu] >= phys[tl]:
        return tu, tl, rates, False
    gap = phys[u] - phys[l]
    fac = [weight[u] / weight[l] * (math.exp(-gap * C2 / t) if gap else 1.0) for t in grid]
    return tl, tu, tuple(v*x for v,x in zip(rates,fac)), True


def orient_serialized_tie(u, l, rates, phys, serialized, weight):
    if serialized[u] >= serialized[l]:
        return u, l, rates, False
    if phys[u] != phys[l]:
        raise ValueError(f"collision row not physically downward: {u} {l}")
    return l, u, tuple(v*weight[u]/weight[l] for v in rates), True


def main():
    ap = argparse.ArgumentParser(description=(__doc__ or 'Build the frozen L33T v3 molecular products.').split("\n")[0])
    ap.add_argument("--inputs", type=Path, required=True)
    ap.add_argument("--spectroscopy", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--species", nargs="+", default=["A", "E"])
    ap.add_argument("--isotope", choices=("12", "13"), default="12")
    ap.add_argument("--main-spectroscopy", type=Path, default=None,
                    help="12CH3OH stage-1 directory (CO calibration, HITRAN/V7 binding); defaults to --spectroscopy")
    ap.add_argument("--perpendicular", choices=("direct-c",), default="direct-c")
    ap.add_argument("--vco-c-spectroscopy", type=Path, default=None, help="stage-1 directory holding A{iso}/VCO_C (direct-c)")
    a = ap.parse_args()
    if a.perpendicular == "direct-c" and a.vco_c_spectroscopy is None:
        ap.error("--perpendicular direct-c needs --vco-c-spectroscopy")
    t0 = time.monotonic()
    ROOT = a.inputs.resolve()
    spec = a.spectroscopy.resolve()
    iso = a.isotope
    main_spec = (a.main_spectroscopy or a.spectroscopy).resolve()
    out = a.output.resolve()
    if out.exists():
        raise FileExistsError("choose a new output directory")
    for p in (PACKAGE / "src", PACKAGE / "scripts"):
        sys.path.insert(0, str(p))
    import candidate_export as ce
    import reapply_source_frequency_policy as rfp
    import highj_scheme_v2 as scheme
    import hybrid_e_lowj as hy
    from methanol_lamda.build import MasterLevel, _calculated_collisions
    from methanol_lamda.quantum import LevelKey
    from methanol_lamda.models import Level, LamdaData, RadiativeTransition, CollisionPartner, CollisionTransition
    from methanol_lamda.lamda import write_lamda, validate_lamda
    from methanol_lamda.physics import einstein_a_from_smu2, interpolate_or_sqrt_extrapolate_rate
    from methanol_lamda.co_stretch import COStretchKey, _co_levels, _main_co_terms, _main_ground_energies, _load_c13_model, MAIN_HITRAN_RELATIVE_PATH, CM1_TO_MHZ
    from methanol_lamda.isotopologues.ch3oh13.collisions import scale_rate_tuple, scaling_factor
    from methanol_lamda.isotopologues.ch3oh13.profile import MOLECULAR_MASS_U as MASS13
    from methanol_lamda.ch3_rock import RockKey, _c13_origins
    from methanol_lamda.isotopologues.ch3oh13.build import _load_main_collision_master
    from methanol_lamda import vibrational_revision as vr

    out.mkdir(parents=True)
    levels_csv = spec / f"A{iso}/PERMANENT/levels.csv"
    main_levels_csv = main_spec / "A12/PERMANENT/levels.csv"
    mass = 32.04186 if iso == "12" else MASS13
    constants_name = "Xu 2008" if iso == "12" else "Xu et al. 2014 Table 2"
    profile = json.loads((ROOT / "config/vibrational-revision-2026-10-16/sources.json").read_text())
    rows_all = list(csv.DictReader(levels_csv.open()))
    main_rows = rows_all if iso == "12" else list(csv.DictReader(main_levels_csv.open()))
    raw = {r["level_key"]: float(r["energy_raw_cm1"]) for r in rows_all}
    delta0 = raw["E:vt0:J0:K+0"] - raw["A:vt0:J0:K0:+"]  # ground A/E torsional splitting
    # Shared sources. HITRAN V7/V8 rows are bound to ground identities only up to J=22,
    # where the owner's fail-closed 2.5e-4 cm-1 lower-term binding holds (beyond it the
    # extrapolated BELGI ground terms differ by up to ~1 cm-1, e.g. A J35 K9). All J<=22
    # vibrational terms are therefore identical to the pinned v4; J>22 use labelled proxies.
    levels_j22 = out / "ground-binding-levels-J22.csv"
    with levels_j22.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(main_rows[0].keys()))
        w.writeheader()
        w.writerows(r for r in main_rows if int(r["j"]) <= 22)
    co_terms, co_audit = _main_co_terms(ROOT / MAIN_HITRAN_RELATIVE_PATH, levels_j22)
    with ce.new_main_ground(levels_j22):
        rock_terms, direct, rock_source_ledger, rock_audit = vr._corrected_rock_source(ROOT, profile)
    ladder_last = {}  # (species, K, parity) -> (J, offset) of the highest measured CO term
    ground_e = {LevelKey.decode(r["level_key"]): float(r["energy_cm1"]) for r in main_rows if r["vt"] == "0"}
    for key, term in co_terms.items():
        lad = (key.species, key.k, key.parity)
        if lad not in ladder_last or key.j > ladder_last[lad][0]:
            ladder_last[lad] = (key.j, term - ground_e[key])

    def co_extend(levels):
        # Replace the species-median offset of uncovered identities by the constant offset of
        # the highest measured J in the same (species, K, parity) ladder: the owner's rock rule.
        out_levels, counts = [], Counter()
        for x in levels:
            if x.source == "extrapolated_species_median_vibrational_offset":
                g = LevelKey(x.key.species, 0, x.key.j, x.key.k, x.key.parity)
                lad = ladder_last.get((x.key.species, x.key.k, x.key.parity))
                if lad is not None and x.key.j > lad[0]:
                    x = replace(x, energy_cm1=ground_e[g] + lad[1],
                                source=f"constant_offset_above_last_measured_J{lad[0]}_same_K_ladder")
            counts[x.source.split("_J")[0] if x.source.startswith("constant_offset") else x.source] += 1
            out_levels.append(x)
        return out_levels, counts
    c13, _ = _c13_origins(ROOT / ce.C13_ROCK_MODEL)
    main_origins = {sp: {int(k): v for k, v in x.items()} for sp, x in profile["main_parallel_rock_origins_cm1"].items()}
    permanent = ce.read_strengths(spec / f"A{iso}/PERMANENT/lines.csv")
    vco_a = ce.read_strengths(spec / f"A{iso}/VCO_A/native-lines.csv")
    perp_path = a.vco_c_spectroscopy.resolve() / f"A{iso}/VCO_C/native-lines.csv"
    vco_c = read_c_strengths(perp_path, LevelKey)
    vco_a_main = vco_a if iso == "12" else ce.read_strengths(main_spec / "A12/VCO_A/native-lines.csv")
    c13co = _load_c13_model(ROOT) if iso == "13" else None
    observations = rfp.frequency_maps(ROOT, ROOT / "observations/repaired-observed-frequency-map.json")[iso]
    master = _load_main_collision_master(ROOT)
    hfac, hinfo = hy.build_factors(hy.load_master(ROOT / "data/normalized/master/collisions.csv"), GRID,
                                   interpolate_or_sqrt_extrapolate_rate)
    newground = _main_ground_energies(main_levels_csv)
    co_all = {}
    for sp in ("A", "E"):
        g0 = [MasterLevel(LevelKey.decode(r["level_key"]), float(r["energy_cm1"]), float(r["weight"]), "l33t", r["level_key"], "")
              for r in main_rows if r["species"] == sp and r["vt"] == "0" and float(r["energy_cm1"]) < ECUT]
        co_all[sp] = {x.key: x for x in co_extend(_co_levels(g0, isotopologue="12CH3OH", main_terms=co_terms, c13_model=None)[0])[0]}
    # The 13-feature empirical CO calibration is a 12CH3OH quantity, transferred to 13CH3OH (D43).
    calibration = ce.calibrate_repaired_co(ROOT, profile, newground, vco_a_main, co_all)
    products = []

    for sp in a.species:
        tag = f"{sp.lower()}-{iso}ch3oh"
        # ---------------- levels ----------------
        tors = [MasterLevel(LevelKey.decode(r["level_key"]), float(r["energy_cm1"]), float(r["weight"]),
                            f"l33t-parity-full9-A{iso}", r["level_key"],
                            "global_parity_overlap_candidate" if sp == "A" else "legacy_ENCAL2_full9")
                for r in rows_all if r["species"] == sp and float(r["energy_cm1"]) < ECUT]
        ground = [x for x in tors if x.key.vt == 0]
        phys, basis, weight = {}, {}, {}
        for x in tors:
            lab = x.key.encode(); phys[lab] = x.energy_cm1; weight[lab] = x.weight; basis[lab] = f"BELGI parity/full-vector ({constants_name} constants)"
        if iso == "12":
            co, _ = _co_levels(ground, isotopologue="12CH3OH", main_terms=co_terms, c13_model=None)
            co, co_counts = co_extend(co)
        else:  # Mukhopadhyay (1995) subband origins and K=0 Delta B / Delta D on the isotope ground terms (D36)
            co, co_counts = _co_levels(ground, isotopologue="13CH3OH", main_terms=None, c13_model=c13co)
        for x in co:
            if x.energy_cm1 < ECUT:
                lab = x.key.encode(); phys[lab] = x.energy_cm1; weight[lab] = x.weight; basis[lab] = x.source
        base = LamdaData(tag, mass, [Level(i, x.energy_cm1, x.weight, x.key.encode()) for i, x in enumerate(ground, 1)])
        rock, rock_ledger = vr._levels(base, iso, rock_terms, main_origins, c13)
        rock_basis = {r["quantum_label"]: r["basis"] + "; origin " + r["origin_basis"] for r in rock_ledger}
        for x in rock:
            if x.energy_cm1 < ECUT:
                phys[x.quantum_label] = x.energy_cm1; weight[x.quantum_label] = x.weight; basis[x.quantum_label] = rock_basis[x.quantum_label]
        offsets = {}
        for vtag, band in NEW_BANDS.items():
            off = band["origin_A"] if sp == "A" else band["origin_E"] - delta0
            offsets[vtag] = off
            for x in ground:
                e = x.energy_cm1 + off
                if e < ECUT:
                    lab = make_label(sp, vtag, x.key.j, x.key.k, x.key.parity)
                    phys[lab] = e; weight[lab] = x.weight
                    basis[lab] = (f"ground rotational term + species band origin {off:.4f} cm-1 (Sunaga 2026 Table III experimental"
                                  + (")" if iso == "12" else "; 12CH3OH origin, no 13C shift: no gas-phase 13C source)"))
        # Serialization: deterministic order, ties nudged by nextafter (never used physically).
        order = sorted(phys, key=lambda lab: (phys[lab], lab))
        ser, index = {}, {}
        prev = -math.inf
        ties = 0
        for i, lab in enumerate(order, 1):
            e = phys[lab]
            if e <= prev:
                e = math.nextafter(prev, math.inf); ties += 1
            ser[lab] = e; prev = e; index[lab] = i
        manifold_of = {lab: split(lab)[1] for lab in phys}
        # ---------------- radiation ----------------
        rad, rad_ledger, zero_gap = {}, [], 0
        def add(x, y, s, origin, freq=None, replace_existing=False):
            nonlocal zero_gap
            if x not in phys or y not in phys:
                return
            nu10_crossband = (manifold_of[x] == "vnu10") != (manifold_of[y] == "vnu10")
            if s < (NU10_FLOOR if nu10_crossband else FLOOR):
                return
            if sp == "A" and phys_parity(x) == phys_parity(y):
                raise ValueError(f"parity-forbidden A edge {x} {y} ({origin})")
            u, l = (x, y) if (phys[x], x) > (phys[y], y) else (y, x)
            gap = phys[u] - phys[l]
            if gap <= 0:
                zero_gap += 1
                return
            pair = frozenset((u, l))
            if pair in rad and not replace_existing:
                raise ValueError(f"duplicate radiative endpoint pair {u} {l} ({origin})")
            nu = freq if freq is not None else gap * CM1_TO_MHZ
            rad[pair] = (u, l, s, nu, einstein_a_from_smu2(nu, s, weight[u]), origin,
                         "D24 observed" if freq is not None else "model energy gap")
        for x, y, s, _ in permanent:
            if x.species != sp:
                continue
            xl, yl = x.encode(), y.encode()
            obs = observations.get(frozenset((xl, yl)))
            add(xl, yl, s, "repaired_full_permanent_dipole" + ("; " + obs["source"] if obs else ""),
                float(obs["frequency_mhz"]) if obs else None)
            if x.vt == y.vt == 0:
                for vtag in VIB_TAGS:
                    add(make_label(sp, vtag, x.j, x.k, x.parity), make_label(sp, vtag, y.j, y.k, y.parity), s,
                        f"vt0_permanent_strength_to_{vtag}_intraband")
        mu2_rock = profile["rock_strength"]["squared_transition_moment_au"][sp] * AU2D2
        for x, y, s, _ in vco_a:
            if x.species != sp:
                continue
            for ru, gl in ([(x, y)] if x == y else [(x, y), (y, x)]):
                add(COStretchKey.from_ground(ru).encode(), gl.encode(), s * calibration["rate_scale"],
                    "repaired_pure_a_CO_crossband_empirical_normalization")
                add(RockKey.from_ground(ru).encode(), gl.encode(), s * mu2_rock / REF2,
                    "repaired_pure_a_parallel_rock_theoretical_proxy")
                for vtag, band in NEW_BANDS.items():
                    if band["axis"] == "a":
                        add(make_label(sp, vtag, ru.j, ru.k, ru.parity), gl.encode(),
                            s * band["mu2_" + sp] * AU2D2 / REF2, f"pure_a_{vtag}_crossband_Sunaga2026_moment")
        perp_dropped = 0.0
        for x, y, s in vco_c:
            if x.species != sp:
                continue
            for ru, gl in ([(x, y)] if x == y else [(x, y), (y, x)]):
                if ru.vt != 0:
                    continue  # a band level carries a vt0 ground function
                for vtag, band in NEW_BANDS.items():
                    if band["axis"] == "c":
                        add(make_label(sp, vtag, ru.j, ru.k, ru.parity), gl.encode(), s * band["mu2_" + sp] * AU2D2 / REF2,
                            f"perpendicular_{vtag}_crossband_Sunaga2026_moment_direct_c" + (f"_to_vt{gl.vt}" if gl.vt else ""))
        for line in direct:
            if line.upper.species != sp or abs(line.upper.j - line.lower.j) > 1 or line.upper.j == line.lower.j == 0:
                continue
            add(RockKey.from_ground(line.upper).encode(), line.lower.encode(), line.smu2_debye2,
                "corrected_raw_V7_empirical" if iso == "12" else "corrected_main_V7_strength_to_13_proxy", replace_existing=True)
        # ---------------- collisions ----------------
        tkeys = {x.key: i for i, x in enumerate(tors, 1)}
        native = _calculated_collisions(master, sp, tkeys, tors, GRID, allow_sqrt_endpoint_extrapolation=True)
        coll = {}
        hybrid_rows = 0
        def put(c, u, l, rates, source, status, note, prov):
            coll[(c, frozenset((u, l)))] = (u, l, tuple(rates), source, status, note, prov)
        for r in native:
            rates = r.rates
            note = r.transform
            if sp == "E" and r.source == DAG:
                f = hfac[("p" if r.collider_id == 2 else "o", r.source_upper_key, r.source_lower_key)]["f"]
                rates = tuple(v * x for v, x in zip(rates, f))
                note += "; collision-v4 E hybrid (RF E x Dagdigian/RF A-quartet scale)"
                hybrid_rows += 1
            put(r.collider_id, r.upper_key.encode(), r.lower_key.encode(), rates, r.source, r.status, note,
                f"{r.source_upper_key}->{r.source_lower_key}")
        ground_rows = [(c, v) for (c, _), v in coll.items() if split(v[0])[1] == "vt0" and split(v[1])[1] == "vt0"]
        copies = 0
        for c, (u, l, rates, source, status, note, prov) in ground_rows:
            gu, gl = split(u), split(l)
            gap = phys[u] - phys[l]
            for man in ("vt1", "vt2", "vt3", "vt4", *VIB_TAGS):
                tu, tl = make_label(sp, man, gu[2], gu[3], gu[4]), make_label(sp, man, gl[2], gl[3], gl[4])
                if tu not in phys or tl not in phys or (c, frozenset((tu, tl))) in coll:
                    continue  # native rows (vt1/vt2) take precedence; identities outside the cut are absent
                vals, n2 = rates, note + f"; vt0 exact-identity copy to {man}"
                tu, tl, vals, reversed_copy = copy_orientation(u, l, tu, tl, vals, phys, weight, GRID)
                if reversed_copy:
                    n2 += "; source detailed balance for energy reversal"
                put(c, tu, tl, vals, source, f"vt0-identity-copy-to-{man}", n2, prov)
                copies += 1
        # Extrapolation of missing intramanifold p-H2 pairs, |dJ| <= DJMAX.
        by_man = defaultdict(list)
        for lab in phys:
            by_man[manifold_of[lab]].append(lab)
        extrap, ext_prov, unsupported, rules = 0, [], Counter(), Counter()
        for man, labs in by_man.items():
            keyof = {lab: scheme.Key(*(lambda s_: (s_[0], man, s_[2], s_[3], s_[4]))(split(lab))) for lab in labs}
            labof = {v: k for k, v in keyof.items()}
            states = {keyof[lab]: scheme.State(keyof[lab], phys[lab], weight[lab]) for lab in labs if split(lab)[2] <= 15}
            rates_net = {}
            for (c, pair), (u, l, vals, *_r) in coll.items():
                if c == 2 and u in keyof and l in keyof and split(u)[2] <= 15 and split(l)[2] <= 15:
                    rates_net[frozenset((keyof[u], keyof[l]))] = scheme.Rate(keyof[u], keyof[l], vals)
            net = scheme.SourceNetwork(states, rates_net, GRID, "p-H2")
            labs_sorted = sorted(labs, key=lambda lab: split(lab)[2])
            jvals = [split(lab)[2] for lab in labs_sorted]
            for i, x in enumerate(labs_sorted):
                jx = jvals[i]
                hi = bisect.bisect_right(jvals, jx + DJMAX)
                for y in labs_sorted[i + 1:hi]:
                    if (2, frozenset((x, y))) in coll:
                        continue
                    h, l = (x, y) if (phys[x], x) > (phys[y], y) else (y, x)
                    kh, kl = keyof[h], keyof[l]
                    choice, rule_prefix = missing_pair_donor(kh, kl, net)
                    if choice is None:
                        unsupported[man] += 1; continue
                    vals = scheme.mapped_rates(choice, net)
                    rule = rule_prefix + choice.rule
                    rules[rule.split("+")[-1] + (" (after K translation)" if rule_prefix else "")] += 1
                    put(2, h, l, vals, "l33t-extrapolation", "l33t-JK-extrapolated",
                        f"{rule}; scale {choice.scale:.6g}; donor {labof.get(choice.donor_upper, choice.donor_upper.label())}->"
                        f"{labof.get(choice.donor_lower, choice.donor_lower.label())}; source-direction reversed {int(choice.reverse_source_direction)}", "")
                    extrap += 1
        # Final orientation by serialized energy (exact physical ties only).
        tie_swaps = 0
        final = []
        for (c, pair), (u, l, vals, source, status, note, prov) in coll.items():
            u, l, vals, reversed_tie = orient_serialized_tie(u, l, vals, phys, ser, weight)
            if reversed_tie:
                tie_swaps += 1
                note += "; exact-tie serialization reorientation"
            if iso == "13":  # exact-identity 12CH3OH transfer with the reviewed reduced-mass factor (D35/D40)
                vals = tuple(scale_rate_tuple(vals, c))
                note += f"; reduced-mass-scaled 12CH3OH proxy (factor {scaling_factor(c):.12g})"
            final.append((c, u, l, vals, source, status, note, prov))
        # ---------------- write ----------------
        public = [Level(index[lab], ser[lab], weight[lab], lab) for lab in order]
        rad_rows = sorted(rad.values(), key=lambda r: (index[r[0]], index[r[1]]))
        radiative = [RadiativeTransition(i, index[u], index[l], A, nu / 1000.0)
                     for i, (u, l, s, nu, A, origin, fp) in enumerate(rad_rows, 1)]
        partners = []
        descs = ((2, "CH3OH-pH2"), (3, "CH3OH-oH2"), (6, "CH3OH-He")) if iso == "12" else \
            ((2, "13CH3OH-pH2 [reduced-mass-scaled 12CH3OH proxy]"), (3, "13CH3OH-oH2 [reduced-mass-scaled 12CH3OH proxy]"),
             (6, "13CH3OH-He [reduced-mass-scaled 12CH3OH proxy]"))
        for c, desc in descs:
            rows = sorted((r for r in final if r[0] == c), key=lambda r: (index[r[1]], index[r[2]]))
            partners.append(CollisionPartner(c, desc, GRID, [CollisionTransition(i, index[r[1]], index[r[2]], r[3])
                                                             for i, r in enumerate(rows, 1)]))
        molecule = (f"{sp}-{iso}CH3OH l33t: all states below {ECUT:g} cm-1, J<=46, vt0-4 + vCO1 vrock1 vnu4 vnu5 vnu6 vnu10 "
                    f"[parity/full-vector BELGI; collision-v4 rules + signed-K/J extrapolation |dJ|<={DJMAX}] "
                    + f"[v3 direct-c nu10; nu10 cross-band floor {NU10_FLOOR:g} D2] " + "[EXPERIMENTAL]")
        data = LamdaData(molecule, mass, public, radiative, partners)
        errors = validate_lamda(data)
        if errors:
            raise ValueError(errors[:10])
        fname = f"{tag}-l33t-E{int(ECUT)}-J46.dat"
        write_lamda(out / fname, data)
        gz_csv(out / f"{tag}.levels.csv.gz", ["index", "label", "manifold", "physical_energy_cm1", "serialized_energy_cm1", "weight", "basis"],
               [[index[lab], lab, manifold_of[lab], repr(phys[lab]), repr(ser[lab]), weight[lab], basis[lab]] for lab in order])
        gz_csv(out / f"{tag}.radiative-ledger.csv.gz", ["index", "upper", "lower", "smu2_debye2", "frequency_mhz", "a_s", "operator", "frequency_policy"],
               [[i, u, l, repr(s), repr(nu), repr(A), origin, fp] for i, (u, l, s, nu, A, origin, fp) in enumerate(rad_rows, 1)])
        gz_csv(out / f"{tag}.collision-provenance.csv.gz", ["collider", "upper", "lower", "source", "status", "transform", "source_pair"],
               [[c, u, l, source, status, note, prov] for c, u, l, vals, source, status, note, prov in
                sorted(final, key=lambda r: (r[0], index[r[1]], index[r[2]]))])
        level_counts = Counter(manifold_of.values())
        rad_counts = Counter(r[5].split(";")[0] for r in rad_rows)
        coll_counts = Counter((c, status) for c, u, l, v, source, status, n, p in final)
        rec = dict(file=fname, sha256=sha(out / fname), bytes=(out / fname).stat().st_size, levels=len(public),
                   levels_by_manifold=dict(sorted(level_counts.items())), J_max=max(split(x)[2] for x in phys),
                   K_max=max(abs(split(x)[3]) for x in phys), serialization_ties=ties, radiative=len(radiative),
                   radiative_by_operator=dict(rad_counts), zero_gap_strengths_not_emitted=zero_gap,
                   perpendicular_strength_dropped_K0_K0=perp_dropped,
                   collisions={f"{c}:{s}": n for (c, s), n in sorted(coll_counts.items())},
                   e_hybrid_native_rows=hybrid_rows, vt0_identity_copies=copies, extrapolated_pH2_rows=extrap,
                   extrapolation_rules=dict(rules), unsupported_pairs_left_absent=dict(unsupported), tie_reorientations=tie_swaps,
                   co_calibration=dict(rate_scale=calibration["rate_scale"], effective_moment_debye=calibration["effective_moment_debye"]),
                   co_level_basis=dict(co_counts), band_offsets_cm1=offsets, ground_AE_splitting_cm1=delta0)
        products.append(rec)
        print(json.dumps({k: rec[k] for k in ("file", "levels", "levels_by_manifold", "radiative", "extrapolated_pH2_rows",
                                              "unsupported_pairs_left_absent", "collisions")}, default=str), flush=True)
        del rad, coll, final, data, partners

    manifest = dict(product_version="l33t-v3", upstream_builder_sha256=UPSTREAM_BUILDER_SHA256,
                    nu10_crossband_floor_smu2_debye2=NU10_FLOOR,
                    owner_choices="2026-11-01: nu10 floor 1e-7 D2; retain audited collisions",
                    status="built-awaiting-reference-verification", isotope=iso, main_spectroscopy=str(main_spec), created_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                    products=products, energy_cutoff_cm1=ECUT, j_max=46, temperature_grid_K=GRID, donor_edge_J=JEDGE,
                    extrapolated_delta_J_max=DJMAX, product_floor_smu2_debye2=FLOOR, new_bands=NEW_BANDS, excluded_bands=EXCLUDED,
                    hybrid_factor_info=hinfo["info"], co_source_audit=co_audit, rock_source_audit=rock_audit,
                    inputs={"spectroscopy_manifest": str(spec / "SPECTROSCOPY-MANIFEST.json"),
                            "spectroscopy_manifest_sha256": sha(spec / "SPECTROSCOPY-MANIFEST.json"),
                            "levels_csv_sha256": sha(levels_csv), "lines_csv_sha256": sha(spec / f"A{iso}/PERMANENT/lines.csv"),
                            "vco_a_sha256": sha(spec / f"A{iso}/VCO_A/native-lines.csv"),
                            "perpendicular_operator": a.perpendicular, "perpendicular_lines": str(perp_path),
                            "perpendicular_lines_sha256": sha(perp_path),
                            "main_levels_csv_sha256": sha(main_levels_csv), "main_vco_a_sha256": sha(main_spec / "A12/VCO_A/native-lines.csv"),
                            "collision_master_sha256": sha(ROOT / "data/normalized/master/collisions.csv"),
                            "observation_map_sha256": sha(ROOT / "observations/repaired-observed-frequency-map.json"),
                            "vibrational_profile_sha256": sha(ROOT / "config/vibrational-revision-2026-10-16/sources.json")},
                    code={"build_l33t.py": sha(Path(__file__)), "candidate_export.py": sha(Path(ce.__file__)),
                          "reapply_source_frequency_policy.py": sha(Path(rfp.__file__)), "highj_scheme_v2.py": sha(Path(scheme.__file__)),
                          "hybrid_e_lowj.py": sha(Path(hy.__file__))},
                    wall_seconds=time.monotonic() - t0, population_or_grid_runs=0)
    from review_support import portable_paths
    manifest = portable_paths(manifest, {ROOT: "@inputs", spec: "@spectroscopy", main_spec: "@main-spectroscopy", out: "@output", PACKAGE: "@package"})
    (out / "BUILD-MANIFEST.json").write_text(json.dumps(manifest, indent=2, default=str) + "\n")


if __name__ == "__main__":
    main()
