#!/usr/bin/env python3
"""E-type low-J collision hybrid: Rabli-Flower E structure, Dagdigian/RF A-type scale.

Every Dagdigian-sourced E-species row (p-H2 and o-H2, the vt0 natives and their
identity copies into other manifolds) in a low-J parent is replaced by

    H_c(pair, T) = RF_pE(pair, T) * R_c(q, T),
    R_c(q, T)    = sum_q S[D_c, A] / sum_q S[RF_p, A],

where q = {(J,|K|), (J',|K'|)} is the quartet of the pair, S = g_u k_ul
exp(-(E_u - E_l)/2kT) with each source's own energies, c = p-H2 or o-H2, and
all curves are put on the parent temperature grid with the builder's own
log-log interpolation / endpoint sqrt(T) policy. Rationale: Dagdigian's E-type
coefficients reproduce his A-type ones (E/A quartet ratio 0.94-0.99), i.e. they
lack the E-type torsion-rotation structure that Rabli-Flower's retain; the A-type
comparison carries the newer PES/close-coupling scale onto that structure.

A copy row is multiplied by the same factor f_c = H_c / D_c as its vt0 source
pair (factors are defined in the Dagdigian source orientation), so orientation,
detailed-balance and isotope multipliers are preserved. Nothing else changes.

Imported by build_l33t.py; this module has no standalone rewrite workflow.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import statistics
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

C2 = 1.4387768775039338
DAG, RF = "dagdigian2024-h2", "rabli-flower2010-h2-vt012"
NOTE = " [E low-J p/o-H2 hybrid: Rabli-Flower E structure x Dagdigian/RF A-quartet scale, collision-v4]"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def jk(label: str) -> tuple[int, int]:
    p = label.split(":")
    return int(p[2][1:]), abs(int(p[3][1:]))


def g(label: str) -> float:
    return 2.0 * jk(label)[0] + 1.0


def quartet(u: str, l: str) -> tuple:
    return tuple(sorted((jk(u), jk(l))))


def load_master(path: Path) -> dict:
    """(tag, species, frozenset(pair)) -> row dict; tags Dp, Do, RFp; vt0 only."""
    out = {}
    with path.open() as f:
        for r in csv.DictReader(f):
            if ":vt0:" not in r["upper_key"] or r["source"] not in (DAG, RF):
                continue
            tag = ("D" if r["source"] == DAG else "RF") + {"2": "p", "3": "o"}[r["collider_id"]]
            if tag == "RFo":
                continue
            row = dict(u=r["upper_key"], l=r["lower_key"],
                       ts=tuple(float(x) for x in r["temperatures_k"].split(";")),
                       ks=tuple(float(x) for x in r["rates_cm3_s"].split(";")),
                       eu=float(r["source_upper_energy_cm1"]), el=float(r["source_lower_energy_cm1"]))
            key = (tag, r["species"], frozenset((row["u"], row["l"])))
            if key in out:
                raise ValueError(f"duplicate master row {key}")
            out[key] = row
    return out


def oriented(row: dict, u: str, l: str) -> tuple:
    """Source coefficients for the process u -> l (detailed balance if the row is l -> u)."""
    if row["u"] == u:
        return row["ks"]
    gap = row["eu"] - row["el"]
    if gap < 0:
        raise ValueError("source row not downward in its own energies")
    return tuple(k * g(row["u"]) / g(row["l"]) * math.exp(-C2 * gap / t) for k, t in zip(row["ks"], row["ts"]))


def build_factors(master: dict, grid: tuple, interp) -> tuple[dict, dict]:
    on_grid = lambda ts, ks: [interp(ts, ks, t) for t in grid]
    sym = lambda row: [g(row["u"]) * k * math.exp(-C2 * (row["eu"] - row["el"]) / (2 * t))
                       for k, t in zip(on_grid(row["ts"], row["ks"]), grid)]
    ratio = {}
    for c in ("p", "o"):
        num, den = defaultdict(lambda: [0.0] * len(grid)), defaultdict(lambda: [0.0] * len(grid))
        for (tag, sp, pair), row in master.items():
            if tag != "D" + c or sp != "A" or ("RFp", "A", pair) not in master:
                continue
            q = quartet(row["u"], row["l"])
            for i, v in enumerate(sym(row)):
                num[q][i] += v
            for i, v in enumerate(sym(master[("RFp", "A", pair)])):
                den[q][i] += v
        ratio[c] = {q: [n / d if d > 0 else None for n, d in zip(num[q], den[q])] for q in num}
    classes = {}
    for c, table in ratio.items():
        by_class = defaultdict(lambda: [[] for _ in grid])
        for q, vals in table.items():
            cls = (abs(q[0][0] - q[1][0]), abs(q[0][1] - q[1][1]))
            for i, v in enumerate(vals):
                if v is not None:
                    by_class[cls][i].append(v)
        classes[c] = {cls: [statistics.median(v) if v else None for v in cols] for cls, cols in by_class.items()}
        classes[c]["global"] = [statistics.median([v for vals in table.values() for v in [vals[i]] if v is not None]) for i in range(len(grid))]
    factors, info = {}, Counter()
    for (tag, sp, pair), drow in master.items():
        if sp != "E" or tag not in ("Dp", "Do"):
            continue
        c = tag[1]
        u, l = drow["u"], drow["l"]
        rf = master.get(("RFp", "E", pair))
        if rf is None:
            raise ValueError(f"no Rabli-Flower E counterpart for {u} -> {l}")
        rfk = on_grid(rf["ts"], oriented(rf, u, l))
        dk = on_grid(drow["ts"], drow["ks"])
        q = quartet(u, l)
        cls = (abs(q[0][0] - q[1][0]), abs(q[0][1] - q[1][1]))
        rq, src = [], []
        for i in range(len(grid)):
            v = ratio[c].get(q, [None] * len(grid))[i]
            if v is not None:
                src.append("quartet")
            else:
                v = classes[c].get(cls, [None] * len(grid))[i]
                src.append("class-median" if v is not None else "global-median")
                if v is None:
                    v = classes[c]["global"][i]
            rq.append(v)
        if any(d <= 0 for d in dk):
            raise ValueError(f"nonpositive Dagdigian coefficient for {u} -> {l}")
        h = [a * b for a, b in zip(rfk, rq)]
        factors[(c, u, l)] = dict(f=[x / d for x, d in zip(h, dk)], quartet=q, ratio_source=src,
                                  rf_zero_nodes=sum(x == 0 for x in rfk))
        info[f"{tag}:{'all-quartet' if set(src) == {'quartet'} else 'fallback'}"] += 1
        info[f"{tag}:rf_zero_nodes"] += sum(x == 0 for x in rfk)
    return factors, dict(info=dict(info), ratio=ratio)
