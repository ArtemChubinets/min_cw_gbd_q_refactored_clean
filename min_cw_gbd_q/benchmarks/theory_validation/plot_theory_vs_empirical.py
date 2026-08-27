#!/usr/bin/env python
"""Plot article theory vs empirical results for q-ary GBD theory-validation.

Mirrors the q=2 workflow of ``random_test_with_theory``: compare the article's
theoretical predictors (E_min, T_total from ``article_theory.full_prediction``)
against the empirical ensemble summaries written by ``theory_runner.py``.

Input
-----
One or more result directories, each containing a ``summary.json`` (a JSON list
of per-code dicts).  ``n`` and ``k`` are parsed from the directory name
(``..._n40_k6`` -> n=40, k=6; ``..._n32k16`` -> n=32, k=16).  ``q`` is read from
each entry's ``q`` field.  This makes the script work both for the per-k ensemble
dirs (each holding all q in one summary.json) and for the per-q feasibility dirs.

Aggregation per (q, k) point, with ``rate = k / n``:

* empirical quality    = mean over *hit* codes of ``mean_weight_given_hit``
                         (error bar = sample std / sqrt(n_hit))
* empirical attempts   = mean over *hit* codes of ``mean_attempts_given_hit``
* theory quality       = ``article_E_min`` (identical for all codes of a point)
* theory complexity    = ``theory["T_total"]``
* empirical complexity = ``cost_factor * empirical_attempts`` where
                         ``cost_factor = T_total / expected_attempts``

The cost factor reduces EXACTLY to ``3 * q**(k/2)`` when ``collision_depth == 1``
and ``k`` is even (the ensemble config) - verified numerically via
``full_prediction``: for even k and d=1,
``scan_fraction(S, 1) = 1`` so
``T_total = q^(k/2) * (1 + 2*1) / (1 - e^-S) = 3 q^(k/2) * expected_attempts``.
Using ``T_total / expected_attempts`` (instead of hard-coding ``3 q^(k/2)``) keeps
the empirical marker correct for any collision depth (e.g. the d=0.5 feasibility
runs), where the naive ``3 q^(k/2)`` factor would be wrong.

Output
------
A PNG with two side-by-side subplots (quality, complexity; the complexity plot is
log2-scaled), three colour series q = 3, 4, 5.  A cap_rate table (rows q, columns
k) is printed to stdout.
"""
import argparse
import json
import math
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


Q_ORDER = [3, 4, 5]
# One colour per q, fixed so the two subplots share the same legend mapping.
Q_COLOR = {3: "tab:blue", 4: "tab:orange", 5: "tab:green"}
Q_MARKER = {3: "o", 4: "s", 5: "^"}


def parse_nk(dirname):
    """Extract (n, k) from a result directory name like ``..._n40_k6``."""
    m_n = re.search(r"n(\d+)", dirname)
    m_k = re.search(r"k(\d+)", dirname)
    if not m_n or not m_k:
        raise ValueError("cannot parse n/k from directory name: %r" % dirname)
    return int(m_n.group(1)), int(m_k.group(1))


def load_summary(directory):
    path = Path(directory) / "summary.json"
    if not path.is_file():
        raise FileNotFoundError("missing summary.json in %s" % directory)
    with open(path, "r") as fh:
        data = json.load(fh)
    if not isinstance(data, list) or not data:
        raise ValueError("summary.json is not a non-empty list: %s" % path)
    return data


def aggregate(codes):
    """Aggregate one (q, k) point from its list of per-code dicts."""
    hit = [c for c in codes if c.get("mean_weight_given_hit") is not None]
    n_codes = len(codes)
    n_hit = len(hit)
    n_cap = n_codes - n_hit

    weights = np.asarray([c["mean_weight_given_hit"] for c in hit], dtype=float)
    attempts = np.asarray([c["mean_attempts_given_hit"] for c in hit], dtype=float)

    def mean_err(vals):
        if vals.size == 0:
            return None, None
        mean = float(vals.mean())
        std = float(vals.std(ddof=1)) if vals.size > 1 else 0.0
        err = std / math.sqrt(vals.size)
        return mean, err

    w_mean, w_err = mean_err(weights)
    a_mean, _ = mean_err(attempts)

    theory = codes[0]["theory"]
    theory_E_min = codes[0].get("article_E_min")
    theory_T_total = theory.get("T_total")
    theory_expected_attempts = theory.get("expected_attempts")

    # cost_factor = T_total / expected_attempts; equals 3*q^(k/2) for d=1, even k.
    cost_factor = None
    if theory_T_total is not None and theory_expected_attempts:
        cost_factor = theory_T_total / theory_expected_attempts

    emp_complexity = None
    if cost_factor is not None and a_mean is not None:
        emp_complexity = cost_factor * a_mean

    return {
        "n_codes": n_codes,
        "n_hit": n_hit,
        "n_cap": n_cap,
        "cap_rate": (n_cap / n_codes) if n_codes else None,
        "w_mean": w_mean,
        "w_err": w_err,
        "a_mean": a_mean,
        "theory_E_min": theory_E_min,
        "theory_T_total": theory_T_total,
        "cost_factor": cost_factor,
        "emp_complexity": emp_complexity,
    }


def collect_points(directories):
    """Return points [(q, k, n, rate, agg)] and the sorted k values."""
    points = []
    for d in directories:
        d = Path(d)
        n, k = parse_nk(d.name)
        codes = load_summary(d)
        by_q = {}
        for c in codes:
            by_q.setdefault(c["q"], []).append(c)
        for q, group in by_q.items():
            agg = aggregate(group)
            points.append((q, k, n, k / n, agg))
    points.sort(key=lambda t: (t[0], t[1]))
    ks = sorted({p[1] for p in points})
    return points, ks


def print_cap_table(points, ks):
    print("\ncap_rate (fraction of capped codes) - rows q, columns k:")
    header = ["q"] + ["k=%d" % k for k in ks]
    print("  ".join("%10s" % h for h in header))
    by_key = {(p[0], p[1]): p[4] for p in points}
    for q in Q_ORDER:
        row = ["q=%d" % q]
        for k in ks:
            agg = by_key.get((q, k))
            if agg is None:
                row.append("-")
            else:
                row.append("%d/%d" % (agg["n_cap"], agg["n_codes"]))
        print("  ".join("%10s" % c for c in row))
    print()


def build_figure(points):
    # Group per q, sorted by rate for line plotting.
    series = {q: [] for q in Q_ORDER}
    for q, k, n, rate, agg in points:
        if q in series:
            series[q].append((rate, agg))

    fig, (ax_q, ax_c) = plt.subplots(
        1, 2, figsize=(12, 5), constrained_layout=True
    )

    # ---- quality (left): weight vs rate ----
    for q in Q_ORDER:
        pts = series[q]
        if not pts:
            continue
        pts = sorted(pts, key=lambda t: t[0])
        color = Q_COLOR[q]

        theory = [(t[0], t[1]["theory_E_min"]) for t in pts
                  if t[1]["theory_E_min"] is not None]
        if theory:
            ax_q.plot([e[0] for e in theory], [e[1] for e in theory],
                      color=color, linewidth=1.8, alpha=0.85,
                      label="q=%d theory (E_min)" % q)

        emp = [(t[0], t[1]["w_mean"], t[1]["w_err"]) for t in pts
               if t[1]["w_mean"] is not None]
        if emp:
            ax_q.errorbar(
                [e[0] for e in emp], [e[1] for e in emp],
                yerr=[(e[2] if e[2] is not None else 0.0) for e in emp],
                color=color, marker=Q_MARKER[q], linestyle="None", markersize=6,
                capsize=3, elinewidth=1.2, label="q=%d empirical" % q)

    ax_q.set_xlabel("rate (k/n)")
    ax_q.set_ylabel("min weight given hit")
    ax_q.set_title("Quality: theory E_min vs empirical mean weight")
    ax_q.grid(True, alpha=0.3)
    ax_q.legend(fontsize=8)

    # ---- complexity (right): log2(complexity) vs rate ----
    for q in Q_ORDER:
        pts = series[q]
        if not pts:
            continue
        pts = sorted(pts, key=lambda t: t[0])
        color = Q_COLOR[q]

        theory = [(t[0], math.log2(t[1]["theory_T_total"])) for t in pts
                  if t[1]["theory_T_total"] is not None]
        if theory:
            ax_c.plot([e[0] for e in theory], [e[1] for e in theory],
                      color=color, linewidth=1.8, alpha=0.85,
                      label="q=%d theory (log2 T_total)" % q)

        emp = [(t[0], math.log2(t[1]["emp_complexity"])) for t in pts
               if t[1]["emp_complexity"] is not None]
        if emp:
            ax_c.plot([e[0] for e in emp], [e[1] for e in emp],
                      color=color, marker=Q_MARKER[q], linestyle="None",
                      markersize=6, label="q=%d empirical" % q)

    ax_c.set_xlabel("rate (k/n)")
    ax_c.set_ylabel("log2(complexity)")
    ax_c.set_title("Complexity: theory T_total vs empirical analog")
    ax_c.grid(True, alpha=0.3)
    ax_c.legend(fontsize=8)

    return fig


def main():
    ap = argparse.ArgumentParser(
        description="Plot article theory vs empirical q-ary GBD results."
    )
    ap.add_argument(
        "directories", nargs="*",
        help="result directories containing summary.json (default: ensemble glob)",
    )
    ap.add_argument(
        "--out", default="theory_vs_empirical.png",
        help="output PNG path (default: theory_vs_empirical.png)",
    )
    args = ap.parse_args()

    if args.directories:
        dirs = args.directories
    else:
        base = Path(__file__).resolve().parent / "results"
        dirs = sorted(str(d) for d in base.glob("ENSEMBLE_20260820_alpha1_d1_n40_k*"))
        if not dirs:
            print("ERROR: no ensemble result directories found under %s" % base,
                  file=sys.stderr)
            return 2

    points, ks = collect_points(dirs)

    print("Points:")
    for q, k, n, rate, agg in points:
        cost = agg["cost_factor"]
        naive = 3.0 * (q ** (k // 2))
        tag = "  [d=1, even k: cost_factor == 3*q^(k/2)]" if (
            cost is not None and abs(cost - naive) < 1e-6 * naive
        ) else ""
        print(
            "  q=%d n=%d k=%d rate=%.3f  n_hit=%d/%d cap=%d/%d  "
            "w_mean=%.3f+/-%.3f  attempts=%.3f  E_min=%.3f  T_total=%.3g  "
            "cost_factor=%.3g%s"
            % (q, n, k, rate, agg["n_hit"], agg["n_codes"], agg["n_cap"],
               agg["n_codes"],
               (agg["w_mean"] if agg["w_mean"] is not None else float("nan")),
               (agg["w_err"] if agg["w_err"] is not None else float("nan")),
               (agg["a_mean"] if agg["a_mean"] is not None else float("nan")),
               (agg["theory_E_min"] if agg["theory_E_min"] is not None
                else float("nan")),
               (agg["theory_T_total"] if agg["theory_T_total"] is not None
                else float("nan")),
               (cost if cost is not None else float("nan")), tag),
        )

    print_cap_table(points, ks)

    fig = build_figure(points)
    fig.savefig(args.out, dpi=150)
    plt.close(fig)
    print("Wrote PNG: %s" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
