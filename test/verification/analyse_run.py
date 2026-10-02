#!/usr/bin/env python3
"""
analyze_run.py - reduce a MOOSE run to the handful of numbers we care about.

Usage:
    analyze_run.py RUN.log                         # solver health only
    analyze_run.py RUN.log --ref REF.csv --test TEST.csv [--tol 1e-5]
    analyze_run.py RUN.log --ref-log REF.log       # side-by-side Newton stats

What it reports
---------------
1. Did the run finish?             time reached, failed attempts, dt cuts
2. How hard was Newton working?    iterations per step (histogram, median of
                                   plastic steps), convergence order p and
                                   contraction ratio on the last iterations
3. Did anything suspicious fire?   C_ep fallbacks, Newton->Primal switches,
                                   DIVERGED_* reasons, FD-check line
4. Same answer as the reference?   max relative difference per CSV column
                                   on common time points

Convergence order p is estimated from the last three residuals of each
converged step:  p = log(r3/r2) / log(r2/r1).  p ~ 1 linear, p ~ 2 quadratic.
Steps with fewer than 4 residuals (elastic, 1-2 iterations) are excluded.
"""
import argparse
import csv
import math
import re
import statistics
import sys
from collections import Counter

ANSI = re.compile(r"\x1b\[[0-9;]*m")
RE_STEP = re.compile(r"Time Step\s+(\d+),\s*time\s*=\s*([-+0-9.eE]+),\s*dt\s*=\s*([-+0-9.eE]+)")
RE_NL = re.compile(r"^\s*(\d+)\s+Nonlinear \|R\|\s*=\s*([-+0-9.eE]+|nan|inf)", re.I)
RE_NL_FAIL = re.compile(r"Nonlinear solve did not converge due to (\w+)")
RE_L_FAIL = re.compile(r"Linear solve did not converge due to (\w+)")
# an undamped Newton step can invert elements; MOOSE recovers by cutting dt, so
# the only trace is these lines. Without them a run that wrecked the mesh seven
# times looks, in the summary, like a run that took very few iterations.
RE_DEGEN = re.compile(r"Degenerate map|libMesh degeneracy exception")
RE_RASHID = re.compile(r"Rashid approximation for the rotation tensor")


def parse_log(path):
    attempts = []          # dicts: step, time, dt, res[], converged
    cur = None
    diverged = Counter()
    lin_fail = Counter()
    fallbacks = 0
    fallback_lines = []
    primal = 0
    fd_line = None
    mandel_line = None
    degen = 0
    rashid = 0

    with open(path, errors="replace") as f:
        for raw in f:
            line = ANSI.sub("", raw)
            m = RE_STEP.search(line)
            if m:
                cur = dict(step=int(m.group(1)), time=float(m.group(2)),
                           dt=float(m.group(3)), res=[], converged=None)
                attempts.append(cur)
                continue
            m = RE_NL.match(line)
            if m and cur is not None:
                try:
                    cur["res"].append(float(m.group(2)))
                except ValueError:
                    cur["res"].append(float("nan"))
                continue
            if "Solve Converged!" in line and cur is not None:
                cur["converged"] = True
            elif "Solve Did NOT Converge!" in line and cur is not None:
                cur["converged"] = False
            m = RE_NL_FAIL.search(line)
            if m:
                diverged[m.group(1)] += 1
            m = RE_L_FAIL.search(line)
            if m:
                lin_fail[m.group(1)] += 1
            if "C_ep fallback" in line:
                fallbacks += 1
                if len(fallback_lines) < 5:
                    fallback_lines.append(line.strip())
            if "switching to Primal" in line:
                primal += 1
            if fd_line is None and "FD check" in line:
                fd_line = line.strip()
            if mandel_line is None and "Mandel round-trip" in line:
                mandel_line = line.strip()
            if RE_DEGEN.search(line):
                degen += 1
            if RE_RASHID.search(line):
                rashid += 1

    return dict(attempts=attempts, diverged=diverged, lin_fail=lin_fail,
                fallbacks=fallbacks, fallback_lines=fallback_lines,
                primal=primal, fd_line=fd_line, mandel_line=mandel_line,
                degen=degen, rashid=rashid)


def newton_stats(attempts):
    conv = [a for a in attempts if a["converged"]]
    its = [max(len(a["res"]) - 1, 0) for a in conv]
    # iterations burned on attempts that were then thrown away. Reporting only
    # the converged ones made a run that failed 7 times and cut dt by 128x
    # print "total 2", which then read as 11.5x FEWER than a healthy run.
    wasted = sum(max(len(a["res"]) - 1, 0) for a in attempts if not a["converged"])
    orders, ratios = [], []
    for a in conv:
        r = [x for x in a["res"] if x > 0 and math.isfinite(x)]
        if len(r) >= 4:
            r1, r2, r3 = r[-3], r[-2], r[-1]
            if r1 > r2 > r3 > 0 and r2 / r1 < 1.0:
                orders.append(math.log(r3 / r2) / math.log(r2 / r1))
                ratios.append(r3 / r2)
    plastic_its = [i for i in its if i >= 3]
    return dict(n_conv=len(conv), its=its, plastic_its=plastic_its,
                orders=orders, ratios=ratios, wasted=wasted,
                its_all=sum(its) + wasted)


def fmt_hist(its):
    c = Counter(its)
    return "  ".join(f"{k}:{c[k]}" for k in sorted(c))


def report_log(name, P):
    A = P["attempts"]
    S = newton_stats(A)
    failed = [a for a in A if a["converged"] is False]
    dts = [a["dt"] for a in A]
    print(f"--- {name}")
    if not A:
        print("    no time steps found (wrong file?)")
        S["last_ok"] = float("nan")
        return S
    last_ok = max((a["time"] for a in A if a["converged"]), default=float("nan"))
    print(f"    attempts {len(A)}   converged {S['n_conv']}   failed {len(failed)}"
          f"   time reached {last_ok:g}   dt min/max {min(dts):.3g}/{max(dts):.3g}")
    S["last_ok"] = last_ok
    if S["its"]:
        print(f"    Newton its/step  total {sum(S['its'])}   histogram  {fmt_hist(S['its'])}")
    if S["wasted"]:
        print(f"    WASTED: {S['wasted']} more iterations on {len(failed)} failed attempt(s)"
              f"   -> {S['its_all']} total work, not {sum(S['its'])}")
    if S["plastic_its"]:
        print(f"    plastic steps (>=3 its): n={len(S['plastic_its'])}"
              f"  median {statistics.median(S['plastic_its']):g}  max {max(S['plastic_its'])}")
    if S["orders"]:
        print(f"    convergence order p  median {statistics.median(S['orders']):.2f}"
              f"   last-step contraction median {statistics.median(S['ratios']):.3g}"
              f"   ({'quadratic' if statistics.median(S['orders']) > 1.6 else 'linear'})")
    if P["diverged"]:
        print(f"    nonlinear DIVERGED: {dict(P['diverged'])}")
    if P["lin_fail"]:
        print(f"    linear DIVERGED:    {dict(P['lin_fail'])}")
    if P["degen"] or P["rashid"]:
        print(f"    MESH WRECKED: {P['degen']} inverted-element exception(s),"
              f" {P['rashid']} Rashid sqrt(<0) -- the step is overshooting, not"
              f" merely converging slowly")
    if failed:
        a = failed[0]
        rs = ", ".join(f"{x:.3e}" for x in a["res"][:8])
        print(f"    first failure: step {a['step']} t={a['time']:g} dt={a['dt']:.3g}  |R|: {rs}")
    print(f"    C_ep fallbacks {P['fallbacks']}   Newton->Primal {P['primal']}")
    for l in P["fallback_lines"]:
        print(f"      {l}")
    if P["fd_line"]:
        print(f"    {P['fd_line']}")
    if P["mandel_line"]:
        print(f"    {P['mandel_line']}")
    return S


def load_csv(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    out = {}
    for r in rows:
        try:
            t = round(float(r["time"]), 9)
        except (KeyError, ValueError):
            continue
        out[t] = r            # last row wins for duplicate times
    return out, (rows[0].keys() if rows else [])


def compare_csv(ref_path, test_path, tol, atol, tmax=None, cols_re=None):
    R, rcols = load_csv(ref_path)
    T, tcols = load_csv(test_path)
    cols = [c for c in rcols if c in tcols and c != "time"]
    if cols_re:
        cols = [c for c in cols if re.search(cols_re, c)]
    times = sorted(set(R) & set(T))
    if tmax is not None:
        times = [t for t in times if t <= tmax]
        R = {t: v for t, v in R.items() if t <= tmax}
    print(f"--- CSV  ref={ref_path}  test={test_path}")
    print(f"    common time points {len(times)}  (ref {len(R)}, test {len(T)})"
          f"   last common t={times[-1] if times else float('nan'):g}")
    if not times:
        print("    VERDICT: nothing to compare")
        return False
    worst_ok = True
    for c in cols:
        ref_vals = [float(R[t][c]) for t in times]
        scale = max(abs(v) for v in ref_vals) or 1.0
        # t_at starts at the first time so identical columns (all diffs = 0) still print
        maxd, maxabs, t_at, first_bad = 0.0, 0.0, times[0], None
        for t in times:
            ad = abs(float(T[t][c]) - float(R[t][c]))
            d = ad / scale
            if d > maxd:
                maxd, t_at = d, t
            maxabs = max(maxabs, ad)
            # a point passes if it is close relatively OR absolutely; the
            # absolute floor keeps numerically-zero columns (lateral stresses,
            # grain differences in the elastic range) from failing on noise
            if first_bad is None and d > tol and ad > atol:
                first_bad = t
        ok = first_bad is None
        worst_ok &= ok
        tag = "ok  " if ok else "DIFF"
        extra = "" if ok else f"   first > tol at t={first_bad:g}"
        note = "  (abs floor)" if ok and maxd > tol else ""
        print(f"    {tag} {c:26s} max rel {maxd:.2e}  max abs {maxabs:.2e}  at t={t_at:g}{extra}{note}")
    covered = len(times) == len(R)
    print(f"    VERDICT: {'PASS' if worst_ok and covered else 'FAIL'}"
          f"{'' if covered else '  (test run did not cover the whole reference)'}")
    return worst_ok and covered


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("log")
    ap.add_argument("--ref-log")
    ap.add_argument("--ref")
    ap.add_argument("--test")
    ap.add_argument("--tol", type=float, default=1e-5,
                    help="max relative difference per column, scaled by the column's max |value|")
    ap.add_argument("--cols", default=None,
                    help="regex: compare only matching CSV columns (e.g. 'global')")
    ap.add_argument("--tmax", type=float, default=None,
                    help="compare only up to this time (e.g. before a bifurcation)")
    ap.add_argument("--atol", type=float, default=1e-6,
                    help="absolute floor: differences below this always pass")
    a = ap.parse_args()

    S = report_log(a.log, parse_log(a.log))
    if a.ref_log:
        Sr = report_log(a.ref_log + "  [reference]", parse_log(a.ref_log))
        if S["its"] and Sr["its"]:
            t_t, t_r = S.get("last_ok"), Sr.get("last_ok")
            same_t = (t_t is not None and t_r is not None
                      and math.isfinite(t_t) and math.isfinite(t_r)
                      and abs(t_t - t_r) <= 1e-9 * max(abs(t_r), 1.0))
            print(f"--- Newton total: test {S['its_all']} vs reference {Sr['its_all']}"
                  f"  (including iterations on failed attempts)")
            if same_t:
                # a ratio only means anything when both runs did the same work
                print(f"    both reached t={t_t:g}:"
                      f"  {Sr['its_all'] / max(S['its_all'], 1):.1f}x fewer")
            else:
                print(f"    NOT COMPARABLE: test reached t={t_t:g}, reference t={t_r:g}."
                      f" A run that stalls early does fewer iterations BECAUSE it failed;"
                      f" that is not a speedup.")
    if a.ref and a.test:
        compare_csv(a.ref, a.test, a.tol, a.atol, a.tmax, a.cols)


if __name__ == "__main__":
    sys.exit(main())
