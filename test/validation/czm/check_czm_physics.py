#!/usr/bin/env python3
"""
check_czm_physics.py -- judge CZM validation output against the MODEL,
not against an old gold file.

A gold file only says "the same as last time". These four tests were
"the same as last time" for three months while the traction law applied the
softening twice, and then reported "Differs from gold" for a 1000x physics
error. This script asks a different question: does every row of the output
satisfy the constitutive law, given the opening it reports?

Requires the deterministic configuration the regression inputs now use:
quality_std_dev = 0, spatial_quality_std_dev = 0, damage_viscosity = 0.
Then every interface QP is the single contact f = 1 and the law is exact:

    r = delta_eff / delta_0(phi)                     (Wang 2025 mixed-mode delta_0)
    r <  1 :  T = T_peak * e * r^mu  * exp(-r^mu)     (loading branch)
    r >= 1 :  T = T_peak * r * (1 - s(r) * D_law(r))  (secant branch)
              D_law = 1 - e * r^(eta-1) * exp(-r^eta)
              s     = 0.5 * (1 + tanh((r - 1) / eps_act))
    T_peak  = sqrt((sigma_n*rn)^2 + (tau_s*rs)^2 + (tau_t*rt)^2)   (elliptic)

With s = 1 the secant branch equals the exponential envelope exactly,
because 1 - D_law = e*r^(eta-1)*exp(-r^eta). Only the tanh activation near
r = 1 separates them, and the prediction includes it.

Usage (repo root, after running the suite):
    python3 test/validation/czm/check_czm_physics.py            # *_out.csv
    python3 test/validation/czm/check_czm_physics.py --gold     # gold/*.csv
Exit status 0 = every check passed.
"""
import argparse
import csv
import glob
import math
import os
import re
import sys

E = math.e

# Solver-health limits (measured after the yield-stress scaling: mean 3.5 / max 5
# in mixed_mode, 2.2 / 3 in mode_ii; the oscillating runs had 10.9 / 17 and 13.6 / 99).
NEWTON_MEAN = 6.0
NEWTON_MAX = 10.0
HERE = os.path.dirname(os.path.abspath(__file__))


# ---------------------------------------------------------------- parameters
def czm_params(inp):
    """Read the HomogenizedExponentialCZM block of an input file."""
    txt = open(inp).read()
    m = re.search(r"type\s*=\s*HomogenizedExponentialCZM(.*?)\n\s*\[\]", txt, re.S)
    if not m:
        return None
    blk = re.sub(r"#[^\n]*", "", m.group(1))

    def get(name, default=None):
        mm = re.search(r"^\s*" + name + r"\s*=\s*([-+0-9.eE]+)", blk, re.M)
        if mm:
            return float(mm.group(1))
        if default is None:
            raise KeyError(f"{name} not set in {inp}")
        return default

    p = dict(SN=get("normal_strength"), TS=get("shear_strength_s"),
             TT=get("shear_strength_t"), D0N=get("delta_0_normal"),
             D0T=get("delta_0_tangent"), MU=get("mu", 1.0), ETA=get("eta", 1.0),
             QSD=get("quality_std_dev", 0.1),
             SQSD=get("spatial_quality_std_dev", 0.15),
             VISC=get("damage_viscosity", 0.0),
             EPS_ACT=get("damage_activation_width", 0.02))
    p["EPS_MAC"] = get("macaulay_regularisation",
                       1e-3 * min(p["D0N"], p["D0T"]))
    p["DC"] = get("delta_c", 5000.0 * max(p["D0N"], p["D0T"]))   # code's clamped default
    return p


# --------------------------------------------------------------------- model
def delta0_mixed(phi, dn, dt):
    if phi < 1e-12:
        return dn
    if phi > 1e6:
        return dt
    return dn * dt * math.sqrt(1 + phi * phi) / math.sqrt(dt * dt + phi * phi * dn * dn)


def traction(p, d_n, d_s=0.0, d_t=0.0):
    """Predicted (|T_cohesive|, T_normal_total, D) for one deterministic QP."""
    eps = p["EPS_MAC"]
    sr = math.sqrt(d_n * d_n + eps * eps)
    dn_pos, dn_neg = 0.5 * (d_n + sr), 0.5 * (sr - d_n)
    dtm = math.hypot(d_s, d_t)
    deff = math.hypot(dn_pos, dtm)
    if deff <= 1e-12:
        return 0.0, -E * p["SN"] / p["D0N"] * dn_neg, 0.0
    phi = dtm / max(dn_pos, 1e-6 * p["D0N"])
    d0 = delta0_mixed(phi, p["D0N"], p["D0T"])
    dtot = math.sqrt(dn_pos ** 2 + d_s ** 2 + d_t ** 2)
    rn, rs, rt = dn_pos / dtot, abs(d_s) / dtot, abs(d_t) / dtot
    tpk = math.sqrt((p["SN"] * rn) ** 2 + (p["TS"] * rs) ** 2 + (p["TT"] * rt) ** 2)
    r = deff / d0
    if r < 1.0:
        ex = r ** p["MU"]
        T, D = tpk * E * ex * math.exp(-ex), 0.0
    else:
        s = 0.5 * (1 + math.tanh((r - 1) / p["EPS_ACT"]))
        re_ = r ** p["ETA"]
        dlaw = 1 - math.exp(1 - re_) * re_ / r
        D = max(0.0, s * dlaw)
        T = tpk * r * (1 - D)
    t_norm = T * dn_pos / deff - E * p["SN"] / p["D0N"] * dn_neg
    return T, t_norm, D


def peak_factor(eta, eps_act):
    """max over r >= 1 of r*(1 - s(r)*D_law(r)): the secant branch's peak / T_peak."""
    best, r = 1.0, 1.0
    while r < 1.0 + 20 * eps_act:
        s = 0.5 * (1 + math.tanh((r - 1) / eps_act))
        re_ = r ** eta
        best = max(best, r * (1 - s * (1 - math.exp(1 - re_) * re_ / r)))
        r += eps_act * 1e-3
    return best


class QP:
    """One deterministic interface QP WITH history, mirroring the C++:
    kappa = running max of delta_eff (irreversibility), D = max(D_target, D_old),
    secant branch T_peak*delta_eff/d0 once kappa >= d0 (so unloading returns
    linearly to the origin with D frozen), and D = 1 once kappa >= delta_c."""

    def __init__(self, p):
        self.p, self.kappa, self.D = p, 0.0, 0.0

    def step(self, d_n, d_s=0.0, d_t=0.0):
        p = self.p
        eps = p["EPS_MAC"]
        sr = math.sqrt(d_n * d_n + eps * eps)
        dn_pos, dn_neg = 0.5 * (d_n + sr), 0.5 * (sr - d_n)
        dtm = math.hypot(d_s, d_t)
        deff = math.hypot(dn_pos, dtm)
        pen = E * p["SN"] / p["D0N"] * dn_neg
        if deff <= 1e-12:
            return dict(T=0.0, Tn=-pen, Tt=0.0, D=self.D)
        phi = dtm / max(dn_pos, 1e-6 * p["D0N"])
        d0 = delta0_mixed(phi, p["D0N"], p["D0T"])
        dtot = math.sqrt(dn_pos ** 2 + d_s ** 2 + d_t ** 2)
        rn, rs, rt = dn_pos / dtot, abs(d_s) / dtot, abs(d_t) / dtot
        tpk = math.sqrt((p["SN"] * rn) ** 2 + (p["TS"] * rs) ** 2 + (p["TT"] * rt) ** 2)
        r = deff / d0
        ex = r ** (p["MU"] if r < 1.0 else p["ETA"])
        T_und = tpk * E * ex * math.exp(-ex)
        self.kappa = max(self.kappa, deff)
        T_base = tpk * deff / d0 if self.kappa >= d0 else T_und
        if self.kappa >= p["DC"]:
            Dt = 1.0
        else:
            rk = self.kappa / d0
            s = 0.5 * (1 + math.tanh((rk - 1) / p["EPS_ACT"]))
            re_ = rk ** p["ETA"]
            Dt = s * (1 - math.exp(1 - re_) * re_ / rk)
        self.D = 1.0 if self.kappa >= p["DC"] else max(Dt, self.D, 0.0)
        T = T_base * (1 - self.D)
        return dict(T=T, Tn=T * dn_pos / deff - pen, Tt=T * dtm / deff, D=self.D)


# -------------------------------------------------------------------- checks
class Report:
    def __init__(self):
        self.fail = 0
        self.lines = []

    def add(self, test, what, val, tol, ok, note=""):
        self.fail += 0 if ok else 1
        self.lines.append(f"  {test:20s} {what:36s} {val:10.3e}  <= {tol:7.0e}  "
                          f"{'ok  ' if ok else 'FAIL'} {note}")


def rows(path):
    return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(open(path))]


def worst(errs):
    return max(errs) if errs else 0.0


def check_test(name, inp, csvp, rep, tol):
    p = czm_params(inp)
    R = rows(csvp)
    det = p["QSD"] == 0 and p["SQSD"] == 0 and p["VISC"] == 0

    # 1. ceiling: no traction column may exceed the law's own supremum.
    # NOT simply the strength: just above r = 1 the tanh activation s is still
    # < 1, so damage trails D_law while the secant ramp T_peak*r keeps rising.
    # The peak therefore overshoots T_peak by ~0.14*eps_act (+0.28% at the
    # default 0.02), independent of eta. The pre-secant formulation could not.
    pf = peak_factor(p["ETA"], p["EPS_ACT"])
    ceil = max(p["SN"], p["TS"], p["TT"]) * pf * (1 + 1e-6)
    if not det:
        ceil *= 1 + 0.1 * 2.857 + 3 * p["SQSD"]       # stochastic: generous bound
    tcols = [c for c in R[0] if re.search(r"traction|^T_", c)]
    tmax = max(abs(r[c]) for r in R for c in tcols)
    rep.add(name, f"max traction / ceiling ({ceil:.1f} MPa)", tmax / ceil, 1.0,
            tmax <= ceil, f"max = {tmax:.2f} MPa  (activation overshoot {100*(pf-1):+.2f}%)")

    # 1b. solver health. Rows are CONVERGED states, so the law checks below
    # cannot see how Newton got there: mixed_mode once passed them exactly
    # while every step oscillated between |R| 2.4e3 and 3.4e3 for ~10
    # iterations (bulk return map switching on in the iterates). A healthy
    # run of these 48-DOF problems needs a handful of iterations per step.
    if "newton_its" in R[0]:
        its = [r["newton_its"] for r in R if r["time"] > 0]
        if its:
            mean, mx = sum(its) / len(its), max(its)
            rep.add(name, "Newton iterations per step, mean", mean, NEWTON_MEAN,
                    mean <= NEWTON_MEAN, f"{len(its)} steps")
            rep.add(name, "Newton iterations per step, max", mx, NEWTON_MAX,
                    mx <= NEWTON_MAX)

    if not det:
        rep.add(name, "exact law check", 0.0, 0.0, True,
                "SKIPPED: stochastic or viscous input")
        return

    # 2. exact row-wise law WITH history, wherever the CSV gives the kinematics.
    # The tangential strength is ambiguous: which local tangent (s or t) a
    # global y-shear lands on is MOOSE's choice of interface frame. Both are
    # tried and the better one reported, so a mismatch cannot hide behind it.
    errT, errD = [], []

    def run(p_, kin, meas_T, meas_D):
        q, eT, eD = QP(p_), [], []
        for r in R:
            k = kin(r)
            if k is None:
                continue
            out = q.step(*k)
            m = meas_T(r, out)
            if m is not None:
                pred, val = m
                eT.append(abs(val - pred) / max(abs(pred), 1.0))
            if meas_D and meas_D in r:
                eD.append(abs(r[meas_D] - out["D"]))
        return eT, eD

    if "avg_traction" in R[0] and "avg_separation" in R[0]:          # mode I (any path)
        errT, errD = run(p, lambda r: (r["avg_separation"],),
                         lambda r, o: (o["Tn"], r["avg_traction"]), "avg_damage")
        first = next((r for r in R if r["time"] > 0), None)
        if first:
            rep.add(name, "first-step normal jump >= 0 (tension)",
                    max(0.0, -first["avg_separation"]), 1e-12,
                    first["avg_separation"] >= -1e-12,
                    f"avg_separation = {first['avg_separation']:+.3e}"
                    "  (was -4.35e-07 before the factor-e fix)")
    else:
        if "T_mag" in R[0] and "delta_n" in R[0]:                    # mixed mode
            kin = lambda r: (r["delta_n"], r["delta_s"])
            meas = lambda r, o: (math.hypot(o["Tn"], o["Tt"]), r["T_mag"])
        elif "avg_tangent_traction" in R[0]:                         # mode II
            kin = lambda r: (r.get("avg_normal_jump", 0.0), r["avg_tangent_jump"])
            meas = lambda r, o: (o["Tt"], abs(r["avg_tangent_traction"]))
        else:
            kin = None
        if kin:
            best = None
            for lab, ts in (("tau=shear_strength_s", p["TS"]), ("tau=shear_strength_t", p["TT"])):
                pp = dict(p, TS=ts)
                eT, eD = run(pp, kin, meas, "avg_damage")
                if eT and (best is None or worst(eT) < worst(best[1])):
                    best = (lab, eT, eD)
            if best:
                name = f"{name}[{best[0]}]"
                errT, errD = best[1], best[2]

    if errT:
        rep.add(name, "traction vs law, worst row (rel)", worst(errT), tol,
                worst(errT) <= tol, f"{len(errT)} rows")
    if errD:
        rep.add(name, "damage vs law, worst row (abs)", worst(errD), tol,
                worst(errD) <= tol, f"{len(errD)} rows")
    if not errT and not errD:
        rep.add(name, "exact law check", 0.0, 0.0, True,
                "no usable kinematic columns -- ceiling only")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gold", action="store_true", help="check gold/*.csv instead of *_out.csv")
    ap.add_argument("--tol", type=float, default=1e-3)
    a = ap.parse_args()
    rep = Report()
    for inp in sorted(glob.glob(os.path.join(HERE, "*", "*.i"))):
        d, base = os.path.dirname(inp), os.path.splitext(os.path.basename(inp))[0]
        csvp = os.path.join(d, "gold", f"{base}_out.csv") if a.gold else os.path.join(d, f"{base}_out.csv")
        if not os.path.exists(csvp):
            rep.add(base, "output present", 1, 0, False, f"missing {os.path.relpath(csvp)}")
            continue
        check_test(base, inp, csvp, rep, a.tol)
    print("\n".join(rep.lines))
    print(f"\n{'PASS' if rep.fail == 0 else 'FAIL'}: {rep.fail} failed check(s)")
    return 1 if rep.fail else 0


if __name__ == "__main__":
    sys.exit(main())
