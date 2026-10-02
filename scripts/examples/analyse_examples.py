#!/usr/bin/env python3
"""
analyse_examples.py -- check every example CSV against values predicted from the
presets, not from a previous run.

    python3 analyse_examples.py                 # everything
    python3 analyse_examples.py --filter coral
    python3 analyse_examples.py --gold          # also diff against gold/ if present

WHAT IS CHECKED, and why these checks and not others
----------------------------------------------------
1. ELASTIC SLOPE. The pre-yield slope of stress against strain, compared with
   E_i = 1/(C^-1)_ii for the normal cases and G_ij for the shear cases. This is
   an independent check of the elasticity side: it uses no plasticity at all.

2. YIELD-SURFACE CONSISTENCY -- the main quantitative test.
   At EVERY row where the material is plastic, the stress must sit on the
   current yield surface:
       phi(sigma) = sqrt(s:F:s) + f_lin.s  ==  r(kappa)
   Both sides come from the CSV: the six stress components and kappa (the
   effective_plastic_strain column). F, f_lin and r(.) come from the preset
   formulas. Nothing here depends on how far the run went, on the time step, or
   on a stored reference, which is what makes it the strongest of the three.
   A violation means the return map converged to a point off the surface, or
   the preset resolved to the wrong constants.

3. FIRST YIELD and PEAK, reported rather than asserted.
   First yield is where kappa leaves zero, linearly interpolated between the
   bracketing rows, compared with r0 x strength. The raw first plastic row
   overshoots by up to one elastic increment -- about 4% of the yield stress at
   these settings -- so the interpolated value is the honest one and even that
   is a discretisation estimate. The peak is compared with r(kappa) evaluated at
   the row where the peak occurs, since the peak is wherever r happens to be
   largest within the strain range, not a material constant.

   Note that peak != strength for the bone presets. UMAT strengths are ULTIMATE
   strengths and r(0) = RDY = 0.7, and the runs stop at 2% strain, which puts
   r near 0.93-0.94 rather than 1. See MATERIAL_MODEL_THEORY.md section 4.

TOLERANCES are set from a clean baseline, not guessed, and the run prints the
headroom it actually used, so a FAIL is a regression rather than a marginal
default. The floor on check 2 is set by FINITE strain: the CSV reports Cauchy
stress while the yield surface is evaluated on the same measure the return map
uses, and at 2% strain the two differ in the third digit.
"""

import argparse
import csv
import glob
import importlib.util
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# Reuse the preset formulas from verify_presets.py. It is a script, not a
# module, so its top-level output is suppressed during import.
# ---------------------------------------------------------------------------
_vp_path = os.path.join(HERE, "verify_presets.py")
if not os.path.exists(_vp_path):
    sys.exit("verify_presets.py must sit next to this script (it holds the preset formulas)")
_spec = importlib.util.spec_from_file_location("vp", _vp_path)
vp = importlib.util.module_from_spec(_spec)


class _Quiet:
    def write(self, *a):
        pass

    def flush(self):
        pass


_stdout, sys.stdout = sys.stdout, _Quiet()
try:
    _spec.loader.exec_module(vp)
finally:
    sys.stdout = _stdout

BF = {"trabecular_iso": 0, "trabecular_ti": 1, "trabecular_fabric_ti": 2,
      "trabecular_fabric_ortho": 3, "compact_iso": 4, "compact_ti": 5}
MM = np.array([0.85, 0.95, 1.20])
MAJOR = 3
R0, KSLOPE, KMAX, KWIDTH = 0.7, 100.0, 0.015, 8.0
# _F_matrix / stress 6-vector order: xx yy zz yz xz xy
SIX = ["stress_xx", "stress_yy", "stress_zz", "stress_yz", "stress_xz", "stress_xy"]
AXIS_SLOT = {"xx": 0, "yy": 1, "zz": 2, "yz": 3, "xz": 4, "xy": 5}

# Coral interface, from the elliptic interaction and the Wang 2025 mixed-mode
# form. Same numbers as the file headers.
CORAL_IFACE = {
    "coral_czm_mode_I_opening":  dict(peak=626.0, col="normal_traction",  quiet="tangent_traction"),
    "coral_czm_mode_II_shear":   dict(peak=374.0, col="tangent_traction", quiet="normal_traction"),
    "coral_czm_mixed_mode":      dict(peak=515.6, col=None,               quiet=None),
}


def rho_of(flag):
    return 0.90 if flag.startswith("compact") else 0.20


def r_of(flag, k):
    """Post-yield strength ratio for the preset attached to this flag."""
    k = np.asarray(k, dtype=float)
    if flag.startswith("compact"):          # exp_hardening, UMAT PYFL=1
        return R0 + (1.0 - R0) * (1.0 - np.exp(-KSLOPE * k))
    # simple_softening, UMAT PYFL=2
    return R0 + (1.0 - R0) * (np.exp(-((k - KMAX) ** 2) / (KWIDTH * KMAX ** 2))
                              - np.exp(-1.0 / KWIDTH - KSLOPE * k))


def read_csv(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        return {}
    out = {}
    for key in rows[0]:
        try:
            out[key] = np.array([float(r[key]) for r in rows])
        except (ValueError, TypeError):
            pass
    return out


def fmt(v, w=10, p=4):
    return f"{v:{w}.{p}g}"


class Report:
    def __init__(self):
        self.rows, self.n_fail, self.n_pass, self.missing = [], 0, 0, []
        self.notes = []
        self.worst = {}          # label -> (worst rel/residual seen, tolerance)

    def add(self, name, label, value, target, tol, unit="", info=False, bound=False):
        """bound=True marks a RESIDUAL: there is no predicted value, only a
        tolerance it has to stay under. Printing 0 in the 'predicted' column for
        those reads as a prediction that the yield surface is zero, which is not
        what it means."""
        if bound:
            # a residual: compare the value itself against the bound. It must be
            # COUNTED -- marking it info (which a nan 'rel' against a zero target
            # used to do) silently kept 54 yield checks out of the pass tally.
            rel = value
        elif target in (None, 0.0) or not np.isfinite(value):
            rel = float("nan")
        else:
            rel = abs(value - target) / abs(target)
        if info:
            verdict = "info"
        elif not np.isfinite(rel):
            verdict = "info"
        else:
            verdict = "ok" if rel <= tol else "FAIL"
            if verdict == "ok":
                self.n_pass += 1
            else:
                self.n_fail += 1
            w, _ = self.worst.get(label, (0.0, tol))
            self.worst[label] = (max(w, rel), tol)
        self.rows.append((name, label, value, target, rel, unit, verdict,
                          f"<= {tol:g}" if bound else None))

    def print(self):
        if not self.rows:
            print("  nothing to report")
            return
        w = max(len(r[0]) for r in self.rows)
        q = max(len(r[1]) for r in self.rows)
        print(f"  {'case':<{w}}  {'quantity':<{q}} {'measured':>12} {'predicted':>12} "
              f"{'rel':>9}  verdict")
        print("  " + "-" * (w + q + 46))
        for row in self.rows:
            name, label, val, tgt, rel, unit, verdict = row[:7]
            shown = row[7] if len(row) > 7 else None
            t = shown if shown else ("--" if tgt is None else fmt(tgt, 12))
            r = "--" if not np.isfinite(rel) else f"{rel:9.2e}"
            print(f"  {name:<{w}}  {label:<{q}} {fmt(val, 12)} {t:>12} {r:>9}  {verdict}")


def analyse_bone(root, rep, tol_elastic, tol_yield, tol_peak):
    for flag in sorted(BF):
        d = os.path.join(root, "bone", flag)
        if not os.path.isdir(d):
            continue
        rho = rho_of(flag)
        F, flin = vp.moose_yield(BF[flag], MAJOR, rho, MM)
        C = vp.moose_elastic(BF[flag], MAJOR, rho, MM)
        S = np.linalg.inv(C[:3, :3])

        for inp in sorted(glob.glob(os.path.join(d, "*.i"))):
            case = os.path.basename(inp)[:-2]
            csv_path = os.path.join(d, f"{case}_out.csv")
            name = f"{flag}/{case}"
            if not os.path.exists(csv_path):
                rep.missing.append(name)
                continue
            cols = read_csv(csv_path)
            kind, key = case.split("_")
            slot = AXIS_SLOT[key]

            s = np.column_stack([cols[c] for c in SIX])        # n x 6
            kap = cols["plastic_strain"]
            strain = cols.get(f"strain_{key}")

            # ---- 1. elastic slope -------------------------------------------
            el = kap <= 0.0
            if strain is not None and el.sum() >= 3:
                # last few purely elastic points, skipping t = 0
                idx = np.where(el)[0][1:]
                if len(idx) >= 2:
                    num = s[idx, slot]
                    den = strain[idx] * (2.0 if kind == "shear" else 1.0)  # engineering shear
                    good = np.abs(den) > 1e-14
                    if good.sum() >= 2:
                        slope = float(np.polyfit(den[good], num[good], 1)[0])
                        target = C[slot, slot] if kind == "shear" else 1.0 / S[slot, slot]
                        rep.add(name, "elastic modulus", slope, target, tol_elastic, "MPa")

            # ---- 2. yield-surface consistency (the real test) ----------------
            pl = kap > 1e-12
            if pl.sum() >= 2:
                sp = s[pl]
                quad = np.einsum("ij,jk,ik->i", sp, F, sp)
                phi = np.sqrt(np.clip(quad, 0.0, None)) + sp @ flin
                r = r_of(flag, kap[pl])
                resid = np.max(np.abs(phi / r - 1.0))
                # bound=True: add() compares the residual against the bound,
                # counts it, and prints "<= tol" in the predicted column (a
                # literal 0 there read as "the yield surface is predicted to be
                # zero", which is not what a residual target means).
                rep.add(name, "yield consistency |phi/r-1|", resid, 0.0, tol_yield,
                        bound=True)

            # ---- 3. first yield and peak, informational ----------------------
            if pl.sum() and (~pl).sum():
                j = int(np.argmax(pl))                       # first plastic row
                if j > 0:
                    # interpolate where kappa leaves zero
                    k0, k1 = kap[j - 1], kap[j]
                    w = 0.0 if k1 == k0 else (0.0 - k0) / (k1 - k0)
                    sy = s[j - 1, slot] + w * (s[j, slot] - s[j - 1, slot])
                    denom = math.sqrt(F[slot, slot]) + (flin[slot] if kind != "compression"
                                                        else -flin[slot])
                    tgt = R0 / denom if denom > 0 else None
                    rep.add(name, "first yield", abs(sy), tgt, 0.05, "MPa", info=True)

            if pl.sum():
                ip = int(np.argmax(np.abs(s[:, slot])))
                peak = abs(s[ip, slot])
                rk = float(r_of(flag, kap[ip]))
                denom = math.sqrt(F[slot, slot]) + (flin[slot] if kind != "compression"
                                                    else -flin[slot])
                tgt = rk / denom if denom > 0 else None

                # "peak = r x strength" is only true if the state really is the
                # single component this case is meant to probe. Measure that
                # instead of assuming it: compare the slot's own contribution to
                # phi against the total. purity = 1 is a pure state.
                sp1 = s[ip]
                phi_tot = (math.sqrt(max(float(sp1 @ F @ sp1), 0.0)) + float(flin @ sp1))
                slot_only = np.zeros(6)
                slot_only[slot] = sp1[slot]
                phi_slot = (math.sqrt(max(float(slot_only @ F @ slot_only), 0.0))
                            + float(flin @ slot_only))
                purity = phi_slot / phi_tot if phi_tot else float("nan")
                pure = np.isfinite(purity) and abs(purity - 1.0) < 1e-2

                rep.add(name, "peak stress", peak, tgt, tol_peak, "MPa", info=not pure)
                if not pure:
                    rep.add(name, "  purity at peak", purity, 1.0, 0.0, "", info=True)
                    # the components that are carrying the rest, largest first
                    off = [(abs(sp1[k]), SIX[k], sp1[k]) for k in range(6) if k != slot]
                    off.sort(reverse=True)
                    big = ", ".join(f"{nm.replace('stress_','')}={v:+.4g}"
                                    for _, nm, v in off[:3] if abs(v) > 1e-12)
                    rep.notes.append(
                        f"{name}: peak row t-index {ip}, {SIX[slot].replace('stress_','')}"
                        f"={sp1[slot]:+.5g}, kappa={kap[ip]:.5g}, phi={phi_tot:.5g}, "
                        f"r={rk:.5g}; other components: {big or 'none'}")


def analyse_coral(root, rep, tol_peak, tol_quiet, tol_damage=0.2):
    """Checks against the model's own two mixing rules, evaluated at the mode
    mix the run actually took, rather than against the nominal pure-mode
    numbers alone. Both are reported."""
    SN, TS = 626.0, 374.0            # normal and shear strengths
    D0N, D0T = 1.91e-4, 2.17e-4      # characteristic openings
    ETA = 0.25                       # softening exponent, for the damage law

    def mixed_delta0(phi):
        if phi < 1e-12:
            return D0N
        if phi > 1e6:
            return D0T
        return D0N * D0T * math.sqrt(1 + phi * phi) / math.sqrt(D0T ** 2 + phi ** 2 * D0N ** 2)

    d = os.path.join(root, "coral")
    if not os.path.isdir(d):
        return
    for inp in sorted(glob.glob(os.path.join(d, "*.i"))):
        case = os.path.basename(inp)[:-2]
        csv_path = os.path.join(d, f"{case}_out.csv")
        name = f"coral/{case}"
        if not os.path.exists(csv_path):
            rep.missing.append(name)
            continue
        cols = read_csv(csv_path)
        spec = CORAL_IFACE.get(case)
        if spec is None:
            continue
        tn, tt = cols.get("normal_traction"), cols.get("tangent_traction")
        jn, jt = cols.get("normal_jump"), cols.get("tangent_jump")
        if tn is None or tt is None:
            rep.missing.append(name + " (traction columns)")
            continue

        res = np.sqrt(tn ** 2 + tt ** 2)
        ip = int(np.argmax(res))

        # Elliptic interaction at the mode mix the run actually reached. This
        # tests the mixing rule itself, and unlike the nominal number it stays
        # valid when the path is not exactly 45 degrees.
        if jn is not None and jt is not None:
            dtot = math.hypot(float(jn[ip]), float(jt[ip]))
            if dtot > 1e-20:
                rn, rt = abs(float(jn[ip])) / dtot, abs(float(jt[ip])) / dtot
                t_pred = math.sqrt((SN * rn) ** 2 + (TS * rt) ** 2)
                rep.add(name, "peak |T| at observed mix", float(res[ip]), t_pred, tol_peak, "MPa")
                # the traction follows the jump direction, so the components
                # are the resultant resolved along it
                if rn > 1e-3:
                    rep.add(name, "  peak T_normal", float(tn[ip]), t_pred * rn, tol_peak, "MPa")
                if rt > 1e-3:
                    rep.add(name, "  peak T_tangent", float(tt[ip]), t_pred * rt, tol_peak, "MPa")

        rep.add(name, "peak |T| vs nominal", float(res[ip]), spec["peak"], tol_peak, "MPa",
                info=True)

        if spec["col"] is not None:
            off = float(np.max(np.abs(cols[spec["quiet"]])))
            ref = float(np.max(np.abs(cols[spec["col"]]))) or 1.0
            rep.rows.append((name, "off-mode / loaded", off / ref, None,
                             off / ref, "", "ok" if off / ref <= tol_quiet else "FAIL",
                             f"<= {tol_quiet:g}"))
            if off / ref <= tol_quiet:
                rep.n_pass += 1
            else:
                rep.n_fail += 1
        elif jn is not None and jt is not None:
            m = np.abs(jn) > 1e-12
            if m.sum():
                rep.add(name, "tangent/normal jump", float(np.median(jt[m] / jn[m])), 1.0, 0.05)

        dmg = cols.get("interface_damage")
        if dmg is not None and len(dmg) > 1:
            drops = float(np.min(np.diff(dmg)))
            rep.rows.append((name, "damage min increment", drops, None,
                             abs(min(drops, 0.0)), "", "ok" if drops >= -1e-9 else "FAIL",
                             ">= 0"))
            if drops >= -1e-9:
                rep.n_pass += 1
            else:
                rep.n_fail += 1

            # Predict the final damage from the damage law at the opening the
            # run actually reached, instead of assuming it fails completely.
            # D = 1 - r^(eta-1) exp(1 - r^eta),  r = delta_eff / delta_0(mix).
            de = cols.get("interface_delta_eff")
            if de is not None and jn is not None and jt is not None:
                dn_end, dt_end = abs(float(jn[-1])), abs(float(jt[-1]))
                phi_mix = dt_end / max(dn_end, 1e-6 * D0N)
                d0 = mixed_delta0(phi_mix)
                rk = float(de[-1]) / d0
                if rk > 1.0:
                    d_pred = 1.0 - rk ** (ETA - 1.0) * math.exp(1.0 - rk ** ETA)
                    rep.add(name, "final damage", float(dmg[-1]), d_pred, tol_damage)
                    rep.notes.append(
                        f"{name}: delta_eff/delta_0 = {rk:.3f} at the end "
                        f"(delta_0 mix = {d0:.3e} mm), so the interface is NOT expected "
                        f"to fail completely; the damage law gives {d_pred:.3f}.")
                else:
                    rep.add(name, "final damage", float(dmg[-1]), None, 0.0, "", info=True)


def diff_gold(root, rep_filter):
    """Optional regression diff against gold/<case>_out.csv, when present."""
    found = 0
    worst = 0.0
    for gold in sorted(glob.glob(os.path.join(root, "*", "*", "gold", "*_out.csv"))) + \
                sorted(glob.glob(os.path.join(root, "*", "gold", "*_out.csv"))):
        live = os.path.join(os.path.dirname(os.path.dirname(gold)), os.path.basename(gold))
        if rep_filter and rep_filter not in gold:
            continue
        if not os.path.exists(live):
            continue
        g, l = read_csv(gold), read_csv(live)
        found += 1
        for k in g:
            if k not in l or len(g[k]) != len(l[k]):
                print(f"  {live}: column '{k}' missing or different length")
                worst = float("inf")
                continue
            scale = max(np.max(np.abs(g[k])), 1e-30)
            worst = max(worst, float(np.max(np.abs(g[k] - l[k])) / scale))
    if found:
        print(f"\nGold regression: {found} file(s) compared, worst relative difference "
              f"{worst:.2e}  " + ("ok" if worst < 1e-10 else "CHANGED"))
    else:
        print("\nGold regression: no gold/ files found (create them with "
              "GOLD=1 ./run_examples.sh)")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=HERE, help="examples folder (default: alongside this script)")
    ap.add_argument("--filter", default="", help="substring filter on the case path")
    ap.add_argument("--gold", action="store_true", help="also diff against gold/ if present")
    ap.add_argument("--tol-elastic", type=float, default=2e-2,
                    help="relative tolerance on the elastic slope (default 2e-2)")
    ap.add_argument("--tol-yield", type=float, default=2e-3,
                    help="absolute tolerance on |phi/r - 1| (default 2e-3, about "
                         "4x the 4.5e-4 worst case on a clean baseline)")
    ap.add_argument("--tol-peak", type=float, default=3e-2,
                    help="relative tolerance on peak stress / traction (default 3e-2)")
    ap.add_argument("--tol-quiet", type=float, default=2e-2,
                    help="max off-mode traction as a fraction of the loaded one (default 2e-2). "
                         "The smooth Macaulay bracket leaves a small positive normal gap even "
                         "when the normal jump is held at zero, so this cannot be driven to "
                         "round-off; 1e-3 was too tight and failed mode II at 1.2e-2.")
    a = ap.parse_args()

    rep = Report()
    analyse_bone(a.root, rep, a.tol_elastic, a.tol_yield, a.tol_peak)
    analyse_coral(a.root, rep, a.tol_peak, a.tol_quiet)

    if a.filter:
        rep.rows = [r for r in rep.rows if a.filter in r[0]]

    print("=" * 78)
    print("EXAMPLE ANALYSIS -- measured against values predicted from the presets")
    print("=" * 78)
    rep.print()
    print()
    if rep.missing:
        print(f"{len(rep.missing)} case(s) have no CSV yet (run them first):")
        for m in rep.missing[:10]:
            print(f"    {m}")
        if len(rep.missing) > 10:
            print(f"    ... and {len(rep.missing) - 10} more")
        print()
    if rep.notes:
        print("Diagnostics for the rows that could not be asserted:")
        for n in rep.notes:
            print(f"    {n}")
        print()
    print(f"{rep.n_pass} passed, {rep.n_fail} failed"
          + ("   VERDICT: PASS" if rep.n_fail == 0 and rep.n_pass else
             "   VERDICT: FAIL" if rep.n_fail else "   VERDICT: nothing checked"))
    if rep.worst:
        print()
        print("headroom (worst observed vs tolerance):")
        for label, (w, tol) in sorted((k.strip(), v) for k, v in rep.worst.items()):
            print(f"  {label:<34} {w:9.2e}  /  {tol:7.1e}"
                  f"   {'x%.0f' % (tol / w) if w > 0 else 'exact':>8}")

    if a.gold:
        diff_gold(a.root, a.filter)

    return 1 if rep.n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
