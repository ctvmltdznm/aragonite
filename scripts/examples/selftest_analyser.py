#!/usr/bin/env python3
"""
selftest_analyser.py -- verify analyse_examples.py against synthetic CSVs.

No MOOSE involved. For each bone flag and load case it writes a CSV that
satisfies the constitutive model EXACTLY by construction:

    elastic rows   sigma = E * epsilon        (or G * gamma for shear)
    plastic rows   sigma chosen so that phi(sigma) == r(kappa) identically

so a correct analyser must report an elastic-slope error at round-off and
|phi/r - 1| at round-off. Anything larger means the analyser has the stress
6-vector order, the slot mapping, the F/f_lin evaluation or the post-yield law
wrong -- all of which are silent failure modes that would otherwise only show up
as a confusing FAIL on real data.

It then perturbs one case by 5% and confirms the analyser actually catches it.

    python3 selftest_analyser.py

Writes into a temporary directory and removes it afterwards.
"""

import csv
import math
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import analyse_examples as A  # noqa: E402


def synth_case(flag, case, n=100, bad=False):
    """Rows for one case, exactly on the yield surface once plastic."""
    rho = A.rho_of(flag)
    F, flin = A.vp.moose_yield(A.BF[flag], A.MAJOR, rho, A.MM)
    C = A.vp.moose_elastic(A.BF[flag], A.MAJOR, rho, A.MM)
    S = np.linalg.inv(C[:3, :3])
    kind, key = case.split("_")
    slot = A.AXIS_SLOT[key]
    shear = kind == "shear"
    sign = -1.0 if kind == "compression" else 1.0

    E = C[slot, slot] if shear else 1.0 / S[slot, slot]
    # sigma on the surface: |sigma|*sqrt(F_ss) + f_s*sigma = r  (f_s = 0 in shear)
    denom = math.sqrt(F[slot, slot]) + (0.0 if shear else sign * flin[slot])
    sy = A.R0 / denom                       # stress at first yield
    ey = sy / E                             # strain at first yield

    rows = []
    for i in range(n + 1):
        t = i / n
        eps = 0.02 * t
        if eps <= ey:
            sig, kap = E * eps, 0.0
        else:
            # let kappa grow with the strain past yield; the exact rate does not
            # matter, only that sigma then sits on the surface at that kappa
            kap = (eps - ey) * 0.9
            sig = float(A.r_of(flag, kap)) / denom
        s6 = [0.0] * 6
        s6[slot] = sign * sig
        if bad and eps > ey:
            s6[slot] *= 1.05            # 5% off the surface
        rows.append(dict(time=t,
                         stress_xx=s6[0], stress_yy=s6[1], stress_zz=s6[2],
                         stress_yz=s6[3], stress_xz=s6[4], stress_xy=s6[5],
                         **{f"strain_{key}": sign * eps / (2.0 if shear else 1.0)},
                         plastic_strain=kap))
    return rows


D0N, D0T, ETA = 1.91e-4, 2.17e-4, 0.25


def _mixed_delta0(phi):
    if phi < 1e-12:
        return D0N
    if phi > 1e6:
        return D0T
    return D0N * D0T * math.sqrt(1 + phi * phi) / math.sqrt(D0T ** 2 + phi ** 2 * D0N ** 2)


def synth_coral(case, n=100):
    """Synthetic interface response that obeys the model's own rules: the peak
    follows the elliptic interaction at the mode mix, and damage follows the
    damage law at the opening reached. An earlier version ramped damage
    linearly to 1, which the analyser correctly rejected -- the law gives
    0.54 to 0.70 over these openings, not 1."""
    spec = A.CORAL_IFACE[case]
    SN, TS = 626.0, 374.0
    rows = []
    for i in range(n + 1):
        t = i / n
        amp = max(0.0, math.sin(math.pi * t ** 0.7))
        if case.endswith("mixed_mode"):
            jn = jt = 4.3e-4 * t
        elif spec["col"] == "normal_traction":
            jn, jt = 6e-4 * t, 0.0
        else:
            jn, jt = 0.0, 6.5e-4 * t

        dtot = math.hypot(jn, jt)
        if dtot > 1e-20:
            rn, rt = abs(jn) / dtot, abs(jt) / dtot
        else:
            rn, rt = 1.0, 0.0
        tpk = math.sqrt((SN * rn) ** 2 + (TS * rt) ** 2)
        tn, tt = tpk * amp * rn, tpk * amp * rt

        deff = dtot
        phi = jt / max(jn, 1e-6 * D0N)
        d0 = _mixed_delta0(phi)
        rk = deff / d0
        dmg = 0.0 if rk <= 1.0 else 1.0 - rk ** (ETA - 1.0) * math.exp(1.0 - rk ** ETA)

        rows.append(dict(time=t, normal_traction=tn, tangent_traction=tt,
                         normal_jump=jn, tangent_jump=jt,
                         interface_damage=dmg, interface_delta_eff=deff,
                         stress_xx=0.0, stress_xy=0.0, plastic_strain=0.0))
    return rows


def write_csv(path, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def build(root, bad_case=None):
    cases = [f"{k}_{a}" for k in ("tension", "compression") for a in ("xx", "yy", "zz")] + \
            [f"shear_{p}" for p in ("xy", "xz", "yz")]
    for flag in A.BF:
        for case in cases:
            d = os.path.join(root, "bone", flag)
            os.makedirs(d, exist_ok=True)
            open(os.path.join(d, f"{case}.i"), "w").close()      # analyser globs *.i
            write_csv(os.path.join(d, f"{case}_out.csv"),
                      synth_case(flag, case, bad=(bad_case == (flag, case))))
    for case in A.CORAL_IFACE:
        d = os.path.join(root, "coral")
        os.makedirs(d, exist_ok=True)
        open(os.path.join(d, f"{case}.i"), "w").close()
        write_csv(os.path.join(d, f"{case}_out.csv"), synth_coral(case))


def run(root, *extra):
    r = subprocess.run([sys.executable, os.path.join(HERE, "analyse_examples.py"),
                        "--root", root, *extra],
                       capture_output=True, text=True)
    return r.stdout + r.stderr


def main():
    ok = True

    root = tempfile.mkdtemp(prefix="selftest_clean_")
    try:
        build(root)
        out = run(root)
        worst_y = worst_e = 0.0
        for ln in out.splitlines():
            # the trailing headroom table repeats the same labels with the
            # tolerance in the second-to-last field, so stop at it
            if ln.startswith("headroom"):
                break
            p = ln.split()
            if not p or p[-1] not in ("ok", "FAIL", "info"):
                continue
            if "phi/r" in ln:
                worst_y = max(worst_y, float(p[-2]))
            if "elastic" in ln and "modulus" in ln:
                worst_e = max(worst_e, float(p[-2]))
        nfail = int([l for l in out.splitlines() if "passed," in l][0].split()[2])
        print(f"clean synthetic data:")
        print(f"   worst |phi/r - 1|        {worst_y:.3e}   "
              + ("OK" if worst_y < 1e-10 else "TOO LARGE -- analyser is wrong"))
        print(f"   worst elastic-slope rel  {worst_e:.3e}   "
              + ("OK" if worst_e < 1e-10 else "TOO LARGE -- analyser is wrong"))
        print(f"   reported failures        {nfail}   "
              + ("OK" if nfail == 0 else "should be 0"))
        ok &= worst_y < 1e-10 and worst_e < 1e-10 and nfail == 0
    finally:
        shutil.rmtree(root, ignore_errors=True)

    root = tempfile.mkdtemp(prefix="selftest_bad_")
    try:
        bad = ("trabecular_fabric_ortho", "shear_xz")
        build(root, bad_case=bad)
        out = run(root)
        line = [l for l in out.splitlines()
                if f"{bad[0]}/{bad[1]}" in l and "phi/r" in l]
        caught = bool(line) and "FAIL" in line[0]
        print(f"injected 5% error in {bad[0]}/{bad[1]}:")
        print(f"   caught                   {caught}   " + ("OK" if caught else "MISSED"))
        if line:
            print(f"   {line[0].strip()}")
        ok &= caught
    finally:
        shutil.rmtree(root, ignore_errors=True)

    print("\nSELFTEST", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
