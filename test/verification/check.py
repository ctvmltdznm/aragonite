"""Simulate the patched runDebugChecks softening check against every post-yield
mode, plus a deliberately broken derivative, to confirm it passes correct code
and still catches real errors."""
import math

def make(mode, r0=1.0, h_=0.0, res=0.7, kslope=30.0, kmax=0.001, kmin=0.02, kwidth=8.0, bug=False):
    def r(k):
        if mode == "perfect":            return 1.0
        if mode == "exp_hardening":      return r0 + (1-r0+h_)*(1-math.exp(-kslope*k))
        if mode == "linear_hardening":   return 1.0 + kslope*k
        if mode == "simple_softening":
            return r0 + (1-r0)*(math.exp(-((k-kmax)**2)/(kwidth*kmax**2))
                                - math.exp(-1/kwidth - kslope*k))
        if mode == "exp_softening":
            return 1.0 if k < kmax else res + (1-res)*math.exp(-kslope*(k-kmax))
        if mode == "piecewise_softening":
            if k < kmax: return 1.0
            if k < kmin: return 1.0 - (1-res)*(k-kmax)/(kmin-kmax)
            return res
    def dr(k):
        f = 1.7 if bug else 1.0          # deliberate error
        if mode == "perfect":            return 0.0
        if mode == "exp_hardening":      return f*(1-r0+h_)*kslope*math.exp(-kslope*k)
        if mode == "linear_hardening":   return f*kslope
        if mode == "simple_softening":
            return f*(1-r0)*(-2*(k-kmax)/(kwidth*kmax**2)*math.exp(-((k-kmax)**2)/(kwidth*kmax**2))
                             + kslope*math.exp(-1/kwidth - kslope*k))
        if mode == "exp_softening":
            return 0.0 if k < kmax else -f*(1-res)*kslope*math.exp(-kslope*(k-kmax))
        if mode == "piecewise_softening":
            if k < kmax: return 0.0
            if k < kmin: return -f*(1-res)/(kmin-kmax)
            return 0.0
    return r, dr, kmax, kmin

def check(r, dr, kmax, kmin):
    k_ref = kmax if kmax > 0 else 0.01
    rmax, kinks, first = 0.0, 0, 0.0
    for k in (0.2*k_ref, 0.9*k_ref, 1.0*k_ref, 1.5*k_ref, 4.0*k_ref):
        h = 1e-6*k_ref
        r0v = r(k); fwd = (r(k+h)-r0v)/h; bwd = (r0v-r(k-h))/h; an = dr(k)
        sc = max(abs(an), abs(fwd), abs(bwd), 1.0)
        if abs(fwd-bwd) > 1e-4*sc:
            if not kinks: first = k
            kinks += 1; continue
        rmax = max(rmax, abs(0.5*(fwd+bwd)-an)/sc)
    jump = 0.0
    for ks in (kmax, kmin):
        if ks > 0:
            h = 1e-6*max(k_ref, ks)
            jump = max(jump, abs(r(ks+h)-r(ks-h)))
    return rmax, jump, kinks, first

modes = ["perfect","exp_hardening","linear_hardening","simple_softening",
         "exp_softening","piecewise_softening"]
print(f"{'mode':22} {'dr err':>10} {'r jump':>10} {'kinks':>6}  verdict")
for m in modes:
    kw = dict(r0=0.7, h_=0.0) if m in ("exp_hardening","simple_softening") else {}
    rmax, jump, kinks, first = check(*make(m, **kw))
    v = "OK" if rmax < 1e-5 and jump < 1e-4 else "FAILED"
    note = f"  (kink @ {first:g})" if kinks else ""
    print(f"{m:22} {rmax:10.2e} {jump:10.2e} {kinks:6d}  {v}{note}")

print("\nsame, with a deliberately wrong derivative (x1.7) -- all must FAIL:")
for m in modes:
    if m == "perfect": continue
    kw = dict(r0=0.7, h_=0.0) if m in ("exp_hardening","simple_softening") else {}
    rmax, jump, kinks, first = check(*make(m, bug=True, **kw))
    v = "OK" if rmax < 1e-5 and jump < 1e-4 else "FAILED"
    print(f"{m:22} {rmax:10.2e} {jump:10.2e} {kinks:6d}  {v}"
          + ("   <-- MISSED" if v == "OK" else ""))

print("\nand a deliberately DISCONTINUOUS r (exp_softening residual bumped below kmax):")
def r_bad(k): return 1.05 if k < 0.001 else 0.7 + 0.3*math.exp(-30*(k-0.001))
def dr_bad(k): return 0.0 if k < 0.001 else -0.3*30*math.exp(-30*(k-0.001))
rmax, jump, kinks, first = check(r_bad, dr_bad, 0.001, 0.02)
print(f"{'exp_softening(jump)':22} {rmax:10.2e} {jump:10.2e} {kinks:6d}  "
      + ("OK   <-- MISSED" if jump < 1e-4 else "FAILED"))

