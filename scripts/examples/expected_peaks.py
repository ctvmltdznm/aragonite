"""Predicted first-yield stress for every example, so the README can say what to
look for and so we catch any preset that would not yield at 2% strain."""
import numpy as np, importlib.util, sys
spec = importlib.util.spec_from_file_location("vp", "verify_presets.py")
vp = importlib.util.module_from_spec(spec)
_out = sys.stdout
class _N:
    def write(self,*a): pass
    def flush(self): pass
sys.stdout = _N(); spec.loader.exec_module(vp); sys.stdout = _out

MM = np.array([0.85, 0.95, 1.20])
R0 = 0.7
CASES = {"trabecular_iso":(0,0.20,False), "trabecular_ti":(1,0.20,False),
         "trabecular_fabric_ti":(2,0.20,True), "trabecular_fabric_ortho":(3,0.20,True),
         "compact_iso":(4,0.90,False), "compact_ti":(5,0.90,True)}
MAJOR = 3
print(f"{'model':24} {'case':16} {'yield stress':>13} {'yield strain':>13} {'margin @2%':>11}")
print("-"*82)
worst = 1e9
for name,(bf,rho,_) in CASES.items():
    F,f = vp.moose_yield(bf, MAJOR, rho, MM)
    C   = vp.moose_elastic(bf, MAJOR, rho, MM)
    S   = np.linalg.inv(C)
    for i,ax in enumerate(("xx","yy","zz")):
        for sgn,kind in ((+1,"tension"),(-1,"compression")):
            # uniaxial stress s along axis i: phi = |s|*sqrt(F_ii) + f_i*s = r0
            sy = R0/(np.sqrt(F[i,i]) + sgn*f[i]) * sgn
            ey = S[i,i]*sy                       # uniaxial stress -> e = S_ii * s
            print(f"{name:24} {kind+'_'+ax:16} {sy:13.3f} {ey:13.5f} {0.02/abs(ey):10.1f}x")
            worst = min(worst, 0.02/abs(ey))
    for p,(ax,sl) in enumerate([("xy",5),("xz",4),("yz",3)]):
        tau_y = R0/np.sqrt(F[sl,sl])
        g_y = tau_y/C[sl,sl]                     # engineering shear
        print(f"{name:24} {'shear_'+ax:16} {tau_y:13.3f} {g_y:13.5f} {0.02/g_y:10.1f}x")
        worst = min(worst, 0.02/g_y)
print("-"*82)
print(f"smallest margin over all 54 cases: {worst:.1f}x  "
      + ("OK, every case yields well before t=1" if worst > 1.5 else "TOO TIGHT"))
