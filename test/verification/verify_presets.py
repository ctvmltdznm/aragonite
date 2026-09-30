"""Compare UMAT flag blocks (UMAT_QUADRIC_PRIMAL_Major.f) with the MOOSE preset mapping
(MaterialModelPresets.h + applyBonePreset / ComputeFabricElasticityTensor).
Yield: phi = sqrt(s:F:s) + f.s on random stresses. Elasticity: sigma = C:eps on random strains."""
import numpy as np
rng = np.random.default_rng(0)
r2 = np.sqrt(2.0)
def tsfu(rho, ex, d=1.0):
    return rho**ex if rho <= 0.5 else rho**ex + (d-1)*((rho-0.5)/0.5)**ex

P = {  # UMAT constants
 0: dict(E0=8534.64, V0=0.246, KS=1.63, s0p=61.43, s0n=89.24, z0=0.1876, PP=1.686),
 1: dict(E0=6561.94, EAA=18661.3, V0=0.323, VA0=0.32, MUA0=3739.02, KS=1.63,
         s0p=50.63, s0n=68.07, sap=98.68, san=165.80, z0=0.4707, za=0.1809, taua=51.62, PP=1.69),
 2: dict(E0=9994.65, V0=0.2278, MU0=3361.14, KS=1.62, LS=1.1, s0p=66.012, s0n=98.876,
         z0=0.2182, tau0=41.889, PP=1.686, QQ=1.05),
 3: dict(E0=9994.65, V0=0.2278, MU0=3361.14, KS=1.62, LS=1.1, s0p=66.012, s0n=98.876,
         z0=0.218, tau0=41.889, PP=1.69, QQ=1.05),
 4: dict(E0=19327.0, V0=0.3434, KS=1.63, s0p=144.7, s0n=234.2, z0=0.49, PP=1.69),
 5: dict(E0=15079.6, EAA=24578.5, V0=0.4620, VA0=0.354, MUA0=6578.01, KS=1.63,
         s0p=56.54, s0n=201.28, sap=176.40, san=268.00, z0=0.0074, za=1.4045, taua=82.55, PP=1.686),
}
# Abaqus/Mandel order: 11 22 33 12 13 23. Plane index in that order: 3->(0,1) 4->(0,2) 5->(1,2)
APL = {3: (0, 1), 4: (0, 2), 5: (1, 2)}

def umat(bf, major, rho, mm, umat_bug=False):
    c = P[bf]; a = major-1
    CC = np.zeros((6, 6)); SS = None; FF4 = np.zeros((6, 6)); FF = np.zeros(6)
    tk = tsfu(rho, c['KS']); tp = tsfu(rho, c['PP'])
    S = lambda sp, sn: (sp+sn)/2/sp/sn
    if bf in (0, 4, 1, 5):
        E0, V0 = c['E0'], c['V0']; MU0 = E0/2/(1+V0)
        Ev = [c.get('EAA') if (bf in (1, 5) and i == a) else E0 for i in range(3)]
        for i in range(3): CC[i, i] = 1/Ev[i]
        for k, (i, j) in APL.items():
            ax = bf in (1, 5) and a in (i, j)
            CC[k, k] = 0.5/(c['MUA0'] if ax else MU0)
            CC[i, j] = CC[j, i] = (-c['VA0']/c['EAA']) if ax else (-V0/E0)
        CC /= tk; SS = np.linalg.inv(CC)
        S0 = S(c['s0p'], c['s0n']); T0 = np.sqrt(0.5/S0**2/(1+c['z0']))
        if bf in (0, 4):
            for i in range(3): FF4[i, i] = S0**2; FF[i] = -(c['s0p']-c['s0n'])/2/c['s0p']/c['s0n']
            for k, (i, j) in APL.items(): FF4[i, j] = FF4[j, i] = -c['z0']*S0**2; FF4[k, k] = 0.5/T0**2
        else:
            SA = S(c['sap'], c['san'])
            for i in range(3):
                FF4[i, i] = SA**2 if i == a else S0**2
                FF[i] = (-(c['sap']-c['san'])/2/c['sap']/c['san']) if i == a else (-(c['s0p']-c['s0n'])/2/c['s0p']/c['s0n'])
            for k, (i, j) in APL.items():
                ax = a in (i, j)
                FF4[i, j] = FF4[j, i] = (-c['za']*SA**2) if ax else (-c['z0']*S0**2)
                FF4[k, k] = 0.5/(c['taua'] if ax else T0)**2
        FF4 /= tp**2; FF /= tp
    else:
        m = list(mm)
        if bf == 2:
            b, cc = [i for i in range(3) if i != a]
            m[b] = m[cc] = 0.5*(m[b]+m[cc])
        E0, V0, MU0, LS, QQ = c['E0'], c['V0'], c['MU0'], c['LS'], c['QQ']
        Giso = E0/2/(1+V0)
        SS = np.zeros((6, 6))
        lam = E0*V0/(1+V0)/(1-2*V0)
        for i in range(3):
            SS[i, i] = E0*(1-V0)/(1+V0)/(1-2*V0)*m[i]**(2*LS)*tk
        for k, (i, j) in APL.items():
            tr = bf == 2 and a not in (i, j)
            SS[k, k] = 2*(Giso if tr else MU0)*(m[i]*m[j])**LS*tk
            if umat_bug and bf == 2 and major == 3 and k == 5:
                SS[k, k] = 2*MU0*(m[i]*m[j])**LS          # UMAT line: missing TSFU
            SS[i, j] = SS[j, i] = lam*(m[i]*m[j])**LS*tk
        S0 = S(c['s0p'], c['s0n']); T0 = np.sqrt(0.5/S0**2/(1+c['z0']))
        for i in range(3):
            FF4[i, i] = S0**2/(tp*m[i]**(2*QQ))**2
            FF[i] = -(c['s0p']-c['s0n'])/2/c['s0p']/c['s0n']/(tp*m[i]**(2*QQ))
        for k, (i, j) in APL.items():
            FF4[i, j] = FF4[j, i] = -(c['z0']*(m[i]/m[j])**(2*QQ))*S0**2/(tp*m[i]**(2*QQ))**2
            tr = bf == 2 and a not in (i, j)
            FF4[k, k] = 0.5/((T0 if tr else c['tau0'])*tp*(m[i]*m[j])**QQ)**2
    return SS, FF4, FF

# ---------------- MOOSE side (mirror of the C++) ----------------
PL = [(0, 1), (0, 2), (1, 2)]  # planes 12 13 23
def moose_yield(bf, major, rho, mm):
    c = P[bf]; a = major-1; t = tsfu(rho, c['PP'])
    S = lambda sp, sn: (sp+sn)/(2*sp*sn)
    tau_iso = 1/(S(c['s0p'], c['s0n'])*np.sqrt(2*(1+c['z0'])))
    sT = [0]*3; sC = [0]*3; tau = [0]*3; z = [0]*3
    if bf in (0, 4):
        sT = [c['s0p']*t]*3; sC = [c['s0n']*t]*3; tau = [tau_iso*t]*3; z = [c['z0']]*3
    elif bf in (1, 5):
        for i in range(3):
            sT[i] = (c['sap'] if i == a else c['s0p'])*t; sC[i] = (c['san'] if i == a else c['s0n'])*t
        for p, (i, j) in enumerate(PL):
            ax = a in (i, j)
            tau[p] = (c['taua'] if ax else tau_iso)*t
            z[p] = c['za']*S(sT[a], sC[a])**2/S(sT[i], sC[i])**2 if ax else c['z0']
    else:
        m = list(mm); q = c['QQ']
        if bf == 2:
            b, cc = (a+1) % 3, (a+2) % 3; m[b] = m[cc] = 0.5*(m[b]+m[cc])
        for i in range(3):
            sT[i] = c['s0p']*t*m[i]**(2*q); sC[i] = c['s0n']*t*m[i]**(2*q)
        for p, (i, j) in enumerate(PL):
            tr = bf == 2 and i != a and j != a
            tau[p] = (tau_iso if tr else c['tau0'])*t*(m[i]*m[j])**q
            z[p] = c['z0']*(m[i]/m[j])**(2*q)
    # constructor F build, MOOSE order 11 22 33 23 13 12, tensor shear
    F = np.zeros((6, 6)); f = np.zeros(6)
    for i in range(3): F[i, i] = S(sT[i], sC[i])**2; f[i] = (sC[i]-sT[i])/(2*sC[i]*sT[i])
    F[3, 3] = 1/tau[2]**2; F[4, 4] = 1/tau[1]**2; F[5, 5] = 1/tau[0]**2
    F[0, 1] = F[1, 0] = -z[0]*F[0, 0]; F[0, 2] = F[2, 0] = -z[1]*F[0, 0]; F[1, 2] = F[2, 1] = -z[2]*F[1, 1]
    return F, f

def moose_elastic(bf, major, rho, mm):
    c = P[bf]; a = major-1; t = tsfu(rho, c['KS']); E0, nu0 = c['E0'], c['V0']
    Giso = E0/(2*(1+nu0)); E = [0]*3; G = [0]*3; So = [0]*3
    if bf in (0, 4):
        E = [E0*t]*3; G = [Giso*t]*3; So = [-nu0/(E0*t)]*3
    elif bf in (1, 5):
        for i in range(3): E[i] = (c['EAA'] if i == a else E0)*t
        for p, (i, j) in enumerate(PL):
            ax = a in (i, j)
            G[p] = (c['MUA0'] if ax else Giso)*t
            So[p] = -c['VA0']/(c['EAA']*t) if ax else -nu0/(E0*t)
    else:
        m = list(mm); l = c['LS']
        if bf == 2:
            b, cc = (a+1) % 3, (a+2) % 3; m[b] = m[cc] = 0.5*(m[b]+m[cc])
        for i in range(3): E[i] = E0*t*m[i]**(2*l)
        for p, (i, j) in enumerate(PL):
            mmv = (m[i]*m[j])**l; tr = bf == 2 and i != a and j != a
            G[p] = (Giso if tr else c['MU0'])*t*mmv; So[p] = -nu0/(E0*t*mmv)
    nu12, nu13, nu23 = -So[0]*E[0], -So[1]*E[0], -So[2]*E[1]
    E1, E2, E3 = E
    nu21, nu31, nu32 = nu12*E2/E1, nu13*E3/E1, nu23*E3/E2
    d = 1 - nu12*nu21 - nu23*nu32 - nu31*nu13 - 2*nu21*nu32*nu13
    C = np.zeros((6, 6))  # engineering Voigt, order 11 22 33 23 13 12
    C[0, 0] = E1*(1-nu23*nu32)/d; C[1, 1] = E2*(1-nu13*nu31)/d; C[2, 2] = E3*(1-nu12*nu21)/d
    C[0, 1] = C[1, 0] = E1*(nu21+nu31*nu23)/d; C[0, 2] = C[2, 0] = E1*(nu31+nu21*nu32)/d
    C[1, 2] = C[2, 1] = E2*(nu32+nu12*nu31)/d
    C[3, 3], C[4, 4], C[5, 5] = G[2], G[1], G[0]
    return C

def to_mandel_abq(T):   # 11 22 33 12 13 23
    return np.array([T[0, 0], T[1, 1], T[2, 2], r2*T[0, 1], r2*T[0, 2], r2*T[1, 2]])
def to_voigt_moose(T):  # 11 22 33 23 13 12 (tensor shear)
    return np.array([T[0, 0], T[1, 1], T[2, 2], T[1, 2], T[0, 2], T[0, 1]])

worst = {}
for bf in range(6):
    for major in (1, 2, 3):
        for trial in range(20):
            rho = rng.uniform(0.05, 0.45)
            mm = rng.uniform(0.6, 1.4, 3); mm *= 3/mm.sum()
            SS, FF4, FF = umat(bf, major, rho, mm)
            F, f = moose_yield(bf, major, rho, mm)
            C = moose_elastic(bf, major, rho, mm)
            for k in range(10):
                A = rng.normal(size=(3, 3)); sig = (A+A.T)*50*rho
                sm = to_mandel_abq(sig); sv = to_voigt_moose(sig)
                phu = np.sqrt(sm@FF4@sm)+FF@sm; phm = np.sqrt(sv@F@sv)+f@sv
                e_y = abs(phu-phm)/abs(phu)
                B = rng.normal(size=(3, 3)); eps = (B+B.T)*1e-3
                su = SS@to_mandel_abq(eps)                   # Mandel stress, Abaqus order
                ev = to_voigt_moose(eps); ev[3:] *= 2        # engineering shear strain
                sv2 = C@ev
                su_t = np.array([su[0], su[1], su[2], su[5]/r2, su[4]/r2, su[3]/r2])
                e_e = np.max(abs(su_t-sv2))/np.max(abs(su_t))
                key = (bf, major)
                worst[key] = max(worst.get(key, (0, 0))[0], e_y), max(worst.get(key, (0, 0))[1], e_e)
for (bf, mj), (ey, ee) in sorted(worst.items()):
    print(f"BF={bf} MAJOR={mj}:  max rel err yield {ey:.2e}   elastic {ee:.2e}")

# the UMAT SSSS(6,6) omission, BF=2 MAJOR=3
rho = 0.15; mm = np.array([0.9, 0.95, 1.15])
SSb, _, _ = umat(2, 3, rho, mm, umat_bug=True); SSf, _, _ = umat(2, 3, rho, mm)
print(f"\nBF=2 MAJOR=3, rho={rho}: UMAT SSSS(6,6)/corrected = {SSb[5,5]/SSf[5,5]:.2f}  (= 1/rho^KS = {1/tsfu(rho,1.62):.2f})")


# =============================================================================
# 2. ALL 49 elastic_model x plastic_model COMBINATIONS
# =============================================================================
# The two flags are independent, so every pairing must give a usable material:
# a positive-definite stiffness and a convex (PSD) yield surface. Coral and
# legacy stand-ins use the aragonite production values.
MODELS = ["trabecular_iso", "trabecular_ti", "trabecular_fabric_ti",
          "trabecular_fabric_ortho", "compact_iso", "compact_ti", "coral", "legacy"]
BF_OF = {"trabecular_iso": 0, "trabecular_ti": 1, "trabecular_fabric_ti": 2,
         "trabecular_fabric_ortho": 3, "compact_iso": 4, "compact_ti": 5}
RDY = {0: 0.7, 1: 0.7, 2: 0.7, 3: 0.7, 4: 0.7, 5: 0.7}   # bone presets, UMAT RDY
CORAL_C = [171800, 57500, 30200, 106700, 46900, 84200, 42100, 31100, 46600]
CORAL_S = dict(sT=[4980, 4100, 5340], tau=[4510, 5080, 5060], zeta=[0.27, 0.27, 0.27])
RHO, MM, MAJOR = 0.20, np.array([0.85, 0.95, 1.20]), 3
LEGACY_E = dict(E0=9994.65, nu0=0.2278, k=2.0, l=1.0)   # legacy defaults, user E_0/nu_0

def C_of(model):
    """6x6 engineering Voigt stiffness, order 11 22 33 23 13 12."""
    if model in BF_OF:
        return moose_elastic(BF_OF[model], MAJOR, RHO, MM)
    if model == "coral":
        C11, C12, C13, C22, C23, C33, C44, C55, C66 = CORAL_C
        C = np.zeros((6, 6))
        C[:3, :3] = [[C11, C12, C13], [C12, C22, C23], [C13, C23, C33]]
        C[3, 3], C[4, 4], C[5, 5] = C44, C55, C66
        return C
    # legacy: fabric Zysset-Curnier with the legacy defaults
    E0, nu0, k, l = LEGACY_E["E0"], LEGACY_E["nu0"], LEGACY_E["k"], LEGACY_E["l"]
    t = RHO**k; G0 = E0/(2*(1+nu0))
    E = [E0*t*MM[i]**(2*l) for i in range(3)]
    G = [G0*t*(MM[i]*MM[j])**l for i, j in PL]
    So = [-nu0/(E0*t*(MM[i]*MM[j])**l) for i, j in PL]
    nu12, nu13, nu23 = -So[0]*E[0], -So[1]*E[0], -So[2]*E[1]
    E1, E2, E3 = E
    nu21, nu31, nu32 = nu12*E2/E1, nu13*E3/E1, nu23*E3/E2
    d = 1 - nu12*nu21 - nu23*nu32 - nu31*nu13 - 2*nu21*nu32*nu13
    C = np.zeros((6, 6))
    C[0,0]=E1*(1-nu23*nu32)/d; C[1,1]=E2*(1-nu13*nu31)/d; C[2,2]=E3*(1-nu12*nu21)/d
    C[0,1]=C[1,0]=E1*(nu21+nu31*nu23)/d; C[0,2]=C[2,0]=E1*(nu31+nu21*nu32)/d
    C[1,2]=C[2,1]=E2*(nu32+nu12*nu31)/d
    C[3,3], C[4,4], C[5,5] = G[2], G[1], G[0]
    return C

def F_of(model):
    """(_F_matrix, _f_lin_vector, r0) exactly as the constructor builds them."""
    if model in BF_OF:
        F, f = moose_yield(BF_OF[model], MAJOR, RHO, MM)
        return F, f, RDY[BF_OF[model]]
    sT = CORAL_S["sT"]; sC = sT; tau = CORAL_S["tau"]; z = CORAL_S["zeta"]
    F = np.zeros((6, 6)); f = np.zeros(6)
    for i in range(3):
        F[i, i] = ((sT[i]+sC[i])/(2*sT[i]*sC[i]))**2
        f[i] = (sC[i]-sT[i])/(2*sC[i]*sT[i])
    F[3,3] = 1/tau[2]**2; F[4,4] = 1/tau[1]**2; F[5,5] = 1/tau[0]**2
    F[0,1]=F[1,0]=-z[0]*F[0,0]; F[0,2]=F[2,0]=-z[1]*F[0,0]; F[1,2]=F[2,1]=-z[2]*F[1,1]
    return F, f, 1.0   # coral and legacy both keep r(0) = 1

def mandel6(C):
    """engineering-Voigt stiffness -> Mandel matrix (shear rows/cols x2)."""
    M = C.copy().astype(float)
    M[3:, :] *= r2; M[:, 3:] *= r2
    return M

def psd_minors(F):
    A = F[:3, :3]
    sc = abs(A).max()
    m = [A[0,0]*A[1,1]-A[0,1]**2, A[0,0]*A[2,2]-A[0,2]**2, A[1,1]*A[2,2]-A[1,2]**2]
    return all(x >= -1e-10*sc**2 for x in m) and np.linalg.det(A) >= -1e-10*sc**3

def yield_onset(F, f, r0, axis):
    """uniaxial tensile stress at first yield: solve sqrt(F_aa)s + f_a s = r0."""
    return r0/(np.sqrt(F[axis, axis]) + f[axis])

print("\n" + "="*78)
print("2. elastic_model x plastic_model: the 49 preset pairings plus the legacy row and "
      "column\n   = 64 combinations "
      f"(rho={RHO}, m={[float(x) for x in MM]}, main_direction={MAJOR}).\n"
      "   legacy stands in with its own defaults (E_0/nu_0 user-supplied) and, on the "
      "plastic\n   side, with the aragonite explicit strengths.")
print("="*78)
bad = []
for em in MODELS:
    C = C_of(em); Cm = mandel6(C)
    e_ok = np.linalg.eigvalsh((Cm+Cm.T)/2).min() > 0
    for pm in MODELS:
        F, f, r0 = F_of(pm)
        y_ok = psd_minors(F) and np.linalg.eigvalsh(F).min() >= -1e-10*abs(F).max()
        if not (e_ok and y_ok):
            bad.append((em, pm, e_ok, y_ok))
print(f"combinations checked: {len(MODELS)**2}   failures: {len(bad)}")
for em, pm, e_ok, y_ok in bad:
    print(f"    {em:24s} x {pm:24s}  elastic PD {e_ok}  yield PSD {y_ok}")

print("\nPer-model summary (E in MPa, yield onset = r0 x uniaxial tensile strength):")
print(f"{'model':26s} {'E1':>9s} {'E2':>9s} {'E3':>9s} {'G12':>9s} "
      f"{'sy_xx':>8s} {'sy_yy':>8s} {'sy_zz':>8s} {'r0':>5s}")
for m in MODELS:
    C = C_of(m); S = np.linalg.inv(C[:3, :3])
    E = [1/S[i, i] for i in range(3)]
    F, f, r0 = F_of(m)
    sy = [yield_onset(F, f, r0, i) for i in range(3)]
    print(f"{m:26s} {E[0]:9.1f} {E[1]:9.1f} {E[2]:9.1f} {C[5,5]:9.1f} "
          f"{sy[0]:8.2f} {sy[1]:8.2f} {sy[2]:8.2f} {r0:5.2f}")

# =============================================================================
# 3. CONVEXITY GUARD: principal minors vs eigenvalues
# =============================================================================
print("\n" + "="*78)
print("3. checkYieldSurfaceConvexity(): principal-minor test vs eigenvalues")
print("="*78)
mism = 0
for bf in range(6):
    for major in (1, 2, 3):
        F, _ = moose_yield(bf, major, 0.2, MM)
        mism += psd_minors(F) != (np.linalg.eigvalsh(F[:3, :3]).min() >= -1e-12)
for _ in range(20000):
    d = rng.uniform(0.1, 10, 3); z = rng.uniform(-2, 2, 3)
    A = np.diag(d)
    A[0,1]=A[1,0]=-z[0]*d[0]; A[0,2]=A[2,0]=-z[1]*d[0]; A[1,2]=A[2,1]=-z[2]*d[1]
    F = np.zeros((6, 6)); F[:3, :3] = A
    ev = np.linalg.eigvalsh(A).min()
    if abs(ev) > 1e-8:
        mism += psd_minors(F) != (ev >= 0)
print(f"presets (18) + 20000 random zeta triples: mismatches {mism}")

# old Eq.55/Eq.56 form, for the record
def old_check(F):
    F11, F22, F33 = F[0,0], F[1,1], F[2,2]
    z12, z13, z23 = -F[0,1]/F11, -F[0,2]/F11, -F[1,2]/F22
    ok55 = (abs(z12) <= F22/F11+1e-10 and abs(z13) <= F33/F11+1e-10
            and abs(z23) <= F33/F22+1e-10)
    det = (F22**2*F33**2 - F11**2*F33**2*z12**2 - F11**2*F22**2*z13**2
           + 2*F11**2*F22**2*z12*z13*z23 - F22**4*z23**2)
    return ok55 and det >= -1e-10
print("presets the OLD Eq.55/56 guard rejected although convex:")
for bf in range(6):
    for major in (1, 2, 3):
        F, _ = moose_yield(bf, major, 0.2, MM)
        if psd_minors(F) and not old_check(F):
            print(f"    BF={bf} MAJOR={major}")

# =============================================================================
# 4. POST-YIELD LAWS vs UMAT RADK/DRADK, and derivative consistency
# =============================================================================
print("\n" + "="*78)
print("4. Post-yield laws: MOOSE vs UMAT RADK, and dr/dk vs finite differences")
print("="*78)
RDY0, KSLOPE, KMAX, KWIDTH = 0.7, 100.0, 0.015, 8.0

def umat_radk(pyfl, k):
    if pyfl == 1:
        return RDY0 + (1-RDY0)*(1-np.exp(-KSLOPE*k))
    off = 1.0/KWIDTH
    return RDY0 + (1-RDY0)*(np.exp(-((k-KMAX)**2)/(KWIDTH*KMAX**2)) - np.exp(-off-KSLOPE*k))

def umat_dradk(pyfl, k):
    if pyfl == 1:
        return (1-RDY0)*KSLOPE*np.exp(-KSLOPE*k)
    off = 1.0/KWIDTH
    return (1-RDY0)*((-k+KMAX)*np.exp(-((k-KMAX)**2)/(KWIDTH*KMAX**2))/(0.5*KWIDTH*KMAX**2)
                     + KSLOPE*np.exp(-off-KSLOPE*k))

def moose_r(mode, k, r0=RDY0, h=0.0):
    if mode == "exp_hardening":
        return r0 + (1-r0+h)*(1-np.exp(-KSLOPE*k))
    return r0 + (1-r0)*(np.exp(-((k-KMAX)**2)/(KWIDTH*KMAX**2)) - np.exp(-1/KWIDTH-KSLOPE*k))

def moose_dr(mode, k, r0=RDY0, h=0.0):
    if mode == "exp_hardening":
        return (1-r0+h)*KSLOPE*np.exp(-KSLOPE*k)
    return (1-r0)*(-2*(k-KMAX)/(KWIDTH*KMAX**2)*np.exp(-((k-KMAX)**2)/(KWIDTH*KMAX**2))
                   + KSLOPE*np.exp(-1/KWIDTH-KSLOPE*k))

ks = np.array([0.0, 0.003, 0.0135, 0.015, 0.0225, 0.06])
for mode, pyfl in [("exp_hardening", 1), ("simple_softening", 2)]:
    dr_err = max(abs(moose_r(mode, k) - umat_radk(pyfl, k)) for k in ks)
    dd_err = max(abs(moose_dr(mode, k) - umat_dradk(pyfl, k)) for k in ks)
    hh = 1e-7
    fd_err = max(abs((moose_r(mode, k+hh)-moose_r(mode, k-hh))/(2*hh) - moose_dr(mode, k))
                 / max(abs(moose_dr(mode, k)), 1.0) for k in ks[1:])
    print(f"{mode:18s} PYFL={pyfl}  r vs UMAT {dr_err:.2e}   dr/dk vs UMAT {dd_err:.2e}"
          f"   dr/dk vs FD {fd_err:.2e}   r(0) = {moose_r(mode, 0.0):.3f}")
# legacy exp_hardening must be untouched by the generalisation (r0 = 1)
leg = max(abs(moose_r("exp_hardening", k, r0=1.0, h=0.7) - (1 + 0.7*(1-np.exp(-KSLOPE*k))))
          for k in ks)
print(f"legacy exp_hardening (r0=1, h=0.7) unchanged: max diff {leg:.2e}")

print("\n" + "="*78)
print("Representative presets for run_tests.sh (the expensive layer):")
print("  M1  trabecular_fabric_ortho  - bone production path, fabric + simple_softening")
print("  M2  compact_ti, main_direction=3  - TI zeta conversion, exp_hardening, the")
print("      surface the old convexity guard rejected")
print("  M3  coral  - aragonite production path, explicit 9+9 constants")
print("  M4  legacy - regression guard: pre-flag behaviour must not move")
print("="*78)
