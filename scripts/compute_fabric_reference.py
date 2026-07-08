#!/usr/bin/env python3
"""
compute_fabric_reference.py
Compute the explicit symmetric9 stiffness that ComputeFabricElasticityTensor
should produce, for the fabric<->orthotropy equivalence test (Checkpoint 3).

Relations (Zysset & Curnier 1995, matching the class header):
  E_i   = E0 * rho^k * m_i^(2l)
  G_ij  = G0 * rho^k * m_i^l * m_j^l
  nu_ij = nu0 * (m_i/m_j)^l          (compliance is symmetric: nu_ij/E_i = nu_ji/E_j)

Outputs the 9 orthotropic stiffness entries for MOOSE fill_method=symmetric9:
  C11 C12 C13 C22 C23 C33 C44 C55 C66
"""
import numpy as np

def fabric_stiffness(E0, nu0, rho, m1, m2, m3, k=1.8, l=1.0, G0=None):
    if G0 is None:
        G0 = E0 / (2.0 * (1.0 + nu0))
    rk = rho**k
    m = [m1, m2, m3]
    E  = [E0 * rk * mi**(2*l) for mi in m]
    G23 = G0 * rk * m[1]**l * m[2]**l
    G13 = G0 * rk * m[0]**l * m[2]**l
    G12 = G0 * rk * m[0]**l * m[1]**l
    nu12 = nu0 * (m[0]/m[1])**l
    nu13 = nu0 * (m[0]/m[2])**l
    nu23 = nu0 * (m[1]/m[2])**l
    # Orthotropic compliance (6x6, Voigt), symmetric by construction
    S = np.zeros((6,6))
    S[0,0]=1/E[0]; S[1,1]=1/E[1]; S[2,2]=1/E[2]
    S[0,1]=-nu12/E[0]; S[0,2]=-nu13/E[0]; S[1,2]=-nu23/E[1]
    S[1,0]=S[0,1]; S[2,0]=S[0,2]; S[2,1]=S[1,2]
    S[3,3]=1/G23; S[4,4]=1/G13; S[5,5]=1/G12
    C = np.linalg.inv(S)
    entries = [C[0,0],C[0,1],C[0,2],C[1,1],C[1,2],C[2,2],C[3,3],C[4,4],C[5,5]]
    return entries, dict(E=E,G=(G23,G13,G12),nu=(nu12,nu13,nu23))

if __name__ == "__main__":
    # Illustrative trabecular-bone-like case
    E0, nu0, rho = 15000.0, 0.3, 0.3
    m1, m2, m3, k, l = 1.2, 1.0, 0.8, 1.8, 1.0
    C, info = fabric_stiffness(E0,nu0,rho,m1,m2,m3,k,l)
    print("Base: E0=%g nu0=%g rho=%g  fabric m=(%g,%g,%g) k=%g l=%g"%(E0,nu0,rho,m1,m2,m3,k,l))
    print("Directional E:", [round(x,2) for x in info['E']])
    print("G(23,13,12):",  [round(x,2) for x in info['G']])
    print("nu(12,13,23):", [round(x,4) for x in info['nu']])
    print("\nsymmetric9 C_ijkl (C11 C12 C13 C22 C23 C33 C44 C55 C66):")
    print("'" + " ".join("%.4f"%c for c in C) + "'")

def _voigt_stress(C_entries, eps_eng):
    import numpy as np
    C11,C12,C13,C22,C23,C33,C44,C55,C66 = C_entries
    C = np.array([[C11,C12,C13,0,0,0],[C12,C22,C23,0,0,0],[C13,C23,C33,0,0,0],
                  [0,0,0,C44,0,0],[0,0,0,0,C55,0],[0,0,0,0,0,C66]])
    return C @ np.array(eps_eng)

def write_gold(path, C_entries):
    # Uniform strain matching the test BCs (engineering shears = 2*eps_offdiag)
    eps = [1e-3, -5e-4, 2e-4, 2*3e-4, 2*4e-4, 2*5e-4]  # xx yy zz yz xz xy
    sig = _voigt_stress(C_entries, eps)
    # Voigt order here is xx yy zz yz xz xy; map to postprocessor columns
    s_xx,s_yy,s_zz,s_yz,s_xz,s_xy = sig
    import os
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path,'w') as f:
        f.write("time,s_xx,s_yy,s_zz,s_xy,s_xz,s_yz\n")
        f.write("0,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n"%(s_xx,s_yy,s_zz,s_xy,s_xz,s_yz))

if __name__ == "__main__" and len(__import__('sys').argv) > 1 and __import__('sys').argv[1]=='--gold':
    C_o,_ = fabric_stiffness(15000,0.3,0.3,1.2,1.0,0.8,1.8,1.0)
    C_i,_ = fabric_stiffness(15000,0.3,1.0,1.0,1.0,1.0,1.8,1.0)
    write_gold('tv/fabric/fabric_orthotropy/gold/fabric_orthotropy_out.csv', C_o)
    write_gold('tv/fabric/fabric_isotropic/gold/fabric_isotropic_out.csv',   C_i)
    print("golds written")
