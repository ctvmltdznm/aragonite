// MaterialModelPresets.h
//
// Single source of truth for the material-model flags shared by
//   ComputeFabricElasticityTensor     (parameter: elastic_model)
//   OrthotropicPlasticityStressUpdate (parameter: plastic_model)
//
// Flag IDs 0-5 are the UMAT bone flags, PROPS(1) in UMAT_QUADRIC_PRIMAL_Major.f
// (Schwiedrzik 2012/2013). All bone constants below are transcribed from
// that file. Line numbers refer to it.
//
//   0 trabecular_iso           ISORGWTB, PHZ 2013        (UMAT 359-440)
//   1 trabecular_ti            TIRGWTB,  PHZ 2013        (UMAT 441-651)
//   2 trabecular_fabric_ti     FABRGWTB, PHZ 2013        (UMAT 652-914)
//   3 trabecular_fabric_ortho  FABRGWTB, PHZ 2013        (UMAT 915-1015)
//   4 compact_iso              ISORBCT,  JJS 2013        (UMAT 1016-1098)
//   5 compact_ti               TIRBCT,   JJS 2013        (UMAT 1099-1310)
//   6 coral                    full orthotropy, values from our aragonite runs
//                              (fill_and_generate_configs.py, BULK_MATS)
//   7 legacy                   pre-flag behaviour, unchanged (default)
//
// Conventions (all verified against the UMAT):
//  * Density scaling: stiffness x TSFU(rho,k,delta), strengths x TSFU(rho,p,delta).
//    DELTA = 1 in every UMAT preset.
//  * The UMAT works in Mandel notation, order 11 22 33 12 13 23, with
//    FFFF(shear) = 0.5/tau^2. In our tensor-shear Voigt (order 11 22 33 23 13 12)
//    with F(shear) = 1/tau^2 the shear strengths therefore carry over 1:1.
//  * Shear strengths that the UMAT does not take as input (isotropic plane of
//    the iso / TI / fabric-TI models) follow from the isotropy constraint
//        tau = 1 / (S0 sqrt(2 (1 + zeta0))),  S0 = (s+ + s-) / (2 s+ s-).
//  * Post-yield: UMAT strengths are ULTIMATE strengths; yield onsets at
//    r(0) = RDY = 0.7. Hence every bone preset sets initial_yield_ratio = 0.7.
//    PYFL=2 -> simple_softening, PYFL=1 -> exp_hardening (UMAT forms).
//  * Not carried over: densification (DENSFL) and the UMAT damage law
//    (KONSTD/CRITD). Neither exists in the MOOSE implementation.

#pragma once

#include "MooseEnum.h"
#include "MooseTypes.h"

namespace MaterialModelPresets
{

enum class Model : int
{
  TRABECULAR_ISO = 0,
  TRABECULAR_TI = 1,
  TRABECULAR_FABRIC_TI = 2,
  TRABECULAR_FABRIC_ORTHO = 3,
  COMPACT_ISO = 4,
  COMPACT_TI = 5,
  CORAL = 6,
  LEGACY = 7
};

inline MooseEnum
modelEnum()
{
  return MooseEnum("trabecular_iso=0 trabecular_ti=1 trabecular_fabric_ti=2 "
                   "trabecular_fabric_ortho=3 compact_iso=4 compact_ti=5 coral=6 legacy=7",
                   "legacy");
}

inline std::string
modelDoc()
{
  return "Material model flag. IDs 0-5 equal the UMAT bone flags (PROPS(1)): "
         "0 trabecular_iso, 1 trabecular_ti, 2 trabecular_fabric_ti, "
         "3 trabecular_fabric_ortho, 4 compact_iso, 5 compact_ti. "
         "6 coral: full orthotropy, every constant explicit (defaults from our aragonite "
         "runs); coral RVEs normally resolve grain/needle interfaces with "
         "HomogenizedExponentialCZM in addition (see MATERIAL_MODELS.md). "
         "7 legacy: pre-flag behaviour, unchanged. "
         "Every preset value can be overridden by setting the parameter explicitly.";
}

/// Doc for the SHARED `material_model` parameter, which both
/// ComputeFabricElasticityTensor and OrthotropicPlasticityStressUpdate accept.
/// Set it once, normally in [GlobalParams], and both blocks follow it.
inline std::string
sharedModelDoc()
{
  return "Material model for BOTH the elastic and the plastic response. Set it once, "
         "normally in [GlobalParams], and ComputeFabricElasticityTensor and "
         "OrthotropicPlasticityStressUpdate both follow it, so the two cannot drift "
         "apart. A per-block elastic_model or plastic_model overrides it for that block; "
         "mixing the two is supported but warned about. " + modelDoc();
}

inline bool isBone(Model m) { return static_cast<int>(m) <= 5; }
inline bool isTI(Model m)
{
  return m == Model::TRABECULAR_TI || m == Model::COMPACT_TI;
}
inline bool isIso(Model m)
{
  return m == Model::TRABECULAR_ISO || m == Model::COMPACT_ISO;
}
inline bool isFabric(Model m)
{
  return m == Model::TRABECULAR_FABRIC_TI || m == Model::TRABECULAR_FABRIC_ORTHO;
}
/// models that use main_direction
inline bool usesMainDirection(Model m)
{
  return isTI(m) || m == Model::TRABECULAR_FABRIC_TI;
}

// ============================================================================
// ELASTIC PRESETS
// ============================================================================
struct ElasticPreset
{
  Real E_0 = 0;     ///< E0: isotropic / transverse Young's modulus [MPa]
  Real nu_0 = 0;    ///< V0
  Real G_0 = 0;     ///< MU0 (fabric models); 0 -> E0/(2(1+nu0))
  Real E_a = 0;     ///< EAA: axial Young's modulus (TI)
  Real nu_a = 0;    ///< VA0: compliance S_at = -nu_a / E_a (TI)
  Real G_a = 0;     ///< MUA0: shear modulus of planes containing the axis (TI)
  Real k = 0;       ///< KS: density exponent
  Real l = 0;       ///< LS: fabric exponent (fabric models)
  Real delta = 1.0; ///< DELTA
};

inline ElasticPreset
elasticPreset(Model m)
{
  ElasticPreset e;
  switch (m)
  {
    case Model::TRABECULAR_ISO: // UMAT 367-370
      e.E_0 = 8534.64; e.nu_0 = 0.246; e.k = 1.63;
      break;
    case Model::TRABECULAR_TI: // UMAT 449-455
      e.E_0 = 6561.94; e.E_a = 18661.3; e.nu_0 = 0.323; e.nu_a = 0.32;
      e.G_a = 3739.02; e.k = 1.63;
      break;
    case Model::TRABECULAR_FABRIC_TI:    // UMAT 661-666
    case Model::TRABECULAR_FABRIC_ORTHO: // UMAT 924-929
      e.E_0 = 9994.65; e.nu_0 = 0.2278; e.G_0 = 3361.14; e.k = 1.62; e.l = 1.1;
      break;
    case Model::COMPACT_ISO: // UMAT 1025-1028
      e.E_0 = 19327.0; e.nu_0 = 0.3434; e.k = 1.63;
      break;
    case Model::COMPACT_TI: // UMAT 1108-1114
      e.E_0 = 15079.6; e.E_a = 24578.5; e.nu_0 = 0.4620; e.nu_a = 0.354;
      e.G_a = 6578.01; e.k = 1.63;
      break;
    default:
      break;
  }
  return e;
}

/// Coral: symmetric9 C_ijkl [MPa] = C1111 C1122 C1133 C2222 C2233 C3333 C2323 C1313 C1212
/// From fill_and_generate_configs.py (BULK_MATS). PROVISIONAL, to be confirmed.
inline std::vector<Real>
coralCijkl()
{
  return {171800, 57500, 30200, 106700, 46900, 84200, 42100, 31100, 46600};
}

// ============================================================================
// PLASTIC PRESETS
// ============================================================================
struct PlasticPreset
{
  // yield surface
  Real s0p = 0, s0n = 0;   ///< SIGD0P / SIGD0N  (transverse / isotropic / base)
  Real sap = 0, san = 0;   ///< SIGDAP / SIGDAN  (axial, TI)
  Real zeta0 = 0;          ///< ZETA0
  Real zetaa = 0;          ///< ZETAA0 (TI), referenced to the AXIAL F in the UMAT
  Real tau0 = 0;           ///< TAUD0 (fabric models only; derived otherwise)
  Real taua = 0;           ///< TAUDA0 (TI)
  Real p = 0;              ///< PP: density exponent
  Real q = 0;              ///< QQ: fabric exponent
  Real delta = 1.0;        ///< DELTA

  // post-yield + viscosity (UMAT 321-345, flags per model)
  std::string postyield = "perfect";
  Real rdy = 1.0;          ///< RDY -> initial_yield_ratio
  Real kslope = 0, kmax = 0, kmin = 0, kwidth = 8.0, residual = 0.0;
  std::string viscosity = "linear";
  Real eta = 1e-4;         ///< ETA
  Real m = 1.0;            ///< MM
};

inline PlasticPreset
plasticPreset(Model m)
{
  PlasticPreset s;
  // UMAT post-yield constants common to all bone flags (UMAT 321-345)
  auto bonePostYield = [&s](int pyfl)
  {
    s.rdy = 0.7;      // RDY
    s.kslope = 100.0; // KSLOPE
    s.kmax = 0.015;   // KMAX
    s.kwidth = 8.0;   // KWIDTH
    s.residual = 0.0; // exp_hardening saturates at r = 1 (UMAT PYFL=1)
    s.postyield = (pyfl == 2) ? "simple_softening" : "exp_hardening";
    s.viscosity = "linear"; // VISCFL = 1
    s.eta = 1e-4;           // ETA
    s.m = 1.0;              // MM
  };

  switch (m)
  {
    case Model::TRABECULAR_ISO: // UMAT 376-379, PYFL=2
      s.s0p = 61.43; s.s0n = 89.24; s.zeta0 = 0.1876; s.p = 1.686;
      bonePostYield(2);
      break;
    case Model::TRABECULAR_TI: // UMAT 461-468, PYFL=2
      s.s0p = 50.63; s.s0n = 68.07; s.sap = 98.68; s.san = 165.80;
      s.zeta0 = 0.4707; s.zetaa = 0.1809; s.taua = 51.62; s.p = 1.69;
      bonePostYield(2);
      break;
    case Model::TRABECULAR_FABRIC_TI: // UMAT 672-677, PYFL=2
      s.s0p = 66.012; s.s0n = 98.876; s.zeta0 = 0.2182; s.tau0 = 41.889;
      s.p = 1.686; s.q = 1.05;
      bonePostYield(2);
      break;
    case Model::TRABECULAR_FABRIC_ORTHO: // UMAT 935-940, PYFL=2
      s.s0p = 66.012; s.s0n = 98.876; s.zeta0 = 0.218; s.tau0 = 41.889;
      s.p = 1.69; s.q = 1.05;
      bonePostYield(2);
      break;
    case Model::COMPACT_ISO: // UMAT 1034-1037, PYFL=1
      s.s0p = 144.7; s.s0n = 234.2; s.zeta0 = 0.49; s.p = 1.69;
      bonePostYield(1);
      break;
    case Model::COMPACT_TI: // UMAT 1120-1127, PYFL=1
      s.s0p = 56.54; s.s0n = 201.28; s.sap = 176.40; s.san = 268.00;
      s.zeta0 = 0.0074; s.zetaa = 1.4045; s.taua = 82.55; s.p = 1.686;
      bonePostYield(1);
      break;
    case Model::CORAL: // fill_and_generate_configs.py BULK_MATS, PROVISIONAL
      s.postyield = "exp_softening";
      s.rdy = 1.0;
      s.residual = 0.7;
      s.kslope = 30.0;
      s.kmax = 0.001;
      s.kmin = 0.02;
      s.viscosity = "linear";
      s.eta = 0.002;
      s.m = 0.001;
      break;
    default:
      break;
  }
  return s;
}

/// Coral directional strengths [MPa], PROVISIONAL (fill_and_generate_configs.py).
/// Order: xx, yy, zz tension; compression defaults to tension. Shear: xy, xz, yz.
/// zeta12/13/23 have NO default: zero is not acceptable and the coral value is
/// not known, so the user must set them.
struct CoralStrengths
{
  Real sxx = 4980, syy = 4100, szz = 5340;
  Real txy = 4510, txz = 5080, tyz = 5060;
};

/// UMAT TSFU (UMAT 2026-2035)
inline Real
tsfu(Real rho, Real ex, Real delta)
{
  if (rho <= 0.0)
    return 0.0;
  if (rho <= 0.5)
    return std::pow(rho, ex);
  return std::pow(rho, ex) + (delta - 1.0) * std::pow((rho - 0.5) / 0.5, ex);
}

} // namespace MaterialModelPresets
