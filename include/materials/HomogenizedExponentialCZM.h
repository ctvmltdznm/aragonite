#pragma once
#include "CZMComputeLocalTractionTotalBase.h"
#include <random>

class HomogenizedExponentialCZM : public CZMComputeLocalTractionTotalBase
{
public:
  static InputParameters validParams();
  HomogenizedExponentialCZM(const InputParameters & parameters);

protected:
  virtual void initialSetup() override;
  virtual void initQpStatefulProperties() override;
  virtual void computeInterfaceTractionAndDerivatives() override;

  Real computeMixedModeDelta0(Real phi, Real dn, Real dt) const;
  Real computePeakStrengthForContact(Real dn, Real ds, Real dt,
                                     Real sn, Real ss, Real st) const;
  Real solveForFailureDisplacement(Real ratio, Real exponent);

  // ── MD parameters ──────────────────────────────────────────────────────────
  const Real _normal_strength;
  const Real _shear_strength_s;
  const Real _shear_strength_t;
  const Real _delta_0_normal;
  const Real _delta_0_tangent;
  Real       _delta_c;
  const Real _mu;
  const Real _eta;
  const Real _failure_traction_ratio;
  const Real _damage_viscosity;

  // ── Homogenisation ─────────────────────────────────────────────────────────
  // _quality_std_dev        : within-QP contact spread (GH quadrature)
  // _spatial_quality_std_dev: across-QP spatial heterogeneity
  //   Each interface QP draws f_spatial ~ N(1, spatial_quality_std_dev^2) once
  //   (deterministically from element ID + QP index), shifting the mean of its
  //   local GH distribution. This prevents synchronised failure across the interface.
  const Real _quality_std_dev;
  const Real _spatial_quality_std_dev;
  const unsigned int _spatial_random_seed;

  // ── Regularisation widths (were hard-coded 1e-6 and 0.02) ─────────────────
  // _eps_mac: smooth-Macaulay width. Defaults to 1e-3*min(delta_0) so it is a
  //           fixed fraction of the opening scale, not an absolute length.
  // _eps_act: tanh damage-activation width, as a fraction of delta_0. Its
  //           logarithmic gain at onset is 2/_eps_act (~100 by default), which
  //           is what amplified round-off into a 35x-per-step divergence in the
  //           gold comparison. Regression tests should widen it.
  Real       _eps_mac;
  const Real _eps_act;

  // set by solveForFailureDisplacement when the r = 5000 clamp is hit, so the
  // banner can flag that delta_c is not the requested failure_traction_ratio
  bool _delta_c_clamped = false;

  // ── Gauss-Hermite quadrature (5-point, N(0,1)) ────────────────────────────
  static constexpr unsigned int N_QUAD = 5;
  static const Real GH_NODES[N_QUAD];
  static const Real GH_WEIGHTS[N_QUAD];

  // ── Stateful properties per quadrature class ──────────────────────────────
  // _qp_damage : current damage D_g (updated each converged step)
  // _qp_kappa  : history variable = max delta_eff seen by class g
  //              guarantees damage irreversibility under any loading path
  MaterialProperty<std::vector<Real>> & _qp_damage;
  const MaterialProperty<std::vector<Real>> & _qp_damage_old;
  MaterialProperty<std::vector<Real>> & _qp_kappa;
  const MaterialProperty<std::vector<Real>> & _qp_kappa_old;

  // ── Scalar output ──────────────────────────────────────────────────────────
  MaterialProperty<Real> & _damage;
  const MaterialProperty<Real> & _damage_old;
  MaterialProperty<Real> & _delta_eff_prop;
  const MaterialProperty<Real> & _delta_eff_old;

  // ── Diagnostics: the material's own view of the interface state ──────────
  // Not stateful, so they cost one store each per evaluation. They exist so
  // a verification script can compare the code's internal phi, delta_0,
  // kappa and LOCAL jump components with what the postprocessors report,
  // instead of inferring them from side averages.
  MaterialProperty<Real> & _diag_phi;          // czm_mode_ratio
  MaterialProperty<Real> & _diag_delta0;       // czm_delta0      (d0_base, f = 1)
  MaterialProperty<Real> & _diag_kappa;        // czm_kappa       (GH-weighted mean)
  MaterialProperty<Real> & _diag_peak;         // czm_peak_traction (elliptic, f = 1)
  MaterialProperty<Real> & _diag_jump_s;       // czm_jump_s      (LOCAL tangent 1)
  MaterialProperty<Real> & _diag_jump_t;       // czm_jump_t      (LOCAL tangent 2)
};
