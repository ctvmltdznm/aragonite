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
  Real computeInterfaceArea();
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
  const Real _normal_gap_tol;

  // ── Homogenisation ─────────────────────────────────────────────────────────
  // _quality_std_dev        : within-QP contact spread (GH quadrature)
  // _spatial_quality_std_dev: across-QP spatial heterogeneity
  //   Each interface QP draws f_spatial ~ N(1, spatial_quality_std_dev^2) once
  //   (deterministically from element ID + QP index), shifting the mean of its
  //   local GH distribution. This prevents synchronised failure across the interface.
  const Real _quality_std_dev;
  const Real _spatial_quality_std_dev;
  const unsigned int _spatial_random_seed;

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
};