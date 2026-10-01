// ComputeFabricElasticityTensor.h
//
// Anisotropic elasticity with material-model flags (elastic_model).
//
//   legacy (default)          Zysset & Curnier (1995) fabric elasticity, exactly as
//                             before: E_0, nu_0 required; exponent_k=2, exponent_l=1,
//                             G_0 = E_0/(2(1+nu_0)) unless given.
//   trabecular_iso / _ti / _fabric_ti / _fabric_ortho, compact_iso / _ti
//                             UMAT bone flags 0-5; constants from MaterialModelPresets.h,
//                             each overridable.
//   coral                     full orthotropy from C_ijkl (symmetric9); no density or
//                             fabric scaling.
//
// The material-frame tensor is built once in the constructor. If coupled Euler angles
// are given it is rotated per quadrature point with MOOSE RotationTensor +
// RankFourTensor::rotate, the same convention ComputeElasticityTensorCoupled uses and
// that OrthotropicPlasticityStressUpdate::updateState is written against.
//
// Zysset-Curnier directional constants (fabric models):
//   E_i = E_0 rho^k m_i^(2l),  G_ij = G_0 rho^k (m_i m_j)^l,  nu_ij = nu_0 (m_i/m_j)^l

#pragma once

#include "ComputeElasticityTensorBase.h"
#include "MaterialModelPresets.h"

class ComputeFabricElasticityTensor : public ComputeElasticityTensorBase
{
public:
  static InputParameters validParams();
  ComputeFabricElasticityTensor(const InputParameters & parameters);

protected:
  virtual void computeQpElasticityTensor() override;

  /// Orthotropic stiffness from engineering constants, convention S_ij = -nu_ij / E_i.
  RankFourTensor buildOrthotropicStiffness(Real E1, Real E2, Real E3,
                                           Real G12, Real G13, Real G23,
                                           Real nu12, Real nu13, Real nu23) const;

  /// Stiffness from moduli E_i, shear moduli G_(12,13,23) and off-diagonal
  /// compliances S_12, S_13, S_23 (checks positive definiteness).
  RankFourTensor stiffnessFromCompliance(const Real E[3], const Real G[3],
                                         const Real S_off[3]);

  /// value of `name` if set by the user, otherwise the preset
  Real resolve(const std::string & name, Real preset) const;
  /// elastic_model if set, else the shared material_model, else legacy.
  /// Non-virtual and base-state only: safe to call from the initialiser list.
  MaterialModelPresets::Model resolveModel() const;
  /// warn about parameters the chosen model ignores
  void warnIgnored(const std::vector<std::string> & names) const;

  /// TSFU: density scaling with cortical correction (UMAT 2026-2035)
  Real computeTSFU(Real rho, Real exponent, Real delta) const;

  const MaterialModelPresets::Model _model;

  // resolved inputs (for reporting)
  Real _E_0, _nu_0, _G_0, _E_a, _nu_a, _G_a;
  const Real _density_rho;
  Real _fabric_m1, _fabric_m2, _fabric_m3;
  Real _exponent_k, _exponent_l, _delta_cortical;
  const unsigned int _main_direction;

  // resolved material-frame engineering constants (for reporting)
  Real _E1, _E2, _E3;
  Real _G12, _G13, _G23;
  Real _nu12, _nu13, _nu23;

  /// material-frame stiffness, built once
  RankFourTensor _C_material;

  /// coupled per-element Euler angles (optional)
  const bool _has_coupled_angles;
  const VariableValue & _euler_angle_1;
  const VariableValue & _euler_angle_2;
  const VariableValue & _euler_angle_3;
};

