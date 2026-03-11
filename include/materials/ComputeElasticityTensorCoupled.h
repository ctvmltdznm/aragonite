//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* Custom version of ComputeElasticityTensor that supports:
//* - Coupled Euler angles (per-element orientations from exodus files)
//* - Optional density scaling for porous materials

#pragma once

#include "ComputeElasticityTensor.h"

/**
 * ComputeElasticityTensorCoupled computes an elasticity tensor with rotation
 * from coupled variables (supporting per-element Euler angles from exodus files).
 * 
 * Also supports optional density-based stiffness scaling:
 *   C_scaled = C_input × ρ^k
 * 
 * This is useful for porous materials like trabecular bone or coral aragonite
 * where the full-density stiffness is known but simulation is at reduced density.
 */
class ComputeElasticityTensorCoupled : public ComputeElasticityTensor
{
public:
  static InputParameters validParams();

  ComputeElasticityTensorCoupled(const InputParameters & parameters);

protected:
  virtual void computeQpElasticityTensor() override;

  /// Coupled Euler angles (can vary per element!)
  const bool _has_coupled_angles;
  const VariableValue & _Euler_angle_1;
  const VariableValue & _Euler_angle_2;
  const VariableValue & _Euler_angle_3;
  
  /// Optional density scaling for porous materials
  const Real _density_rho;                  ///< Relative density (0-1)
  const Real _elasticity_density_exponent;  ///< Density exponent k
  const bool _apply_density_scaling;        ///< Pre-computed flag for efficiency
  Real _density_scale_factor;               ///< Pre-computed ρ^k
};
