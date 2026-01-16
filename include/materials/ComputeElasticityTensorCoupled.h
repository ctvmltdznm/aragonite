//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* Custom version of ComputeElasticityTensor that supports coupled Euler angles
//* This allows reading per-element orientations from exodus files

#pragma once

#include "ComputeElasticityTensor.h"

/**
 * ComputeElasticityTensorCoupled computes an elasticity tensor with rotation
 * from coupled variables (supporting per-element Euler angles from exodus files)
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
};
