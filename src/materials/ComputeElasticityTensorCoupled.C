//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* Custom version of ComputeElasticityTensor that supports coupled Euler angles

#include "ComputeElasticityTensorCoupled.h"
#include "RotationTensor.h"

registerMooseObject("aragoniteApp", ComputeElasticityTensorCoupled);

InputParameters
ComputeElasticityTensorCoupled::validParams()
{
  InputParameters params = ComputeElasticityTensor::validParams();
  
  params.addClassDescription("Compute an elasticity tensor with rotation from coupled Euler "
                             "angles, allowing per-element orientations from exodus files");
  
  // Add coupled Euler angles as OPTIONAL parameters with different names
  // to avoid conflict with parent class scalar parameters
  params.addCoupledVar("coupled_euler_angle_1", 
                       "First Euler angle (phi1) in degrees - coupled to AuxVariable for per-element angles");
  params.addCoupledVar("coupled_euler_angle_2", 
                       "Second Euler angle (Phi) in degrees - coupled to AuxVariable for per-element angles");
  params.addCoupledVar("coupled_euler_angle_3", 
                       "Third Euler angle (phi2) in degrees - coupled to AuxVariable for per-element angles");
  
  return params;
}

ComputeElasticityTensorCoupled::ComputeElasticityTensorCoupled(
    const InputParameters & parameters)
  : ComputeElasticityTensor(parameters),
    _has_coupled_angles(isCoupled("coupled_euler_angle_1") && 
                        isCoupled("coupled_euler_angle_2") && 
                        isCoupled("coupled_euler_angle_3")),
    _Euler_angle_1(_has_coupled_angles ? coupledValue("coupled_euler_angle_1") : _zero),
    _Euler_angle_2(_has_coupled_angles ? coupledValue("coupled_euler_angle_2") : _zero),
    _Euler_angle_3(_has_coupled_angles ? coupledValue("coupled_euler_angle_3") : _zero)
{
  // If coupled angles are provided, use them for per-element orientations
  // Otherwise, parent class handles scalar euler_angle_1/2/3 parameters
  if (_has_coupled_angles && isParamValid("euler_angle_1"))
    mooseWarning("Both coupled (coupled_euler_angle_X) and scalar (euler_angle_X) "
                 "Euler angles provided. Using coupled variables (per-element angles).");
}

void
ComputeElasticityTensorCoupled::computeQpElasticityTensor()
{
  // First, call parent to compute the base elasticity tensor
  // This handles all the fill_method logic, material constants, etc.
  ComputeElasticityTensor::computeQpElasticityTensor();
  
  // If we have coupled Euler angles, apply rotation using per-element values
  if (_has_coupled_angles)
  {
    // Get Euler angles at current quadrature point (can be different per element!)
    Real phi1 = _Euler_angle_1[_qp];
    Real Phi = _Euler_angle_2[_qp];
    Real phi2 = _Euler_angle_3[_qp];
    
    // Create rotation tensor from Euler angles (Bunge convention)
    // RotationTensor expects RealVectorValue(phi1, Phi, phi2)
    RealVectorValue euler_angles(phi1, Phi, phi2);
    RotationTensor R(euler_angles);
    
    // Rotate the elasticity tensor from material frame to global frame
    _elasticity_tensor[_qp].rotate(R);
  }
  // If no coupled angles, parent class handles rotation from scalar parameters
}