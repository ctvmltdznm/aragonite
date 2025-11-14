//* This file is part of the MOOSE framework
//* https://www.mooseframework.org

#include "OrthotropicPlasticityStressUpdate.h"

registerMooseObject("aragoniteApp", OrthotropicPlasticityStressUpdate);

RankFourTensor OrthotropicPlasticityStressUpdate::_identityTensor(RankFourTensor::initIdentity);

InputParameters
OrthotropicPlasticityStressUpdate::validParams()
{
  InputParameters params = StressUpdateBase::validParams();
  params.addClassDescription("Orthotropic plasticity with direction-dependent yield stresses and "
                            "linear softening behavior for aragonite crystals");

  params.addRequiredParam<Real>("sigma_xx_max", "Yield stress in material xx direction (MPa)");
  params.addRequiredParam<Real>("sigma_yy_max", "Yield stress in material yy direction (MPa)");
  params.addRequiredParam<Real>("sigma_zz_max", "Yield stress in material zz direction (MPa)");
  params.addRequiredParam<Real>("tau_xy_max", "Shear yield stress in material xy plane (MPa)");
  params.addRequiredParam<Real>("tau_xz_max", "Shear yield stress in material xz plane (MPa)");
  params.addRequiredParam<Real>("tau_yz_max", "Shear yield stress in material yz plane (MPa)");
  params.addParam<Real>("hardening_modulus", -40000.0, "Hardening modulus (MPa), negative for softening");
  params.addParam<Real>("max_plastic_strain", 0.02, "Maximum plastic strain for softening");
  params.addParam<Real>("absolute_tolerance", 5.0, "Absolute tolerance (MPa)");
  params.addParam<unsigned int>("max_iterations", 100, "Maximum iterations");

  return params;
}

OrthotropicPlasticityStressUpdate::OrthotropicPlasticityStressUpdate(
    const InputParameters & parameters)
  : StressUpdateBase(parameters),
    _sigma_xx_max(getParam<Real>("sigma_xx_max")),
    _sigma_yy_max(getParam<Real>("sigma_yy_max")),
    _sigma_zz_max(getParam<Real>("sigma_zz_max")),
    _tau_xy_max(getParam<Real>("tau_xy_max")),
    _tau_xz_max(getParam<Real>("tau_xz_max")),
    _tau_yz_max(getParam<Real>("tau_yz_max")),
    _hardening_modulus(getParam<Real>("hardening_modulus")),
    _max_plastic_strain(getParam<Real>("max_plastic_strain")),
    _absolute_tolerance(getParam<Real>("absolute_tolerance")),
    _max_iterations(getParam<unsigned int>("max_iterations")),
    _equivalent_plastic_strain(declareProperty<Real>("equivalent_plastic_strain")),
    _equivalent_plastic_strain_old(getMaterialPropertyOld<Real>("equivalent_plastic_strain")),
    _plastic_strain(declareProperty<RankTwoTensor>("plastic_strain_tensor")),
    _plastic_strain_old(getMaterialPropertyOld<RankTwoTensor>("plastic_strain_tensor")),
    _active_direction(0)
{
}

void
OrthotropicPlasticityStressUpdate::initQpStatefulProperties()
{
  _equivalent_plastic_strain[_qp] = 0.0;
  _plastic_strain[_qp].zero();
}

void
OrthotropicPlasticityStressUpdate::propagateQpStatefulProperties()
{
  _equivalent_plastic_strain[_qp] = _equivalent_plastic_strain_old[_qp];
  _plastic_strain[_qp] = _plastic_strain_old[_qp];
}

void
OrthotropicPlasticityStressUpdate::updateState(
    RankTwoTensor & strain_increment,
    RankTwoTensor & inelastic_strain_increment,
    const RankTwoTensor & /*rotation_increment*/,
    RankTwoTensor & stress_new,
    const RankTwoTensor & stress_old,
    const RankFourTensor & elasticity_tensor,
    const RankTwoTensor & elastic_strain_old,
    bool /*compute_full_tangent_operator*/,
    RankFourTensor & /*tangent_operator*/)
{
  // Get old plastic strain
  Real eps_eq_old = _equivalent_plastic_strain_old[_qp];
  
  // Elastic predictor
  RankTwoTensor elastic_strain_trial = elastic_strain_old + strain_increment;
  stress_new = elasticity_tensor * elastic_strain_trial;

  // Compute softening factor based on OLD plastic strain
  Real eps_for_softening = eps_eq_old;
  if (eps_for_softening > _max_plastic_strain)
    eps_for_softening = _max_plastic_strain;
  
  Real softening_factor = 1.0 + _hardening_modulus * eps_for_softening / 2000.0;
  if (softening_factor < 0.1)
    softening_factor = 0.1;

  // Find active direction
  Real f_xx = std::abs(stress_new(0,0)) / (_sigma_xx_max * softening_factor);
  Real f_yy = std::abs(stress_new(1,1)) / (_sigma_yy_max * softening_factor);
  Real f_zz = std::abs(stress_new(2,2)) / (_sigma_zz_max * softening_factor);
  Real f_xy = std::abs(stress_new(0,1)) / (_tau_xy_max * softening_factor);
  Real f_xz = std::abs(stress_new(0,2)) / (_tau_xz_max * softening_factor);
  Real f_yz = std::abs(stress_new(1,2)) / (_tau_yz_max * softening_factor);

  Real f_max = f_xx;
  _active_direction = 0;
  if (f_yy > f_max) { f_max = f_yy; _active_direction = 1; }
  if (f_zz > f_max) { f_max = f_zz; _active_direction = 2; }
  if (f_xy > f_max) { f_max = f_xy; _active_direction = 3; }
  if (f_xz > f_max) { f_max = f_xz; _active_direction = 4; }
  if (f_yz > f_max) { f_max = f_yz; _active_direction = 5; }

  // Check if elastic
  if (f_max <= 1.0)
  {
    inelastic_strain_increment.zero();
    _equivalent_plastic_strain[_qp] = eps_eq_old;
    _plastic_strain[_qp] = _plastic_strain_old[_qp];
    return;
  }

  // PLASTIC CORRECTION
  
  // Get base yield stress (without softening)
  Real sigma_0 = 0.0;
  switch (_active_direction)
  {
    case 0: sigma_0 = _sigma_xx_max; break;
    case 1: sigma_0 = _sigma_yy_max; break;
    case 2: sigma_0 = _sigma_zz_max; break;
    case 3: sigma_0 = _tau_xy_max; break;
    case 4: sigma_0 = _tau_xz_max; break;
    case 5: sigma_0 = _tau_yz_max; break;
  }
  
  // Current yield stress with OLD softening
  Real sigma_yield_old = sigma_0 * softening_factor;
  
  // Trial stress in active direction
  Real stress_trial = 0.0;
  switch (_active_direction)
  {
    case 0: stress_trial = std::abs(stress_new(0,0)); break;
    case 1: stress_trial = std::abs(stress_new(1,1)); break;
    case 2: stress_trial = std::abs(stress_new(2,2)); break;
    case 3: stress_trial = std::abs(stress_new(0,1)); break;
    case 4: stress_trial = std::abs(stress_new(0,2)); break;
    case 5: stress_trial = std::abs(stress_new(1,2)); break;
  }
  
  // Effective modulus
  Real E_eff = 70000.0;
  switch (_active_direction)
  {
    case 0: E_eff = 140432.0; break;
    case 1: E_eff = 70297.0; break;
    case 2: E_eff = 63413.0; break;
    case 3: E_eff = 46600.0; break;
    case 4: E_eff = 31100.0; break;
    case 5: E_eff = 42100.0; break;
  }
  
  // Plastic strain increment
  Real denominator = E_eff - _hardening_modulus;
  Real delta_eps_plastic = (stress_trial - sigma_yield_old) / denominator;
  
  // Safety limits
  if (delta_eps_plastic < 0.0)
    delta_eps_plastic = 0.0;
  if (delta_eps_plastic > 0.001)
    delta_eps_plastic = 0.001;
  
  // Update equivalent plastic strain
  Real eps_eq_new = eps_eq_old + delta_eps_plastic;
  _equivalent_plastic_strain[_qp] = eps_eq_new;
  
  // CRITICAL: Compute NEW yield stress with UPDATED plastic strain
  Real eps_new_for_softening = eps_eq_new;
  if (eps_new_for_softening > _max_plastic_strain)
    eps_new_for_softening = _max_plastic_strain;
  
  Real softening_factor_new = 1.0 + _hardening_modulus * eps_new_for_softening / 2000.0;
  if (softening_factor_new < 0.1)
    softening_factor_new = 0.1;
  
  Real sigma_yield_new = sigma_0 * softening_factor_new;
  
  // Flow direction
  RankTwoTensor flow_direction;
  flow_direction.zero();
  
  switch (_active_direction)
  {
    case 0:
      flow_direction(0,0) = (stress_new(0,0) > 0.0) ? 1.0 : -1.0;
      break;
    case 1:
      flow_direction(1,1) = (stress_new(1,1) > 0.0) ? 1.0 : -1.0;
      break;
    case 2:
      flow_direction(2,2) = (stress_new(2,2) > 0.0) ? 1.0 : -1.0;
      break;
    case 3:
      flow_direction(0,1) = (stress_new(0,1) > 0.0) ? 1.0 : -1.0;
      flow_direction(1,0) = flow_direction(0,1);
      break;
    case 4:
      flow_direction(0,2) = (stress_new(0,2) > 0.0) ? 1.0 : -1.0;
      flow_direction(2,0) = flow_direction(0,2);
      break;
    case 5:
      flow_direction(1,2) = (stress_new(1,2) > 0.0) ? 1.0 : -1.0;
      flow_direction(2,1) = flow_direction(1,2);
      break;
  }
  
  // Update inelastic strain
  inelastic_strain_increment = delta_eps_plastic * flow_direction;
  _plastic_strain[_qp] = _plastic_strain_old[_qp] + inelastic_strain_increment;
  
  // Update stress - set it to NEW yield stress (with softening!)
  // This is the KEY: stress goes to the NEW (lower) yield surface
  Real sign = (stress_trial > 0.0) ? 1.0 : -1.0;
  
  switch (_active_direction)
  {
    case 0:
      stress_new(0,0) = sigma_yield_new * sign;
      break;
    case 1:
      stress_new(1,1) = sigma_yield_new * sign;
      break;
    case 2:
      stress_new(2,2) = sigma_yield_new * sign;
      break;
    case 3:
      stress_new(0,1) = sigma_yield_new * sign;
      stress_new(1,0) = stress_new(0,1);
      break;
    case 4:
      stress_new(0,2) = sigma_yield_new * sign;
      stress_new(2,0) = stress_new(0,2);
      break;
    case 5:
      stress_new(1,2) = sigma_yield_new * sign;
      stress_new(2,1) = stress_new(1,2);
      break;
  }
}

Real
OrthotropicPlasticityStressUpdate::computeYieldFunction(const RankTwoTensor & stress)
{
  Real eps_eq = _equivalent_plastic_strain[_qp];
  
  if (eps_eq > _max_plastic_strain)
    eps_eq = _max_plastic_strain;
  
  Real softening_factor = 1.0 + _hardening_modulus * eps_eq / 2000.0;
  if (softening_factor < 0.1)
    softening_factor = 0.1;

  Real f_xx = std::abs(stress(0,0)) / (_sigma_xx_max * softening_factor);
  Real f_yy = std::abs(stress(1,1)) / (_sigma_yy_max * softening_factor);
  Real f_zz = std::abs(stress(2,2)) / (_sigma_zz_max * softening_factor);
  Real f_xy = std::abs(stress(0,1)) / (_tau_xy_max * softening_factor);
  Real f_xz = std::abs(stress(0,2)) / (_tau_xz_max * softening_factor);
  Real f_yz = std::abs(stress(1,2)) / (_tau_yz_max * softening_factor);

  Real f_max = f_xx;
  _active_direction = 0;
  if (f_yy > f_max) { f_max = f_yy; _active_direction = 1; }
  if (f_zz > f_max) { f_max = f_zz; _active_direction = 2; }
  if (f_xy > f_max) { f_max = f_xy; _active_direction = 3; }
  if (f_xz > f_max) { f_max = f_xz; _active_direction = 4; }
  if (f_yz > f_max) { f_max = f_yz; _active_direction = 5; }

  Real sigma_yield_active = getCurrentYieldStress(_active_direction);
  return (f_max - 1.0) * sigma_yield_active;
}

Real
OrthotropicPlasticityStressUpdate::getCurrentYieldStress(unsigned int direction)
{
  Real sigma_0 = 0.0;
  switch (direction)
  {
    case 0: sigma_0 = _sigma_xx_max; break;
    case 1: sigma_0 = _sigma_yy_max; break;
    case 2: sigma_0 = _sigma_zz_max; break;
    case 3: sigma_0 = _tau_xy_max; break;
    case 4: sigma_0 = _tau_xz_max; break;
    case 5: sigma_0 = _tau_yz_max; break;
  }

  Real eps_eq = _equivalent_plastic_strain[_qp];
  if (eps_eq > _max_plastic_strain)
    eps_eq = _max_plastic_strain;
  
  Real softening_factor = 1.0 + _hardening_modulus * eps_eq / 2000.0;
  Real sigma_y = sigma_0 * softening_factor;

  Real min_stress = sigma_0 * 0.1;
  if (sigma_y < min_stress)
    sigma_y = min_stress;

  return sigma_y;
}

RankTwoTensor
OrthotropicPlasticityStressUpdate::computeFlowDirection(const RankTwoTensor & stress)
{
  // Not used in current implementation, but required by interface
  RankTwoTensor flow;
  flow.zero();
  return flow;
}
