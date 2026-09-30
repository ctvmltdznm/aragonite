// ComputeFabricElasticityTensor.C
// Anisotropic elasticity with material-model flags. See header.

#include "ComputeFabricElasticityTensor.h"
#include "RotationTensor.h"

registerMooseObject("aragoniteApp", ComputeFabricElasticityTensor);

using MaterialModelPresets::Model;

InputParameters
ComputeFabricElasticityTensor::validParams()
{
  InputParameters params = ComputeElasticityTensorBase::validParams();

  params.addClassDescription(
      "Anisotropic elasticity selected by 'elastic_model': UMAT bone flags 0-5 "
      "(isotropic, transversely isotropic, fabric-based Zysset-Curnier 1995), coral "
      "(full orthotropy from C_ijkl) or legacy fabric elasticity. Optional per-element "
      "rotation from coupled Euler angles.");

  params.addParam<MooseEnum>("elastic_model", MaterialModelPresets::modelEnum(),
                             MaterialModelPresets::modelDoc());

  // --------------------------------------------------------------------------
  // Engineering constants (preset defaults depend on elastic_model)
  // --------------------------------------------------------------------------
  params.addParam<Real>("E_0", "Isotropic / transverse / base Young's modulus E0 [MPa]. "
                               "Required for legacy; preset default otherwise.");
  params.addParam<Real>("nu_0", "Isotropic / transverse / base Poisson's ratio nu0. "
                                "Required for legacy; preset default otherwise.");
  params.addParam<Real>("G_0", "Base shear modulus mu0 [MPa] (fabric models). "
                               "Legacy default: E_0/(2(1+nu_0)).");
  params.addParam<Real>("E_a", "Axial Young's modulus E_a [MPa] (TI models)");
  params.addParam<Real>("nu_a", "Axial Poisson's ratio nu_a (TI models), S_at = -nu_a/E_a: "
                                "transverse strain under axial load");
  params.addParam<Real>("G_a", "Shear modulus of the planes containing the axis [MPa] (TI models)");

  params.addParam<std::vector<Real>>(
      "C_ijkl", "coral: stiffness in symmetric9 order "
                "(C1111 C1122 C1133 C2222 C2233 C3333 C2323 C1313 C1212) [MPa]");

  // --------------------------------------------------------------------------
  // Fabric, density, exponents (defaults are the legacy defaults; presets
  // replace the exponents unless set explicitly)
  // --------------------------------------------------------------------------
  params.addParam<Real>("fabric_m1", 1.0, "Fabric eigenvalue m1 (normalised: m1+m2+m3=3)");
  params.addParam<Real>("fabric_m2", 1.0, "Fabric eigenvalue m2");
  params.addParam<Real>("fabric_m3", 1.0, "Fabric eigenvalue m3");
  params.addParam<Real>("density_rho", 1.0,
                        "Relative density rho = BV/TV, (0,1]. Must be set for bone models.");
  params.addParam<Real>("exponent_k", 2.0, "Density exponent k (UMAT KS). Legacy default 2.0.");
  params.addParam<Real>("exponent_l", 1.0, "Fabric exponent l (UMAT LS). Legacy default 1.0.");
  params.addParam<Real>("delta_cortical", 1.0,
                        "Cortical density correction (UMAT DELTA), only affects rho > 0.5");
  params.addRangeCheckedParam<unsigned int>(
      "main_direction", 1, "main_direction>=1 & main_direction<=3",
      "Material axis of transverse isotropy (1,2,3), UMAT PROPS(6). "
      "Used by trabecular_ti, compact_ti, trabecular_fabric_ti.");

  // --------------------------------------------------------------------------
  // Orientation
  // --------------------------------------------------------------------------
  params.addCoupledVar("coupled_euler_angle_1", "Bunge phi1 [deg], per element (optional)");
  params.addCoupledVar("coupled_euler_angle_2", "Bunge Phi [deg], per element (optional)");
  params.addCoupledVar("coupled_euler_angle_3", "Bunge phi2 [deg], per element (optional)");

  return params;
}

ComputeFabricElasticityTensor::ComputeFabricElasticityTensor(const InputParameters & parameters)
  : ComputeElasticityTensorBase(parameters),
    _model(static_cast<Model>(static_cast<int>(getParam<MooseEnum>("elastic_model")))),
    _E_0(0), _nu_0(0), _G_0(0), _E_a(0), _nu_a(0), _G_a(0),
    _density_rho(getParam<Real>("density_rho")),
    _fabric_m1(getParam<Real>("fabric_m1")),
    _fabric_m2(getParam<Real>("fabric_m2")),
    _fabric_m3(getParam<Real>("fabric_m3")),
    _exponent_k(getParam<Real>("exponent_k")),
    _exponent_l(getParam<Real>("exponent_l")),
    _delta_cortical(getParam<Real>("delta_cortical")),
    _main_direction(getParam<unsigned int>("main_direction")),
    _E1(0), _E2(0), _E3(0), _G12(0), _G13(0), _G23(0), _nu12(0), _nu13(0), _nu23(0),
    _has_coupled_angles(isCoupled("coupled_euler_angle_1") &&
                        isCoupled("coupled_euler_angle_2") &&
                        isCoupled("coupled_euler_angle_3")),
    _euler_angle_1(_has_coupled_angles ? coupledValue("coupled_euler_angle_1") : _zero),
    _euler_angle_2(_has_coupled_angles ? coupledValue("coupled_euler_angle_2") : _zero),
    _euler_angle_3(_has_coupled_angles ? coupledValue("coupled_euler_angle_3") : _zero)
{
  if ((isCoupled("coupled_euler_angle_1") || isCoupled("coupled_euler_angle_2") ||
       isCoupled("coupled_euler_angle_3")) && !_has_coupled_angles)
    paramError("coupled_euler_angle_1", "Provide all three coupled Euler angles or none.");

  if (_density_rho <= 0.0 || _density_rho > 1.0)
    paramError("density_rho", "density_rho must be in (0, 1]");
  if (_fabric_m1 <= 0.0 || _fabric_m2 <= 0.0 || _fabric_m3 <= 0.0)
    mooseError("Fabric eigenvalues must be positive");

  const MooseEnum model_name = getParam<MooseEnum>("elastic_model");
  Moose::out << "\n=== ELASTICITY: elastic_model = " << model_name << " ===\n";

  // ==========================================================================
  // LEGACY: exactly the previous implementation
  // ==========================================================================
  if (_model == Model::LEGACY)
  {
    warnIgnored({"E_a", "nu_a", "G_a", "C_ijkl", "main_direction"});
    if (!isParamValid("E_0"))
      paramError("E_0", "elastic_model = legacy requires E_0");
    if (!isParamValid("nu_0"))
      paramError("nu_0", "elastic_model = legacy requires nu_0");
    _E_0 = getParam<Real>("E_0");
    _nu_0 = getParam<Real>("nu_0");
    _G_0 = isParamValid("G_0") ? getParam<Real>("G_0") : _E_0 / (2.0 * (1.0 + _nu_0));

    if (_E_0 <= 0.0)
      mooseError("E_0 must be positive");
    if (_nu_0 < -1.0 || _nu_0 > 0.5)
      mooseError("nu_0 must be in range (-1, 0.5) for stability");

    const Real fabric_sum = _fabric_m1 + _fabric_m2 + _fabric_m3;
    if (std::abs(fabric_sum - 3.0) > 0.01)
      mooseWarning("Fabric eigenvalues should sum to 3.0, got ", fabric_sum,
                   ". Consider normalizing.");

    const Real rho_k = computeTSFU(_density_rho, _exponent_k, _delta_cortical);
    const Real m1_l = std::pow(_fabric_m1, _exponent_l);
    const Real m2_l = std::pow(_fabric_m2, _exponent_l);
    const Real m3_l = std::pow(_fabric_m3, _exponent_l);

    _E1 = _E_0 * rho_k * std::pow(_fabric_m1, 2.0 * _exponent_l);
    _E2 = _E_0 * rho_k * std::pow(_fabric_m2, 2.0 * _exponent_l);
    _E3 = _E_0 * rho_k * std::pow(_fabric_m3, 2.0 * _exponent_l);
    _G12 = _G_0 * rho_k * m1_l * m2_l;
    _G13 = _G_0 * rho_k * m1_l * m3_l;
    _G23 = _G_0 * rho_k * m2_l * m3_l;
    _nu12 = _nu_0 * std::pow(_fabric_m1 / _fabric_m2, _exponent_l);
    _nu13 = _nu_0 * std::pow(_fabric_m1 / _fabric_m3, _exponent_l);
    _nu23 = _nu_0 * std::pow(_fabric_m2 / _fabric_m3, _exponent_l);
    const Real nu21 = _nu_0 * std::pow(_fabric_m2 / _fabric_m1, _exponent_l);
    const Real nu31 = _nu_0 * std::pow(_fabric_m3 / _fabric_m1, _exponent_l);
    const Real nu32 = _nu_0 * std::pow(_fabric_m3 / _fabric_m2, _exponent_l);

    const Real delta = 1.0 - _nu12 * nu21 - _nu23 * nu32 - nu31 * _nu13 - 2.0 * nu21 * nu32 * _nu13;
    if (delta <= 0.0)
      mooseError("Computed orthotropic constants violate positive definiteness. delta = ", delta,
                 ". Try reducing nu_0 or adjusting fabric values.");

    _C_material = buildOrthotropicStiffness(_E1, _E2, _E3, _G12, _G13, _G23, _nu12, _nu13, _nu23);
  }
  // ==========================================================================
  // CORAL: full orthotropy, every constant explicit
  // ==========================================================================
  else if (_model == Model::CORAL)
  {
    warnIgnored({"E_0", "nu_0", "G_0", "E_a", "nu_a", "G_a", "density_rho", "fabric_m1",
                 "fabric_m2", "fabric_m3", "exponent_k", "exponent_l", "delta_cortical",
                 "main_direction"});
    const std::vector<Real> C = isParamSetByUser("C_ijkl")
                                    ? getParam<std::vector<Real>>("C_ijkl")
                                    : MaterialModelPresets::coralCijkl();
    if (C.size() != 9)
      paramError("C_ijkl", "coral expects 9 constants (symmetric9), got ", C.size());
    _C_material.fillFromInputVector(C, RankFourTensor::symmetric9);

    Moose::out << "C_ijkl (symmetric9, MPa)" << (isParamSetByUser("C_ijkl") ? "" : " [coral default]")
               << ":";
    for (const auto c : C)
      Moose::out << " " << c;
    Moose::out << "\nNote: coral RVEs normally also carry HomogenizedExponentialCZM on grain/needle "
                  "interfaces (see MATERIAL_MODELS.md).\n";
  }
  // ==========================================================================
  // BONE MODELS (UMAT flags 0-5)
  // ==========================================================================
  else
  {
    using namespace MaterialModelPresets;
    const ElasticPreset pre = elasticPreset(_model);

    if (!isParamSetByUser("density_rho"))
      paramError("density_rho", "Bone elastic models require density_rho (BV/TV)");
    if (isFabric(_model) && !(isParamSetByUser("fabric_m1") && isParamSetByUser("fabric_m2") &&
                              isParamSetByUser("fabric_m3")))
      mooseError("elastic_model = ", model_name, " requires fabric_m1, fabric_m2, fabric_m3");

    if (isIso(_model))
      warnIgnored({"E_a", "nu_a", "G_a", "G_0", "fabric_m1", "fabric_m2", "fabric_m3",
                   "exponent_l", "main_direction", "C_ijkl"});
    else if (isTI(_model))
      warnIgnored({"G_0", "fabric_m1", "fabric_m2", "fabric_m3", "exponent_l", "C_ijkl"});
    else if (_model == Model::TRABECULAR_FABRIC_ORTHO)
      warnIgnored({"E_a", "nu_a", "G_a", "main_direction", "C_ijkl"});
    else // fabric TI
      warnIgnored({"E_a", "nu_a", "G_a", "C_ijkl"});

    _E_0 = resolve("E_0", pre.E_0);
    _nu_0 = resolve("nu_0", pre.nu_0);
    _exponent_k = resolve("exponent_k", pre.k);
    _delta_cortical = resolve("delta_cortical", pre.delta);
    const Real G_iso = _E_0 / (2.0 * (1.0 + _nu_0));
    const Real t = computeTSFU(_density_rho, _exponent_k, _delta_cortical);
    const unsigned int a = _main_direction - 1; // axis index
    // plane order 0:(1,2) 1:(1,3) 2:(2,3)
    const unsigned int pl[3][2] = {{0, 1}, {0, 2}, {1, 2}};

    Real E[3], G[3], S[3];

    if (isIso(_model))
    {
      for (unsigned int i = 0; i < 3; ++i)
      {
        E[i] = _E_0 * t;
        G[i] = G_iso * t;
        S[i] = -_nu_0 / (_E_0 * t);
      }
    }
    else if (isTI(_model))
    {
      _E_a = resolve("E_a", pre.E_a);
      _nu_a = resolve("nu_a", pre.nu_a);
      _G_a = resolve("G_a", pre.G_a);
      for (unsigned int i = 0; i < 3; ++i)
        E[i] = (i == a ? _E_a : _E_0) * t;
      for (unsigned int p = 0; p < 3; ++p)
      {
        const bool axial_plane = (pl[p][0] == a || pl[p][1] == a);
        G[p] = (axial_plane ? _G_a : G_iso) * t;
        S[p] = axial_plane ? -_nu_a / (_E_a * t) : -_nu_0 / (_E_0 * t);
      }
    }
    else // fabric models
    {
      _G_0 = resolve("G_0", pre.G_0);
      _exponent_l = resolve("exponent_l", pre.l);
      Real m[3] = {_fabric_m1, _fabric_m2, _fabric_m3};

      const Real fabric_sum = m[0] + m[1] + m[2];
      if (std::abs(fabric_sum - 3.0) > 0.01)
        mooseWarning("Fabric eigenvalues should sum to 3.0, got ", fabric_sum, ".");

      const bool ti = (_model == Model::TRABECULAR_FABRIC_TI);
      if (ti)
      {
        // UMAT: average the two transverse eigenvalues
        const unsigned int b = (a + 1) % 3, c = (a + 2) % 3;
        m[b] = m[c] = 0.5 * (m[b] + m[c]);
        _fabric_m1 = m[0]; _fabric_m2 = m[1]; _fabric_m3 = m[2];
      }
      const Real l = _exponent_l;
      for (unsigned int i = 0; i < 3; ++i)
        E[i] = _E_0 * t * std::pow(m[i], 2.0 * l);
      for (unsigned int p = 0; p < 3; ++p)
      {
        const unsigned int i = pl[p][0], j = pl[p][1];
        const Real mm = std::pow(m[i] * m[j], l);
        // fabric TI: transverse plane uses E0/(2(1+nu0)) so the plane is isotropic
        const bool transverse = ti && i != a && j != a;
        G[p] = (transverse ? G_iso : _G_0) * t * mm;
        S[p] = -_nu_0 / (_E_0 * t * mm);
      }
    }

    _C_material = stiffnessFromCompliance(E, G, S);

    Moose::out << "rho=" << _density_rho << " k=" << _exponent_k << " delta=" << _delta_cortical
               << " TSFU=" << t;
    if (isTI(_model) || _model == Model::TRABECULAR_FABRIC_TI)
      Moose::out << " main_direction=" << _main_direction;
    if (isFabric(_model))
      Moose::out << "\nm (used)=" << _fabric_m1 << " " << _fabric_m2 << " " << _fabric_m3
                 << " l=" << _exponent_l << " G_0=" << _G_0;
    Moose::out << "\nE_0=" << _E_0 << " nu_0=" << _nu_0;
    if (isTI(_model))
      Moose::out << " E_a=" << _E_a << " nu_a=" << _nu_a << " G_a=" << _G_a;
    Moose::out << "\n";
  }

  if (_model != Model::CORAL)
  {
    Moose::out << "E1=" << _E1 << " E2=" << _E2 << " E3=" << _E3 << " MPa\n"
               << "G12=" << _G12 << " G13=" << _G13 << " G23=" << _G23 << " MPa\n"
               << "nu12=" << _nu12 << " nu13=" << _nu13 << " nu23=" << _nu23
               << "  (S_ij = -nu_ij/E_i)\n";
  }
  Moose::out << "Rotation: " << (_has_coupled_angles ? "coupled Euler angles" : "none")
             << "\n=====================================================\n\n";
}

void
ComputeFabricElasticityTensor::computeQpElasticityTensor()
{
  _elasticity_tensor[_qp] = _C_material;

  if (_has_coupled_angles)
  {
    const RealVectorValue euler(_euler_angle_1[_qp], _euler_angle_2[_qp], _euler_angle_3[_qp]);
    const RotationTensor R(euler);
    _elasticity_tensor[_qp].rotate(R);
  }
}

RankFourTensor
ComputeFabricElasticityTensor::stiffnessFromCompliance(const Real E[3], const Real G[3],
                                                        const Real S_off[3])
{
  _E1 = E[0]; _E2 = E[1]; _E3 = E[2];
  _G12 = G[0]; _G13 = G[1]; _G23 = G[2];
  // S_ij = -nu_ij / E_i
  _nu12 = -S_off[0] * E[0];
  _nu13 = -S_off[1] * E[0];
  _nu23 = -S_off[2] * E[1];

  const Real nu21 = _nu12 * _E2 / _E1, nu31 = _nu13 * _E3 / _E1, nu32 = _nu23 * _E3 / _E2;
  const Real delta = 1.0 - _nu12 * nu21 - _nu23 * nu32 - nu31 * _nu13 - 2.0 * nu21 * nu32 * _nu13;
  if (delta <= 0.0 || _E1 <= 0.0 || _E2 <= 0.0 || _E3 <= 0.0 || _G12 <= 0.0 || _G13 <= 0.0 ||
      _G23 <= 0.0)
    mooseError("Elastic constants are not positive definite (delta = ", delta,
               "). Check the overridden parameters.");

  return buildOrthotropicStiffness(_E1, _E2, _E3, _G12, _G13, _G23, _nu12, _nu13, _nu23);
}

RankFourTensor
ComputeFabricElasticityTensor::buildOrthotropicStiffness(Real E1, Real E2, Real E3,
                                                         Real G12, Real G13, Real G23,
                                                         Real nu12, Real nu13, Real nu23) const
{
  // conjugate Poisson's ratios: nu_ij/E_i = nu_ji/E_j
  const Real nu21 = nu12 * E2 / E1;
  const Real nu31 = nu13 * E3 / E1;
  const Real nu32 = nu23 * E3 / E2;
  const Real delta = 1.0 - nu12 * nu21 - nu23 * nu32 - nu31 * nu13 - 2.0 * nu21 * nu32 * nu13;

  const Real C11 = E1 * (1.0 - nu23 * nu32) / delta;
  const Real C22 = E2 * (1.0 - nu13 * nu31) / delta;
  const Real C33 = E3 * (1.0 - nu12 * nu21) / delta;
  const Real C12 = E1 * (nu21 + nu31 * nu23) / delta;
  const Real C13 = E1 * (nu31 + nu21 * nu32) / delta;
  const Real C23 = E2 * (nu32 + nu12 * nu31) / delta;

  // symmetric9: C1111 C1122 C1133 C2222 C2233 C3333 C2323 C1313 C1212
  RankFourTensor C;
  C.fillFromInputVector({C11, C12, C13, C22, C23, C33, G23, G13, G12}, RankFourTensor::symmetric9);
  return C;
}

Real
ComputeFabricElasticityTensor::resolve(const std::string & name, Real preset) const
{
  return isParamSetByUser(name) ? getParam<Real>(name) : preset;
}

void
ComputeFabricElasticityTensor::warnIgnored(const std::vector<std::string> & names) const
{
  for (const auto & n : names)
    if (isParamSetByUser(n))
      mooseWarning("Parameter '", n, "' is ignored for elastic_model = ",
                   getParam<MooseEnum>("elastic_model"));
}

Real
ComputeFabricElasticityTensor::computeTSFU(Real rho, Real exponent, Real delta) const
{
  // UMAT 2026-2035. delta = 1 disables the cortical correction.
  const Real rho_p = std::pow(rho, exponent);
  if (rho > 0.5 && std::abs(delta - 1.0) > 1e-12)
    return rho_p + (delta - 1.0) * std::pow((rho - 0.5) / 0.5, exponent);
  return rho_p;
}
