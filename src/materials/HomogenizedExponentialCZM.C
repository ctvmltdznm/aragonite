// HomogenizedExponentialCZM.C
// PROOF OF CONCEPT: Area-based homogenization with identical contacts

#include "HomogenizedExponentialCZM.h"
#include "Assembly.h"
#include "MooseRandom.h"

registerMooseObject("aragoniteApp", HomogenizedExponentialCZM);

InputParameters
HomogenizedExponentialCZM::validParams()
{
  InputParameters params = CZMComputeLocalTractionTotalBase::validParams();
  params.addClassDescription(
    "Homogenized exponential CZM representing continuum interface as "
    "collection of MD-scale contacts.");
  
  // ============================================================================
  // MD-DERIVED PARAMETERS (keep MD scale!)
  // ============================================================================
  params.addRequiredParam<Real>("normal_strength", 
                                "Peak normal strength from MD (MPa)");
  params.addRequiredParam<Real>("shear_strength_s", 
                                "Peak shear strength in s-direction from MD (MPa)");
  params.addRequiredParam<Real>("shear_strength_t",
                                "Peak shear strength in t-direction from MD (MPa)");
  params.addParam<Real>("delta_0", "Damage initiation displacement (mm). Used as fallback if delta_0_normal/tangent not specified.");
  params.addParam<Real>("delta_0_normal", "Mode I (normal) separation at peak traction (overrides delta_0)");
  params.addParam<Real>("delta_0_tangent", "Mode II (tangential) separation at peak traction (overrides delta_0)");

  
  // Delta_c options
  params.addParam<Real>("delta_c",
                        "Failure displacement (mm). If not specified, auto-computed.");
  params.addParam<Real>("failure_traction_ratio", 0.01,
                        "Traction ratio for failure threshold (default 0.01)");
  
  // Exponential coefficients
  params.addParam<Real>("mu", 1.0,
                       "Exponential coefficient for loading (default 1.0)");
  params.addParam<Real>("eta", 1.0,
                       "Exponential coefficient for softening (default 1.0)");
  
  // Viscoplastic regularization
  params.addParam<Real>("damage_viscosity", 0.0,
                       "Viscoplastic time scale for damage (seconds). "
                       "0 = instantaneous, >0 = rate-dependent. Typical: 1e-4 to 1e-2 s");
  
  // Other
  params.addParam<Real>("normal_gap_tol", 1e-8, 
                       "Tolerance for normal gap opening");
  params.addParam<unsigned int>("random_seed", 0, 
                                "Seed for random contact distributions (0 = unseeded)");

  // ============================================================================
  // HOMOGENIZATION PARAMETERS
  // ============================================================================
  params.addRequiredParam<Real>("md_contact_area",
                                "Area of single MD contact (mm²). "
                                "Example: 1e-10 mm² = (10 nm)²");
  
  params.addParam<unsigned int>("max_contacts", 10000,
                                "Maximum number of contacts to simulate per QP. "
                                "Caps n_contacts for computational efficiency. "
                                "Default: 10000");
  
  // ============================================================================
  // DISTRIBUTION PARAMETERS
  // ============================================================================
  params.addParam<Real>("strength_std_dev", 0.1,
                       "Standard deviation of strength distribution (as fraction of mean). "
                       "Default 0.1 = 10% coefficient of variation. "
                       "Set to 0 for deterministic (all contacts identical).");

  params.addParam<Real>("initial_damage_max", 0.0,
                       "Maximum initial damage (0 to 1). "
                       "Each contact assigned random damage U(0, max). "
                       "Default 0.0 = no pre-existing defects. "
                       "Typical: 0.2 for 0-20% initial damage.");

  params.addParam<Real>("delta0_std_dev", 0.0,
                       "Standard deviation of δ_0 distribution (as fraction of mean). "
                       "Default 0.0 = all contacts have same δ_0. "
                       "Typical: 0.1 for 10% variation.");
  
return params;
}

HomogenizedExponentialCZM::HomogenizedExponentialCZM(const InputParameters & parameters)
  : CZMComputeLocalTractionTotalBase(parameters),
    _normal_strength(getParam<Real>("normal_strength")),
    _shear_strength_s(getParam<Real>("shear_strength_s")),
    _shear_strength_t(getParam<Real>("shear_strength_t")),
    _delta_0(isParamValid("delta_0") ? getParam<Real>("delta_0") : 0.0),  // 0.0 = not provided
    _delta_0_normal(isParamValid("delta_0_normal") ? getParam<Real>("delta_0_normal") :
                     (isParamValid("delta_0") ? getParam<Real>("delta_0") : -1.0)),
    _delta_0_tangent(isParamValid("delta_0_tangent") ? getParam<Real>("delta_0_tangent") :
                      (isParamValid("delta_0") ? getParam<Real>("delta_0") : -1.0)),
    //_delta_0(getParam<Real>("delta_0")),
    //_delta_0_normal(isParamValid("delta_0_normal") ? getParam<Real>("delta_0_normal") : getParam<Real>("delta_0")),
    //_delta_0_tangent(isParamValid("delta_0_tangent") ? getParam<Real>("delta_0_tangent") : getParam<Real>("delta_0")),
    // _delta_c(getParam<Real>("delta_c")),
    _delta_c(0.0),  // Will be computed
    _mu(getParam<Real>("mu")),
    _eta(getParam<Real>("eta")),
    _failure_traction_ratio(getParam<Real>("failure_traction_ratio")),
    _damage_viscosity(getParam<Real>("damage_viscosity")),
    _normal_gap_tol(getParam<Real>("normal_gap_tol")),
    _md_contact_area(getParam<Real>("md_contact_area")),
    _max_contacts(getParam<unsigned int>("max_contacts")),
    _n_contacts(declareProperty<Real>("n_contacts")),
    _n_contacts_old(getMaterialPropertyOld<Real>("n_contacts")),
    _contact_damage(declareProperty<std::vector<Real>>("contact_damage")),
    _contact_damage_old(getMaterialPropertyOld<std::vector<Real>>("contact_damage")),
    _damage(declareProperty<Real>("damage")),
    _damage_old(getMaterialPropertyOld<Real>("damage")),
    _delta_eff(declareProperty<Real>("delta_eff")),
    _delta_eff_old(getMaterialPropertyOld<Real>("delta_eff")),
    _interface_area(declareProperty<Real>("interface_area")),
    _n_contacts_prop(declareProperty<Real>("n_contacts")),
    _strength_std_dev(getParam<Real>("strength_std_dev")),
    _initial_damage_max(getParam<Real>("initial_damage_max")),
    _delta0_std_dev(getParam<Real>("delta0_std_dev")),
    _random_seed(getParam<unsigned int>("random_seed"))
{
  // Validate MD parameters
  if (!isParamValid("delta_0") && (!isParamValid("delta_0_normal") || !isParamValid("delta_0_tangent")))
  mooseError("Must specify either 'delta_0' OR both 'delta_0_normal' and 'delta_0_tangent'");

  if (_delta_0_normal <= 0.0)
    mooseError("delta_0_normal must be positive");
  if (_delta_0_tangent <= 0.0)
    mooseError("delta_0_tangent must be positive");

  if (_normal_strength <= 0.0)
    mooseError("normal_strength must be positive");
  if (_shear_strength_s <= 0.0 || _shear_strength_t <= 0.0)
    mooseError("shear strengths must be positive");
  if (_mu <= 0.0 || _eta <= 0.0)
    mooseError("mu and eta must be positive");
  if (_damage_viscosity < 0.0)
    mooseError("damage_viscosity must be non-negative");
  
  // Validate homogenization parameters
  if (_md_contact_area <= 0.0)
    mooseError("md_contact_area must be positive (got ", _md_contact_area, ")");
  if (_max_contacts == 0)
    mooseError("max_contacts must be at least 1");
}

void
HomogenizedExponentialCZM::initialSetup()
{
  CZMComputeLocalTractionTotalBase::initialSetup();
  
  // ============================================================================
  // COMPUTE DELTA_C (once, after all parameters known)
  // ============================================================================
  if (isParamValid("delta_c"))
  {
    _delta_c = getParam<Real>("delta_c");
  }
  else
  {
    Real r_failure = solveForFailureDisplacement(_failure_traction_ratio, _eta);
    Real delta_0_for_deltac = std::max({_delta_0_normal, _delta_0_tangent, _delta_0});
    _delta_c = r_failure * delta_0_for_deltac;
  }
  
  // Validate
  Real delta_0_max = std::max({_delta_0_normal, _delta_0_tangent, _delta_0});
  if (_delta_c <= delta_0_max)
    mooseError("delta_c (", _delta_c, ") must be greater than max(delta_0_normal, delta_0_tangent, delta_0) (", delta_0_max, ")");
  
  // Print info
  Moose::out << "\n==========================================================\n"
             << "HomogenizedExponentialCZM Setup\n"
             << "==========================================================\n"
             << "MD-scale parameters:\n"
             << "  delta_0 = " << _delta_0 << " mm (" << _delta_0*1e6 << " nm)\n"
             << "  delta_0_normal = " << _delta_0_normal << " mm (" << _delta_0_normal*1e6 << " nm)\n"
             << "  delta_0_tangent = " << _delta_0_tangent << " mm (" << _delta_0_tangent*1e6 << " nm)\n"
             << "  delta_c = " << _delta_c << " mm (" << _delta_c*1e6 << " nm)\n"
             << "  sigma_n = " << _normal_strength << " MPa\n"
             << "\nHomogenization:\n"
             << "  MD contact area = " << _md_contact_area << " mm² ("
             << std::sqrt(_md_contact_area)*1e6 << " nm × "
             << std::sqrt(_md_contact_area)*1e6 << " nm)\n"
             << "  max_contacts = " << _max_contacts << "\n"
             << "==========================================================\n\n";
}

Real
HomogenizedExponentialCZM::getEffectiveDelta0() const
{
  if (_delta_0 > 0)
    return _delta_0;
  return std::max(_delta_0_normal, _delta_0_tangent);
}

void
HomogenizedExponentialCZM::initQpStatefulProperties()
{
  // ============================================================================
  // COMPUTE INTERFACE AREA AND NUMBER OF CONTACTS
  // ============================================================================
  
  Real area = computeInterfaceArea();
  _interface_area[_qp] = area;
  
  // Compute number of MD contacts
  Real n_real = area / _md_contact_area;
  
  // Apply bounds
  n_real = std::max(1.0, n_real);  // At least 1 contact
  n_real = std::min(n_real, static_cast<Real>(_max_contacts));  // Cap at max
  
  _n_contacts[_qp] = n_real;
  unsigned int n = static_cast<unsigned int>(n_real);
  
  // ============================================================================
  // RESIZE VECTORS AND INITIALIZE
  // ============================================================================
  
  // Resize vectors for this QP
  _contact_damage[_qp].resize(n);
  _contact_strength_n.resize(n);
  _contact_strength_s.resize(n);
  _contact_strength_t.resize(n);
  _contact_delta0.resize(n);
  
  // Initialize average damage and delta_eff
  _damage[_qp] = 0.0;
  _delta_eff[_qp] = 0.0;
  
  // Store n_contacts for output
  _n_contacts_prop[_qp] = n;
  
  // Print info for first QP only
/*  if (_qp == 0)
  {
    Moose::out << "First QP initialization:\n"
               << "  Interface area = " << area << " mm² (" 
               << std::sqrt(area)*1e3 << " μm × " << std::sqrt(area)*1e3 << " μm)\n"
               << "  n_contacts = " << n << " MD contacts\n"
               << "  Represents " << n << " × (" << std::sqrt(_md_contact_area)*1e6 
               << " nm)² sub-interfaces\n\n";
  }
*/
  
  static bool seeded = false;
  if (!seeded && _random_seed != 0) {
    MooseRandom::seed(_random_seed);
    seeded = true;
  }

  // ============================================================================
  // INITIALIZE EACH CONTACT WITH DISTRIBUTIONS
  // ============================================================================
  
  for (unsigned int i = 0; i < n; ++i)
  {
    // ========================================================================
    // INITIALIZE DAMAGE STATE
    // ========================================================================
    _contact_damage[_qp][i] = 0.0;
    
    // ========================================================================
    // ASSIGN STRENGTH FROM GAUSSIAN DISTRIBUTION
    // ========================================================================
    if (_strength_std_dev > 1e-12)
    {
      // Gaussian distribution: N(μ, σ²)
      // Use Box-Muller transform for normal distribution
      Real random_normal = MooseRandom::randNormal(0.0, 1.0);  // Standard normal N(0,1)
      
      // Scale to desired distribution
      Real factor_n = 1.0 + _strength_std_dev * random_normal;
      Real factor_s = 1.0 + _strength_std_dev * MooseRandom::randNormal(0.0, 1.0);
      Real factor_t = 1.0 + _strength_std_dev * MooseRandom::randNormal(0.0, 1.0);
      
      // Clamp to reasonable range [μ - 3σ, μ + 3σ] (99.7% of distribution)
      Real min_factor = 1.0 - 3.0 * _strength_std_dev;
      Real max_factor = 1.0 + 3.0 * _strength_std_dev;
      factor_n = std::max(min_factor, std::min(max_factor, factor_n));
      factor_s = std::max(min_factor, std::min(max_factor, factor_s));
      factor_t = std::max(min_factor, std::min(max_factor, factor_t));
      
      _contact_strength_n[i] = _normal_strength * factor_n;
      _contact_strength_s[i] = _shear_strength_s * factor_s;
      _contact_strength_t[i] = _shear_strength_t * factor_t;
    }
    else
    {
      // Deterministic (all identical)
      _contact_strength_n[i] = _normal_strength;
      _contact_strength_s[i] = _shear_strength_s;
      _contact_strength_t[i] = _shear_strength_t;
    }

    // ========================================================================
    // ASSIGN δ_0 FROM GAUSSIAN DISTRIBUTION (NEW!)
    // ========================================================================
    Real base_delta0 = (_delta_0 > 0) ? _delta_0 : std::max(_delta_0_normal, _delta_0_tangent);

    if (_delta0_std_dev > 1e-12)
    {
      Real factor = 1.0 + _delta0_std_dev * MooseRandom::randNormal(0.0, 1.0);
    
      // Clamp to [0.5, 1.5] - reasonable physical range
      factor = std::max(0.5, std::min(1.5, factor));
    
      _contact_delta0[i] = base_delta0 * factor;
    }
    else
    {
      _contact_delta0[i] = base_delta0;
    }

    // ========================================================================
    // ASSIGN PRE-EXISTING DAMAGE
    // ========================================================================
    if (_initial_damage_max > 1e-12)
    {
      // Uniform distribution [0, max]
      Real D_0_i = _initial_damage_max * MooseRandom::rand();
      _contact_damage[_qp][i] = D_0_i;
    }
    else
    {
      // No pre-existing damage
      _contact_damage[_qp][i] = 0.0;
    }
  }
}

Real
HomogenizedExponentialCZM::computeInterfaceArea()
{
  // For interface elements, the "volume" is actually the area
  // This works for both 2D (line) and 3D (surface) interfaces
  
  // MOOSE provides element volume through Assembly
  Real elem_volume = _assembly.elemVolume();
  
  // Divide by number of QPs to get area per QP
  // (This is an approximation - assumes uniform distribution)
  Real n_qp = _qrule->n_points();
  Real area_per_qp = elem_volume / n_qp;
  
  return area_per_qp;
}

Real
HomogenizedExponentialCZM::computePeakStrength(const Real delta_n, 
                                               const Real delta_s, 
                                               const Real delta_t)
{
  Real delta_total = std::sqrt(delta_n*delta_n + delta_s*delta_s + delta_t*delta_t);
  
  if (delta_total < 1e-12)
    return _normal_strength;
  
  // Mohr-Coulomb mixed-mode
  Real ratio_n = std::abs(delta_n) / delta_total;
  Real ratio_s = std::abs(delta_s) / delta_total;
  Real ratio_t = std::abs(delta_t) / delta_total;
  
  Real T_peak_sq = std::pow(_normal_strength * ratio_n, 2.0) +
                   std::pow(_shear_strength_s * ratio_s, 2.0) +
                   std::pow(_shear_strength_t * ratio_t, 2.0);
  
  return std::sqrt(T_peak_sq);
}

Real
HomogenizedExponentialCZM::computeDamage(Real delta_eff)
{
  if (delta_eff < getEffectiveDelta0())
  //if (delta_eff < _delta_0)
    return 0.0;
  
  if (delta_eff >= _delta_c)
    return 1.0;
  
  Real r = delta_eff / getEffectiveDelta0(); //_delta_0;
  
  if (_mu == 1.0 && _eta == 1.0)
  {
    return 1.0 - std::exp(1.0 - r);
  }
  else
  {
    Real exponent = std::pow(r, _eta);
    return 1.0 - std::exp(1.0 - exponent) * exponent / r;
  }
}

Real
HomogenizedExponentialCZM::computeMixedModeDelta0(Real phi, Real delta_0_n, Real delta_0_t) const
{
  // Wang 2025 Eq. 25: Mixed-mode damage initiation threshold
  // δ₀m = (δ₀ⁿ δ₀ᵗ √(1+φ²)) / √((δ₀ᵗ)² + φ²(δ₀ⁿ)²)
  //
  // Special cases:
  //   φ=0 (pure tension): δ₀m = δ₀ⁿ
  //   φ→∞ (pure shear):   δ₀m = δ₀ᵗ
  //   δ₀ⁿ = δ₀ᵗ (isotropic): δ₀m = δ₀ⁿ = δ₀ᵗ
  
  // Handle pure modes
  if (phi < 1e-12)  // Pure normal
    return delta_0_n;
  
  if (phi > 1e6)  // Pure shear
    return delta_0_t;
  
  // General mixed mode
  Real phi_sq = phi * phi;
  Real delta_0_t_sq = delta_0_t * delta_0_t;
  Real delta_0_n_sq = delta_0_n * delta_0_n;
  
  Real numerator = delta_0_n * delta_0_t * std::sqrt(1.0 + phi_sq);
  Real denominator = std::sqrt(delta_0_t_sq + phi_sq * delta_0_n_sq);
  
  return numerator / denominator;
}

Real
HomogenizedExponentialCZM::computeDamageDeriv(Real delta_eff, Real delta_0_use, Real phi)
{
  // Use provided delta_0, or fall back to member variable if not provided
  Real delta_0_actual = (delta_0_use > 1e-12) ? delta_0_use : getEffectiveDelta0(); //_delta_0;
  
  // Always use Wang 2025 Eq. 25
  // If isotropic (delta_0_normal == delta_0_tangent), this automatically
  // gives δ₀m = δ₀ (constant for all φ)
  Real delta_0_threshold = computeMixedModeDelta0(phi, _delta_0_normal, _delta_0_tangent);
  
  // For damage evolution rate r = δ_eff / δ₀_mixed
  // Use mode-dependent scaling
  Real phi_sq = phi * phi;
  Real delta_0_mixed = std::sqrt(_delta_0_tangent * _delta_0_tangent + 
                                 phi_sq * _delta_0_normal * _delta_0_normal) / 
                       std::sqrt(1.0 + phi_sq);

  Real delta_c_scaled = _delta_c;

  // Check bounds using Wang Eq. 25 threshold
  if (delta_eff < delta_0_threshold || delta_eff >= delta_c_scaled)
    return 0.0;

  // Compute r using mode-dependent delta_0
  Real r = delta_eff / delta_0_mixed;
  
  if (_mu == 1.0 && _eta == 1.0)
  {
    return (1.0 / delta_0_mixed) * std::exp(1.0 - r);
  }
  else
  {
    Real exponent = std::pow(r, _eta);
    Real d_exponent = _eta * std::pow(r, _eta - 1.0) / delta_0_mixed;  // Use delta_0_mixed insteaf of _actual for mixed-mode
    
    Real term1 = -std::exp(1.0 - exponent) * d_exponent * exponent / r;
    Real term2 = std::exp(1.0 - exponent) * d_exponent / r;
    Real term3 = -std::exp(1.0 - exponent) * exponent / (r*r) / delta_0_mixed;  // Use delta_0_mixed
    
    return term1 + term2 + term3;
  }
}
/*
void
HomogenizedExponentialCZM::computeDamageGradient(Real delta_eff,
                                                 Real delta_0_contact,
                                                 Real delta_n_pos,
                                                 Real delta_s,
                                                 Real delta_t,
                                                 Real phi,
                                                 Real dD_grad[3])
{
  // Compute full gradient dD/d(δₙ, δₛ, δₜ) including phi dependence
  // 
  // dD/dδᵢ = ∂D/∂δ_eff · ∂δ_eff/∂δᵢ + ∂D/∂δ₀_mixed · ∂δ₀_mixed/∂φ · ∂φ/∂δᵢ
  
  // Initialize to zero
  dD_grad[0] = 0.0;
  dD_grad[1] = 0.0;
  dD_grad[2] = 0.0;
  
  if (delta_eff < 1e-12)
    return;

  Real delta_0_actual = (delta_0_contact > 1e-12) ? delta_0_contact : getEffectiveDelta0();
  
  // Use Wang Eq. 25 for threshold
  Real delta_0_threshold = computeMixedModeDelta0(phi, _delta_0_normal, _delta_0_tangent);
  
  if (delta_eff < delta_0_threshold || delta_eff >= _delta_c)
    return RealVectorValue(0, 0, 0);
  
  // For evolution rate
  Real phi_sq = phi * phi;
  Real delta_0_mixed = std::sqrt(_delta_0_tangent * _delta_0_tangent + 
                                 phi_sq * _delta_0_normal * _delta_0_normal) / 
                       std::sqrt(1.0 + phi_sq);
  
  // --- TERM 1: ∂D/∂δ_eff · ∂δ_eff/∂δᵢ (current implementation) ---
  Real dD_ddelta_eff = computeDamageDeriv(delta_eff, delta_0_contact, phi);
  
  // ∂δ_eff/∂δₙ = δₙ/δ_eff
  // ∂δ_eff/∂δₛ = δₛ/δ_eff
  // ∂δ_eff/∂δₜ = δₜ/δ_eff
  Real ddelta_eff_dn = delta_n_pos / delta_eff;
  Real ddelta_eff_ds = delta_s / delta_eff;
  Real ddelta_eff_dt = delta_t / delta_eff;
  
  dD_grad[0] = dD_ddelta_eff * ddelta_eff_dn;
  dD_grad[1] = dD_ddelta_eff * ddelta_eff_ds;
  dD_grad[2] = dD_ddelta_eff * ddelta_eff_dt;
  
  // --- TERM 2: ∂D/∂δ₀_mixed · ∂δ₀_mixed/∂φ · ∂φ/∂δᵢ (NEW!) ---
  
  if (delta_n_pos < 1e-12)
    return;  // Can't compute phi derivatives when δₙ=0
  
  // ∂D/∂δ₀_mixed
  Real dD_ddelta0 = computeDamageDerivWrtDelta0(delta_eff, delta_0_mixed);
  
  // ∂δ₀_mixed/∂φ = -δ₀ · φ / (1 + φ²)^(3/2)
  Real one_plus_phi2 = 1.0 + phi*phi;
  Real ddelta0_mixed_dphi = -delta_0_actual * phi / 
                            std::pow(one_plus_phi2, 1.5);
  
  // ∂φ/∂δₙ = -δₜ/δₙ² = -φ/δₙ
  // ∂φ/∂δₛ = 0  (phi only depends on tangential magnitude, not components)
  // ∂φ/∂δₜ = 1/δₙ  (assuming δₜ = sqrt(δₛ² + δₜ²) approximation)
  
  // More precisely, if δₜ_mag = sqrt(δₛ² + δₜ²), then:
  // ∂φ/∂δₛ = ∂(δₜ_mag/δₙ)/∂δₛ = (δₛ/δₜ_mag)/δₙ
  // ∂φ/∂δₜ = ∂(δₜ_mag/δₙ)/∂δₜ = (δₜ/δₜ_mag)/δₙ
  
  Real delta_t_mag = std::sqrt(delta_s*delta_s + delta_t*delta_t);
  
  Real dphi_dn = -phi / delta_n_pos;
  Real dphi_ds = (delta_t_mag > 1e-12) ? (delta_s / delta_t_mag) / delta_n_pos : 0.0;
  Real dphi_dt = (delta_t_mag > 1e-12) ? (delta_t / delta_t_mag) / delta_n_pos : 0.0;
  
  // Add contribution from phi dependence
  Real factor = dD_ddelta0 * ddelta0_mixed_dphi;
  
  dD_grad[0] += factor * dphi_dn;
  dD_grad[1] += factor * dphi_ds;
  dD_grad[2] += factor * dphi_dt;
}
*/

Real
HomogenizedExponentialCZM::computeDamageDerivWrtDelta0(Real delta_eff, Real delta_0_mixed, Real eta)
{
  // Compute ∂D/∂δ₀_mixed holding δ_eff constant
  // D = 1 - r^(η-1) · exp(1 - r^η)  where r = δ_eff/δ₀
  // ∂D/∂δ₀ = (∂D/∂r) · (∂r/∂δ₀) = (dD/dr) · (-δ_eff/δ₀²)
  
  // Note: This function receives delta_0_mixed already computed
  // We should check against original delta_0, but we only have delta_0_mixed here
  // To maintain consistency, we need to back-calculate delta_0_actual
  // For now, assume delta_0_actual ≈ delta_0_mixed * sqrt(1+phi²) for rough check
  // Better solution: pass delta_0_actual as separate parameter

  if (delta_eff >= _delta_c)  // Only check upper bound
//  if (delta_eff < delta_0_mixed || delta_eff >= _delta_c)
    return 0.0;

  Real r = delta_eff / delta_0_mixed;
  
  // Compute dD/dr
  Real dD_dr;
  if (_mu == 1.0 && eta == 1.0)
  {
    dD_dr = std::exp(1.0 - r);
  }
  else
  {
    Real r_eta = std::pow(r, eta);
    Real exp_term = std::exp(1.0 - r_eta);
    Real r_eta_minus_2 = (r > 1e-12) ? std::pow(r, eta - 2.0) : 0.0;
    
    dD_dr = r_eta_minus_2 * exp_term * ((eta - 1.0) - eta * r_eta);
  }
  
  // ∂r/∂δ₀ = -δ_eff/δ₀²
  Real dr_ddelta0 = -delta_eff / (delta_0_mixed * delta_0_mixed);
  
  return dD_dr * dr_ddelta0;
}

Real
HomogenizedExponentialCZM::solveForFailureDisplacement(Real traction_ratio, Real exponent)
{
  Real r = 5.0;  // Initial guess
  const Real e = std::exp(1.0);
  const int max_iter = 50;
  const Real tol = 1e-8;
  
  for (int iter = 0; iter < max_iter; ++iter)
  {
    Real r_exp = std::pow(r, exponent);
    Real f = e * r_exp * std::exp(-r_exp) - traction_ratio;
    
    if (std::abs(f) < tol)
      break;
    
    Real dr_exp = exponent * std::pow(r, exponent - 1.0);
    Real df = e * dr_exp * std::exp(-r_exp) * (1.0 - r_exp);
    
    if (std::abs(df) < 1e-12)
      mooseError("Cannot solve for failure displacement");
    
    Real r_new = r - f / df;
    
    if (r_new < 1.0)
      r_new = (r + 1.0) / 2.0;
    if (r_new > 20.0)
      r_new = 20.0;
    
    r = r_new;
  }
  
  if (r <= 1.0)
    mooseError("Failed to solve for failure displacement");
  
  return r;
}

void
HomogenizedExponentialCZM::computeInterfaceTractionAndDerivatives()
{
  // Get displacement jumps in LOCAL coordinates
  const Real delta_n = _interface_displacement_jump[_qp](0);
  const Real delta_s = _interface_displacement_jump[_qp](1);
  const Real delta_t = _interface_displacement_jump[_qp](2);
  
  // Apply Macaulay bracket
  Real delta_n_pos = (delta_n > _normal_gap_tol) ? delta_n : 0.0;

  // Compute tangential magnitude
  Real delta_t_mag = std::sqrt(delta_s*delta_s + delta_t*delta_t);
  
  // Effective opening (SAME for all contacts)
  Real delta_eff_val = std::sqrt(delta_n_pos*delta_n_pos + delta_t_mag*delta_t_mag);
  _delta_eff[_qp] = delta_eff_val;
  
  // Mode ratio
  Real phi = 0.0;
  if (delta_eff_val > 1e-8)
  {
    if (delta_n_pos > 1e-8)
      phi = delta_t_mag / delta_n_pos;
    else
      phi = 100.0;  // Pure shear

    // Cap phi to prevent numerical issues
    if (phi > 100.0) phi = 100.0;
  }

  // Peak strength (SAME for all contacts in this POC)
  Real T_peak = computePeakStrength(delta_n_pos, delta_s, delta_t);
  
  // ============================================================================
  // LOOP OVER ALL MD CONTACTS - COMPUTE AVERAGE
  // ============================================================================
  const unsigned int n = static_cast<unsigned int>(_n_contacts[_qp]);
  Real T_undamaged_sum = 0.0;
  Real dT_dd_damaged_sum[3] = {0.0, 0.0, 0.0};  // Direction-dependent!
  Real D_sum = 0.0;
  
  const Real e = std::exp(1.0);
  // Real r = delta_eff_val / _delta_0;
  
  for (unsigned int i = 0; i < n; ++i)
  {
    Real r = delta_eff_val / _contact_delta0[i];
    // ========================================================================
    // COMPUTE PEAK STRENGTH FOR THIS CONTACT (using its specific strengths)
    // ========================================================================
    Real T_peak_i = computePeakStrengthForContact(delta_n_pos, delta_s, delta_t, 
                                                   _contact_strength_n[i],
                                                   _contact_strength_s[i],
                                                   _contact_strength_t[i]);

    // ========================================================================
    // COMPUTE UNDAMAGED TRACTION FOR THIS CONTACT
    // ========================================================================
    Real T_undamaged_i = 0.0;
    Real dT_dd_undamaged_i = 0.0;
    
    if (delta_eff_val > 1e-12)
    {
      Real exponent = 0.0;
      
      if (r < 1.0)
      {
        // Loading phase
        exponent = std::pow(r, _mu);
        T_undamaged_i = T_peak_i * e * exponent * std::exp(-exponent);
        dT_dd_undamaged_i = T_peak_i * e / _contact_delta0[i] *   // Use contact-specific delta_0
                           _mu * std::pow(r, _mu - 1.0) * (1.0 - exponent) * std::exp(-exponent);
      }
      else
      {
        // Softening phase
        exponent = std::pow(r, _eta);
        T_undamaged_i = T_peak_i * e * exponent * std::exp(-exponent);
        dT_dd_undamaged_i = T_peak_i * e / _contact_delta0[i] *   // Use contact-specific delta_0
                           _eta * std::pow(r, _eta - 1.0) * (1.0 - exponent) * std::exp(-exponent);
      }
    }

    // ========================================================================
    // COMPUTE DAMAGE FOR THIS CONTACT
    // ========================================================================
    Real D_target = computeDamageForContact(delta_eff_val, _contact_delta0[i], phi);

    /* OLD 1D INCOMPLETE JACOBIAN -- NON DIRECTIONAL
    // Compute base derivative (includes phi in delta_0_mixed but not derivative of phi)
    Real dD_dd_base = computeDamageDeriv(delta_eff_val, _contact_delta0[i], phi);
    
    // Add correction for changing phi
    Real dD_dd_total = dD_dd_base;
    
    if (delta_n_pos > 1e-12 && delta_eff_val > 1e-12)
    // Compute delta_0 threshold and mixed for current contact
    {
      Real delta_0_actual = _contact_delta0[i];
      
      // Use Wang Eq. 25 for threshold
      Real delta_0_threshold = computeMixedModeDelta0(phi, _delta_0_normal, _delta_0_tangent);
      
      // For evolution
      Real phi_sq = phi * phi;
      Real delta_0_mixed = std::sqrt(_delta_0_tangent * _delta_0_tangent + 
                                     phi_sq * _delta_0_normal * _delta_0_normal) / 
                           std::sqrt(1.0 + phi_sq);
      
      if (delta_eff_val >= delta_0_threshold && delta_eff_val < _delta_c)      
      {
        // ∂D/∂δ₀_mixed
        Real dD_ddelta0 = computeDamageDerivWrtDelta0(delta_eff_val, delta_0_mixed, _eta);
        
        // ∂δ₀_mixed/∂φ = -δ₀·φ / (1+φ²)^(3/2)
        Real one_plus_phi2 = 1.0 + phi*phi;
        Real ddelta0_dphi = -delta_0_actual * phi / std::pow(one_plus_phi2, 1.5);
        
        // ∂φ/∂δ along loading direction (scalar approximation)
        // For radial loading: ∂(δₜ/δₙ)/∂δ = 0
        // For non-radial (curved interface): ∂φ/∂δ ≈ |δₜ|/(δₙ·δ)
        Real delta_t_mag_val = std::sqrt(delta_s*delta_s + delta_t*delta_t);
        Real dphi_approx = delta_t_mag_val / (delta_n_pos * delta_eff_val);
        
        // Chain rule: ∂D/∂δ₀ · ∂δ₀/∂φ · ∂φ/∂δ
        Real dD_correction = dD_ddelta0 * ddelta0_dphi * dphi_approx;
        
        dD_dd_total = dD_dd_base + dD_correction;
      }
    }
    */
    
    // ========================================================================
    // COMPUTE FULL 3D DAMAGE DERIVATIVE (including phi-dependence)
    // ========================================================================
    
    Real dD_dd[3] = {0.0, 0.0, 0.0};  // [normal, shear_s, shear_t]
    
    if (delta_eff_val > 1e-12)
    {
      // ------------------------------------------------------------------------
      // PART 1: Base derivative (from delta_eff changing)
      // ∂D/∂δ_eff · ∂δ_eff/∂δᵢ
      // ------------------------------------------------------------------------
      
      Real dD_ddelta_eff = computeDamageDeriv(delta_eff_val, _contact_delta0[i], phi);
      
      // ∂δ_eff/∂δᵢ = δᵢ/δ_eff
      Real ddelta_eff_dd[3];
      ddelta_eff_dd[0] = delta_n_pos / delta_eff_val;  // ∂δ_eff/∂δₙ
      ddelta_eff_dd[1] = delta_s / delta_eff_val;      // ∂δ_eff/∂δₛ
      ddelta_eff_dd[2] = delta_t / delta_eff_val;      // ∂δ_eff/∂δₜ
      
      for (int dir = 0; dir < 3; dir++)
        dD_dd[dir] = dD_ddelta_eff * ddelta_eff_dd[dir];
      
      // ------------------------------------------------------------------------
      // PART 2: Phi correction (from mode ratio changing)
      // ∂D/∂δ₀_mixed · ∂δ₀_mixed/∂φ · ∂φ/∂δᵢ
      // ------------------------------------------------------------------------
      
      if (delta_n_pos > 1e-8 && delta_eff_val > 1e-12)
      {
        Real delta_0_actual = _contact_delta0[i];
        
        // Always use Wang Eq. 25
        Real delta_0_threshold = computeMixedModeDelta0(phi, _delta_0_normal, _delta_0_tangent);
        
        // For evolution
        Real phi_sq = phi * phi;
        Real delta_0_mixed = std::sqrt(_delta_0_tangent * _delta_0_tangent + 
                                       phi_sq * _delta_0_normal * _delta_0_normal) / 
                             std::sqrt(1.0 + phi_sq);
        
        // Only add correction if in damage regime
        if (delta_eff_val >= delta_0_threshold && delta_eff_val < _delta_c)
        {
          // ∂D/∂δ₀_mixed
          Real dD_ddelta0 = computeDamageDerivWrtDelta0(delta_eff_val, delta_0_mixed, _eta);
          
          // ∂δ₀_mixed/∂φ = -δ₀·φ / (1+φ²)^(3/2)
          Real one_plus_phi2 = 1.0 + phi*phi;
          Real ddelta0_dphi = -delta_0_actual * phi / std::pow(one_plus_phi2, 1.5);
          
          // Chain: ∂D/∂δ₀ · ∂δ₀/∂φ
          Real dD_dphi = dD_ddelta0 * ddelta0_dphi;
          
          // ∂φ/∂δᵢ where φ = |δₜ|/δₙ
          Real delta_t_mag_val = std::sqrt(delta_s*delta_s + delta_t*delta_t);
          
          if (delta_t_mag_val > 1e-12)
          {
            // ∂φ/∂δₙ = -φ/δₙ
            Real dphi_ddn = -phi / delta_n_pos;
            
            // ∂φ/∂δₛ = (δₛ/|δₜ|)/δₙ
            Real dphi_dds = (delta_s / delta_t_mag_val) / delta_n_pos;
            
            // ∂φ/∂δₜ = (δₜ/|δₜ|)/δₙ
            Real dphi_ddt = (delta_t / delta_t_mag_val) / delta_n_pos;
            
            // Add corrections
            dD_dd[0] += dD_dphi * dphi_ddn;
            dD_dd[1] += dD_dphi * dphi_dds;
            dD_dd[2] += dD_dphi * dphi_ddt;
          }
          // else: pure normal loading, no phi correction needed
        }
      }
    }

    // Apply viscoplasticity
    Real D_i;
    Real dD_dd_i[3];  // Now an array
    // Real dD_dd_i;  // Scalar
    
    if (_damage_viscosity > 1e-12)
    {
      Real D_rate = (D_target - _contact_damage_old[_qp][i]) / _damage_viscosity;
      D_i = _contact_damage_old[_qp][i] + D_rate * _dt;
      D_i = std::max(0.0, std::min(1.0, D_i));

      // Scale all derivatives by viscosity factor
      Real visc_factor = _dt / _damage_viscosity;
      for (int dir = 0; dir < 3; dir++)
        dD_dd_i[dir] = dD_dd[dir] * visc_factor;
      
      // dD_dd_i = dD_dd_total * _dt / _damage_viscosity; // Old 1D
    }
    else
    {
      D_i = D_target;
      for (int dir = 0; dir < 3; dir++)
        dD_dd_i[dir] = dD_dd[dir];
      // dD_dd_i = dD_dd_total; // Old 1D
    }
    
    _contact_damage[_qp][i] = D_i;
    
    // ========================================================================
    // DAMAGED TRACTION DERIVATIVE (Old 1D)
    // ========================================================================
/*    Real dT_dd_damaged_i = (1.0 - D_i) * dT_dd_undamaged_i - T_undamaged_i * dD_dd_i;
    
    // Accumulate
    T_undamaged_sum += T_undamaged_i;
    dT_dd_damaged_sum += dT_dd_damaged_i;
    D_sum += D_i;
    */
    
    // ========================================================================
    // DAMAGED TRACTION DERIVATIVE (direction-dependent!)
    // ========================================================================
    
    // For each direction: dTⁱ/dδⱼ (diagonal approximation)
    // Normal direction gets normal stiffness, shear directions get shear stiffness
    Real K_elastic[3];
    K_elastic[0] = _normal_strength / _contact_delta0[i];
    K_elastic[1] = _shear_strength_s / _contact_delta0[i];
    K_elastic[2] = _shear_strength_t / _contact_delta0[i];
    
    Real dT_dd_damaged_i[3];
    for (int dir = 0; dir < 3; dir++)
    {
      // Damaged traction derivative: mix of elastic stiffness and damage evolution
      Real dT_dd_undamaged_dir = (delta_eff_val > 1e-12) ? dT_dd_undamaged_i : K_elastic[dir];
      dT_dd_damaged_i[dir] = (1.0 - D_i) * dT_dd_undamaged_dir - T_undamaged_i * dD_dd_i[dir];
    }
    
    // Accumulate
    T_undamaged_sum += T_undamaged_i;
    for (int dir = 0; dir < 3; dir++)
      dT_dd_damaged_sum[dir] += dT_dd_damaged_i[dir];
    D_sum += D_i;
  }  // End loop over contacts
  
  // ============================================================================
  // COMPUTE AVERAGES
  // ============================================================================
  Real T_undamaged_avg = T_undamaged_sum / static_cast<Real>(n);
  //Real dT_dd_undamaged_avg = dT_dd_undamaged_sum / static_cast<Real>(n);
  // Real dT_dd_total = dT_dd_damaged_sum / static_cast<Real>(n);  // Already includes damage effects -- Old 1D
  Real dT_dd_total[3];  // Now direction-dependent!
  for (int dir = 0; dir < 3; dir++)
    dT_dd_total[dir] = dT_dd_damaged_sum[dir] / static_cast<Real>(n);

  Real D_avg = D_sum / static_cast<Real>(n);
  
  // Store average damage (for AuxVariable output)
  _damage[_qp] = D_avg;
  
  // Store n_contacts as material property for output
  _n_contacts_prop[_qp] = n;
  
  // ============================================================================
  // APPLY AVERAGE DAMAGE TO AVERAGE TRACTION
  // ============================================================================
  Real T_damaged = T_undamaged_avg * (1.0 - D_avg);
  
  // Compute traction vector
  if (delta_eff_val > 1e-12)
  {
    Real ratio_n = delta_n_pos / delta_eff_val;
    Real ratio_s = delta_s / delta_eff_val;
    Real ratio_t = delta_t / delta_eff_val;
    
    _interface_traction[_qp](0) = T_damaged * ratio_n;
    _interface_traction[_qp](1) = T_damaged * ratio_s;
    _interface_traction[_qp](2) = T_damaged * ratio_t;
  }
  else
  {
    _interface_traction[_qp].zero();
  }
  
  // ============================================================================
  // COMPUTE JACOBIAN (using direction-dependent derivatives)
  // ============================================================================
  _dinterface_traction_djump[_qp].zero();
  
  if (delta_eff_val > 1e-12)
  {
    // Use direction-dependent stiffness for better accuracy
    // Diagonal approximation: each direction uses its own derivative
    // (No need to apply damage formula again - it's already in dT_dd_total)
    
    Real coef1 = T_damaged / delta_eff_val;

    for (unsigned int i = 0; i < 3; ++i)
    {
      for (unsigned int j = 0; j < 3; ++j)
      {
        Real jump_i = (i == 0) ? delta_n_pos : _interface_displacement_jump[_qp](i);
        Real jump_j = (j == 0) ? delta_n_pos : _interface_displacement_jump[_qp](j);
        
        // Use direction-dependent derivative for diagonal terms
        Real coef2_j = (dT_dd_total[j] - coef1) / (delta_eff_val * delta_eff_val);
        
        _dinterface_traction_djump[_qp](i, j) = 
          coef1 * (i == j ? 1.0 : 0.0) + coef2_j * jump_i * jump_j;
      }
    }
  }
  else
  {
    // Elastic stiffness (mode-dependent)
    _dinterface_traction_djump[_qp](0, 0) = _normal_strength / _delta_0_normal;      // Normal
    _dinterface_traction_djump[_qp](1, 1) = _shear_strength_s / _delta_0_tangent;    // Shear-s
    _dinterface_traction_djump[_qp](2, 2) = _shear_strength_t / _delta_0_tangent;    // Shear-t

    //Real K_init = T_peak / _delta_0;
    //for (unsigned int i = 0; i < 3; ++i)
    //  _dinterface_traction_djump[_qp](i, i) = K_init;
  }
}

Real
HomogenizedExponentialCZM::computePeakStrengthForContact(const Real delta_n,
                                                         const Real delta_s,
                                                         const Real delta_t,
                                                         const Real sigma_n,
                                                         const Real tau_s,
                                                         const Real tau_t)
{
  Real delta_total = std::sqrt(delta_n*delta_n + delta_s*delta_s + delta_t*delta_t);
  
  if (delta_total < 1e-12)
    return sigma_n;
  
  // Mode ratios
  Real ratio_n = std::abs(delta_n) / delta_total;
  Real ratio_s = std::abs(delta_s) / delta_total;
  Real ratio_t = std::abs(delta_t) / delta_total;
  
  // Mixed-mode strength
  Real T_peak_sq = std::pow(sigma_n * ratio_n, 2.0) +
                   std::pow(tau_s * ratio_s, 2.0) +
                   std::pow(tau_t * ratio_t, 2.0);
  
  return std::sqrt(T_peak_sq);
}

Real
HomogenizedExponentialCZM::computeDamageForContact(Real delta_eff,
                                                   Real delta_0_contact,
                                                   Real phi)
{
  // Always use Wang 2025 Eq. 25
  // For isotropic case (δ₀ⁿ = δ₀ᵗ), this automatically gives constant threshold
  Real delta_0_mixed = computeMixedModeDelta0(phi, _delta_0_normal, _delta_0_tangent);
  
  // Check initiation using Wang Eq. 25 threshold (original delta_0 for isotropic)
  if (delta_eff < delta_0_mixed)
    return 0.0;
  
  if (delta_eff >= _delta_c)
    return 1.0;
  
  // Damage evolution uses mode-dependent scaling
  Real r = delta_eff / delta_0_mixed;

  // USE THE SAME FORMULA AS computeDamage()!
  if (_mu == 1.0 && _eta == 1.0)
  {
    return 1.0 - std::exp(1.0 - r);
  }
  else
  {
    Real exponent = std::pow(r, _eta);

  // DEBUG - print only when damage actually starts
  Real  damage_value = 1.0 - std::exp(1.0 - exponent) * exponent / r;
  static bool first_real_damage = false;
  if (!first_real_damage && damage_value > 1e-6) {
    std::cout << "FIRST DAMAGE: delta_eff=" << delta_eff*1e6 << " nm, "
              << "delta_0_mixed/threshold=" << delta_0_mixed*1e6 << " nm, "
              << "r=" << r << ", damage=" << damage_value << ", phi=" << phi << std::endl;
    first_real_damage = true;
  }

    return 1.0 - std::exp(1.0 - exponent) * exponent / r;
  }
}