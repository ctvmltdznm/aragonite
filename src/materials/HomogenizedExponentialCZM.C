// HomogenizedExponentialCZM.C
// Probabilistic CZM homogenisation via 5-point Gauss-Hermite quadrature.
//
// A quality factor f ~ N(1, quality_std_dev^2) scales both T_peak and delta_0
// proportionally, so weaker contacts fail at smaller openings (progressive
// failure). The homogenised traction is E_f[T(delta, f)], computed exactly for
// smooth integrands with the GH rule. The result is independent of the number
// of contacts and of the interface area.
//
// Numerical improvements over earlier versions:
//   1. Smooth Macaulay bracket everywhere — traction and Jacobian are derived
//      from the same smooth function, eliminating the Newton 2-cycle.
//   2. History variable kappa = max(delta_eff) for irreversibility — damage
//      is strictly monotone; no healing under unloading or partial closure.
//   3. Smooth tanh damage activation — removes the kink in dD/ddelta at
//      delta_eff = delta_0 that caused Newton oscillation near initiation.
//   4. Consistent viscous tangent — visc_factor is set to zero whenever the
//      update is clamped (kappa old >= delta_eff), matching the actual response.
//   5. Chain-rule Jacobian for column 0 — d(delta_n_pos)/d(delta_n) multiplies
//      every dT_i/d(delta_n) entry, fully consistent with the smooth bracket.

#include "HomogenizedExponentialCZM.h"
#include "Assembly.h"

#include <algorithm>
#include <cmath>

registerMooseObject("aragoniteApp", HomogenizedExponentialCZM);

// 5-point Gauss-Hermite nodes and weights for N(0,1)
// z_i = sqrt(2)*x_i^GH,  w_i = w_i^GH/sqrt(pi).   Sum(w_i) = 1.
const Real HomogenizedExponentialCZM::GH_NODES[N_QUAD]   =
    { -2.856970, -1.355630, 0.0, 1.355630, 2.856970 };
const Real HomogenizedExponentialCZM::GH_WEIGHTS[N_QUAD] =
    {  0.011257,  0.222076, 0.533333, 0.222076, 0.011257 };

InputParameters
HomogenizedExponentialCZM::validParams()
{
  InputParameters params = CZMComputeLocalTractionTotalBase::validParams();
  params.addClassDescription(
    "Homogenised exponential CZM with probabilistic progressive failure "
    "and smooth numerics (smooth Macaulay bracket, kappa history variable, "
    "tanh damage activation).");

  params.addRequiredParam<Real>("normal_strength",  "Peak normal traction (MPa)");
  params.addRequiredParam<Real>("shear_strength_s", "Peak shear traction s (MPa)");
  params.addRequiredParam<Real>("shear_strength_t", "Peak shear traction t (MPa)");
  // LENGTH UNIT: all openings are in MODEL length units, i.e. whatever the mesh
  // uses. The shipped inputs are micron-scaled (a 2 x 1 x 1 two-grain aragonite
  // domain), so delta_0_normal = 1.2e-4 means 0.12 nm -- an MD-scale separation.
  // The banner below converts with *1e3 (um -> nm) accordingly. These strings
  // used to say "(mm)", which contradicted both the banner and the inputs.
  params.addRequiredParam<Real>("delta_0_normal",
    "Mode-I characteristic opening, in model length units (um in the shipped inputs)");
  params.addRequiredParam<Real>("delta_0_tangent",
    "Mode-II characteristic opening, in model length units (um in the shipped inputs)");
  params.addParam<Real>("delta_c",
    "Complete failure opening, in model length units. Auto-computed if omitted -- but see "
    "the warning from solveForFailureDisplacement: for small eta the solve clamps and "
    "failure_traction_ratio is silently ignored. Prefer setting this explicitly.");
  params.addParam<Real>("failure_traction_ratio", 0.01,
    "Traction/T_peak at which delta_c is defined (default 0.01).");
  params.addParam<Real>("mu",  1.0, "Loading branch exponent (>= 1 recommended)");
  params.addParam<Real>("eta", 1.0, "Softening branch exponent (>= 1 recommended)");
  params.addParam<Real>("damage_viscosity", 0.0,
    "Viscous regularisation time (simulation time units). 0 = rate-independent.");
  params.addParam<Real>("quality_std_dev", 0.1,
    "CV of quality factor f within a QP (GH quadrature spread). 0 = deterministic.");
  params.addParam<Real>("spatial_quality_std_dev", 0.15,
    "Std dev of spatial quality factor across QPs. Each QP draws f_spatial ~ N(1, sigma^2) "
    "once from its element ID, shifting the mean of its local GH distribution. "
    "Prevents synchronised interface failure and spurious damage localisation. "
    "0 = all QPs identical (not recommended for polycrystal simulations).");
  params.addParam<unsigned int>("spatial_random_seed", 1234,
    "Seed for spatial quality factor field. Change to generate different realisations.");
  params.addParam<Real>("macaulay_regularisation",
    "Width eps_mac of the smooth max(delta_n, 0). Defaults to 1e-3*min(delta_0), which keeps "
    "the residual traction at zero gap scale-free; a hard-coded absolute value makes it depend "
    "on which material is being run. The spurious traction at zero gap is "
    "O(eps_mac^2/delta_0^2) once K_contact carries the factor e (see below).");
  params.addParam<Real>("damage_activation_width", 0.02,
    "Width eps_act of the tanh damage activation, as a fraction of delta_0. The logarithmic "
    "gain of the activation envelope in its pre-onset tail is 2/eps_act, i.e. ~100 at the "
    "default: a relative perturbation in delta_eff is amplified ~100x in D at initiation. "
    "Widen it for reproducible regression tests, narrow it for a sharper onset.");
  return params;
}

HomogenizedExponentialCZM::HomogenizedExponentialCZM(const InputParameters & parameters)
  : CZMComputeLocalTractionTotalBase(parameters),
    _normal_strength(getParam<Real>("normal_strength")),
    _shear_strength_s(getParam<Real>("shear_strength_s")),
    _shear_strength_t(getParam<Real>("shear_strength_t")),
    _delta_0_normal(getParam<Real>("delta_0_normal")),
    _delta_0_tangent(getParam<Real>("delta_0_tangent")),
    _delta_c(0.0),
    _mu(getParam<Real>("mu")),
    _eta(getParam<Real>("eta")),
    _failure_traction_ratio(getParam<Real>("failure_traction_ratio")),
    _damage_viscosity(getParam<Real>("damage_viscosity")),
    _quality_std_dev(getParam<Real>("quality_std_dev")),
    _spatial_quality_std_dev(getParam<Real>("spatial_quality_std_dev")),
    _spatial_random_seed(getParam<unsigned int>("spatial_random_seed")),
    _eps_act(getParam<Real>("damage_activation_width")),
    _qp_damage(declareProperty<std::vector<Real>>("qp_damage")),
    _qp_damage_old(getMaterialPropertyOld<std::vector<Real>>("qp_damage")),
    _qp_kappa(declareProperty<std::vector<Real>>("qp_kappa")),
    _qp_kappa_old(getMaterialPropertyOld<std::vector<Real>>("qp_kappa")),
    _damage(declareProperty<Real>("damage")),
    _damage_old(getMaterialPropertyOld<Real>("damage")),
    _delta_eff_prop(declareProperty<Real>("delta_eff")),
    _delta_eff_old(getMaterialPropertyOld<Real>("delta_eff")),
    _diag_phi(declareProperty<Real>("czm_mode_ratio")),
    _diag_delta0(declareProperty<Real>("czm_delta0")),
    _diag_kappa(declareProperty<Real>("czm_kappa")),
    _diag_peak(declareProperty<Real>("czm_peak_traction")),
    _diag_jump_s(declareProperty<Real>("czm_jump_s")),
    _diag_jump_t(declareProperty<Real>("czm_jump_t"))
{
  if (_delta_0_normal  <= 0.0) mooseError("delta_0_normal must be positive");
  if (_delta_0_tangent <= 0.0) mooseError("delta_0_tangent must be positive");
  if (_mu  <= 0.0 || _eta <= 0.0) mooseError("mu and eta must be positive");
  if (_damage_viscosity < 0.0)    mooseError("damage_viscosity must be non-negative");
  if (_quality_std_dev  < 0.0 || _quality_std_dev > 1.0)
    mooseError("quality_std_dev must be in [0, 1]");
  if (_eps_act <= 0.0) mooseError("damage_activation_width must be positive");

  _eps_mac = isParamValid("macaulay_regularisation")
               ? getParam<Real>("macaulay_regularisation")
               : 1.0e-3 * std::min(_delta_0_normal, _delta_0_tangent);
  if (_eps_mac <= 0.0) mooseError("macaulay_regularisation must be positive");
}

void
HomogenizedExponentialCZM::initialSetup()
{
  CZMComputeLocalTractionTotalBase::initialSetup();

  if (isParamValid("delta_c"))
    _delta_c = getParam<Real>("delta_c");
  else
  {
    Real r = solveForFailureDisplacement(_failure_traction_ratio, _eta);
    _delta_c = r * std::max(_delta_0_normal, _delta_0_tangent);
  }

  if (_delta_c <= std::max(_delta_0_normal, _delta_0_tangent))
    mooseError("delta_c must exceed delta_0. Got delta_c = ", _delta_c);

  Moose::out
    << "\n====================================================\n"
    << "HomogenizedExponentialCZM\n"
    << "====================================================\n"
    << "MD parameters:\n"
    << "  normal_strength   = " << _normal_strength   << " MPa\n"
    << "  shear_strength_s  = " << _shear_strength_s  << " MPa\n"
    << "  delta_0_normal    = " << _delta_0_normal  << " model units = "
                                  << _delta_0_normal  * 1e3 << " nm (um-scaled mesh)\n"
    << "  delta_0_tangent   = " << _delta_0_tangent << " model units = "
                                  << _delta_0_tangent * 1e3 << " nm\n"
    << "  delta_c           = " << _delta_c        << " model units = "
                                  << _delta_c         * 1e3 << " nm"
                                  << (_delta_c_clamped ? "   [CLAMPED - see warning]" : "") << "\n"
    << "  delta_c / delta_0 = " << _delta_c / std::max(_delta_0_normal, _delta_0_tangent) << "\n"
    << "  mu=" << _mu << ", eta=" << _eta << "\n"
    << "Homogenisation:\n"
    << "  quality_std_dev         = " << _quality_std_dev * 100 << "% (within-QP GH spread)\n"
    << "  spatial_quality_std_dev = " << _spatial_quality_std_dev * 100 << "% (across-QP)\n"
    << "  spatial_random_seed     = " << _spatial_random_seed << "\n"
    << "  quadrature points       = " << N_QUAD << " (Gauss-Hermite)\n"
    << "  damage_viscosity        = " << _damage_viscosity
    << (_damage_viscosity > 1e-12 ? "" : "   (rate-independent)") << "\n"
    << "Regularisation:\n"
    << "  macaulay_regularisation = " << _eps_mac << " ("
       << 100.0 * _eps_mac / std::min(_delta_0_normal, _delta_0_tangent) << "% of delta_0)\n"
    << "  damage_activation_width = " << _eps_act
       << "  (tanh gain at onset ~ " << 2.0 / _eps_act << "x)\n"
    << "====================================================\n\n";
}

void
HomogenizedExponentialCZM::initQpStatefulProperties()
{
  _qp_damage[_qp].assign(N_QUAD, 0.0);
  _qp_kappa[_qp].assign(N_QUAD, 0.0);
  _damage[_qp]        = 0.0;
  _delta_eff_prop[_qp] = 0.0;
}

// =============================================================================
// CORE: traction and Jacobian
// =============================================================================
void
HomogenizedExponentialCZM::computeInterfaceTractionAndDerivatives()
{
  const Real delta_n = _interface_displacement_jump[_qp](0);
  const Real delta_s = _interface_displacement_jump[_qp](1);
  const Real delta_t = _interface_displacement_jump[_qp](2);

  // ── 1. Smooth Macaulay bracket ────────────────────────────────────────────
  // Replace the hard switch max(0, delta_n) with a smooth approximation.
  // _eps_mac defaults to 1e-3*min(delta_0) so the regularisation is a fixed
  // FRACTION of the physical opening scale rather than an absolute length: a
  // hard-coded 1e-6 was 0.83% of delta_0 for the bone CZM test and 0.52% for
  // coral, so the size of the artefact depended on the material.
  // The derivative factor d_dn enters the Jacobian via the chain rule.
  const Real sqrt_reg  = std::sqrt(delta_n * delta_n + _eps_mac * _eps_mac);
  const Real delta_n_pos = 0.5 * (delta_n + sqrt_reg);       // smooth max(delta_n, 0)
  const Real d_dn        = 0.5 * (1.0 + delta_n / sqrt_reg); // d(delta_n_pos)/d(delta_n)
  // d_dn → 0 deep in compression, → 1 in tension, smooth and differentiable everywhere.

  const Real delta_t_mag = std::sqrt(delta_s * delta_s + delta_t * delta_t);

  const Real delta_eff = std::sqrt(delta_n_pos * delta_n_pos + delta_t_mag * delta_t_mag);
  _delta_eff_prop[_qp] = delta_eff;

  // Mode ratio (phi = 0 pure tension, phi → ∞ pure shear).
  // phi_den regularises against the smooth Macaulay bracket giving a tiny but
  // nonzero delta_n_pos for strongly compressive QPs. No hard clamp needed:
  // computeMixedModeDelta0 returns delta_0_tangent for phi → ∞.
  Real phi = 0.0;
  if (delta_eff > 1e-12)
  {
    const Real phi_den = std::max(delta_n_pos, 1.0e-6 * _delta_0_normal);
    phi = delta_t_mag / phi_den;
  }

  // Mixed-mode delta_0 base (linear in delta_0_n and delta_0_t, so scales with f)
  const Real d0_base = computeMixedModeDelta0(phi, _delta_0_normal, _delta_0_tangent);

  // diagnostics (see header)
  _diag_phi[_qp]    = phi;
  _diag_delta0[_qp] = d0_base;
  _diag_jump_s[_qp] = delta_s;
  _diag_jump_t[_qp] = delta_t;
  _diag_peak[_qp]   = computePeakStrengthForContact(delta_n_pos, delta_s, delta_t,
                        _normal_strength, _shear_strength_s, _shear_strength_t);

  // Jump vector for Jacobian (uses smooth delta_n_pos)
  const Real jump[3] = { delta_n_pos, delta_s, delta_t };

  // ── 2. Gauss-Hermite quadrature ───────────────────────────────────────────
  const Real e_const  = std::exp(1.0);
  const Real eps_act  = _eps_act;  // tanh activation width (fraction of delta_0)

  // ── Spatial quality factor ────────────────────────────────────────────────
  // Each CZM QP draws f_spatial ~ N(1, spatial_quality_std_dev^2) once,
  // deterministically from element ID + QP + seed. This shifts the mean of
  // the local GH distribution for this QP so neighbouring QPs have different
  // strengths and fail at different loading levels — preventing synchronised
  // failure across the interface face (damage localisation into one element).
  //
  // The GH quadrature provides within-QP contact spread around f_spatial.
  // Set spatial_quality_std_dev = 0 to recover the deterministic (all-QPs-
  // identical) behaviour for testing and single-interface validation.
  Real f_spatial = 1.0;
  if (_spatial_quality_std_dev > 0.0)
  {
    const unsigned int qp_seed = _spatial_random_seed
                                 + _current_elem->id() * 100u
                                 + static_cast<unsigned int>(_qp);
    std::mt19937 rng(qp_seed);
    std::normal_distribution<Real> dist(1.0, _spatial_quality_std_dev);
    f_spatial = std::max(0.1, dist(rng));
  }

  Real T_sum          = 0.0;
  Real dT_sum         = 0.0;
  Real D_sum          = 0.0;
  Real dT_und_eff_sum = 0.0;  // sum_g w_g*(1-D_g)*dT_und_g: elastic-only dT, for phi Jacobian
  Real K_sum          = 0.0;  // sum_g w_g*kappa_g, diagnostic only

  for (unsigned int g = 0; g < N_QUAD; g++)
  {
    // Quality factor = spatial mean (f_spatial) + within-QP GH spread.
    // f_spatial varies across QPs (spatial heterogeneity).
    // GH spread provides contact-level variability within each QP.
    Real f_g = f_spatial + _quality_std_dev * GH_NODES[g];
    f_g = std::max(0.1, f_g);

    const Real d0_g  = d0_base * f_g;
    const Real dc_g  = _delta_c * f_g;

    // Elliptic interaction peak strength
    const Real T_peak_g = computePeakStrengthForContact(
        delta_n_pos, delta_s, delta_t,
        _normal_strength  * f_g,
        _shear_strength_s * f_g,
        _shear_strength_t * f_g);

    // ── Undamaged traction and derivative ──
    Real T_und_g  = 0.0;
    Real dT_und_g = 0.0;

    if (delta_eff > 1e-12 && d0_g > 1e-20)
    {
      const Real r_g = delta_eff / d0_g;
      if (r_g < 1.0)
      {
        const Real ex = std::pow(r_g, _mu);
        T_und_g  = T_peak_g * e_const * ex * std::exp(-ex);
        dT_und_g = T_peak_g * e_const / d0_g
                 * _mu * std::pow(r_g, _mu - 1.0) * (1.0 - ex) * std::exp(-ex);
      }
      else
      {
        const Real ex = std::pow(r_g, _eta);
        T_und_g  = T_peak_g * e_const * ex * std::exp(-ex);
        dT_und_g = T_peak_g * e_const / d0_g
                 * _eta * std::pow(r_g, _eta - 1.0) * (1.0 - ex) * std::exp(-ex);
      }
    }

    // ── History variable kappa: max delta_eff seen by this class ──────────
    // Enforces irreversibility: damage can only grow, never heal.
    const Real kappa_old_g = _qp_kappa_old[_qp][g];
    const Real kappa_g     = std::max(kappa_old_g, delta_eff);
    _qp_kappa[_qp][g]      = kappa_g;

    // Damage acts on the secant ray through the peak once this class has passed it.
    // (1 - D) * T_peak * r reproduces e r^eta exp(-r^eta) exactly on monotonic loading,
    // and unloading/reloading becomes linear to the origin (secant).
    Real T_base_g  = T_und_g;
    Real dT_base_g = dT_und_g;
    if (kappa_g >= d0_g && d0_g > 1e-20)
    {
      T_base_g  = T_peak_g * delta_eff / d0_g;
      dT_base_g = T_peak_g / d0_g;
    }

    // dkappa/d(delta_eff): 1 if currently at the maximum, 0 during unloading
    const Real dkappa_deff = (delta_eff >= kappa_old_g - 1e-14) ? 1.0 : 0.0;

    // ── Smooth damage activation and law (in terms of kappa) ──────────────
    // Hard if (kappa >= d0_g) creates a Jacobian kink at initiation.
    // Replace with a tanh envelope: s(r_k) → 0 before damage, 1 after.
    Real D_target_g  = 0.0;
    Real dD_dkappa_g = 0.0;

    if (d0_g > 1e-20 && kappa_g > 1e-20)
    {
      const Real r_k = kappa_g / d0_g;

      if (kappa_g < dc_g)
      {
        // Smooth activation envelope (continuous, no kink at r_k = 1)
        const Real cosh_arg = std::cosh((r_k - 1.0) / eps_act);
        const Real s        = 0.5 * (1.0 + std::tanh((r_k - 1.0) / eps_act));
        const Real ds_dr    = 0.5 / (eps_act * cosh_arg * cosh_arg);

        // Damage law: D = 1 - r^(eta-1) exp(1 - r^eta), continuous at r=1 (D=0)
        const Real r_eta   = std::pow(r_k, _eta);
        const Real expterm = std::exp(1.0 - r_eta);

        Real D_law = 0.0;
        Real dD_dr = 0.0;
        if (_eta == 1.0)
        {
          D_law = 1.0 - expterm;
          dD_dr = expterm;
        }
        else
        {
          D_law = 1.0 - expterm * r_eta / r_k;
          dD_dr = std::pow(r_k, _eta - 2.0) * expterm * (_eta * r_eta - _eta + 1.0);
        }

        // D_target = s * D_law  (smooth activation × damage law)
        D_target_g   = s * D_law;
        dD_dkappa_g  = (ds_dr * D_law + s * dD_dr) / d0_g;
      }
      else
      {
        D_target_g  = 1.0;
        dD_dkappa_g = 0.0;
      }
    }

    // ── Viscous update with correct tangent ──────────────────────────────
    // D_new = D_old + (D_target - D_old) * dt / viscosity
    // Tangent: dD_new/d(delta_eff) = (dt/viscosity) * dD_target/d(delta_eff)
    //          only when the update is not clamped.
    const Real D_old_g = _qp_damage_old[_qp][g];
    Real D_g;
    Real visc_factor;

    if (_damage_viscosity > 1e-12)
    {
      // D_target is already monotone via kappa; clamp to [D_old, 1] for safety
      D_target_g = std::max(D_target_g, D_old_g);
      D_g        = D_old_g + (D_target_g - D_old_g) * _dt / _damage_viscosity;
      D_g        = std::max(0.0, std::min(1.0, D_g));
      // If clipped at 1, derivative is 0; if D_target <= D_old, derivative is 0
      visc_factor = (D_target_g > D_old_g + 1e-14) ? _dt / _damage_viscosity : 0.0;
    }
    else
    {
      // Rate-independent: D follows D_target immediately
      D_target_g = std::max(D_target_g, D_old_g);
      D_g        = D_target_g;
      visc_factor = (D_target_g > D_old_g + 1e-14) ? 1.0 : 0.0;
    }
    _qp_damage[_qp][g] = D_g;

    // ── Past delta_c: force complete failure, bypass viscous delay ─────────
    // Once kappa >= dc_g the contact is fully separated.
    if (kappa_g >= dc_g)
    {
      D_g                = 1.0;
      _qp_damage[_qp][g] = 1.0;
      visc_factor        = 0.0;
    }

    // dD/d(delta_eff) — full chain: D(D_target(kappa(delta_eff)))
    const Real dD_deff_g = visc_factor * dD_dkappa_g * dkappa_deff;

    // ── Accumulate ──
    const Real w_g = GH_WEIGHTS[g];
    T_sum          += w_g * T_base_g * (1.0 - D_g);
    dT_sum         += w_g * ((1.0 - D_g) * dT_base_g - T_base_g * dD_deff_g);
    D_sum          += w_g * D_g;
    K_sum          += w_g * kappa_g;
    dT_und_eff_sum += w_g * (1.0 - D_g) * dT_base_g;
  }

  _damage[_qp] = D_sum;
  _diag_kappa[_qp] = K_sum;

  // ── 3. Traction vector ────────────────────────────────────────────────────
  // T[0] uses smooth delta_n_pos, so it smoothly → 0 in compression.
  _interface_traction[_qp].zero();
  _dinterface_traction_djump[_qp].zero();

  if (delta_eff > 1e-12)
  {
    _interface_traction[_qp](0) = T_sum * delta_n_pos / delta_eff;
    _interface_traction[_qp](1) = T_sum * delta_s     / delta_eff;
    _interface_traction[_qp](2) = T_sum * delta_t     / delta_eff;

    // ── 4. Symmetric Jacobian (delta_eff term) ────────────────────────────
    const Real coef1 = T_sum / delta_eff;
    const Real coef2 = (dT_sum - coef1) / (delta_eff * delta_eff);

    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j)
        _dinterface_traction_djump[_qp](i, j) = coef1*(i==j) + coef2*jump[i]*jump[j];

    // ── 4b. Phi correction — rank-1 update ───────────────────────────────
    // T = T(delta_eff, phi) because phi enters delta_0_mixed and T_peak.
    // The symmetric formula above omits the dT/dphi * dphi/ddelta_j terms.
    // Missing term: J(i,j) += (dT_sum_dphi/delta_eff) * jump[i] * B[j]
    // where B[j] = dphi/d(delta_j_pos) (column 0 gets d_dn applied below).
    //
    // dT_sum/dphi = -(delta_eff/d0_base) * dT_und_eff_sum * dd0_dphi
    // This follows from the chain rule through r = delta_eff/d0_g = delta_eff/(f_g*d0_base).
    // dd0_dphi computed by central finite difference (computeMixedModeDelta0 is cheap).
    {
      Real dT_sum_dphi = 0.0;
      Real B[3] = {0.0, 0.0, 0.0};

      if (d0_base > 1e-20 && delta_t_mag > 1e-20 && delta_n_pos > 1e-20)
      {
        const Real dphi_fd  = 1.0e-5 * phi + 1.0e-4;
        const Real d0_plus  = computeMixedModeDelta0(phi + dphi_fd, _delta_0_normal, _delta_0_tangent);
        const Real d0_minus = computeMixedModeDelta0(phi - dphi_fd, _delta_0_normal, _delta_0_tangent);
        const Real dd0_dphi = (d0_plus - d0_minus) / (2.0 * dphi_fd);

        dT_sum_dphi = -(delta_eff / d0_base) * dT_und_eff_sum * dd0_dphi;

        // B[0] w.r.t. delta_n_pos (chain rule to delta_n applied by d_dn below)
        B[0] = -delta_t_mag / (delta_n_pos * delta_n_pos);
        B[1] =  delta_s     / (delta_t_mag * delta_n_pos);
        B[2] =  delta_t     / (delta_t_mag * delta_n_pos);
      }

      if (dT_sum_dphi != 0.0)
      {
        const Real phi_factor = dT_sum_dphi / delta_eff;
        for (unsigned int i = 0; i < 3; ++i)
          for (unsigned int j = 0; j < 3; ++j)
            _dinterface_traction_djump[_qp](i, j) += phi_factor * jump[i] * B[j];
      }
    }

    // ── 5. Chain rule for column 0 ────────────────────────────────────────
    // Applies to BOTH the delta_eff term and the phi correction above.
    // d_dn converts dphi/d(delta_n_pos) [stored in B[0]] to dphi/d(delta_n).
    for (unsigned int i = 0; i < 3; ++i)
      _dinterface_traction_djump[_qp](i, 0) *= d_dn;
  }
  else
  {
    // Near-zero opening: mode-specific initial elastic stiffness.
    // Same factor e as K_contact above -- the initial slope of
    // T_peak*e*r*exp(-r) is e*T_peak/delta_0. Without it the near-zero
    // Jacobian is 2.71828x too soft and disagrees with the r>0 branch.
    // Column 0 scaled by d_dn (= 0.5 at delta_n = 0, continuous).
    _dinterface_traction_djump[_qp](0, 0) = e_const * _normal_strength  / _delta_0_normal  * d_dn;
    _dinterface_traction_djump[_qp](1, 1) = e_const * _shear_strength_s / _delta_0_tangent;
    _dinterface_traction_djump[_qp](2, 2) = e_const * _shear_strength_t / _delta_0_tangent;
  }

  // ── 6. Compressive penalty contact ───────────────────────────────────────
  // Without this, grain surfaces in compression exert zero restoring force
  // from the CZM, leaving bulk elements — especially those at triple junctions
  // — as the sole resistance to grain interpenetration. Under multi-grain
  // loading these elements become pathologically over-stressed.
  //
  // A smooth compressive penalty:
  //   T_penalty[0] = -K_contact × delta_n_neg
  // where delta_n_neg = smooth_max(0, -delta_n) = 0.5*(sqrt_reg - delta_n).
  //
  // The Jacobian contribution is the direct derivative (NOT going through
  // the d_dn chain rule, which was already applied above):
  //   dT_penalty[0]/d(delta_n) = K_contact × (1 - d_dn)
  // → full stiffness in compression (d_dn→0), zero in tension (d_dn→1), smooth.
  //
  // K_contact must equal the INITIAL COHESIVE STIFFNESS so the Jacobian is
  // continuous at delta_n = 0. For T = T_peak*e*r*exp(-r) that slope is
  //     dT/d(delta) |_{r->0} = e * T_peak / delta_0
  // NOT T_peak/delta_0. The factor e was missing, leaving the compressive
  // branch 2.71828x too soft, so the two eps_mac-driven terms at delta_n = 0
  // did not cancel and the interface carried a spurious TENSILE prestress
  //     (e-1) * T_peak * eps_mac / (2 * delta_0)
  // = 5.0 MPa at T_peak = 700, delta_0 = 1.2e-4, eps_mac = 1e-6. Observable
  // as a NEGATIVE normal jump early in a pure-tension test: mode_i.i showed
  // avg_separation < 0 until the bulk stress passed 4.98 MPa.
  {
    const Real delta_n_neg  = 0.5 * (sqrt_reg - delta_n);  // smooth max(0, -delta_n)
    const Real K_contact    = e_const * _normal_strength / _delta_0_normal;
    _interface_traction[_qp](0)           -= K_contact * delta_n_neg;
    _dinterface_traction_djump[_qp](0, 0) += K_contact * (1.0 - d_dn);
  }
}

// =============================================================================
// Helpers
// =============================================================================

Real
HomogenizedExponentialCZM::computeMixedModeDelta0(Real phi, Real dn, Real dt) const
{
  // Wang 2025 Eq. 25
  if (phi < 1e-12) return dn;
  if (phi > 1e6)   return dt;
  const Real phi2 = phi * phi;
  return dn * dt * std::sqrt(1.0 + phi2) / std::sqrt(dt * dt + phi2 * dn * dn);
}

Real
HomogenizedExponentialCZM::computePeakStrengthForContact(
    Real delta_n, Real delta_s, Real delta_t,
    Real sigma_n, Real tau_s, Real tau_t) const
{
  // Elliptic interaction criterion (not Mohr-Coulomb; comment corrected)
  const Real dtot = std::sqrt(delta_n*delta_n + delta_s*delta_s + delta_t*delta_t);
  if (dtot < 1e-12) return sigma_n;
  const Real rn = std::abs(delta_n) / dtot;
  const Real rs = std::abs(delta_s) / dtot;
  const Real rt = std::abs(delta_t) / dtot;
  return std::sqrt(std::pow(sigma_n*rn, 2.0) + std::pow(tau_s*rs, 2.0)
                   + std::pow(tau_t*rt, 2.0));
}

Real
HomogenizedExponentialCZM::solveForFailureDisplacement(Real ratio, Real exponent)
{
  // Solve  e * r^exponent * exp(-r^exponent) = ratio  for r > 1.
  // Substituting x = r^exponent gives e*x*exp(-x) = ratio, whose large root is
  // ~ -log(ratio/e) + log(-log(ratio/e)); then r = x^(1/exponent). For small
  // exponent that exponentiation is brutal: ratio = 0.01, exponent = 0.15 needs
  // r ~ 7.7e5, far above the 5000 clamp below. The clamp used to be applied
  // SILENTLY, so delta_c came out as 5000*delta_0, failure_traction_ratio was
  // ignored, and the "kappa >= delta_c -> D = 1" complete-failure branch never
  // fired. That was the missing safety net when the traction ran to 3727 MPa
  // against a model ceiling of 900 MPa. Now it is reported.
  Real r = 5.0;
  const Real e = std::exp(1.0);
  bool converged = false;
  for (int iter = 0; iter < 50; ++iter)
  {
    const Real rexp = std::pow(r, exponent);
    const Real f    = e * rexp * std::exp(-rexp) - ratio;
    if (std::abs(f) < 1e-8) { converged = true; break; }
    const Real df = e * exponent * std::pow(r, exponent-1.0) * std::exp(-rexp) * (1.0 - rexp);
    if (std::abs(df) < 1e-12) mooseError("Cannot solve for failure displacement");
    Real rn = r - f / df;
    if (rn < 1.0)  rn = 0.5 * (r + 1.0);
    if (rn > 5000.0) { rn = 5000.0; _delta_c_clamped = true; }
    r = rn;
  }
  if (r <= 1.0) mooseError("Failure displacement solution failed");

  if (_delta_c_clamped || !converged)
  {
    const Real x = -std::log(ratio / e);
    mooseWarning(
      "HomogenizedExponentialCZM: could not honour failure_traction_ratio = ", ratio,
      " at eta = ", exponent, ".\n"
      "  The solve for delta_c hit the r = 5000 clamp, so delta_c = 5000*delta_0 and the\n"
      "  requested residual traction is NOT where the interface fails. The\n"
      "  'kappa >= delta_c -> D = 1' branch will not fire for any realistic opening;\n"
      "  complete failure then relies entirely on D -> 1 asymptotically.\n"
      "  Needed r is about ", std::pow(x + std::log(x), 1.0 / exponent), ".\n"
      "  FIX: set 'delta_c' explicitly (it is a physical parameter), or raise\n"
      "  'failure_traction_ratio' until this warning disappears.");
  }
  return r;
}
