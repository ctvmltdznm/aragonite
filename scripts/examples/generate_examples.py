#!/usr/bin/env python3
"""
generate_examples.py -- build the examples/ folder for the material-model flags.

Produces 54 bone inputs (6 flags x 9 load cases) plus one coral two-element
cohesive-interface input, each self-contained and runnable with no arguments:

    aragonite-opt -i examples/bone/trabecular_iso/tension_xx.i

Regenerate after any change to the flags or presets:

    python3 generate_examples.py --out examples

Load cases, on a single HEX8 element:
  tension_xx/yy/zz, compression_xx/yy/zz
      Symmetry BCs on the three minus faces and a driven plus face; the lateral
      plus faces stay traction free, so the state is UNIAXIAL STRESS and the
      peak maps directly onto the corresponding directional strength.
  shear_xy/xz/yz
      Affine (Taylor) BCs on every face: u_i = gamma * x_j * t. In the material
      frame both the orthotropic stiffness and the quadric yield surface are
      block diagonal, so an affine simple shear produces a PURE shear stress and
      the peak maps onto the corresponding tau. (At finite gamma there is a
      second-order normal stress; at gamma = 2% it is ~1e-4 of the shear stress.)

Why yields land BELOW the nominal strength: every bone preset carries the UMAT's
initial_yield_ratio = RDY = 0.7, so the input strengths are ULTIMATE strengths
and first yield occurs at 0.7 x strength. See MATERIAL_MODEL_THEORY.md section 4.
"""

import argparse
import os
import textwrap

# --------------------------------------------------------------------------
# Models. density_rho is chosen per family: trabecular bone sits at a typical
# BV/TV, compact (cortical) bone near solid. Running a compact preset at
# BV/TV = 0.2 is not a meaningful material.
# --------------------------------------------------------------------------
FABRIC = dict(fabric_m1=0.85, fabric_m2=0.95, fabric_m3=1.20)  # sums to 3

MODELS = {
    "trabecular_iso": dict(
        rho=0.20, fabric=False, ti=False,
        blurb="UMAT flag 0. Isotropic trabecular bone (ISORGWTB, PHZ 2013).\n"
              "All three axes identical; the shear strength is NOT an input, it\n"
              "follows from the isotropy constraint tau = 1/(S0*sqrt(2(1+zeta0)))."),
    "trabecular_ti": dict(
        rho=0.20, fabric=False, ti=True,
        blurb="UMAT flag 1. Transversely isotropic trabecular bone (TIRGWTB).\n"
              "main_direction = 3, so zz is the axial direction and xx/yy are\n"
              "transverse. Exercises the zeta_a -> zeta_ij reference-index\n"
              "conversion (see MATERIAL_MODEL_THEORY.md section 2.3)."),
    "trabecular_fabric_ti": dict(
        rho=0.20, fabric=True, ti=True,
        blurb="UMAT flag 2. Fabric-based transversely isotropic trabecular bone.\n"
              "The two transverse fabric eigenvalues are averaged internally, and\n"
              "the transverse-plane shear modulus and strength use the isotropic\n"
              "forms so that plane stays isotropic."),
    "trabecular_fabric_ortho": dict(
        rho=0.20, fabric=True, ti=False,
        blurb="UMAT flag 3. Fabric-based orthotropic trabecular bone (FABRGWTB).\n"
              "Full orthotropy generated from one base strength plus the fabric\n"
              "eigenvalues: sigma_i = sigma_0 * rho^p * m_i^(2q)."),
    "compact_iso": dict(
        rho=0.90, fabric=False, ti=False,
        blurb="UMAT flag 4. Isotropic compact (cortical) bone (ISORBCT, JJS 2013).\n"
              "Post-yield is exp_hardening (UMAT PYFL=1), not softening."),
    "compact_ti": dict(
        rho=0.90, fabric=False, ti=True,
        blurb="UMAT flag 5. Transversely isotropic compact bone (TIRBCT).\n"
              "main_direction = 3. This is the surface the old Eq. 55/56 convexity\n"
              "guard wrongly rejected for main_direction 2 and 3."),
}

STRAIN = 0.02      # nominal strain / engineering shear at t = 1
AXIS = {"xx": ("disp_x", "right"), "yy": ("disp_y", "top"), "zz": ("disp_z", "front")}
# shear plane -> (driven component, coordinate it varies with)
SHEAR = {"xy": ("x", "y"), "xz": ("x", "z"), "yz": ("y", "z")}

ALL_BOUNDARIES = "'left right bottom top back front'"

HEADER = """\
# ============================================================================
# {title}
#
{blurb}
#
# LOAD CASE: {case_doc}
#
# Single HEX8 element, 1 x 1 x 1 mm. Run as is, no arguments:
#     aragonite-opt -i {fname}
#
# The elastic and plastic responses are BOTH driven by one flag,
# material_model, set once in [GlobalParams]. Setting elastic_model or
# plastic_model in a single block would override it there and warn.
# ============================================================================
"""


def materials_block(name, cfg):
    """[GlobalParams] fragment and the [Materials] block for one model."""
    shared = [f"  material_model = {name}",
              f"  density_rho    = {cfg['rho']}"]
    if cfg["fabric"]:
        shared += [f"  fabric_m1      = {FABRIC['fabric_m1']}",
                   f"  fabric_m2      = {FABRIC['fabric_m2']}",
                   f"  fabric_m3      = {FABRIC['fabric_m3']}"]
    if cfg["ti"]:
        shared += ["  main_direction = 3"]
    return "\n".join(shared)


def bcs_normal(case):
    """Uniaxial stress: symmetry on the minus faces, one driven plus face."""
    key = case.split("_")[1]
    var, face = AXIS[key]
    sign = "-" if case.startswith("compression") else ""
    return textwrap.dedent(f"""\
        [BCs]
          # Symmetry on the three minus faces: removes rigid-body motion and lets
          # the material contract freely on the opposite faces.
          [sym_x]
            type = DirichletBC
            variable = disp_x
            boundary = left
            value = 0
          []
          [sym_y]
            type = DirichletBC
            variable = disp_y
            boundary = bottom
            value = 0
          []
          [sym_z]
            type = DirichletBC
            variable = disp_z
            boundary = back
            value = 0
          []
          # Driven face. The other two plus faces are traction free, so the state
          # is uniaxial stress and the peak is the directional strength x r(0).
          [drive]
            type = FunctionDirichletBC
            variable = {var}
            boundary = {face}
            function = '{sign}{STRAIN}*t'
          []
        []""")


def bcs_shear(case):
    """Affine (Taylor) BCs: every displacement component prescribed on every face.

    This fully determines the deformation gradient, F = I + gamma*e_i(x)e_j and
    nothing else, so the field is HOMOGENEOUS and every quadrature point sees the
    same state. That is what makes the element-average CSV columns checkable at
    all: phi(sigma) == r(kappa) holds pointwise, and with a uniform field it also
    holds for the averages, to 1e-5.

    The state is NOT pure shear, and cannot be. The quadric surface carries a
    linear term, so it is pressure sensitive, and the associated flow for a pure
    shear stress has normal components (N_11, N_22, N_33 run from 0.08 to 0.44 of
    the flow direction). Affine BCs forbid that plastic dilatation, so a normal
    reaction stress develops -- up to about -107 MPa for compact_ti -- and
    sigma_xy ends up above tau * r. analyse_examples.py reports the resulting
    state purity and treats the peak comparison as informational because of it.

    DO NOT try to free the lateral directions to allow the dilatation. It was
    tried: prescribing only the driven component affinely and holding the other
    two on one face each. Removing those constraints does not just free the
    normal strains, it also frees d(u_y)/d(x) and the whole u_z field, so the
    element escapes into other shear modes -- a shear_xy run came back with an
    xz stress 18% of the driven component -- and the field stops being
    homogeneous. Once it is non-uniform, phi(average sigma) != r(average kappa)
    and the consistency check collapses from 1e-5 to 0.4, not because the model
    is wrong but because the averages no longer describe any single material
    point. Homogeneity is worth more here than mode purity.

    Getting both would need the normal faces traction-free AND d(u_y)/d(x) tied,
    i.e. periodic or multi-point constraints. That is the right answer for a real
    RVE and the wrong amount of machinery for a single-element teaching example.
    """
    key = case.split("_")[1]
    comp, coord = SHEAR[key]
    others = [c for c in "xyz" if c != comp]

    fn = textwrap.dedent(f"""\
        [Functions]
          [shear_fn]
            type = ParsedFunction
            expression = '{STRAIN} * {coord} * t'
          []
        []""")

    out = [textwrap.dedent(f"""\
        [BCs]
          # Affine (Taylor) boundary conditions on EVERY face:
          #     u_{comp} = gamma * {coord} * t,   the other two components zero.
          # Engineering shear gamma = {STRAIN} at t = 1. This prescribes the whole
          # deformation gradient, so the field is homogeneous and every quadrature
          # point sees the same state -- which is what makes the element-average
          # CSV columns mean anything.
          #
          # The stress is NOT pure shear: the yield surface is pressure sensitive,
          # so shear flow is dilatant, these BCs forbid the dilatation, and a
          # normal reaction stress appears. The peak therefore sits above
          # tau * r(kappa) and analyse_examples.py reports it as informational
          # with the state purity alongside. Freeing the lateral directions to
          # fix that breaks homogeneity and is much worse; see the docstring in
          # generate_examples.py.
          [u_{comp}]
            type = FunctionDirichletBC
            variable = disp_{comp}
            boundary = {ALL_BOUNDARIES}
            function = shear_fn
          []""")]
    for o in others:
        out.append(f"  [u_{o}]\n"
                   f"    type = DirichletBC\n"
                   f"    variable = disp_{o}\n"
                   f"    boundary = {ALL_BOUNDARIES}\n"
                   f"    value = 0\n"
                   f"  []")
    out.append("[]")
    return fn + "\n\n" + "\n".join(out)


BODY = """
[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 1
    ny = 1
    nz = 1
    xmax = 1
    ymax = 1
    zmax = 1
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'

  # ---- the material model, and the microstructure it needs -----------------
  # One place. ComputeFabricElasticityTensor and
  # OrthotropicPlasticityStressUpdate both read these, so the elastic and
  # plastic responses cannot drift apart.
{shared}
[]

[Physics/SolidMechanics/QuasiStatic]
  [all]
    strain = FINITE
    incremental = true
    add_variables = true
    # REQUIRED with FINITE strain and solve_type = NEWTON. Without it MOOSE
    # drops the rotation-increment terms and convergence stays linear however
    # exact the material tangent is.
    use_finite_deform_jacobian = true
    generate_output = 'stress_xx stress_yy stress_zz stress_xy stress_xz stress_yz
                       strain_xx strain_yy strain_zz strain_xy strain_xz strain_yz'
  []
[]

[Materials]
  [elasticity]
    type = ComputeFabricElasticityTensor
    # model and microstructure come from [GlobalParams]
  []
  [stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = 'plasticity'
  []
  [plasticity]
    type = OrthotropicPlasticityStressUpdate
    # model and microstructure come from [GlobalParams]
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
    absolute_tolerance = 1e-8
    max_iterations_newton = 50
    # Startup self-checks: finite-difference test of the yield gradient, Hessian
    # and softening derivative at stress states with shear, and a Mandel
    # round-trip of the elasticity tensor. Cheap; turn off for production.
    debug_checks = true
  []
[]

{bcs}

[Postprocessors]
  [stress_xx]
    type = ElementAverageValue
    variable = stress_xx
  []
  [stress_yy]
    type = ElementAverageValue
    variable = stress_yy
  []
  [stress_zz]
    type = ElementAverageValue
    variable = stress_zz
  []
  [stress_xy]
    type = ElementAverageValue
    variable = stress_xy
  []
  [stress_xz]
    type = ElementAverageValue
    variable = stress_xz
  []
  [stress_yz]
    type = ElementAverageValue
    variable = stress_yz
  []
  [strain_{strain_pp}]
    type = ElementAverageValue
    variable = strain_{strain_pp}
  []
  [plastic_strain]
    type = ElementAverageMaterialProperty
    mat_prop = effective_plastic_strain
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  line_search = bt
  dt = 0.01
  end_time = 1.0
  nl_abs_tol = 1e-11
  nl_rel_tol = 1e-8
  nl_max_its = 50
[]

[Outputs]
  csv = true
[]
"""

CASE_DOC = {
    "tension": "uniaxial TENSION along {a}, {s}% nominal strain at t = 1.\n#            Peak stress_{aa} = sigma_{aa}_tension x r(0).",
    "compression": "uniaxial COMPRESSION along {a}, {s}% nominal strain at t = 1.\n#            Peak |stress_{aa}| = sigma_{aa}_compression x r(0).",
    "shear": "simple SHEAR in the {a} plane, engineering gamma = {s}% at t = 1.\n"
             "#            The field is homogeneous, so phi(sigma) = r(kappa) holds in the\n"
             "#            element averages to ~1e-5. The stress is NOT pure shear though:\n"
             "#            the yield surface is pressure sensitive, affine BCs forbid the\n"
             "#            dilatant part of the flow, and the normal reaction pushes the\n"
             "#            peak above tau_{aa} x r. That peak is reported, not asserted.",
}


def write_bone(out_dir):
    n = 0
    for name, cfg in MODELS.items():
        d = os.path.join(out_dir, "bone", name)
        os.makedirs(d, exist_ok=True)
        for kind in ("tension", "compression", "shear"):
            keys = ("xx", "yy", "zz") if kind != "shear" else ("xy", "xz", "yz")
            for k in keys:
                case = f"{kind}_{k}"
                fname = f"{case}.i"
                bcs = bcs_shear(case) if kind == "shear" else bcs_normal(case)
                doc = CASE_DOC[kind].format(a=k if kind == "shear" else k[0],
                                            aa=k, s=STRAIN * 100)
                head = HEADER.format(
                    title=f"{name} -- {case}",
                    blurb="\n".join("# " + ln for ln in cfg["blurb"].split("\n")),
                    case_doc=doc, fname=fname)
                body = BODY.format(shared=materials_block(name, cfg), bcs=bcs,
                                   strain_pp=k)
                with open(os.path.join(d, fname), "w") as f:
                    f.write(head + body)
                n += 1
    return n


# --------------------------------------------------------------------------
# Coral: two elements separated by a cohesive interface
# --------------------------------------------------------------------------
# Mesh scale matters here and is the single most confusing parameter in the
# model. delta_0 is a CHARACTERISTIC OPENING, and its ratio to the element size
# governs the numerics: with delta_0 << h the interface is already in softening
# at the first increment and Newton stalls. The MD values (delta_0 ~ 0.19 nm)
# are three orders below any mesh that resolves a coral grain, so production
# runs regularise delta_0 up to roughly the element size. Here the element edge
# is 1e-4 mm = 0.1 um and delta_0 = 1.91e-4 mm, a ratio near 2.
# Common text for the three coral interface tests. They differ only in the
# header, the BCs and the drive amplitude, so they diff cleanly against each
# other -- which is the point: one interface, three loading modes.
CORAL_HEAD = """\
# ============================================================================
# coral -- two elements separated by a cohesive interface
# LOADING: {mode}
#
# Two HEX8 grains with different crystal orientations, bonded by a
# HomogenizedExponentialCZM interface. This is the smallest complete version of
# the aragonite RVE setup: orthotropic elasticity and quadric plasticity in the
# grains, a cohesive law at the boundary between them.
#
# Run as is, no arguments:
#     aragonite-opt -i {fname}
#
{what}
#
# MODE MIXITY. The interface normal is x, so jump component 0 is the opening and
# components 1 and 2 are the sliding. The model mixes the modes twice over:
#   peak traction, by an elliptic interaction on the direction cosines of the
#     jump vector,  T_peak = sqrt((sigma_n*rn)^2 + (tau_s*rs)^2 + (tau_t*rt)^2)
#   characteristic opening, by the Wang 2025 mixed-mode form, which runs from
#     delta_0_normal at pure opening to delta_0_tangent at pure sliding.
# With sigma_n = 626, tau_s = tau_t = 374, delta_0_normal = 1.91e-4 and
# delta_0_tangent = 2.17e-4 mm, that predicts:
#     mode I      T_peak = 626.0 MPa   delta_0_eff = 1.91e-4 mm
#     mode II     T_peak = 374.0 MPa   delta_0_eff = 2.17e-4 mm
#     45 deg mix  T_peak = 515.6 MPa   delta_0_eff = 2.03e-4 mm
# Those are the numbers to check this run against.
#
# MESH SCALE AND delta_0 -- read this before changing either.
#   delta_0 is a characteristic opening, and what governs the numerics is its
#   ratio to the element size. The MD values are delta_0 ~ 0.19 nm, three
#   orders of magnitude below any mesh that resolves a coral grain; used
#   directly on this mesh the interface is in softening at the first increment
#   and Newton stalls. Production runs therefore REGULARISE delta_0 up to
#   roughly the element size. Here the element edge is 1e-4 mm (0.1 um) and
#   delta_0 is near 2e-4 mm, a ratio near 2. The peak traction, which sets
#   failure initiation, is physical; the fracture energy is not.
#
# HEX8 is required. MOOSE hex-to-tet splitting leaves interface tractions
# frozen at the regularisation floor.
# ============================================================================
"""

CORAL_BODY = """
[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 2
    ny = 1
    nz = 1
    xmax = 2e-4      # mm, so each element is 1e-4 mm = 0.1 um
    ymax = 1e-4
    zmax = 1e-4
  []
  [split]
    type = SubdomainBoundingBoxGenerator
    input = gen
    block_id = 1
    block_name = grain_2
    bottom_left = '1e-4 0 0'
    top_right = '2e-4 1e-4 1e-4'
  []
  [rename]
    type = RenameBlockGenerator
    input = split
    old_block = '0'
    new_block = 'grain_1'
  []
  [interface]
    # Duplicates the nodes on the block boundary and creates the sideset the
    # cohesive model acts on. With split_interface the sideset is named after
    # the two blocks: grain_1_grain_2.
    type = BreakMeshByBlockGenerator
    input = rename
    split_interface = true
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'

  # One flag for both the elastic and the plastic response of the grains.
  # coral needs no density or fabric: all nine elastic constants and all nine
  # strengths are explicit.
  material_model = coral
[]

[Physics/SolidMechanics]
  [QuasiStatic]
    [all]
      strain = FINITE
      incremental = true
      add_variables = true
      use_finite_deform_jacobian = true
      generate_output = 'stress_xx stress_yy stress_zz stress_xy strain_xx'
    []
  []
  # The cohesive block is a SIBLING of QuasiStatic, not nested inside it.
  [CohesiveZone]
    [czm]
      boundary = 'grain_1_grain_2'
      strain = FINITE
      generate_output = 'traction_x traction_y traction_z
                         normal_traction tangent_traction
                         jump_x jump_y jump_z normal_jump tangent_jump'
    []
  []
[]

[AuxVariables]
  # Per-element crystal orientation. Both the elasticity and the plasticity
  # read these, so the two stay in the same material frame.
  [euler_phi1]
    order = CONSTANT
    family = MONOMIAL
  []
  [euler_Phi]
    order = CONSTANT
    family = MONOMIAL
  []
  [euler_phi2]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[ICs]
  [phi1_g1]
    type = ConstantIC
    variable = euler_phi1
    value = 0
    block = grain_1
  []
  [phi1_g2]
    type = ConstantIC
    variable = euler_phi1
    value = 35
    block = grain_2
  []
  [Phi_g1]
    type = ConstantIC
    variable = euler_Phi
    value = 0
    block = grain_1
  []
  [Phi_g2]
    type = ConstantIC
    variable = euler_Phi
    value = 20
    block = grain_2
  []
  [phi2_g1]
    type = ConstantIC
    variable = euler_phi2
    value = 0
    block = grain_1
  []
  [phi2_g2]
    type = ConstantIC
    variable = euler_phi2
    value = 0
    block = grain_2
  []
[]

[Materials]
  [elasticity]
    type = ComputeFabricElasticityTensor
    # C_ijkl defaults to the aragonite constants; override C_ijkl to change them.
    coupled_euler_angle_1 = euler_phi1
    coupled_euler_angle_2 = euler_Phi
    coupled_euler_angle_3 = euler_phi2
  []
  [stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = 'plasticity'
  []
  [plasticity]
    type = OrthotropicPlasticityStressUpdate
    # Strengths default to the aragonite values. zeta12/13/23 have NO default:
    # zero is not an acceptable placeholder and no calibrated coral value
    # exists yet, so they must be set explicitly. These are provisional.
    zeta12 = 0.27
    zeta13 = 0.27
    zeta23 = 0.27
    euler_angle_1 = euler_phi1
    euler_angle_2 = euler_Phi
    euler_angle_3 = euler_phi2
    absolute_tolerance = 1e-8
    max_iterations_newton = 50
    debug_checks = true
  []
  [czm_interface]
    type = HomogenizedExponentialCZM
    boundary = 'grain_1_grain_2'
    # Peak tractions, from MD. These are physical and set failure initiation.
    normal_strength  = 626
    shear_strength_s = 374
    shear_strength_t = 374
    # Characteristic openings, REGULARISED to the element size (see header).
    delta_0_normal   = 1.91e-4
    delta_0_tangent  = 2.17e-4
    # Loading and softening branch exponents, fitted to the MD traction curves.
    mu  = 0.92
    eta = 0.25
    # Two levels of scatter in the quality factor: within a quadrature point
    # (integrated exactly by a 5-point Gauss-Hermite rule, so the response
    # depends on neither the contact count nor the interface area) and across
    # quadrature points, which stops a whole interface failing in one step.
    quality_std_dev         = 0.1
    spatial_quality_std_dev = 0.1
    spatial_random_seed     = 1234
    damage_viscosity        = 0.1
  []
[]

{bcs}

[Postprocessors]
  # The same set in all three modes, so the files diff cleanly and so mode
  # purity is visible: in mode I the tangential columns stay at zero, in
  # mode II the normal ones do, and in the mixed case both are active.
  [normal_traction]
    type = SideAverageValue
    variable = normal_traction
    boundary = 'grain_1_grain_2'
  []
  [tangent_traction]
    type = SideAverageValue
    variable = tangent_traction
    boundary = 'grain_1_grain_2'
  []
  [normal_jump]
    type = SideAverageValue
    variable = normal_jump
    boundary = 'grain_1_grain_2'
  []
  [tangent_jump]
    type = SideAverageValue
    variable = tangent_jump
    boundary = 'grain_1_grain_2'
  []
  [interface_damage]
    type = SideAverageMaterialProperty
    property = damage
    boundary = 'grain_1_grain_2'
  []
  [interface_delta_eff]
    type = SideAverageMaterialProperty
    property = delta_eff
    boundary = 'grain_1_grain_2'
  []
  [stress_xx]
    type = ElementAverageValue
    variable = stress_xx
  []
  [stress_xy]
    type = ElementAverageValue
    variable = stress_xy
  []
  [plastic_strain]
    type = ElementAverageMaterialProperty
    mat_prop = effective_plastic_strain
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  line_search = bt
  dt = 0.005
  end_time = 1.0
  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-8
  nl_max_its = 50
[]

[Outputs]
  csv = true
[]
"""

# BCs. Mode I keeps symmetry conditions, which let the bulk contract freely.
# Modes II and mixed prescribe ALL THREE components on BOTH end faces, because
# any unconstrained normal motion would let the interface open and contaminate
# the mode. The bulk deforms by ~1e-6 mm before the interface reaches its peak,
# three orders below the applied displacement, so the jump is essentially the
# applied displacement in every case.
CORAL_BCS_I = """\
[BCs]
  [sym_x]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
  [sym_y]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  []
  [sym_z]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0
  []
  [open_x]
    # ~3 x delta_0_normal, through the peak and well into softening.
    type = FunctionDirichletBC
    variable = disp_x
    boundary = right
    function = '6e-4*t'
  []
[]"""

CORAL_BCS_II = """\
[BCs]
  # grain_1 fully clamped
  [fix_x]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
  [fix_y]
    type = DirichletBC
    variable = disp_y
    boundary = left
    value = 0
  []
  [fix_z]
    type = DirichletBC
    variable = disp_z
    boundary = left
    value = 0
  []
  # grain_2 slides in y with NO x motion, so the interface cannot open and the
  # loading stays pure mode II.
  [no_open_x]
    type = DirichletBC
    variable = disp_x
    boundary = right
    value = 0
  []
  [slide_y]
    # ~3 x delta_0_tangent, through the peak and well into softening.
    type = FunctionDirichletBC
    variable = disp_y
    boundary = right
    function = '6.5e-4*t'
  []
  [no_z]
    type = DirichletBC
    variable = disp_z
    boundary = right
    value = 0
  []
[]"""

CORAL_BCS_MIX = """\
[BCs]
  # grain_1 fully clamped
  [fix_x]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
  [fix_y]
    type = DirichletBC
    variable = disp_y
    boundary = left
    value = 0
  []
  [fix_z]
    type = DirichletBC
    variable = disp_z
    boundary = left
    value = 0
  []
  # grain_2 moves at 45 degrees in the x-y plane: equal opening and sliding, so
  # the jump direction cosines are rn = rs = 1/sqrt(2) and the elliptic
  # interaction predicts T_peak = 515.6 MPa.
  [open_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = right
    function = '4.3e-4*t'
  []
  [slide_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = right
    function = '4.3e-4*t'
  []
  [no_z]
    type = DirichletBC
    variable = disp_z
    boundary = right
    value = 0
  []
[]"""

CORAL_CASES = {
    "coral_czm_mode_I_opening.i": dict(
        mode="MODE I, pure opening. grain_2 is pulled along x, the interface normal.",
        bcs=CORAL_BCS_I,
        what="""\
# WHAT TO LOOK FOR, in the CSV
#   normal_traction rises to 626 MPa, then softens.
#   tangent_traction stays at zero: this is a pure mode.
#   normal_jump is the opening; interface_damage goes 0 -> 1 monotonically.
#   The grains stay almost entirely elastic. The interface is far more
#   compliant than the bulk, so nearly all the applied displacement becomes
#   opening rather than element stretch."""),
    "coral_czm_mode_II_shear.i": dict(
        mode="MODE II, pure sliding. grain_2 slides along y with x held, so the "
             "interface\n#          shears without opening.",
        bcs=CORAL_BCS_II,
        what="""\
# WHAT TO LOOK FOR, in the CSV
#   tangent_traction rises to 374 MPa (the SHEAR strength, lower than the 626
#     of mode I), then softens.
#   normal_traction and normal_jump stay at zero: x is held on both faces.
#   tangent_jump is the slip; interface_damage goes 0 -> 1 monotonically.
#   Compare the peak against coral_czm_mode_I_opening.i: same interface, same
#   damage law, different strength purely because of the loading direction."""),
    "coral_czm_mixed_mode.i": dict(
        mode="MIXED MODE, 45 degrees. grain_2 moves equally in x (opening) and "
             "y (sliding).",
        bcs=CORAL_BCS_MIX,
        what="""\
# WHAT TO LOOK FOR, in the CSV
#   normal_traction and tangent_traction both rise and both soften; neither
#     column is zero, which is what distinguishes this from the two pure modes.
#   The resultant peak should be 515.6 MPa, between the mode I 626 and the
#     mode II 374, exactly as the elliptic interaction predicts for
#     rn = rs = 1/sqrt(2).
#   normal_jump and tangent_jump stay equal to each other throughout, which is
#     the check that the 45 degree path held."""),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="examples", help="output folder (default: examples)")
    a = ap.parse_args()

    n = write_bone(a.out)
    d = os.path.join(a.out, "coral")
    os.makedirs(d, exist_ok=True)
    for fname, cfg in CORAL_CASES.items():
        head = CORAL_HEAD.format(mode=cfg["mode"], what=cfg["what"], fname=fname)
        with open(os.path.join(d, fname), "w") as f:
            f.write(head + CORAL_BODY.format(bcs=cfg["bcs"]))
    print(f"wrote {n} bone inputs and {len(CORAL_CASES)} coral inputs under {a.out}/")


if __name__ == "__main__":
    main()
