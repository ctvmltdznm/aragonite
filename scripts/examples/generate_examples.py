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
# Gold file:
#     aragonite-opt -i {fname} --generate-gold
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
    """Affine simple shear on every face -> homogeneous pure shear stress."""
    key = case.split("_")[1]
    comp, coord = SHEAR[key]
    others = [c for c in "xyz" if c != comp]
    blocks = [textwrap.dedent(f"""\
        [BCs]
          # Affine (Taylor) boundary conditions on EVERY face:
          #     u_{comp} = gamma * {coord} * t,   the other two components zero.
          # The deformation is homogeneous simple shear with engineering shear
          # gamma = {STRAIN} at t = 1. Both the orthotropic stiffness and the quadric
          # yield surface are block diagonal in the material frame, so the stress
          # stays pure shear and the peak is tau_{key} x r(0).
          [u_{comp}]
            type = FunctionDirichletBC
            variable = disp_{comp}
            boundary = {ALL_BOUNDARIES}
            function = shear_fn
          []""")]
    for o in others:
        blocks.append(f"  [u_{o}]\n"
                      f"    type = DirichletBC\n"
                      f"    variable = disp_{o}\n"
                      f"    boundary = {ALL_BOUNDARIES}\n"
                      f"    value = 0\n"
                      f"  []")
    blocks.append("[]")
    fn = textwrap.dedent(f"""\
        [Functions]
          [shear_fn]
            type = ParsedFunction
            expression = '{STRAIN} * {coord} * t'
          []
        []""")
    return fn + "\n\n" + "\n".join(blocks)


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
    "shear": "simple SHEAR in the {a} plane, engineering gamma = {s}% at t = 1.\n#            Peak stress_{aa} = tau_{aa}_max x r(0).",
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
CORAL = """\
# ============================================================================
# coral -- two elements separated by a cohesive interface
#
# Two HEX8 grains with different crystal orientations, bonded by a
# HomogenizedExponentialCZM interface, pulled apart along x. This is the
# smallest complete version of the aragonite RVE setup: orthotropic elasticity
# and quadric plasticity in the grains, a cohesive law at the boundary between
# them.
#
# Run as is, no arguments:
#     aragonite-opt -i coral_czm_two_element.i
#
# WHAT TO LOOK FOR, in coral_czm_two_element_out.csv
#   normal_traction rises to normal_strength (626 MPa), then softens.
#   normal_jump is the interface opening; damage goes 0 -> 1 monotonically.
#   The grains stay almost entirely elastic: the interface is far more
#   compliant than the bulk, so nearly all the applied displacement goes into
#   the opening, not into stretching the elements.
#
# MESH SCALE AND delta_0 -- read this before changing either.
#   delta_0 is a characteristic opening, and what governs the numerics is its
#   ratio to the element size. The MD values are delta_0 ~ 0.19 nm, three
#   orders of magnitude below any mesh that resolves a coral grain; used
#   directly on this mesh the interface is in softening at the first increment
#   and Newton stalls. Production runs therefore REGULARISE delta_0 up to
#   roughly the element size. Here the element edge is 1e-4 mm (0.1 um) and
#   delta_0 = 1.91e-4 mm, a ratio near 2. The peak traction, which is what
#   sets failure initiation, is physical; the fracture energy is not.
#
# HEX8 is required. MOOSE hex-to-tet splitting leaves interface tractions
# frozen at the regularisation floor.
# ============================================================================

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
      generate_output = 'stress_xx stress_yy stress_zz strain_xx'
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
  [pull_x]
    # ~3 x delta_0, enough to take the interface through its peak and well into
    # softening. Almost all of this is interface opening: the bulk stretches
    # only about 7e-7 mm before the interface reaches 626 MPa.
    type = FunctionDirichletBC
    variable = disp_x
    boundary = right
    function = '6e-4*t'
  []
[]

[Postprocessors]
  [stress_xx]
    type = ElementAverageValue
    variable = stress_xx
  []
  [strain_xx]
    type = ElementAverageValue
    variable = strain_xx
  []
  [normal_traction]
    type = SideAverageValue
    variable = normal_traction
    boundary = 'grain_1_grain_2'
  []
  [normal_jump]
    type = SideAverageValue
    variable = normal_jump
    boundary = 'grain_1_grain_2'
  []
  [interface_damage]
    type = SideAverageMaterialProperty
    property = damage
    boundary = 'grain_1_grain_2'
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="examples", help="output folder (default: examples)")
    a = ap.parse_args()

    n = write_bone(a.out)
    d = os.path.join(a.out, "coral")
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "coral_czm_two_element.i"), "w") as f:
        f.write(CORAL)
    print(f"wrote {n} bone inputs and 1 coral input under {a.out}/")


if __name__ == "__main__":
    main()
