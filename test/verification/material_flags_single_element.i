# Single-element driver for the material-model flags.
# Uniaxial compression along z (2 % at t = 1), Poisson free.
#
# This input exists because a MOOSE CLI override can create a block but cannot
# set its TYPE. Overriding Materials/elasticity/elastic_model on an input whose
# elasticity block is a ComputeElasticityTensorCoupled (or is named something
# else) makes MOOSE build an empty [elasticity] block and abort with
#   missing required parameter 'Materials/elasticity/type'
# So the flags are exercised here, on a block that is already a
# ComputeFabricElasticityTensor, not by overriding a production input.
#
# DEFAULTS ARE CORAL, because coral is the only model that needs nothing but
# the strengths it already defaults to. Every other model is selected from the
# command line together with the inputs it requires:
#
#   # bone, fabric-based orthotropic
#   aragonite-opt -i material_flags_single_element.i \
#     Materials/elasticity/elastic_model=trabecular_fabric_ortho \
#     Materials/plasticity/plastic_model=trabecular_fabric_ortho \
#     Materials/elasticity/density_rho=0.2 Materials/plasticity/density_rho=0.2 \
#     Materials/elasticity/fabric_m1=0.85 Materials/plasticity/fabric_m1=0.85 \
#     Materials/elasticity/fabric_m2=0.95 Materials/plasticity/fabric_m2=0.95 \
#     Materials/elasticity/fabric_m3=1.20 Materials/plasticity/fabric_m3=1.20
#
#   # bone, transversely isotropic compact, axis 3
#   aragonite-opt -i material_flags_single_element.i \
#     Materials/elasticity/elastic_model=compact_ti \
#     Materials/plasticity/plastic_model=compact_ti \
#     Materials/elasticity/main_direction=3 Materials/plasticity/main_direction=3 \
#     Materials/elasticity/density_rho=0.2 Materials/plasticity/density_rho=0.2
#
# EXPECTED WARNINGS: selecting a non-coral model leaves zeta12/13/23 set below,
# which that model does not use, so it prints three
#   Parameter 'zeta12' is ignored for plastic_model = ...
# lines. That is the flag machinery reporting correctly, not a failure. No
# single static input can avoid this: any parameter written here is ignored by
# some model.
#
# The two flags are independent: elastic_model and plastic_model may differ.
#
# Gold files: --generate-gold, one per combination, with a non-colliding
# Outputs/file_base.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Physics/SolidMechanics/QuasiStatic/all]
  strain = FINITE
  incremental = true
  add_variables = true
  # REQUIRED with strain = FINITE and solve_type = NEWTON: without it MOOSE
  # omits the rotation-increment terms and convergence stays linear however
  # exact the material tangent is (PRODUCTION_SETTINGS.md).
  use_finite_deform_jacobian = true
  generate_output = 'stress_xx stress_yy stress_zz strain_zz'
[]

[Materials]
  [elasticity]
    type = ComputeFabricElasticityTensor
    elastic_model = coral
    # coral uses the built-in C_ijkl default; override C_ijkl to change it.
    # Bone models additionally need density_rho, and the fabric ones
    # fabric_m1/2/3 -- supply those from the command line.
  []
  [stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = 'plasticity'
  []
  [plasticity]
    type = OrthotropicPlasticityStressUpdate
    plastic_model = coral
    # coral has no calibrated zeta yet; these are placeholders, and they are
    # the source of the "ignored" warnings under any non-coral model.
    zeta12 = 0.27
    zeta13 = 0.27
    zeta23 = 0.27
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
    absolute_tolerance = 1e-8
    max_iterations_newton = 50
    debug_checks = true   # FD + Mandel self-checks; turn off for production
  []
[]

[BCs]
  [x]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
  [y]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  []
  [z0]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0
  []
  [z1]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = front
    function = '-0.02*t'
  []
[]

[Postprocessors]
  [s_zz]
    type = ElementAverageValue
    variable = stress_zz
  []
  [e_zz]
    type = ElementAverageValue
    variable = strain_zz
  []
  [kappa]
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
  nl_abs_tol = 1e-8
  nl_max_its = 50
[]

[Outputs]
  csv = true
  # set a non-colliding name per combination: Outputs/file_base=<elastic>_<plastic>_md<k>
[]
