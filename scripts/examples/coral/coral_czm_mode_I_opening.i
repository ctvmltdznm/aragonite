# ============================================================================
# coral -- two elements separated by a cohesive interface
# LOADING: MODE I, pure opening. grain_2 is pulled along x, the interface normal.
#
# Two HEX8 grains with different crystal orientations, bonded by a
# HomogenizedExponentialCZM interface. This is the smallest complete version of
# the aragonite RVE setup: orthotropic elasticity and quadric plasticity in the
# grains, a cohesive law at the boundary between them.
#
# Run as is, no arguments:
#     aragonite-opt -i coral_czm_mode_I_opening.i
#
# WHAT TO LOOK FOR, in the CSV
#   normal_traction rises to 626 MPa, then softens.
#   tangent_traction stays at zero: this is a pure mode.
#   normal_jump is the opening; interface_damage goes 0 -> 1 monotonically.
#   The grains stay almost entirely elastic. The interface is far more
#   compliant than the bulk, so nearly all the applied displacement becomes
#   opening rather than element stretch.
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
[]

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
