# ============================================================================
# RVE VALIDATION TEST - SHORT RUN FOR FRAMEWORK VERIFICATION
# ============================================================================
#
# RVE WITH GRAIN-LEVEL INTERFACES
# Cohesive zones between grains only (coarser scale)
#
# VALIDATION MODE (current settings):
#   end_time = 5.0        Stops before yielding
#   Purpose: Verify mesh reading, Euler angles, CZM setup work correctly
#   Output: Elastic-only response, sufficient for validation
#
# PRODUCTION MODE (modify for research):
#   Change end_time to 100.0 or higher
#   Output: Full plastic deformation, interface failure, realistic polycrystal
# ============================================================================

[Mesh]
  [file]
    type = FileMeshGenerator
    file = 'tiny_aragonite_blocks.e'
    use_for_exodus_restart = true
  []

  # Create external boundaries
  [left_boundary]
    type = ParsedGenerateSideset
    input = file
    combinatorial_geometry = 'x < 0.01'
    new_sideset_name = 'left'
  []
  
  [right_boundary]
    type = ParsedGenerateSideset
    input = left_boundary
    combinatorial_geometry = 'x > 9.99'
    new_sideset_name = 'right'
  []
  
  [bottom_boundary]
    type = ParsedGenerateSideset
    input = right_boundary
    combinatorial_geometry = 'y < 0.01'
    new_sideset_name = 'bottom'
  []
  
  [back_boundary]
    type = ParsedGenerateSideset
    input = bottom_boundary
    combinatorial_geometry = 'z < 0.01'
    new_sideset_name = 'back'
  []
  # CRITICAL STEP: Convert grain_id to blocks
  # This requires that tiny_aragonite.e has subdomain IDs set to grain_id
  # If not, you need to preprocess the mesh (see instructions below)
  
  # Break mesh at grain boundaries
  [break_grains]
    type = BreakMeshByBlockGenerator
    input = back_boundary # file
    split_interface = true
    # This creates boundaries named: block1_block2, block2_block3, etc.
    # where block numbers correspond to grain_id values
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [bulk]
        strain = FINITE
        incremental = true
        add_variables = true
        generate_output = 'stress_xx stress_yy stress_zz strain_xx strain_yy strain_zz vonmises_stress'
      []
    []
    
    # ADD COHESIVE ZONES AT GRAIN BOUNDARIES
    [CohesiveZone]
      [./grain_interfaces]
        # All 6 grain boundaries (found by --mesh-only)
        boundary = 'Block0_Block1 Block0_Block2 Block0_Block3 Block1_Block2 Block1_Block3 Block2_Block3'
        strain = FINITE
        generate_output = 'traction_x traction_y traction_z normal_traction tangent_traction jump_x jump_y jump_z normal_jump tangent_jump'
      [../]
    []
  []
[]

[AuxVariables]
  [grain_id]
    order = CONSTANT
    family = MONOMIAL
    initial_from_file_var = grain_id
  []
  [needle_id]
    order = CONSTANT
    family = MONOMIAL
    initial_from_file_var = needle_id
  []
  [euler_phi1]
    order = CONSTANT
    family = MONOMIAL
    initial_from_file_var = euler_phi1
  []
  [euler_Phi]
    order = CONSTANT
    family = MONOMIAL
    initial_from_file_var = euler_Phi
  []
  [euler_phi2]
    order = CONSTANT
    family = MONOMIAL
    initial_from_file_var = euler_phi2
  []
  [plastic_strain]
    order = CONSTANT
    family = MONOMIAL
  []
  
  # Interface variables (only exist on boundaries)
  [damage]
    order = CONSTANT
    family = MONOMIAL
  []
  [delta_eff]
    order = CONSTANT
    family = MONOMIAL
  []
  [n_contacts]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [plastic_strain]
    type = MaterialRealAux
    variable = plastic_strain
    property = effective_plastic_strain
    execute_on = 'TIMESTEP_END'
  []
  
  # Interface properties - only on CZM boundaries
  [damage_aux]
    type = MaterialRealAux
    variable = damage
    property = damage
    boundary = 'Block0_Block1 Block0_Block2 Block0_Block3 Block1_Block2 Block1_Block3 Block2_Block3'  # Same as CohesiveZone boundary
    check_boundary_restricted = false  # Allow triple junctions
    execute_on = 'TIMESTEP_END'
  []
  
  [delta_eff_aux]
    type = MaterialRealAux
    variable = delta_eff
    property = delta_eff
    boundary = 'Block0_Block1 Block0_Block2 Block0_Block3 Block1_Block2 Block1_Block3 Block2_Block3'
    check_boundary_restricted = false  # Allow triple junctions
    execute_on = 'TIMESTEP_END'
  []
  
  [n_contacts_aux]
    type = MaterialRealAux
    variable = n_contacts
    property = n_contacts
    boundary = 'Block0_Block1 Block0_Block2 Block0_Block3 Block1_Block2 Block1_Block3 Block2_Block3'
    check_boundary_restricted = false  # Allow triple junctions
    execute_on = 'TIMESTEP_END'
  []
[]

[Materials]
  # ========== ELASTICITY (reads per-element angles) ==========
  [elasticity]
    type = ComputeElasticityTensorCoupled
    fill_method = symmetric9
    C_ijkl = '171800 57500 30200 106700 46900 84200 42100 31100 46600'
    coupled_euler_angle_1 = euler_phi1
    coupled_euler_angle_2 = euler_Phi
    coupled_euler_angle_3 = euler_phi2
  []
  
  [stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = 'plasticity'
  []
  
  # ========== PLASTICITY (reads same per-element angles) ==========
  [plasticity]
    type = OrthotropicPlasticityStressUpdate

    # Yield strengths (MPa)
    sigma_xx_tension = 4980
    sigma_yy_tension = 4100
    sigma_zz_tension = 5340
    tau_xy_max = 4510
    tau_xz_max = 5080
    tau_yz_max = 5060

    # Post-yield behavior
    postyield_mode = 'exp_softening'
    residual_strength = 0.7
    kmax = 0.001
    kmin = 0.02
    kslope = 30.0

    # Viscoplasticity (for stability)
    viscosity_mode = 'linear'
    eta = 0.002
    m = 0.001

    # Numerical parameters
    use_damage = false
    absolute_tolerance = 1.0e-5
    max_iterations_newton = 20
    use_primal_cpp = true
    max_iterations_primal = 1000
    line_search_beta = 1e-4
    line_search_eta = 0.25
    
    # Couple to SAME Euler angles (read from exodus)
    # Note: NO "coupled_" prefix for plasticity material
    euler_angle_1 = euler_phi1
    euler_angle_2 = euler_Phi
    euler_angle_3 = euler_phi2
  []
  
  # HOMOGENIZED INTERFACE - Applied to grain boundaries
  [czm_grain_interfaces]
    type = HomogenizedExponentialCZM
    boundary = 'Block0_Block1 Block0_Block2 Block0_Block3 Block1_Block2 Block1_Block3 Block2_Block3'  # Same as CohesiveZone boundary

    normal_strength = 700.0
    shear_strength_s = 400.0
    shear_strength_t = 400.0
    
    delta_0 = 0.1
    
    md_contact_area = 1.0e-4
    max_contacts = 10000
    
    mu = 0.95
    eta = 0.1
    
    damage_viscosity = 2.5
    
    strength_std_dev = 0.0
    initial_damage_max = 0.0
    delta0_std_dev = 0.0
  []
[]

[BCs]
  # Uniaxial tension test
  [fix_left_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0
  []
  
  [load_right_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 'right'
    function = '0.01*t'
  []
  
  [fix_bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0
  []
  
  [fix_back_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'back'
    value = 0
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
  [plastic_strain]
    type = ElementAverageValue
    variable = plastic_strain
  []
  
  # Interface monitoring (average across all 6 boundaries)
  [avg_damage]
    type = SideAverageValue
    variable = damage
    boundary = 'Block0_Block1 Block0_Block2 Block0_Block3 Block1_Block2 Block1_Block3 Block2_Block3'
  []
  [max_damage]
    type = ElementExtremeValue
    variable = damage
  []
  [avg_normal_traction]
    type = SideAverageValue
    variable = normal_traction
    boundary = 'Block0_Block1 Block0_Block2 Block0_Block3 Block1_Block2 Block1_Block3 Block2_Block3'
  []
  [avg_normal_jump]
    type = SideAverageValue
    variable = normal_jump
    boundary = 'Block0_Block1 Block0_Block2 Block0_Block3 Block1_Block2 Block1_Block3 Block2_Block3'
  []
  [avg_n_contacts]
    type = SideAverageValue
    variable = n_contacts
    boundary = 'Block0_Block1 Block0_Block2 Block0_Block3 Block1_Block2 Block1_Block3 Block2_Block3'
  []
  # Per-interface monitoring (Block0_Block1 is largest)
  [damage_01]
    type = SideAverageValue
    variable = damage
    boundary = 'Block0_Block1'
  []
  [traction_01]
    type = SideAverageValue
    variable = normal_traction
    boundary = 'Block0_Block1'
  []
  [damage_06]
    type = SideAverageValue
    variable = damage
    boundary = 'Block0_Block3'
  []
  [traction_06]
    type = SideAverageValue
    variable = normal_traction
    boundary = 'Block0_Block3'
  []
[]

[Executioner]
  type = Transient
  dtmax = 0.5
  dtmin = 1e-10
  dt = 0.1
  end_time = 5
  
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  l_abs_tol = 1e-8
  nl_max_its = 200
  
  line_search = 'l2'
  
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1
    optimal_iterations = 22
    iteration_window = 8
    growth_factor = 1.5
    cutback_factor = 0.5
  []
[]

[Outputs]
  exodus = true
  csv = true
  print_linear_residuals = false
[]

#====================================================#
# 4-GRAIN RVE WITH COHESIVE ZONE INTERFACES         #
#====================================================#
# 
# INTERFACES:
#   Block0_Block1: 1467 sides (largest)
#   Block0_Block2: 230 sides
#   Block0_Block3: 56 sides
#   Block1_Block2: 715 sides
#   Block1_Block3: 189 sides
#   Block2_Block3: 64 sides
#
#====================================================#
