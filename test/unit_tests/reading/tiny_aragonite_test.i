# Polycrystal RVE Test: tiny_aragonite.e
# Reads per-element Euler angles from exodus file
# Tests complete framework: exodus reading + coupled elasticity + plasticity

[Mesh]
  [file]
    type = FileMeshGenerator
    file = 'tiny_aragonite.e'  # Your exodus file
    use_for_exodus_restart = true
  []
  
  # Create boundaries using coordinate expressions
  # ADJUST the coordinate thresholds based on your mesh size!
  # Run: ./check_mesh_bounds.sh tiny_aragonite.e to get correct values
  
  [left_boundary]
    type = ParsedGenerateSideset
    input = file
    combinatorial_geometry = 'x < 0.01'  # Left face (x≈0)
    new_sideset_name = 'left'
  []
  
  [right_boundary]
    type = ParsedGenerateSideset
    input = left_boundary
    combinatorial_geometry = 'x > 9.99'  # Right face - ADJUST THIS!
    new_sideset_name = 'right'
  []
  
  [bottom_boundary]
    type = ParsedGenerateSideset
    input = right_boundary
    combinatorial_geometry = 'y < 0.01'  # Bottom face (y≈0)
    new_sideset_name = 'bottom'
  []
  
  [back_boundary]
    type = ParsedGenerateSideset
    input = bottom_boundary
    combinatorial_geometry = 'z < 0.01'  # Back face (z≈0)
    new_sideset_name = 'back'
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Physics/SolidMechanics/QuasiStatic]
  [all]
    strain = SMALL
    incremental = true
    add_variables = true
    generate_output = 'stress_xx stress_yy stress_zz strain_xx strain_yy strain_zz vonmises_stress'
  []
[]

[AuxVariables]
  # ========== READ FROM EXODUS FILE ==========
  # These read element variables from tiny_aragonite.e
  
  [grain_id]
    order = CONSTANT
    family = MONOMIAL
    initial_from_file_var = grain_id  # Reads "grain_id" from exodus
  []
  
  [needle_id]
    order = CONSTANT
    family = MONOMIAL
    initial_from_file_var = needle_id  # Reads "needle_id" from exodus
  []
  
  [euler_phi1]
    order = CONSTANT
    family = MONOMIAL
    initial_from_file_var = euler_phi1  # Reads "euler_phi1" from exodus
  []
  
  [euler_Phi]
    order = CONSTANT
    family = MONOMIAL
    initial_from_file_var = euler_Phi  # Reads "euler_Phi" from exodus
  []
  
  [euler_phi2]
    order = CONSTANT
    family = MONOMIAL
    initial_from_file_var = euler_phi2  # Reads "euler_phi2" from exodus
  []
  
  # ========== OUTPUT VARIABLES ==========
  [plastic_strain_equiv]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [plastic_strain_equiv]
    type = MaterialRealAux
    property = effective_plastic_strain
    variable = plastic_strain_equiv
    execute_on = 'TIMESTEP_END'
  []
[]

[Materials]
  # ========== ELASTICITY (reads per-element angles) ==========
  [elasticity]
    type = ComputeElasticityTensorCoupled
    fill_method = symmetric9
    # Aragonite elastic constants (MPa)
    C_ijkl = '171800 57500 30200 106700 46900 84200 42100 31100 46600'
    
    # Couple to Euler angles read from exodus
    # Note: "coupled_" prefix for elasticity material
    coupled_euler_angle_1 = euler_phi1
    coupled_euler_angle_2 = euler_Phi
    coupled_euler_angle_3 = euler_phi2
  []
  
  # ========== STRESS ==========
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
[]

[BCs]
  # Fix one corner to prevent rigid body motion
  [fix_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0
  []
  [fix_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0
  []
  [fix_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'back'
    value = 0
  []
  
  # Apply displacement loading in X-direction
  [pull_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 'right'
    function = '0.001*t'  # 0.1% strain per time unit
  []
[]

[Postprocessors]
  # ========== GLOBAL AVERAGES ==========
  [stress_xx_global]
    type = ElementAverageValue
    variable = stress_xx
  []
  [strain_xx_global]
    type = ElementAverageValue
    variable = strain_xx
  []
  [vonmises_global]
    type = ElementAverageValue
    variable = vonmises_stress
  []
  [plastic_strain_global]
    type = ElementAverageValue
    variable = plastic_strain_equiv
  []
  
  # ========== STRESS/STRAIN VARIATION (check if per-element works) ==========
  [stress_xx_max]
    type = ElementExtremeValue
    variable = stress_xx
    value_type = max
  []
  [stress_xx_min]
    type = ElementExtremeValue
    variable = stress_xx
    value_type = min
  []
  [stress_xx_range]
    type = ParsedPostprocessor
    pp_names = 'stress_xx_max stress_xx_min'
    expression = 'stress_xx_max - stress_xx_min'
  []
  
  # ========== EULER ANGLE CHECK (verify exodus reading) ==========
  [euler_phi1_max]
    type = ElementExtremeValue
    variable = euler_phi1
    value_type = max
  []
  [euler_phi1_min]
    type = ElementExtremeValue
    variable = euler_phi1
    value_type = min
  []
  [euler_phi1_avg]
    type = ElementAverageValue
    variable = euler_phi1
  []
  
  # ========== GRAIN COUNT CHECK ==========
  [num_grains]
    type = ElementExtremeValue
    variable = grain_id
    value_type = max
  []
  
  # ========== PLASTIC BEHAVIOR ==========
  [plastic_volume_fraction]
    type = ElementAverageValue
    variable = plastic_strain_equiv
    execute_on = 'TIMESTEP_END'
  []
  [max_plastic_strain]
    type = ElementExtremeValue
    variable = plastic_strain_equiv
    value_type = max
  []
[]

[Executioner]
  type = Transient
  dt = 0.1
  end_time = 1000  # Go through elastic → plastic transition
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 30
  line_search = 'basic'
  
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1
    optimal_iterations = 10
    iteration_window = 2
    growth_factor = 1.2
    cutback_factor = 0.5
  []
[]

[Outputs]
  csv = true
  console = true
  print_linear_residuals = false
  
  [exodus]
    type = Exodus
    file_base = tiny_aragonite_out
    time_step_interval = 10  # Output every 10 timesteps
  []
[]
