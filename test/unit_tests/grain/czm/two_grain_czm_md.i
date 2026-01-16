# Two-Grain Test with Cohesive Zone Model (CZM)
# Based on two_grain_coupled_final.i + MD-derived interface properties
# 
# WHAT'S NEW:
# 1. BreakMeshByBlockGenerator splits nodes at grain boundary
# 2. CohesiveZoneMaster adds interface physics
# 3. BiLinearMixedModeTraction material with MD parameters
# 4. Interface monitoring postprocessors

[Mesh]
  [generated]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 10
    ny = 5
    nz = 5
    xmin = 0
    xmax = 2.0
    ymin = 0
    ymax = 1.0
    zmin = 0
    zmax = 1.0
  []
  
  # Create two blocks for the two grains
  [split_blocks]
    type = SubdomainBoundingBoxGenerator
    input = generated
    block_id = 2
    block_name = 'grain_2'
    bottom_left = '1.0 0 0'
    top_right = '2.0 1.0 1.0'
  []
  
  [rename]
    type = RenameBlockGenerator
    input = split_blocks
    old_block = '0'
    new_block = 'grain_1'
  []
  
  # This duplicates nodes at the interface between grain_1 and grain_2
  # Creates boundary named "grain_1_grain_2"
  [break_mesh]
    type = BreakMeshByBlockGenerator
    input = rename
    split_interface = true  # Default is true
#    add_interface_boundaries = true  # Default is true
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      # Bulk material behavior (same as before)
      [bulk]
        strain = SMALL
        incremental = true
        add_variables = true
        generate_output = 'stress_xx stress_yy stress_zz strain_xx strain_yy strain_zz vonmises_stress'
      []
    []
    
    # Cohesive Zone at grain boundary
    [CohesiveZone]
      [./czm]
        # Boundary created by BreakMeshByBlockGenerator
        boundary = 'grain_1_grain_2'
        
        # Use same strain formulation as bulk
        strain = SMALL
        
        # Output interface variables
        generate_output = 'traction_x traction_y traction_z normal_traction tangent_traction jump_x jump_y jump_z normal_jump tangent_jump'
      [../]
    []
  []
[]

[AuxVariables]
  [plastic_strain_equiv]
    order = CONSTANT
    family = MONOMIAL
  []
  [grain_id]
    order = CONSTANT
    family = MONOMIAL
  []
  
  # Euler angles - per element
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

[AuxKernels]
  [plastic_strain_equiv]
    type = MaterialRealAux
    property = effective_plastic_strain
    variable = plastic_strain_equiv
    execute_on = 'TIMESTEP_END'
  []
[]

[ICs]
  # Grain IDs
  [grain_1_ic]
    type = ConstantIC
    variable = grain_id
    value = 1
    block = 'grain_1'
  []
  [grain_2_ic]
    type = ConstantIC
    variable = grain_id
    value = 2
    block = 'grain_2'
  []
  
  # Euler angles - Grain 1: c-axis along Z (Phi1=0)
  [euler_phi1_grain1]
    type = ConstantIC
    variable = euler_phi1
    value = 0
    block = 'grain_1'
  []
  [euler_Phi_grain1]
    type = ConstantIC
    variable = euler_Phi
    value = 0
    block = 'grain_1'
  []
  [euler_phi2_grain1]
    type = ConstantIC
    variable = euler_phi2
    value = 0
    block = 'grain_1'
  []
  
  # Euler angles - Grain 2: c-axis along X (Phi1=90)
  [euler_phi1_grain2]
    type = ConstantIC
    variable = euler_phi1
    value = 90
    block = 'grain_2'
  []
  [euler_Phi_grain2]
    type = ConstantIC
    variable = euler_Phi
    value = 0
    block = 'grain_2'
  []
  [euler_phi2_grain2]
    type = ConstantIC
    variable = euler_phi2
    value = 0
    block = 'grain_2'
  []
[]

[Materials]
  #========================================
  # BULK MATERIALS (Grains 1 & 2)
  # Same as before - coupled to Euler angles
  #========================================
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
  
  [plasticity]
    type = OrthotropicPlasticityStressUpdate
    
    # Your calibrated yield stresses
    sigma_xx_tension = 4980
    sigma_yy_tension = 4100
    sigma_zz_tension = 5340
    tau_xy_max = 4510
    tau_xz_max = 5080
    tau_yz_max = 5060
    
    # Softening parameters
    postyield_mode = 'exp_softening'
    residual_strength = 0.7
    kmax = 0.001
    kmin = 0.02
    kslope = 30.0
    
    # Viscoplastic regularization
    viscosity_mode = 'linear'
    eta = 0.002
    m = 0.001
    
    use_damage = false
    absolute_tolerance = 1.0e-5
    max_iterations_newton = 20
    use_primal_cpp = true
    max_iterations_primal = 1000
    line_search_beta = 1e-4
    line_search_eta = 0.25
    
    euler_angle_1 = euler_phi1
    euler_angle_2 = euler_Phi
    euler_angle_3 = euler_phi2
  []
  
  #========================================
  # COHESIVE ZONE MATERIAL
  # MD-derived traction-separation law
  #========================================
  [czm_material]
    type = BiLinearMixedModeTraction
    boundary = 'grain_1_grain_2'
    
    # ========================================
    # PENALTY STIFFNESS
    # ========================================
    # In penalty method, this controls pre-failure stiffness
    # 
    # MD gives VERY high stiffness (~1e14 MPa/mm at atomic scale)
    # But for FEM, we use ~10-100x elastic modulus to avoid:
    #   1. Ill-conditioning of stiffness matrix
    #   2. Numerical instabilities
    # 
    # Rule of thumb: K_penalty ≈ 10-100 × max(E_ii)
    #                          ≈ 10-100 × 171800 MPa
    #                          ≈ 1e6 - 1e7 MPa/mm
    # 
    # Start conservative (lower value = more compliant interface)
    penalty_stiffness = 1.0e6  # MPa/mm
    
    # ========================================
    # MODE I: NORMAL (TENSION)
    # From MD: czm_parameters_z.dat
    # ========================================
    normal_strength = 183.79  # MPa (σc from MD)
    GI_c = 9.850e-05          # N/mm (fracture energy from MD)
    
    # ========================================
    # MODE II/III: SHEAR
    # Averaged from YZ and XZ modes
    # ========================================
    # YZ: τc = 115.05 MPa, G_II = 1.282e-04 N/mm
    # XZ: τc = 103.61 MPa, G_II = 1.423e-04 N/mm
    # Average:
    shear_strength = 109.33   # MPa (τc,avg)
    GII_c = 1.353e-04         # N/mm (G_II,avg)
    
    # ========================================
    # MIXED-MODE PARAMETERS
    # ========================================
    # BK criterion exponent for mixed-mode fracture
    # G_c = G_I + (G_II - G_I) * (G_shear/G_total)^eta
    eta = 2.2  # Standard value for composites
    
    # ========================================
    # VISCOUS REGULARIZATION
    # ========================================
    # Helps convergence during softening
    # Too small - convergence issues
    # Too large - artificial stiffening
    viscosity = 1e-3  # Start here, increase if convergence problems
    
    # Optional: Lag separation state for better convergence
    # lag_separation_state = true
  []
[]

[BCs]
  # Same boundary conditions as before
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
  
  # Pull in X direction (perpendicular to interface)
  [pull_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 'right'
    function = '0.001*t'  # 0.001 mm/s displacement rate
  []
[]

[Postprocessors]
  # ========================================
  # BULK BEHAVIOR (same as before)
  # ========================================
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
  
  # Per-grain averages
  [stress_xx_grain1]
    type = ElementAverageValue
    variable = stress_xx
    block = 'grain_1'
  []
  [strain_xx_grain1]
    type = ElementAverageValue
    variable = strain_xx
    block = 'grain_1'
  []
  [vonmises_grain1]
    type = ElementAverageValue
    variable = vonmises_stress
    block = 'grain_1'
  []
  [plastic_strain_grain1]
    type = ElementAverageValue
    variable = plastic_strain_equiv
    block = 'grain_1'
  []
  [stress_xx_grain2]
    type = ElementAverageValue
    variable = stress_xx
    block = 'grain_2'
  []
  [strain_xx_grain2]
    type = ElementAverageValue
    variable = strain_xx
    block = 'grain_2'
  []
  [vonmises_grain2]
    type = ElementAverageValue
    variable = vonmises_stress
    block = 'grain_2'
  []
  [plastic_strain_grain2]
    type = ElementAverageValue
    variable = plastic_strain_equiv
    block = 'grain_2'
  []
  
  # ========================================
  # INTERFACE BEHAVIOR
  # ========================================
  
  # Maximum tractions at interface
  [max_normal_traction]
    type = ElementExtremeValue  # Changed from NodalExtremeValue
    variable = normal_traction
#    boundary = 'grain_1_grain_2'  # should exist only on the interface
  []
  [max_tangent_traction]
    type = ElementExtremeValue  # Changed from NodalExtremeValue
    variable = tangent_traction
#    boundary = 'grain_1_grain_2'
  []
  
  # Maximum displacement jumps
  [max_normal_jump]
    type = ElementExtremeValue  # Changed from NodalExtremeValue
    variable = normal_jump
#    boundary = 'grain_1_grain_2'
  []
  [max_tangent_jump]
    type = ElementExtremeValue  # Changed from NodalExtremeValue
    variable = tangent_jump
#    boundary = 'grain_1_grain_2'
  []
  
  # Average interface values
  [avg_normal_traction]
    type = SideAverageValue
    variable = normal_traction
    boundary = 'grain_1_grain_2'
  []
  [avg_tangent_traction]
    type = SideAverageValue
    variable = tangent_traction
    boundary = 'grain_1_grain_2'
  []
  [avg_normal_jump]
    type = SideAverageValue
    variable = normal_jump
    boundary = 'grain_1_grain_2'
  []
  [avg_tangent_jump]
    type = SideAverageValue
    variable = tangent_jump
    boundary = 'grain_1_grain_2'
  []
  
  # Interface energy dissipation (approximate)
  # This is cumulative energy ≈ ∫ T·dδ
  [interface_work]
    type = ParsedPostprocessor
    pp_names = 'avg_normal_traction avg_normal_jump avg_tangent_traction avg_tangent_jump'
    expression = 'avg_normal_traction * avg_normal_jump + avg_tangent_traction * avg_tangent_jump'
  []
  
  # ========================================
  # VALIDATION: Elastic moduli
  # ========================================
  [E_grain1]
    type = ParsedPostprocessor
    pp_names = 'stress_xx_grain1 strain_xx_grain1'
    expression = 'if(strain_xx_grain1 > 1e-10, stress_xx_grain1 / strain_xx_grain1, 0)'
  []
  [E_grain2]
    type = ParsedPostprocessor
    pp_names = 'stress_xx_grain2 strain_xx_grain2'
    expression = 'if(strain_xx_grain2 > 1e-10, stress_xx_grain2 / strain_xx_grain2, 0)'
  []
[]

[Executioner]
  type = Transient
  
  # Start with smaller time steps for initial equilibration
  dt = 0.001
  end_time = 150  # Pull 150 mm total
  
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist'
  
  # Convergence tolerances
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
  nl_max_its = 100  # Increase max iterations (CZM can be tricky)
  
  # Line search helps with CZM convergence
  line_search = 'l2'
  
  # Adaptive time stepping
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.001
    optimal_iterations = 10
    iteration_window = 2
    growth_factor = 1.1
    cutback_factor = 0.25
    cutback_factor_at_failure = 0.25  # Aggressive cutback if diverged
  []
[]

[Outputs]
  csv = true
  console = true
  print_linear_residuals = false
  
  [exodus]
    type = Exodus
    file_base = two_grain_czm_out
    time_step_interval = 10  # Save every 10 steps (reduce file size)
  []
[]

# ========================================
# USAGE INSTRUCTIONS:
# ========================================
# 
# 1. Run the simulation:
#    $ ~/projects/aragonite/aragonite-opt -i two_grain_czm_md.i
# 
# 2. Monitor output:
#    - Watch for "max_normal_traction" and "max_tangent_traction"
#    - Should reach ~183.79 MPa (normal) and ~109.33 MPa (shear)
#    - Then soften as interface fails
# 
# 3. Check convergence:
#    - If "PETSC ERROR" or "Solve failed":
#      → Reduce dt (try 0.001)
#      → Increase viscosity (try 1e-2)
#      → Enable lag_separation_state = true in czm_material
#    
#    - If "Negative Jacobian":
#      → Reduce penalty_stiffness (try 1e5)
#      → Check mesh quality
# 
# 4. Visualize in ParaView:
#    - Load two_grain_czm_out.e
#    - Plot "normal_traction" on interface → should see stress concentration
#    - Plot "plastic_strain_equiv" in bulk → see which grain yields first
#    - Animate to see failure progression
# 
# 5. Post-process with Python:
#    $ python analyze_czm_results.py two_grain_czm_out.csv
# 
# ========================================
# WHAT TO EXPECT:
# ========================================
# 
# PHASE 1: Elastic loading (t ~ 0-20)
#   - Both grains deform elastically
#   - Interface traction increases linearly
#   - Small displacement jump (~1e-6 mm)
# 
# PHASE 2: Grain yielding (t ~ 20-50)
#   - Grain 1 (weak orientation) yields at ~4100 MPa
#   - Grain 2 (strong orientation) yields at ~4980 MPa
#   - Plastic strain accumulates
# 
# PHASE 3: Interface failure initiation (t ~ ??)
#   - Depends on relative strength: grains vs interface
#   - If interface is WEAKER than grains:
#     * Interface reaches 183.79 MPa → starts softening
#     * Occurs BEFORE grain yielding
#   - If interface is STRONGER than grains:
#     * Grains yield first
#     * Interface fails later under high stress
# 
# PHASE 4: Interface softening/failure (t ~ ??+)
#   - Traction decreases
#   - Displacement jump increases rapidly
#   - Energy dissipated = GI or GII (fracture energy)
#   - Grains start to separate
# 
# ========================================
# COMPARING TO YOUR TWIN BOUNDARY TEST:
# ========================================
# 
# WITHOUT CZM (your original test):
#   - Grains perfectly bonded
#   - Load capacity ~ 4100 MPa (weak grain yield)
#   - Failure mode: Bulk plasticity
# 
# WITH CZM (this test):
#   - If interface weaker: Load capacity ~ 183.79 MPa
#   - If interface stronger: Load capacity ~ 4100 MPa (then interface fails)
#   - Failure mode: Interface debonding OR bulk plasticity
# 
# This reveals competition between:
#   * Intergranular fracture (interface failure)
#   * Transgranular fracture (bulk yielding)
# 
# ========================================
