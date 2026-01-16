# Unit Test: XY-Shear
# Expected: Yield at 4510 MPa, soften to 3157 MPa

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 1
  ny = 1
  nz = 1
  xmax = 1
  ymax = 1
  zmax = 1
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Physics/SolidMechanics/QuasiStatic]
  [all]
    strain = SMALL
    incremental = true
    add_variables = true
    generate_output = 'stress_xx stress_yy stress_zz stress_xy stress_xz stress_yz strain_xx strain_yy strain_zz strain_xy strain_yz strain_xz'
  []
[]

[AuxVariables]
  [plastic_strain]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [plastic_strain]
    type = MaterialRealAux
    property = effective_plastic_strain
    variable = plastic_strain
    execute_on = 'TIMESTEP_END'
  []
[]

[Materials]
  [elasticity]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    # C11=171800, C12=57500, C13=30200, C22=106700, C23=46900, C33=84200
    # C44=42100 (yz), C55=31100 (xz), C66=46600 (xy)
    C_ijkl = '171800 57500 30200 106700 46900 84200 42100 31100 46600'
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
  []
  
  [stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = 'plasticity'
  []
  
  [plasticity]
    type = OrthotropicPlasticityStressUpdate
    
    sigma_xx_tension = 4980
    sigma_yy_tension = 4100
    sigma_zz_tension = 5340
    tau_xy_max = 4510
    tau_xz_max = 5080
    tau_yz_max = 5060
    
    # ========== NEW ENUM-BASED POST-YIELD SYSTEM ==========
    # Options: perfect, exp_hardening, linear_hardening, 
    #          exp_softening, piecewise_softening
    postyield_mode = 'linear_hardening'
    residual_strength = 0.7
    kmax = 0.001
    kmin = 0.02
    kslope = 30.0
    
    viscosity_mode = 'linear'
    eta = 0.002
    m = 0.001
    
    use_damage = false
    absolute_tolerance = 1.0e-5
    max_iterations_newton = 20
    use_primal_cpp = true
    max_iterations_primal = 1000     # Primal CPPA max iterations
    line_search_beta = 1e-4          # Armero line search parameter
    line_search_eta = 0.25           # Line search backtracking factor
    
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
  []
[]

[BCs]
  # ============================================
  # PURE SHEAR BCs: Fix bottom, shear top
  # ============================================
  
  # Bottom face (y=0): completely fixed
  [fix_bottom_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'bottom'
    value = 0
  []
  [fix_bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0
  []
  [fix_bottom_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'bottom'
    value = 0
  []
  
  # Top face (y=1): apply shear in x, fix y
  [shear_top_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 'top'
    # γ_xy = u_x / h = 0.001*t, so σ_xy = G*γ = 46600*0.001*t = 46.6*t MPa
    # Yield at σ_xy = 4510: t = 4510/46.6 ≈ 96.8 s
    function = '0.001*t'
  []
  [fix_top_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'top'
    value = 0
  []
  # Note: disp_z on top is FREE to allow plane strain condition
  
  # Front/back faces: constrain z for plane strain
  [fix_front_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'front'
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
  [strain_xy]
    type = ElementAverageValue
    variable = strain_xy
  []
  [plastic_strain]
    type = ElementAverageValue
    variable = plastic_strain
  []
[]

[Executioner]
  type = Transient
  dt = 0.1
  end_time = 170.0
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
  file_base = xy_shear_out_linear_linear_hardening
  print_linear_residuals = false
[]
