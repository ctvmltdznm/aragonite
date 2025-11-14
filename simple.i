[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 2
    ny = 2
    nz = 8
    xmax = 1
    ymax = 1
    zmax = 10
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Physics/SolidMechanics/QuasiStatic]
  [all]
    strain = FINITE
    add_variables = true
    generate_output = 'stress_xx stress_yy stress_zz strain_xx strain_yy strain_zz'
  []
[]

[AuxVariables]
  [plastic_strain_equiv]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [plastic_strain_equiv]
    type = MaterialRealAux
    property = equivalent_plastic_strain
    variable = plastic_strain_equiv
  []
[]

[Materials]
  [elasticity]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '171800 57500 30200 106700 46900 84200 42100 31100 46600'
  []
  
  [plasticity]
    type = OrthotropicPlasticityStressUpdate
    sigma_xx_max = 4980
    sigma_yy_max = 4100
    sigma_zz_max = 5340
    tau_xy_max = 4510
    tau_xz_max = 5080
    tau_yz_max = 5060
    hardening_modulus = -20000
    max_plastic_strain = 0.02
    absolute_tolerance = 5.0
  []
  
  [stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = 'plasticity'
  []
[]

[BCs]
  # Fix ONLY the bottom surface in all directions
  [fix_bottom_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'back'
    value = 0
  []
  [fix_bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'back'
    value = 0
  []
  [fix_bottom_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'back'
    value = 0
  []
  
  # Pull top surface in z only - x and y are FREE (Poisson contraction)
  [pull_top]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = 'front'
    function = '1.5*t'
  []
  
  # Sides (left/right/top/bottom in xy) are COMPLETELY FREE - no BCs!
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  
  dt = 0.005
  end_time = 1.0
  
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-7
  nl_max_its = 50
  
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.005
    optimal_iterations = 10
    growth_factor = 1.1
    cutback_factor = 0.5
  []
[]

[Outputs]
  exodus = true
  csv = true
  print_linear_residuals = false
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
  [strain_zz]
    type = ElementAverageValue
    variable = strain_zz
  []
  [plastic_strain]
    type = ElementAverageValue
    variable = plastic_strain_equiv
  []
[]
