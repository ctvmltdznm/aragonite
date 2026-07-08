# EXPLICIT orthotropy + orthotropic plasticity (equivalent to fabric_plasticity)
# Single element, uniaxial z-tension ramped past yield (~9% strain).
# Exercises the ELASTIC SLOPE and the PLASTIC YIELD. Compare fabric_plasticity
# vs explicit_plasticity: identical output => fabric scaling == explicit orthotropy.
[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]
[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3  nx = 1  ny = 1  nz = 1
    xmin = 0 xmax = 1  ymin = 0 ymax = 1  zmin = 0 zmax = 1
  []
[]
[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [all]
        strain = SMALL  incremental = true  add_variables = true
        generate_output = 'stress_zz strain_zz stress_xx strain_xx vonmises_stress'
      []
    []
  []
[]
[AuxVariables]
  [euler_phi1] order = CONSTANT family = MONOMIAL []
  [euler_Phi]  order = CONSTANT family = MONOMIAL []
  [euler_phi2] order = CONSTANT family = MONOMIAL []
  [plastic_strain] order = CONSTANT family = MONOMIAL []
[]
[ICs]
  [phi1] type = ConstantIC variable = euler_phi1 value = 0 []
  [Phi]  type = ConstantIC variable = euler_Phi  value = 0 []
  [phi2] type = ConstantIC variable = euler_phi2 value = 0 []
[]
[AuxKernels]
  [plastic_strain]
    type = MaterialRealAux variable = plastic_strain
    property = effective_plastic_strain execute_on = TIMESTEP_END
  []
[]
[Materials]
  [elasticity]
    type = ComputeElasticityTensorCoupled
    fill_method = symmetric9
    # = fabric(E0=100000, nu0=0.3, rho=1, m=(1.15,1.0,0.9), k=l=1) symmetric9
    # so this run is the EXPLICIT equivalent of fabric_plasticity.
    C_ijkl = '178028.85 66346.15 59711.54 134615.38 51923.08 109038.46 34615.38 39807.69 44230.77'
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
    sigma_xx_tension = 4980  sigma_yy_tension = 4100  sigma_zz_tension = 5340
    tau_xy_max = 4510  tau_xz_max = 5080  tau_yz_max = 5060
    postyield_mode = 'exp_softening'  residual_strength = 0.5
    kmax = 0.001  kslope = 30.0
    viscosity_mode = 'linear'  eta = 0.002  m = 0.001
    use_damage = false  absolute_tolerance = 1.0e-5  max_iterations_newton = 100
    use_primal_cpp = true  max_iterations_primal = 1000
    line_search_beta = 1e-4  line_search_eta = 0.25
    euler_angle_1 = euler_phi1  euler_angle_2 = euler_Phi  euler_angle_3 = euler_phi2
  []
[]
[BCs]
  [fix_x] type = DirichletBC variable = disp_x boundary = left   value = 0 []
  [fix_y] type = DirichletBC variable = disp_y boundary = bottom value = 0 []
  [fix_z] type = DirichletBC variable = disp_z boundary = back   value = 0 []
  [pull_z]
    type = FunctionDirichletBC variable = disp_z boundary = front
    function = '0.001 * t'   # tension to 10% strain at t=100 (yield ~6.6%)
  []
[]
[Postprocessors]
  [stress_zz]  type = ElementAverageValue variable = stress_zz []
  [strain_zz]  type = ElementAverageValue variable = strain_zz []
  [stress_xx]  type = ElementAverageValue variable = stress_xx []
  [strain_xx]  type = ElementAverageValue variable = strain_xx []
  [vonmises]   type = ElementAverageValue variable = vonmises_stress []
  [plastic_strain] type = ElementAverageValue variable = plastic_strain []
[]
[Executioner]
  type = Transient  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-6  nl_abs_tol = 1e-8
  nl_max_its = 50  l_max_its = 100
  dt = 0.1  end_time = 120.0
  line_search = 'l2'
[]
[Outputs]
  csv = true
  print_linear_residuals = false
  perf_graph = false
[]
