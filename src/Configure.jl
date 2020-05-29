abstract type JlBoxConfig end
struct GasConfig <: JlBoxConfig
    file#"MCM_test.eqn.txt"MCM_APINENE.eqn.txt
    temp# Kelvin
    RH# RH/100% [0 - 0.99]
    start_time# seconds, used as t0 in solver
    simulation_time# seconds
    batch_step# seconds
    #temp_celsius
    #Psat# Saturation VP of water vapour, to get concentration of H20
    #Pw
    #Wconc#kg/cm3
    H2O#Convert from kg to molecules/cc
    Cfactor#ppb-to-molecules/cc
    reactants_initial_dict::Dict#ppm ["O3"=>18.0,"APINENE"=>30.0])
    constantdict::Dict
end

struct AerosolConfig <: JlBoxConfig
    file#"MCM_test.eqn.txt"MCM_APINENE.eqn.txt
    temp# Kelvin
    RH# RH/100% [0 - 0.99]
    start_time# seconds, used as t0 in solver
    simulation_time# seconds
    batch_step# seconds
    #temp_celsius
    #Psat_w# Saturation VP of water vapour, to get concentration of H20
    #Pw
    #Wconc#kg/cm3
    H2O#Convert from kg to molecules/cc
    Cfactor#ppb-to-molecules/cc
    reactants_initial_dict::Dict#ppb BUT1ENE APINENE
    constantdict::Dict
    
    num_bins
    #Lognormal Distribution
    total_conc#Total particles per cc
    size_std#Standard Deviation
    lowersize#microns
    uppersize#microns
    meansize#microns
    # - Specify the core material. 
    # This code is currently setup to consider *ammonium sulphate* as the core
    y_core_init#Will hold concentration of core material, only initialise here [molecules/cc] 
    core_density_array#[kg/m3] - need to make sure this matches core definition above
    core_mw#[g/mol]
    core_dissociation#Define this according to choice of core type. Please note this value might change
    
    vp_cutoff
    sigma# Assume surface tension of water (mN/m) ???
    property_methods::Dict
end

struct AdjointConfig
    use_cache::Bool
    diff_method::String
    adjoint_solver
    reltol::Real
    abstol::Real
end

struct SolverConfig
    solver
    sparse::Bool
    reltol::Real
    abstol::Real
    dtinit::Real
    dtmax::Real
    positiveness::Bool
    diff_method::String
end

function SolverConfig(solver, sparse, reltol, abstol, dtinit, dtmax, positiveness)
   return SolverConfig(solver,sparse,reltol,abstol,dtinit,dtmax,positiveness,"gas") 
end

function showconfig(config::SolverConfig)
    println("===================Solver Config===================")
    println("Using solver: $(config.solver)")
    println("Using sparse jacobian: $(config.sparse)")
    println("Reltol: $(config.reltol), Abstol: $(config.abstol)")
    println("DtInit: $(config.dtinit) s, DtMax: $(config.dtmax) s")
    println("Positiveness detection: $(config.positiveness)")
    println("Jacobian method: $(config.diff_method)")
    println("===================================================")
end

function showconfig(config::GasConfig)
    println("===============Gas Simulation Config===============")
    println("Mechanism file: $(config.file)")
    println("Start time t0: $(config.start_time) s, Simulation time: $(config.simulation_time) s, Saving interval: $(config.batch_step) s")
    println("Temperature: $(config.temp) K, Relative Humidity: $(config.RH*100) %")
    println("Initial Condition (ppm): $(config.reactants_initial_dict)")
    #println("===================================================")
end

function showconfig(config::AerosolConfig)
    println("=============Aerosol Simulation Config=============")
    println("Mechanism file: $(config.file)")
    println("Start time t0: $(config.start_time) s, Simulation time: $(config.simulation_time) s, Saving interval: $(config.batch_step) s")
    println("Temperature: $(config.temp) K, Relative Humidity: $(config.RH*100) %")
    println("Initial Condition (ppm): $(config.reactants_initial_dict)")
    println("Num_bins: $(config.num_bins), Vp_cutoff: $(config.vp_cutoff) log10(Pa)")
    println("Property methods: $(config.property_methods)")
    #println("===================================================")
end