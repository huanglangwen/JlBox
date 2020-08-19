abstract type JlBoxConfig end
abstract type PhotolysisConfig end
struct DiurnalPhotolysisConfig{T1<:Real, T2<:Real} <: PhotolysisConfig
    declination::T1
    latitude::T2
end
struct FixedPhotolysisConfig{T1<:Real} <: PhotolysisConfig
    cos_zenith::T1
end
struct GasConfig <: JlBoxConfig
    file::String#"MCM_test.eqn.txt"MCM_APINENE.eqn.txt
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
    constant_dict::Dict
    photolysis_config::PhotolysisConfig
    io::Base.IO
end

function GasConfig(file, temp, RH, start_time, simulation_time, batch_step, H2O, Cfactor,
                   reactants_initial_dict, constant_dict, photolysis_config)
    GasConfig(file, temp, RH, start_time, simulation_time, batch_step, H2O, Cfactor,
              reactants_initial_dict, constant_dict, photolysis_config, Base.stdout)
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
    constant_dict::Dict
    photolysis_config::PhotolysisConfig
    
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
    io::Base.IO
end

function AerosolConfig(file, temp, RH, start_time, simulation_time, batch_step, H2O, Cfactor,
                       reactants_initial_dict, constant_dict, photolysis_config, num_bins, total_conc, size_std,
                       lowersize, uppersize, meansize, y_core_init, core_density_array, core_mw,
                       core_dissociation, vp_cutoff, sigma, property_methods)
    AerosolConfig(file, temp, RH, start_time, simulation_time, batch_step, H2O, Cfactor,
                  reactants_initial_dict, constant_dict, photolysis_config, num_bins, total_conc, size_std,
                  lowersize, uppersize, meansize, y_core_init, core_density_array, core_mw,
                  core_dissociation, vp_cutoff, sigma, property_methods, Base.stdout)
end

struct AdjointConfig
    use_cache::Bool
    diff_method::String
    adjoint_solver
    reltol::Real
    abstol::Real
    io::Base.IO
end

function AdjointConfig(use_cache, diff_method, adjoint_solver, reltol, abstol)
    AdjointConfig(use_cache, diff_method, adjoint_solver, reltol, abstol, Base.stdout)
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

function showconfig(config::SolverConfig, io::Base.IO = Base.stdout)
    println(io, "===================Solver Config===================")
    println(io, "Using solver: $(config.solver)")
    println(io, "Using sparse jacobian: $(config.sparse)")
    println(io, "Reltol: $(config.reltol), Abstol: $(config.abstol)")
    println(io, "DtInit: $(config.dtinit) s, DtMax: $(config.dtmax) s")
    println(io, "Positiveness detection: $(config.positiveness)")
    println(io, "Jacobian method: $(config.diff_method)")
    println(io, "===================================================")
    flush(io)
end

function showconfig(config::GasConfig)
    io = config.io
    println(io, "===============Gas Simulation Config===============")
    println(io, "Mechanism file: $(config.file)")
    println(io, "Start time t0: $(config.start_time) s, Simulation time: $(config.simulation_time) s, Saving interval: $(config.batch_step) s")
    println(io, "Temperature: $(config.temp) K, Relative Humidity: $(config.RH*100) %")
    println(io, "Initial Condition (ppm): $(config.reactants_initial_dict)")
    #println(io, "===================================================")
    flush(io)
end

function showconfig(config::AerosolConfig)
    io = config.io
    println(io, "=============Aerosol Simulation Config=============")
    println(io, "Mechanism file: $(config.file)")
    println(io, "Start time t0: $(config.start_time) s, Simulation time: $(config.simulation_time) s, Saving interval: $(config.batch_step) s")
    println(io, "Temperature: $(config.temp) K, Relative Humidity: $(config.RH*100) %")
    println(io, "Initial Condition (ppm): $(config.reactants_initial_dict)")
    println(io, "Num_bins: $(config.num_bins), Vp_cutoff: $(config.vp_cutoff) log10(Pa)")
    println(io, "Property methods: $(config.property_methods)")
    #println(io, "===================================================")
    flush(io)
end