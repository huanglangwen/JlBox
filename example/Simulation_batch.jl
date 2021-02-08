using ArgParse
using Sundials
using JlBox

function configure_gas_solver_sparse()
    prec = JlBox.default_prec()
    psetup = JlBox.default_psetup_gas()
    ndim=1000
    solver=Sundials.CVODE_BDF(linear_solver=:FGMRES,prec=prec,psetup=psetup,prec_side=2,krylov_dim=ndim)
    sparse=true
    reltol=1e-6
    abstol=1.0e-3
    dtinit=1e-6
    dtmax=100.0
    positiveness=false
    solverconfig=JlBox.SolverConfig(solver,sparse,reltol,abstol,dtinit,dtmax,positiveness)
    solverconfig
end

function configure_aerosol_solver_sparse()
    prec = JlBox.default_prec()
    psetup = JlBox.default_psetup("fine_seeding","fine_analytical", 200)
    ndim=5000
    solver=Sundials.CVODE_BDF(linear_solver=:FGMRES,prec=prec,psetup=psetup,prec_side=2,krylov_dim=ndim)
    sparse=true
    reltol=1e-4
    abstol=1.0e-2
    dtinit=1e-6
    dtmax=100.0
    positiveness=false
    diff_method="fine_seeding"
    solverconfig=JlBox.SolverConfig(solver,sparse,reltol,abstol,dtinit,dtmax,positiveness,diff_method)
    solverconfig
end

function configure_gas_solver_dense()
    solver=Sundials.CVODE_BDF()#OrdinaryDiffEq.TRBDF2(autodiff=false)
    sparse=false
    reltol=1e-6
    abstol=1.0e-3
    dtinit=1e-6
    dtmax=100.0
    positiveness=false
    solverconfig=JlBox.SolverConfig(solver,sparse,reltol,abstol,dtinit,dtmax,positiveness)
    solverconfig
end

function configure_aerosol_solver_dense()
    solver=Sundials.CVODE_BDF()#OrdinaryDiffEq.TRBDF2(autodiff=false)
    sparse=false
    reltol=1e-4
    abstol=1.0e-2
    dtinit=1e-6
    dtmax=100.0
    positiveness=false
    diff_method="fine_analytical"
    solverconfig=JlBox.SolverConfig(solver,sparse,reltol,abstol,dtinit,dtmax,positiveness,diff_method)
    solverconfig
end

function configure_gas(file, init_reac, time, io)
    temp=288.15 # Kelvin
    RH=0.5 # RH/100% [0 - 0.99]
    hour_of_day=12.0 # Define a start time  24 hr format
    start_time=hour_of_day*60*60 # seconds, used as t0 in solver
    simulation_time=time # seconds
    batch_step=time # seconds
    temp_celsius=temp-273.15
    Psat=610.78*exp((temp_celsius/(temp_celsius+238.3))*17.2694)# Saturation VP of water vapour, to get concentration of H20
    Pw=RH*Psat
    Wconc=0.002166*(Pw/(temp_celsius+273.16))*1.0e-6 #kg/cm3
    H2O=Wconc*(1.0/(18.0e-3))*6.0221409e+23#Convert from kg to molecules/cc
    Cfactor= 2.55e+10 #ppb-to-molecules/cc
    reactants_initial_dict=Dict(["O3"=>18.0,init_reac=>30.0])#ppm ["O3"=>18.0,"APINENE"=>30.0])BUT1ENE
    constant_dict=Dict([(:temp,temp),(:H2O,H2O)])
    dec=23.79
    lat=50.0
    photolysis_config=JlBox.DiurnalPhotolysisConfig(dec, lat)
    config=JlBox.GasConfig(file,temp,RH,start_time,simulation_time,batch_step,
                       H2O,Cfactor,reactants_initial_dict,constant_dict,photolysis_config, io)
    config
end

function configure_aerosol(file, num_bins, init_reac, time, io)
    temp=288.15 # Kelvin
    RH=0.5 # RH/100% [0 - 0.99]
    hour_of_day=12.0 # Define a start time  24 hr format
    start_time=hour_of_day*60*60 # seconds, used as t0 in solver
    simulation_time=time # seconds
    batch_step=time # seconds
    temp_celsius=temp-273.15
    Psat_w=610.78*exp((temp_celsius/(temp_celsius+238.3))*17.2694)# Saturation VP of water vapour, to get concentration of H20
    Pw=RH*Psat_w
    Wconc=0.002166*(Pw/(temp_celsius+273.16))*1.0e-6 #kg/cm3
    H2O=Wconc*(1.0/(18.0e-3))*6.0221409e+23#Convert from kg to molecules/cc
    Cfactor= 2.55e+10 #ppb-to-molecules/cc
    reactants_initial_dict=Dict(["O3"=>18.0,init_reac=>30.0,"H2O"=>H2O/Cfactor])#ppb BUT1ENE APINENE
    constant_dict=Dict([(:temp,temp)])
    dec=23.79
    lat=50.0
    photolysis_config=JlBox.DiurnalPhotolysisConfig(dec, lat)
    num_bins=num_bins

    #Lognormal Distribution
    total_conc=100 #Total particles per cc
    size_std=2.2 #Standard Deviation
    lowersize=0.01 #microns
    uppersize=1.0 #microns
    meansize=0.2 #microns

    # - Specify the core material. 
    # This code is currently setup to consider *ammonium sulphate* as the core
    y_core_init=1.0e-3.+zeros(Float64,num_bins) #Will hold concentration of core material, only initialise here [molecules/cc] 
    core_density_array=1770.0.+zeros(Float64,num_bins) #[kg/m3] - need to make sure this matches core definition above
    core_mw=132.14.+zeros(Float64,num_bins) #[g/mol]
    core_dissociation=3.0 #Define this according to choice of core type. Please note this value might change

    vp_cutoff=-6.0
    sigma=72.0e-3 # Assume surface tension of water (mN/m) ???
    property_methods=Dict("bp"=>"joback_and_reid","vp"=>"nannoolal","critical"=>"nannoolal","density"=>"girolami")
    config=JlBox.AerosolConfig(file,temp,RH,start_time,simulation_time,batch_step,
                           H2O,Cfactor,reactants_initial_dict,constant_dict,photolysis_config,num_bins,
                           total_conc,size_std,lowersize,uppersize,meansize,y_core_init,
                           core_density_array,core_mw,core_dissociation,vp_cutoff,
                           sigma,property_methods, io)
    config
end

function generate_config(simulation_type, mechanism, num_bins, time, io)
    mechanism_prefix = mechanism == "FULL" ? "mixed_test" : mechanism
    init_reac = mechanism == "FULL" ? "APINENE" : mechanism
    mechanism_file = joinpath(@__DIR__,"../data/MCM_$(mechanism_prefix).eqn.txt")
    if simulation_type == "GAS"
        return configure_gas(mechanism_file, init_reac, time, io)
    elseif simulation_type == "MIXED"
        return configure_aerosol(mechanism_file, num_bins, init_reac, time, io)
    else
        exit(-1)
    end
end

function generate_solver(simulation_type, jacobian_type)
    if jacobian_type == "DENSE"
        if simulation_type == "GAS"
            return configure_gas_solver_dense()
        elseif simulation_type == "MIXED"
            return configure_aerosol_solver_dense()
        else
            exit(-1)
        end
    elseif jacobian_type == "SPARSE"
        if simulation_type == "GAS"
            return configure_gas_solver_sparse()
        elseif simulation_type == "MIXED"
            return configure_aerosol_solver_sparse()
        else
            exit(-1)
        end
    else
        exit(-1)
    end
end

function parse_cmd()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--mechanism", "-m"
            help = "Specify mechanism: APINENE, BCARY, LIMONENE or FULL"
            arg_type = String
        "--type", "-s"
            help = "Specify simulation type: GAS, MIXED"
            arg_type = String
        "--jacobian-type", "-j"
            help = "Specify jacobian type: DENSE, SPARSE"
            arg_type = String
        "--num-bins", "-b"
            help = "Specify number of bins"
            arg_type = Int
            default = 0
        "--verbose-output"
            help = "Specify verbose output path, default to stdout"
            arg_type = String
            default = "stdout"
        "--time", "-t"
            help = "Simulation period in seconds"
            arg_type = Float64
            default = 3600.0
    end
    return parse_args(s)
end

parsed_args = parse_cmd()
iostr = parsed_args["verbose-output"]
simulation_type = parsed_args["type"]
mechanism = parsed_args["mechanism"]
jacobian_type = parsed_args["jacobian-type"]
num_bins = parsed_args["num-bins"]
time = parsed_args["time"]
io = iostr == "stdout" ? Base.stdout : open(iostr, "a")
println("$(simulation_type), $(jacobian_type), $(mechanism), $(num_bins), $(time)")
config = generate_config(simulation_type, mechanism, num_bins, time, io)
solverconfig = generate_solver(simulation_type, jacobian_type)
JlBox.run_simulation(config, solverconfig)

# julia ~/Code/JlBox/example/Simulation_batch.jl -m APINENE -s MIXED -j SPARSE -b 4  --verbose-output batch.log