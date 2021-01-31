using Test, JlBox, DataFrames, Sundials
using OrdinaryDiffEq
using PyCall
using LinearAlgebra
using IncompleteLU

function configure_gas()
    file=joinpath(@__DIR__,"../data/MCM_APINENE.eqn.txt")#"MCM_test.eqn.txt"MCM_APINENE.eqn.txt"MCM_mixed_test.eqn.txt
    temp=298.15 # Kelvin
    RH=0.5 # RH/100% [0 - 0.99]
    hour_of_day=12.0 # Define a start time  24 hr format
    start_time=hour_of_day*60*60 # seconds, used as t0 in solver
    simulation_time= 10800.0 # seconds
    batch_step=100.0 # seconds
    temp_celsius=temp-273.15
    Psat=610.78*exp((temp_celsius/(temp_celsius+238.3))*17.2694)# Saturation VP of water vapour, to get concentration of H20
    Pw=RH*Psat
    Wconc=0.002166*(Pw/(temp_celsius+273.16))*1.0e-6 #kg/cm3
    H2O=Wconc*(1.0/(18.0e-3))*6.0221409e+23#Convert from kg to molecules/cc
    Cfactor= 2.55e+10 #ppb-to-molecules/cc
    reactants_initial_dict=Dict(["O3"=>18.0,"APINENE"=>30.0])#ppm ["O3"=>18.0,"APINENE"=>30.0])BUT1ENE
    constant_dict=Dict([(:temp,temp),(:H2O,H2O)])
    dec=23.79
    lat=50.0
    photolysis_config=JlBox.DiurnalPhotolysisConfig(dec, lat)
    config=JlBox.GasConfig(file,temp,RH,start_time,simulation_time,batch_step,
                       H2O,Cfactor,reactants_initial_dict,constant_dict,photolysis_config)
    config
end

function configure_gas_solver_dense()
    solver=Sundials.CVODE_BDF()
    sparse=false
    reltol=1e-6
    abstol=1.0e-3
    dtinit=1e-6
    dtmax=100.0
    positiveness=false
    solverconfig=JlBox.SolverConfig(solver,sparse,reltol,abstol,dtinit,dtmax,positiveness)
    solverconfig
end

function configure_gas_solver_sparse()
    prec = JlBox.default_prec()
    psetup = JlBox.default_psetup_gas()
    ndim=100
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

function configure_gas_solver_sparse_v()#With matrix free operator
    prec = JlBox.default_prec()
    psetup = JlBox.default_psetup_gas()
    ndim=100
    solver=Sundials.CVODE_BDF(linear_solver=:FGMRES,prec=prec,psetup=psetup,prec_side=2,krylov_dim=ndim)
    sparse=true
    reltol=1e-6
    abstol=1.0e-3
    dtinit=1e-6
    dtmax=100.0
    positiveness=false
    solverconfig=JlBox.SolverConfig(solver,sparse,reltol,abstol,dtinit,dtmax,positiveness,"gas_v")
    solverconfig
end
@testset "Gas Phase" begin
    @testset "Gas Phase Dense" begin 
        config = configure_gas()
        solverconfig = configure_gas_solver_dense()
        sol, reactants2ind, _ = JlBox.run_simulation(config, solverconfig)
        df = JlBox.postprocess_gas(sol, reactants2ind, config)
        # @test sum(df[1,2:]) ≈ sum(df[end,2:]) rtol=1e-2
        @test df[end,:APINENE] ≈ 4.101378622128155e11 rtol=1e-4
    end

    @testset "Gas Phase Sparse" begin 
        config = configure_gas()
        solverconfig = configure_gas_solver_sparse()
        sol, reactants2ind, _ = JlBox.run_simulation(config, solverconfig)
        df = JlBox.postprocess_gas(sol, reactants2ind, config)
        # @test sum(df[1,2:]) ≈ sum(df[end,2:]) rtol=1e-2
        @test df[end,:APINENE] ≈ 4.101378622128155e11 rtol=1e-4
    end

    @testset "Gas Phase Sparse with Matrix-free operator" begin 
        config = configure_gas()
        solverconfig = configure_gas_solver_sparse_v()
        sol, reactants2ind, _ = JlBox.run_simulation(config, solverconfig)
        df = JlBox.postprocess_gas(sol, reactants2ind, config)
        # @test sum(df[1,2:]) ≈ sum(df[end,2:]) rtol=1e-2
        @test df[end,:APINENE] ≈ 4.101378622128155e11 rtol=1e-4
    end
end

function configure_aerosol()
    file=joinpath(@__DIR__,"../data/MCM_APINENE.eqn.txt")#"MCM_test.eqn.txt"MCM_APINENE.eqn.txt
    temp=288.15 # Kelvin
    RH=0.5 # RH/100% [0 - 0.99]
    hour_of_day=12.0 # Define a start time  24 hr format
    start_time=hour_of_day*60*60 # seconds, used as t0 in solver
    simulation_time= 3600.0 # seconds
    batch_step=300.0 # seconds
    temp_celsius=temp-273.15
    Psat_w=610.78*exp((temp_celsius/(temp_celsius+238.3))*17.2694)# Saturation VP of water vapour, to get concentration of H20
    Pw=RH*Psat_w
    Wconc=0.002166*(Pw/(temp_celsius+273.16))*1.0e-6 #kg/cm3
    H2O=Wconc*(1.0/(18.0e-3))*6.0221409e+23#Convert from kg to molecules/cc
    Cfactor= 2.55e+10 #ppb-to-molecules/cc
    reactants_initial_dict=Dict(["O3"=>18.0,"APINENE"=>30.0,"H2O"=>H2O/Cfactor])#ppb BUT1ENE APINENE
    constant_dict=Dict([(:temp,temp)])
    dec=23.79
    lat=50.0
    photolysis_config=JlBox.DiurnalPhotolysisConfig(dec, lat)
    num_bins=16

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
                           sigma,property_methods)
    config
end

function configure_aerosol_solver_dense()
    solver=TRBDF2()
    sparse=false
    reltol=1e-4
    abstol=1.0e-2
    dtinit=1e-6
    dtmax=100.0
    positiveness=false
    diff_method="fine_seeding"
    solverconfig=JlBox.SolverConfig(solver,sparse,reltol,abstol,dtinit,dtmax,positiveness,diff_method)
    solverconfig
end

function configure_aerosol_solver_sparse()
    prec = JlBox.default_prec()
    psetup = JlBox.default_psetup("fine_seeding","fine_analytical", 20)
    ndim=500
    solver=Sundials.CVODE_BDF(linear_solver=:FGMRES,prec=prec,psetup=psetup,prec_side=2,krylov_dim=ndim)
    sparse=true
    reltol=1e-6
    abstol=1.0e-3
    dtinit=1e-6
    dtmax=100.0
    positiveness=false
    diff_method="fine_seeding"
    solverconfig=JlBox.SolverConfig(solver,sparse,reltol,abstol,dtinit,dtmax,positiveness,diff_method)
    solverconfig
end

function configure_aerosol_solver_sparse_v()
    prec = JlBox.default_prec()
    psetup = JlBox.default_psetup("fine_seeding","fine_analytical", 20)
    ndim=500
    solver=Sundials.CVODE_BDF(linear_solver=:FGMRES,prec=prec,psetup=psetup,prec_side=2,krylov_dim=ndim)
    sparse=true
    reltol=1e-6
    abstol=1.0e-3
    dtinit=1e-6
    dtmax=100.0
    positiveness=false
    diff_method="fine_seeding_v"
    solverconfig=JlBox.SolverConfig(solver,sparse,reltol,abstol,dtinit,dtmax,positiveness,diff_method)
    solverconfig
end

@testset "Mixed Phase" begin
    @testset "Mixed Phase Dense" begin
        config = configure_aerosol()
        solverconfig = configure_aerosol_solver_dense()
        @time sol, reactants2ind, param_dict = JlBox.run_simulation(config, solverconfig)
        df = JlBox.postprocess_gas(sol, reactants2ind, config)
        df_SOA = JlBox.postprocess_aerosol(sol, param_dict, config)
        # @test sum(df[1,2:]) ≈ sum(df[end,2:]) rtol=1e-6
        @test df_SOA[13,:SOA] ≈ 9.165607888151833 rtol=1e-4
    end

    @testset "Mixed Phase Sparse" begin
        config = configure_aerosol()
        solverconfig = configure_aerosol_solver_sparse()
        @time sol, reactants2ind, param_dict = JlBox.run_simulation(config, solverconfig)
        df = JlBox.postprocess_gas(sol, reactants2ind, config)
        df_SOA = JlBox.postprocess_aerosol(sol, param_dict, config)
        # @test sum(df[1,2:]) ≈ sum(df[end,2:]) rtol=1e-6
        @test df_SOA[13,:SOA] ≈ 9.165607888151833 rtol=1e-4
    end

    @testset "Mixed Phase Sparse with Matrix-free operator" begin
        config = configure_aerosol()
        solverconfig = configure_aerosol_solver_sparse_v()
        @time sol, reactants2ind, param_dict = JlBox.run_simulation(config, solverconfig)
        df = JlBox.postprocess_gas(sol, reactants2ind, config)
        df_SOA = JlBox.postprocess_aerosol(sol, param_dict, config)
        # @test sum(df[1,2:]) ≈ sum(df[end,2:]) rtol=1e-6
        @test df_SOA[13,:SOA] ≈ 9.165607888151833 rtol=1e-4
    end

end
