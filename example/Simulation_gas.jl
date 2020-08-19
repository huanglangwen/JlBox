using JlBox
using DataFrames
using OrdinaryDiffEq
using Sundials
using CuArrays

function configure_gas()
    file="../data/MCM_test.eqn.txt"#"MCM_test.eqn.txt"MCM_APINENE.eqn.txt"MCM_mixed_test.eqn.txt"MCM_ISOPRENE.eqn.txt
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
    Cfactor=2.55e+10 #ppb-to-molecules/cc
    reactants_initial_dict=Dict(["O3"=>18.0,"BUT1ENE"=>30.0])#ppm ["O3"=>18.0,"APINENE"=>30.0])BUT1ENE"C5H8
    constant_dict=Dict([(:temp,temp),(:H2O,H2O)])
    dec=23.79
    lat=50.0
    photolysis_config=JlBox.DiurnalPhotolysisConfig(dec, lat)
    solver=TRBDF2(autodiff=false, linsolve=LinSolveGPUFactorize())#TRBDF2(autodiff=false)
    #solver=Sundials.CVODE_BDF()#:FGMRES
    reltol=1e-6
    abstol=1.0e-4
    dtinit=1e-6
    dtmax=100.0
    positiveness=false
    use_jacobian=true
    sparse=false
    config=JlBox.GasConfig(file,temp,RH,start_time,simulation_time,batch_step,
                       H2O,Cfactor,reactants_initial_dict,constant_dict,photolysis_config)
    solverconfig=JlBox.SolverConfig(solver,sparse,reltol,abstol,dtinit,dtmax,positiveness)
    config,solverconfig
end

#Profile.init(n = 10^7, delay = 5.)
config, solverconfig = configure_gas()
sol, reactants2ind, _ = JlBox.run_simulation(config, solverconfig)
df = JlBox.postprocess_gas(sol, reactants2ind)
df
#CSV.write("data/results_gas_jac.csv",df)
#@profile run_simulation()
#open("prof.txt", "w") do s
#    Profile.print(IOContext(s, :displaysize => (1000, 500)))
#end
#using ProfileView
#ProfileView.view()
#Profile.clear()
#using Plots
#plot(log10.(sol[reactants2ind["BUT1ENE"],:]))
#plot!(log10.(sol[reactants2ind["C2H5O2"],:]))
