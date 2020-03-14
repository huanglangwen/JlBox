using JlBox
using DataFrames
using OrdinaryDiffEq
using Sundials
using IncompleteLU
using LinearAlgebra
#using AlgebraicMultigrid
#using DiffEqBase
#using Pardiso
#using CSV
function configure_gas()
    file="../data/MCM_mixed_test.eqn.txt"#"MCM_test.eqn.txt"MCM_APINENE.eqn.txt"MCM_mixed_test.eqn.txt
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
    tspan=(0,simulation_time)
    Cfactor= 2.55e+10 #ppb-to-molecules/cc
    reactants_initial_dict=Dict(["O3"=>18.0,"APINENE"=>30.0])#ppm ["O3"=>18.0,"APINENE"=>30.0])BUT1ENE
    constantdict=Dict([(:temp,temp),(:H2O,H2O)])
    solver=nothing
    reltol=1e-6
    abstol=1.0e-4
    positiveness=false
    use_jacobian=true
    JlBox.GasConfigure(file,temp,RH,hour_of_day,start_time,simulation_time,batch_step,
                       H2O,tspan,Cfactor,reactants_initial_dict,constantdict,solver,reltol,abstol,positiveness,use_jacobian)
end

config=configure_gas()
param_dict,reactants2ind=JlBox.prepare_gas(config)
jac_prototype=JlBox.get_sparsity(param_dict,reactants2ind)
#pc = aspreconditioner(ruge_stuben(jac_prototype))
pc = ilu(jac_prototype;τ=1e-6);
#solver=TRBDF2(linsolve = LinSolveGMRES(Pl=pc))#TRBDF2(linsolve = DiffEqBase.PardisoFactorize())#, linsolve=LinSolveGMRES(Pl=pc)
prec=(z,r,p,t,y,fy,gamma,delta,lr)->LinearAlgebra.ldiv!(z,pc,r)
solver=Sundials.CVODE_BDF(linear_solver=:GMRES,prec=prec,prec_side=1,krylov_dim=100)#:FGMRES
sol,reactants2ind=JlBox.run_simulation_gas_sparse(solver,config,jac_prototype,param_dict,reactants2ind)
num_reactants=length(reactants2ind)
ind2reactants=Dict(reactants2ind[key]=>key for key in keys(reactants2ind))
reactants=[ind2reactants[ind] for ind in 1:num_reactants]
df=DataFrames.DataFrame(transpose(sol))
DataFrames.rename!(df,[Symbol(reac) for reac in reactants])
#CSV.write("data/results_gas_jac.csv",df)
df
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
