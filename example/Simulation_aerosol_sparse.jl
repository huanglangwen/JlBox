using JlBox
using DataFrames
using OrdinaryDiffEq
using Sundials
using IncompleteLU
using LinearAlgebra
#using CuArrays
#using CSV

function configure_aerosol()
    file="../data/MCM_mixed_test.eqn.txt"#"MCM_test.eqn.txt"MCM_APINENE.eqn.txt
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
    tspan=(0.,simulation_time)
    Cfactor= 2.55e+10 #ppb-to-molecules/cc
    reactants_initial_dict=Dict(["O3"=>18.0,"APINENE"=>30.0,"H2O"=>H2O/Cfactor])#ppb BUT1ENE APINENE
    constantdict=Dict([(:temp,temp)])
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
    R_gas=8.3144598 #Ideal gas constant [kg m2 s-2 K-1 mol-1]
    NA=6.0221409e+23 #Avogadros number
    sigma=72.0e-3 # Assume surface tension of water (mN/m) ???
    property_methods=Dict("bp"=>"joback_and_reid","vp"=>"nannoolal","critical"=>"nannoolal","density"=>"girolami")
    diff_method="analytical"
    solver=nothing
    reltol=1e-4
    abstol=1.0e-2
    positiveness=false
    use_jacobian=true
    JlBox.AerosolConfigure(file,temp,RH,hour_of_day,start_time,simulation_time,batch_step,
                           H2O,tspan,Cfactor,reactants_initial_dict,constantdict,num_bins,
                           total_conc,size_std,lowersize,uppersize,meansize,y_core_init,
                           core_density_array,core_mw,core_dissociation,vp_cutoff,R_gas,
                           NA,sigma,property_methods,diff_method,solver,reltol,abstol,positiveness,use_jacobian)
end

config=configure_aerosol()
param_dict,reactants2ind,y_cond=JlBox.prepare_aerosol(config)
jac_prototype=JlBox.get_sparsity_aerosol(param_dict,reactants2ind,y_cond)
pc = ilu(jac_prototype;τ=1e-6);
#solver=TRBDF2(autodiff=false)#linsolve=LinSolveGPUFactorize()#Sundials.CVODE_BDF()
prec=(z,r,p,t,y,fy,gamma,delta,lr)->LinearAlgebra.ldiv!(z,pc,r)
solver=Sundials.CVODE_BDF(linear_solver=:FGMRES,prec=prec,prec_side=1,krylov_dim=100)#:FGMRES
@time sol,reactants2ind,SOA_array,num_reactants,_=JlBox.run_simulation_aerosol_sparse(solver,config,jac_prototype,param_dict,reactants2ind,y_cond)
#num_reactants=length(reactants2ind)
sol_mtx=transpose(sol)
ind2reactants=Dict(reactants2ind[key]=>key for key in keys(reactants2ind))
reactants=[Symbol(ind2reactants[ind]) for ind in 1:num_reactants]
t_length=size(sol_mtx)[1]
t_index=range(0,stop=config.simulation_time,length=t_length)
df_SOA=DataFrames.DataFrame(Time=t_index,SOA=SOA_array)[:,[:Time,:SOA]]
df=DataFrames.DataFrame(sol_mtx[1:end,1:num_reactants])
DataFrames.rename!(df,reactants)
#CSV.write("/data/jlbox_results.csv",df)
#CSV.write("/data/jlbox_SOA.csv",df_SOA)
df_SOA
