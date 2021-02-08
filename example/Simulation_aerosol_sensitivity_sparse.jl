using JlBox
using DataFrames
using Sundials
using OrdinaryDiffEq
#using CSV

function configure_aerosol()
    file=joinpath(@__DIR__,"../data/MCM_mixed_test.eqn.txt")#"MCM_test.eqn.txt"MCM_APINENE.eqn.txt"MCM_mixed_test
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
    diff_method="fine_seeding"
    aerosolconfig=JlBox.AerosolConfig(file,temp,RH,start_time,simulation_time,batch_step,
                           H2O,Cfactor,reactants_initial_dict,constant_dict,photolysis_config,num_bins,
                           total_conc,size_std,lowersize,uppersize,meansize,y_core_init,
                           core_density_array,core_mw,core_dissociation,vp_cutoff,
                           sigma,property_methods)
    prec = JlBox.default_prec()
    default_psetup = JlBox.default_psetup("fine_seeding","fine_analytical", 200)
    jac1 = JlBox.select_jacobian("fine_seeding", nothing)
    jac2 = JlBox.select_jacobian("fine_analytical", nothing)
    psetup = function (p,t,u,du,jok,jcurPtr,gamma)#factorize(I-gamma*J)
        if !haskey(p, "preconditioner")
            n,_ = size(p["jac_mtx"])
            jac2(p["jac_mtx"],u,p,t)
            p["trimat"]=JlBox.sp2tridiag(p["jac_mtx"],gamma)
            p["preconditioner"]=JlBox.TridiagonalLUCache(p["trimat"],zeros(n))
            JlBox.factorize!(p["preconditioner"],p["trimat"])
        else
            jac2(p["jac_mtx"],u,p,t)
            JlBox.sp2tridiag!(p["trimat"],p["jac_mtx"],gamma)
            JlBox.factorize!(p["preconditioner"],p["trimat"])
            @debug "PSETUP CALLED"
        end
        jcurPtr[]=false
    end
    ndim=5000
    solver=Sundials.CVODE_BDF(linear_solver=:FGMRES,prec=prec,psetup=default_psetup,prec_side=2,krylov_dim=ndim)
    sparse=true
    reltol=1e-6
    abstol=1.0e-3
    dtinit=1e-8
    dtmax=100.0
    positiveness=false
    diff_method="fine_seeding"
    solverconfig=JlBox.SolverConfig(solver,sparse,reltol,abstol,dtinit,dtmax,positiveness,diff_method)
    use_cache=true
    adjoint_solver=Sundials.CVODE_BDF(linear_solver=:FGMRES,prec=prec,psetup=psetup,prec_side=2,krylov_dim=ndim)
    reltol=1e-8
    abstol=1e-6
    adjointconfig=JlBox.AdjointConfig(use_cache,diff_method,adjoint_solver,reltol,abstol)
    aerosolconfig,solverconfig,adjointconfig
end

aerosolconfig,solverconfig,adjointconfig=configure_aerosol()
@time dSOA_mass_drate,dSOA_mass_percentk=JlBox.run_simulation_aerosol_adjoint(aerosolconfig,solverconfig,adjointconfig)
df=DataFrames.DataFrame(transpose(dSOA_mass_drate))
df2=DataFrames.DataFrame(transpose(dSOA_mass_percentk))
#CSV.write("/data/jlbox_sensitivity_dSOAdrate_results.csv",df)
#CSV.write("/data/jlbox_sensitivity_dSOApercentk_results.csv",df2)
df2
