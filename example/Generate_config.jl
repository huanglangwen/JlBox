
#https://www.atmos-chem-phys.net/18/15743/2018/acp-18-15743-2018.pdf
#https://www.sciencedirect.com/science/article/pii/S1352231010008459
#https://www.sciencedirect.com/science/article/pii/S0187623613710622
function configure_aerosol(time_hour, temp, RH, reactants_initial_dict, io)
    file=joinpath(@__DIR__,"../data/MCM_mixed_test.eqn.txt")#"MCM_test.eqn.txt"MCM_APINENE.eqn.txtMCM_BCARY.eqn.txt
    hour_of_day=6.0 # Define a start time  24 hr format
    start_time=hour_of_day*60*60 # seconds, used as t0 in solver
    simulation_time= time_hour*3600.0 # seconds
    batch_step=300.0 # seconds
    temp_celsius=temp-273.15
    Psat_w=610.78*exp((temp_celsius/(temp_celsius+238.3))*17.2694)# Saturation VP of water vapour, to get concentration of H20
    Pw=RH*Psat_w
    Wconc=0.002166*(Pw/(temp_celsius+273.16))*1.0e-6 #kg/cm3
    H2O=Wconc*(1.0/(18.0e-3))*6.0221409e+23#Convert from kg to molecules/cc
    Cfactor= 2.55e+10 #ppb-to-molecules/cc
    constant_dict=Dict([(:temp,temp)])
    dec=23.79 #Summer Solstice
    lat=39.466667 #Valencia, Spain
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
                           sigma,property_methods, io)
    config
end

function experiment_b(time_hour, i, io::Base.IO = open("results/log_Experiment_B_rank"*string(i)*".log", "w"), RH::Union{Nothing,Float64} = nothing)
    isoprene = [107, 92, 122, 0, 99, 87, 55] #C5H8
    apinene = [66, 50, 71, 63, 59, 50, 79] #APINENE
    limonene = [58, 50, 40, 65, 53, 51, 76] #LIMONENE
    NO = [34, 48, 41, 32, 150, 244, 198] #NO
    NO2 = [128, 0, 0, 0, 0, 89, 0] #NO2
    HONO = [99, 87, 53, 101, 307, 40, 165] #HONO
    SO2 = [0, 0, 0, 0, 0, 513, 461] #SO2
    Tlo = [302, 298, 297, 294, 295, 295, 302]
    Thi = [307, 300, 300, 298, 297, 300, 305]
    RHlo = [0.5, 26, 19, 8, 8, 15, 20] ./ 100.0
    RHhi = [3, 30, 22, 13, 11, 19, 30] ./ 100.0
    RHval = RH isa Nothing ? (RHlo[i]+RHhi[i])/2 : RH
    configure_aerosol(time_hour, (Tlo[i]+Thi[i])/2, RHval, 
        Dict("C5H8"=>isoprene[i], "APINENE"=>apinene[i],
            "LIMONENE"=>limonene[i], "NO"=>NO[i],
            "NO2"=>NO2[i], "HONO"=>HONO[i], "SO2"=>SO2[i]),
        io) 
end
function experiment_b(time_hour)
    [experiment_b(time_hour, i) for i in 1:7]
end

function experiment_a(time_hour, i, io::Base.IO = open("results/log_Experiment_A_rank$(i).log", "w"), RH::Union{Nothing,Float64} = nothing)
    toluene = [102, 200, 48, 98, 97, 93, 107, 116, 81] #TOLUENE
    oxylene = [22, 49, 11, 24, 21, 22, 26, 29, 21] #OXYL
    TMB = [153, 300, 106, 160, 146, 146, 160, 19, 118] #TM135B
    Octane = [85, 155, 42, 79, 81, 78, 89, 10, 65] #NC8H18
    NO = [19, 23, 23, 37, 4, 21, 21, 57, 31] #NO
    NO2 = [0, 0, 0, 0, 8, 0, 0, 0, 0] #NO2
    HONO = [99, 75, 71, 156, 52, 94, 89, 119, 90] #HONO
    Tlo = [299, 302, 302, 297, 297, 300, 306, 302, 299]
    Thi = [305, 305, 307, 307, 308, 308, 309, 305, 303]
    RHlo = [10, 9, 6, 6, 7, 0.4, 7, 15, 28] ./ 100.0
    RHhi = [16, 18, 14, 13, 14, 0.4, 10, 18, 37] ./ 100.0
    RHval = RH isa Nothing ? (RHlo[i]+RHhi[i])/2 : RH
    configure_aerosol(time_hour, (Tlo[i]+Thi[i])/2, RHval, 
        Dict("TOLUENE"=>toluene[i], "OXYL"=>oxylene[i],
            "TM135B"=>TMB[i], "NO"=>NO[i],
            "NO2"=>NO2[i], "HONO"=>HONO[i], "NC8H18"=>Octane[i]),
        io) 
end
function experiment_a(time_hour)
    [experiment_a(time_hour, i) for i in 1:9]
end

function configure_aerosol_solver_sparse()
    prec = JlBox.default_prec()
    psetup = JlBox.default_psetup("fine_seeding","fine_analytical", 200)
    ndim=5000
    solver=Sundials.CVODE_BDF(linear_solver=:FGMRES,prec=prec,psetup=psetup,prec_side=2,krylov_dim=ndim)
    sparse=true
    reltol=1e-6
    abstol=1.0e-3
    dtinit=1e-8
    dtmax=100.0
    positiveness=false
    diff_method="fine_seeding"
    solverconfig=JlBox.SolverConfig(solver,sparse,reltol,abstol,dtinit,dtmax,positiveness,diff_method)
    solverconfig
end
