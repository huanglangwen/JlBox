struct GasConfigure
    file#"MCM_test.eqn.txt"MCM_APINENE.eqn.txt
    temp# Kelvin
    RH# RH/100% [0 - 0.99]
    hour_of_day# Define a start time  24 hr format
    start_time# seconds, used as t0 in solver
    simulation_time# seconds
    batch_step# seconds
    temp_celsius
    Psat# Saturation VP of water vapour, to get concentration of H20
    Pw
    Wconc#kg/cm3
    H2O#Convert from kg to molecules/cc
    tspan
    Cfactor#ppb-to-molecules/cc
    reactants_initial_dict#ppm ["O3"=>18.0,"APINENE"=>30.0])
    constantdict
end

struct AerosolConfigure
    file#"MCM_test.eqn.txt"MCM_APINENE.eqn.txt
    temp# Kelvin
    RH# RH/100% [0 - 0.99]
    hour_of_day# Define a start time  24 hr format
    start_time# seconds, used as t0 in solver
    simulation_time# seconds
    batch_step# seconds
    temp_celsius
    Psat_w# Saturation VP of water vapour, to get concentration of H20
    Pw
    Wconc#kg/cm3
    H2O#Convert from kg to molecules/cc
    tspan
    Cfactor#ppb-to-molecules/cc
    reactants_initial_dict#ppb BUT1ENE APINENE
    constantdict
    
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
    R_gas#Ideal gas constant [kg m2 s-2 K-1 mol-1]
    NA#Avogadros number
    sigma# Assume surface tension of water (mN/m) ???
    property_methods
end

