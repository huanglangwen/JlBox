const file="data/MCM_APINENE.eqn.txt"#"MCM_test.eqn.txt"
const temp=288.15 # Kelvin
const RH=0.5 # RH/100% [0 - 0.99]
const hour_of_day=12.0 # Define a start time  24 hr format
const start_time=hour_of_day*60*60 # seconds, used as t0 in solver
const simulation_time= 3600.0 # seconds
const batch_step=300.0 # seconds
const temp_celsius=temp-273.15
const Psat_w=610.78*exp((temp_celsius/(temp_celsius+238.3))*17.2694)# Saturation VP of water vapour, to get concentration of H20
const Pw=RH*Psat_w
const Wconc=0.002166*(Pw/(temp_celsius+273.16))*1.0e-6 #kg/cm3
const H2O=Wconc*(1.0/(18.0e-3))*6.0221409e+23#Convert from kg to molecules/cc
const tspan=(0,simulation_time)
const Cfactor= 2.55e+10 #ppb-to-molecules/cc
const reactants_initial_dict=Dict(["O3"=>18.0,"APINENE"=>30.0,"H2O"=>H2O/Cfactor])#ppb 
const constantdict=Dict([(:temp,temp)])

const num_bins=16
#Lognormal Distribution
const total_conc=100 #Total particles per cc
const size_std=2.2 #Standard Deviation
const lowersize=0.01 #microns
const uppersize=1.0 #microns
const meansize=0.2 #microns

# - Specify the core material. 
# This code is currently setup to consider *ammonium sulphate* as the core
const y_core_init=1.0e-3+zeros(Float64,num_bins) #Will hold concentration of core material, only initialise here [molecules/cc] 
const core_density_array=1770.0+zeros(Float64,num_bins) #[kg/m3] - need to make sure this matches core definition above
const core_mw=132.14+zeros(Float64,num_bins) #[g/mol]
const core_dissociation=3.0 #Define this according to choice of core type. Please note this value might change

const vp_cutoff=-6.0
const R_gas=8.3144598 #Ideal gas constant [kg m2 s-2 K-1 mol-1]
const NA=6.0221409e+23 #Avogadros number
const sigma=72.0e-3 # Assume surface tension of water (mN/m) ???
const property_methods=Dict("bp"=>"joback_and_reid","vp"=>"nannoolal","critical"=>"nannoolal","density"=>"girolami")
