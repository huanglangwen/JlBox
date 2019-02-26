const file="data/MCM_APINENE.eqn.txt"#"MCM_test.eqn.txt"MCM_APINENE.eqn.txt
const temp=298.15 # Kelvin
const RH=0.5 # RH/100% [0 - 0.99]
const hour_of_day=12.0 # Define a start time  24 hr format
const start_time=hour_of_day*60*60 # seconds, used as t0 in solver
const simulation_time= 10800.0 # seconds
const batch_step=100.0 # seconds
const temp_celsius=temp-273.15
const Psat=610.78*exp((temp_celsius/(temp_celsius+238.3))*17.2694)# Saturation VP of water vapour, to get concentration of H20
const Pw=RH*Psat
const Wconc=0.002166*(Pw/(temp_celsius+273.16))*1.0e-6 #kg/cm3
const H2O=Wconc*(1.0/(18.0e-3))*6.0221409e+23#Convert from kg to molecules/cc
const tspan=(0,simulation_time)
const Cfactor= 2.55e+10 #ppb-to-molecules/cc
const reactants_initial_dict=Dict(["O3"=>18.0,"APINENE"=>30.0])#ppm ["O3"=>18.0,"BUT1ENE"=>30.0])
const constantdict=Dict([(:temp,temp),(:H2O,H2O)])