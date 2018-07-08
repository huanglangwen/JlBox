include("JlBoxModule.jl")
using Parse:parse_reactants,gen_evaluate_rates,constant_folding!,extract_constants!
using DifferentialEquations

#global start_time
const file="MCM_test.eqn.txt"#"MCM_APINENE.eqn.txt"
const temp=288.15 # Kelvin
const RH=0.5 # RH/100% [0 - 0.99]
#Define a start time 
const hour_of_day=12.0 # 24 hr format
const start_time=hour_of_day*60*60 # seconds, used as t0 in solver
const simulation_time= 7200.0 # seconds
const batch_step=100.0 # seconds
#2)Generate constants used in rate of reaction calculations
#Convert RH to concentration of water vapour molecules [this will change when in Parcel model mode]
const temp_celsius=temp-273.15
# Saturation VP of water vapour, to get concentration of H20
const Psat=610.78*exp((temp_celsius/(temp_celsius+238.3))*17.2694)
const Pw=RH*Psat
const Wconc=0.002166*(Pw/(temp_celsius+273.16))*1.0e-6 #kg/cm3
#Convert from kg to molecules/cc
const H2O=Wconc*(1.0/(18.0e-3))*6.0221409e+23
const tspan=(0,simulation_time)
const Cfactor= 2.55e+10 #ppb-to-molecules/cc

function loss_gain!(num_reactants::Int,num_eqns::Int,
                   reactants::Array{Float64,1},#num_reactants
                   stoich_mtx::SparseMatrixCSC{Float64,Int64},#num_reactants*num_eqns
                   rate_values::Array{Float64,1},#num_eqns
                   dydt::Array{Float64,1}#num_reactants
                   )
    lossgain_mtx=spzeros(num_reactants,num_eqns)
    for eqn_ind in 1:num_eqns
        prod=rate_values[eqn_ind]
        reactant_inds=findn(stoich_mtx[:,eqn_ind])
        for reactant_ind in reactant_inds
            stoich=stoich_mtx[reactant_ind,eqn_ind]
            if stoich>0
                prod*=reactants[reactant_ind]^stoich
            end
        end
        #lossgain_mtx[eqn_ind,:]=stoich_mtx[eqn_ind,:]*prod
        for reactant_ind in reactant_inds
            lossgain_mtx[reactant_ind,eqn_ind]=stoich_mtx[reactant_ind,eqn_ind]*prod
        end
    end
    lossgain_mtx=transpose(lossgain_mtx)#num_eqns*num_reactants
    for reactant_ind in 1:num_reactants
        dydt[reactant_ind]=sum(nonzeros(lossgain_mtx[:,reactant_ind]))*(-1)#dydt negative for reactants, positive for products 
    end #*reactants[reactant_ind]=>wrong!!!
    return dydt
end

function dydt!(reactants::Array{Float64,1},p,t)::Array{Float64,1}
    dy,rate_values,J,stoich_mtx,RO2_inds,num_eqns,num_reactants=p
    time_of_day_seconds=start_time+t
    RO2=sum(reactants[RO2_inds])
    evaluate_rates!(time_of_day_seconds,RO2,H2O,temp,rate_values,J)# =>ratevalues
    loss_gain!(num_reactants,num_eqns,reactants,stoich_mtx,rate_values,dy)
    return dy
end

function run_simulation()
    println("Parsing Reactants")
    stoich_mtx,RO2_inds,num_eqns,num_reactants,reactants2ind=parse_reactants(file)
    reactants_initial_dict=Dict(["O3"=>18.0,"BUT1ENE"=>30.0])#ppm 
    reactants_initial=zeros(Float64,num_reactants)
    @printf("num_eqns: %d, num_reactants: %d\n",num_eqns,num_reactants)
    for (k,v) in reactants_initial_dict
        reactants_initial[reactants2ind[k]]=v*Cfactor#pbb to molcules/cc
    end
    println("Generating evaluate_rates()")
    evaluate_rates_expr=gen_evaluate_rates(file)
    println("Done Generation")
    constantdict=Dict([(:temp,temp),(:H2O,H2O)])
    rate_values=zeros(Float64,num_eqns)
    J=zeros(Float64,62)
    dy=zeros(Float64,num_reactants)
    println("Performing constant folding")
    constant_folding!(evaluate_rates_expr,constantdict,rate_values);
    extract_constants!(evaluate_rates_expr);
    println("Evaluating evaluate_rates codes")
    eval(evaluate_rates_expr)#function evaluate_rates!(ttime::Float64,
                             #RO2::Float64,H2O::Float64,temp::Float64,
                             #rate_values::Array{Float64,1},J::Array{Float64,1})
    println("Solving ODE")
    prob = ODEProblem{false}(dydt!,reactants_initial,tspan,(dy,rate_values,J,stoich_mtx,RO2_inds,num_eqns,num_reactants))
    sol = solve(prob,CVODE_BDF(linear_solver=:Dense),reltol=1e-6,abstol=1.0e-3,
                tstops=0:batch_step:simulation_time,saveat=batch_step,# save_everystep=true,
                dt=1.0e-6, #Initial step-size
                dtmax=100.0,
                max_order = 5,
                max_convergence_failures = 1000,
                progress=true
                )
    return sol,reactants2ind
end

@time sol,reactants2ind=run_simulation()
#Profile.init(n = 10^7, delay = 5.)
#@profile run_simulation()
#open("prof.txt", "w") do s
#    Profile.print(IOContext(s, :displaysize => (1000, 500)))
#end
#using ProfileView
#ProfileView.view()
#Profile.clear()
using Plots
plot(log10.(sol[reactants2ind["BUT1ENE"],:]))
plot!(log10.(sol[reactants2ind["C2H5O2"],:]))
