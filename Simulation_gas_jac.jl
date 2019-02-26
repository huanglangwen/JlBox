include("JlBoxModule.jl")
using .Compute:prepare_gas,read_configure!,loss_gain!
using .Jacobian:gas_jac!
using DifferentialEquations
using SparseArrays
using DataFrames
using CSV
using Printf

include("Configure_gas.jl")
read_configure!("Configure_gas.jl")
param_dict,reactants2ind,evaluate_rates_expr=prepare_gas()

dydt_expr=quote function dydt!(dydt::Array{Float64,1},reactants::Array{Float64,1},p::Dict,t::Real)
    #dy,rate_values,J,stoich_mtx,stoich_list,reactants_list,RO2_inds,num_eqns,num_reactants=p
    rate_values,J,stoich_mtx,stoich_list,reactants_list,RO2_inds,num_eqns,num_reactants=
        [p[ind] for ind in 
            ["rate_values","J","stoich_mtx","stoich_list","reactants_list","RO2_inds",
             "num_eqns","num_reactants"]
        ]
    #dy,rate_values,rate_prods,J,RO2_inds,num_eqns,num_reactants=p
    time_of_day_seconds=start_time+t
    RO2=sum(reactants[RO2_inds])
    evaluate_rates!(time_of_day_seconds,RO2,H2O,temp,rate_values,J)# =>ratevalues
    loss_gain!(num_reactants,num_eqns,reactants,stoich_mtx,stoich_list,reactants_list,rate_values,dydt)
    #loss_gain_static!(num_reactants,num_eqns,reactants,rate_values,rate_prods,dy)
    p["iter"]+=1
    if p["iter"]%200==0
        @printf("time:%e,APINENE:%e\n",t,reactants[1])
    end
    nothing
end
end
open("generated_code.jl", "w") do f
    write(f,repr(evaluate_rates_expr.args[2])[3:end-1])# repr(expr)->:(code)--str[3:end-1]->code
    write(f,"\n")
    write(f,repr(dydt_expr.args[2])[3:end-1])
end
include("generated_code.jl")

num_reactants=param_dict["num_reactants"]
reactants_initial=zeros(Float64,num_reactants)
for (k,v) in reactants_initial_dict
    reactants_initial[reactants2ind[k]]=v*Cfactor#pbb to molcules/cc
end
#lossgain_jac_mtx=zeros(num_reactants,num_reactants)#num_output(dydt)*num_input(y)
println("Solving ODE with Jacobian")
odefun=ODEFunction(dydt!; jac=gas_jac!)#, jac_prototype=lossgain_jac_mtx
param_dict["iter"]=0
prob = ODEProblem(dydt!,reactants_initial,tspan,param_dict)
@time sol = solve(prob,CVODE_BDF(linear_solver=:Dense),reltol=1e-6,abstol=1.0e-3,#,Rodas5(autodiff=false)
            tstops=0:batch_step:simulation_time,saveat=batch_step,# save_everystep=true,
            dt=1.0e-6, #Initial step-size
            dtmax=100.0,
            max_order = 5,
            max_convergence_failures = 1000,
            #progress=true
            )
num_reactants=length(reactants2ind)
ind2reactants=Dict(reactants2ind[key]=>key for key in keys(reactants2ind))
reactants=[ind2reactants[ind] for ind in 1:num_reactants]
df=DataFrame(transpose(sol))
names!(df,[Symbol(reac) for reac in reactants])
CSV.write("data/results_jac.csv",df)
df
