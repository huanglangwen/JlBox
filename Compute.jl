#include("JlBoxModule.jl")
using Parse:parse_reactants,gen_evaluate_rates
using Optimize:constant_folding!,extract_constants!
using DifferentialEquations

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
    reactants_initial=zeros(Float64,num_reactants)
    @printf("num_eqns: %d, num_reactants: %d\n",num_eqns,num_reactants)
    for (k,v) in reactants_initial_dict
        reactants_initial[reactants2ind[k]]=v*Cfactor#pbb to molcules/cc
    end
    println("Generating evaluate_rates()")
    evaluate_rates_expr=gen_evaluate_rates(file)
    println("Done Generation")
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
