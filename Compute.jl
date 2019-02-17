#include("JlBoxModule.jl")
using Parse_eqn:parse_reactants,gen_evaluate_rates
using Optimize:constant_folding!,extract_constants!,generate_loss_gain,mk_reactants_list
using DifferentialEquations
using StaticArrays

function loss_gain!(num_reactants::Int,num_eqns::Int,
                   reactants::Array{Float64,1},#num_reactants
                   stoich_mtx::SparseMatrixCSC{Float64,Int64},#num_reactants*num_eqns
                   stoich_list::Array{Tuple{Int8,SVector{15,Int8},SVector{16,Int64}},1},#num_eqns, both reac and prod
                   reactants_list::Array{Tuple{Int8,SVector{15,Int8},SVector{16,Int64}},1},#num_eqns, only reac
                   rate_values::Array{Float64,1},#num_eqns
                   dydt::Array{Float64,1}#num_reactants
                   )
    lossgain_mtx=spzeros(num_reactants,num_eqns)
    for eqn_ind in 1:num_eqns
        prod=rate_values[eqn_ind]
        #reactant_inds=findn(stoich_mtx[:,eqn_ind])
        num_reacs,stoichvec,indvec=reactants_list[eqn_ind]
        num_stoichs,_,stoich_indvec=stoich_list[eqn_ind]
        for i in 1:num_reacs
            reactant_ind=indvec[i]
            stoich=stoichvec[i]
            #stoich=stoich_mtx[reactant_ind,eqn_ind]
            prod*=reactants[reactant_ind]^stoich#reactants_list come from reactants_mtx (for catalyse A+B=A+C)
            #if stoich>0
            #    prod*=reactants[reactant_ind]^stoich
            #end
        end
        #lossgain_mtx[:,eqn_ind]=stoich_mtx[:,eqn_ind]*prod #!!!!! HUGE PERFORMANCE COST
        #for reactant_ind in reactant_inds
        for i in 1:num_stoichs
            reactant_ind=stoich_indvec[i]
            lossgain_mtx[reactant_ind,eqn_ind]=stoich_mtx[reactant_ind,eqn_ind]*prod
        end
    end
    lossgain_mtx=transpose(lossgain_mtx)#num_eqns*num_reactants
    for reactant_ind in 1:num_reactants
        dydt[reactant_ind]=sum(nonzeros(lossgain_mtx[:,reactant_ind]))*(-1)#dydt negative for reactants, positive for products 
    end #*reactants[reactant_ind]=>wrong!!!
    return dydt
end

function dydt!(reactants::Array{Float64,1},p::Dict,t::Real)::Array{Float64,1}
    #dy,rate_values,J,stoich_mtx,stoich_list,reactants_list,RO2_inds,num_eqns,num_reactants=p
    dydt,rate_values,J,stoich_mtx,stoich_list,reactants_list,RO2_inds,num_eqns,num_reactants=
        [p[ind] for ind in 
            ["dydt","rate_values","J","stoich_mtx","stoich_list","reactants_list","RO2_inds",
             "num_eqns","num_reactants"]
        ]
    #dy,rate_values,rate_prods,J,RO2_inds,num_eqns,num_reactants=p
    time_of_day_seconds=start_time+t
    RO2=sum(reactants[RO2_inds])
    evaluate_rates!(time_of_day_seconds,RO2,H2O,temp,rate_values,J)# =>ratevalues
    loss_gain!(num_reactants,num_eqns,reactants,stoich_mtx,stoich_list,reactants_list,rate_values,dydt)
    #loss_gain_static!(num_reactants,num_eqns,reactants,rate_values,rate_prods,dy)
    return dydt
end

function dydt_aerosol!(y::Array{Float64,1},p::Dict,t::Real)::Array{Float64,1}
    num_reactants=p["num_reactants"]
    y_gas=y[1:num_reactants]#view(xs,lo:hi) passes ref instead of copy
    dy_dt_gas=dydt!(y_gas,p,t)
    C_g_i_t=y[include_inds]

end

function prepare_gas()
    println("Parsing Reactants")
    stoich_mtx,reactants_mtx,RO2_inds,num_eqns,num_reactants,reactants2ind=parse_reactants(file)
    reactants_list=mk_reactants_list(num_reactants,num_eqns,reactants_mtx)
    stoich_list=mk_reactants_list(num_reactants,num_eqns,stoich_mtx)

    
    @printf("num_eqns: %d, num_reactants: %d\n",num_eqns,num_reactants)

    println("Generating evaluate_rates()")
    evaluate_rates_expr=gen_evaluate_rates(file)
    println("Done Generation")
    rate_values=zeros(Float64,num_eqns)
    rate_prods=zeros(Float64,num_eqns)
    J=zeros(Float64,62)
    dydt=zeros(Float64,num_reactants)
    println("Performing constant folding")
    constant_folding!(evaluate_rates_expr,constantdict,rate_values);
    extract_constants!(evaluate_rates_expr);
    println("Evaluating evaluate_rates&loss_gain codes")
    eval(evaluate_rates_expr)
    println("Solving ODE")
    param_dict=Dict("dydt"=>dydt,"rate_values"=>rate_values,"J"=>J,"stoich_mtx"=>stoich_mtx,
                    "stoich_list"=>stoich_list,"reactants_list"=>reactants_list,"RO2_inds"=>RO2_inds,
                    "num_eqns"=>num_eqns,"num_reactants"=>num_reactants)
    return param_dict,reactants2ind
end

function run_simulation_aerosol()
    nothing
end


function run_simulation_gas()
    param_dict,reactants2ind=prepare_gas()
    num_reactants=param_dict["num_reactants"]
    reactants_initial=zeros(Float64,num_reactants)
    for (k,v) in reactants_initial_dict
        reactants_initial[reactants2ind[k]]=v*Cfactor#pbb to molcules/cc
    end
    prob = ODEProblem{false}(dydt!,reactants_initial,tspan,
                            param_dict,
                            #(dy,rate_values,J,stoich_mtx,stoich_list,reactants_list,RO2_inds,num_eqns,num_reactants)
                            )
    sol = solve(prob,CVODE_BDF(linear_solver=:Dense),reltol=1e-6,abstol=1.0e-3,
                tstops=0:batch_step:simulation_time,saveat=batch_step,# save_everystep=true,
                dt=1.0e-6, #Initial step-size
                dtmax=100.0,
                max_order = 5,
                max_convergence_failures = 1000,
                #progress=true
                )
    return sol,reactants2ind
end
