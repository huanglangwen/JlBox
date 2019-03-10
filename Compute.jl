using ..Parse_eqn:parse_reactants,parse_init
using ..Optimize:mk_reactants_list
using ..Jacobian:loss_gain_jac!
using DifferentialEquations
using StaticArrays
using SparseArrays
using Printf
using TimerOutputs
const to = TimerOutput()
function loss_gain!(num_reactants::Int,num_eqns::Int,
                   reactants::Array{Float64,1},#num_reactants
                   stoich_mtx::SparseMatrixCSC{Float64,Int64},#num_reactants*num_eqns
                   stoich_list::Array{Tuple{Int8,SVector{15,Int8},SVector{16,Int64}},1},#num_eqns, both reac and prod
                   reactants_list::Array{Tuple{Int8,SVector{15,Int8},SVector{16,Int64}},1},#num_eqns, only reac
                   rate_values::Array{Float64,1},#num_eqns
                   dydt::Array{Float64,1}#num_reactants
                   )
    lossgain_mtx=spzeros(num_reactants,num_eqns)
    @timeit to "Dydt part1" begin
        for eqn_ind in 1:num_eqns
            prod=rate_values[eqn_ind]
            num_reacs,stoichvec,indvec=reactants_list[eqn_ind]
            num_stoichs,_,stoich_indvec=stoich_list[eqn_ind]
            for i in 1:num_reacs
                reactant_ind=indvec[i]
                stoich=stoichvec[i]
                prod*=reactants[reactant_ind]^stoich#reactants_list come from reactants_mtx (for catalyse A+B=A+C)
            end
            for i in 1:num_stoichs
                reactant_ind=stoich_indvec[i]
                lossgain_mtx[reactant_ind,eqn_ind]=stoich_mtx[reactant_ind,eqn_ind]*prod
            end
        end
    end
    
    @timeit to "Dydt part2" begin
        is,js,vs=findnz(lossgain_mtx)
        lossgain_mtx_T=sparse(js,is,vs,num_eqns,num_reactants)
    end#lossgain_mtx_T=transpose(lossgain_mtx)#num_eqns*num_reactants
    @timeit to "Dydt part3" begin
        for reactant_ind in 1:num_reactants
            dydt[reactant_ind]=sum(nonzeros(lossgain_mtx_T[:,reactant_ind]))*(-1)#dydt negative for reactants, positive for products 
        end
    end
    nothing
end

function dydt!(dydt,reactants::Array{Float64,1},p::Dict,t::Real)
    rate_values,stoich_mtx,stoich_list,reactants_list,num_eqns,num_reactants=
        [p[ind] for ind in 
            ["rate_values","stoich_mtx","stoich_list","reactants_list",
             "num_eqns","num_reactants"]
        ]
    @timeit to "Dydt eval" loss_gain!(num_reactants,num_eqns,reactants,stoich_mtx,stoich_list,reactants_list,rate_values,dydt)
    p["Current_iter"]+=1
    citer=p["Current_iter"]
    if citer%(p["ShowIterPeriod"])==0
        @printf("Current Iteration: %d, time_step: %e\n",citer,t)
    end
    nothing
end


function gas_jac!(jac_mtx,reactants::Array{Float64,1},p::Dict,t::Real)
    rate_values,stoich_mtx,stoich_list,reactants_list,num_eqns,num_reactants=
        [p[ind] for ind in 
            ["rate_values","stoich_mtx","stoich_list","reactants_list",
             "num_eqns","num_reactants"]
        ]
    #Probably have to re-eval rate_values again
    @timeit to "Jacobian eval" loss_gain_jac!(num_reactants,num_eqns,reactants,stoich_mtx,stoich_list,reactants_list,rate_values,jac_mtx)
end

function prepare_gas()
    println("Parsing Reactants")
    rate_values,stoich_mtx,reactants_mtx,num_eqns,num_reactants,reactants2ind=parse_reactants(mechfile)
    reactants_init=parse_init(initfile,reactants2ind,num_reactants)
    reactants_list=mk_reactants_list(num_reactants,num_eqns,reactants_mtx)
    stoich_list=mk_reactants_list(num_reactants,num_eqns,stoich_mtx)
    @printf("num_eqns: %d, num_reactants: %d\n",num_eqns,num_reactants)

    param_dict=Dict("rate_values"=>rate_values,"stoich_mtx"=>stoich_mtx,#"dydt"=>dydt,
                    "stoich_list"=>stoich_list,"reactants_list"=>reactants_list,
                    "num_eqns"=>num_eqns,"num_reactants"=>num_reactants)
    return reactants_init,param_dict,reactants2ind
end

function read_configure!(filename::String)
    @printf("Reading Config file %s\n",filename)
    open(filename) do f
        for s in readlines(f)
            if (length(s)>2)
                if s[1]!='#'
                    eval(Meta.parse(s))#eval runs in Module scope while include runs in global scope
                end
            end
        end
    end
end

function gen_simulation_gas()
    read_configure!("Configure_gas.jl")
    reactants_init,param_dict,reactants2ind=prepare_gas()
    num_reactants=param_dict["num_reactants"]
    println("Solving ODE")
    param_dict["Current_iter"]=0
    param_dict["ShowIterPeriod"]=100
    odefun=ODEFunction(dydt!; jac=gas_jac!)
    prob = ODEProblem{true}(odefun,reactants_init,tspan,param_dict)
    #sol = solve(prob,solver,dense=false)#,reltol=1e-8,abstol=1.0e-8
                #callback=PositiveDomain(reactants_init,abstol=1.0e-3)
                #isoutofdomain=(u,p,t) -> any(x -> x < 0, u)
    return prob
end
