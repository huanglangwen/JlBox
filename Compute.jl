using ..Parse_eqn:parse_reactants,parse_init
using ..Optimize:mk_reactants_list
using ..Jacobian:loss_gain_jac!,mk_sparsity_pattern
using SparseDiffTools
using DiffEqBase
using StaticArrays
using SparseArrays
using Printf
using TimerOutputs
const to = TimerOutput()
function loss_gain!(num_reactants::Int,num_eqns::Int,
                   reactants,
                   lossgain_mtx::SparseMatrixCSC{<:Real,Int64},#num_reactants*num_eqns
                   stoich_mtx::SparseMatrixCSC{<:Real,Int64},#num_reactants*num_eqns
                   stoich_list::Array{Tuple{Int8,SVector{15,Int8},SVector{16,Int64}},1},#num_eqns, both reac and prod
                   reactants_list::Array{Tuple{Int8,SVector{15,Int8},SVector{16,Int64}},1},#num_eqns, only reac
                   rate_values,
                   dydt#num_reactants
                   )
    #lossgain_mtx=spzeros(num_reactants,num_eqns)
    fill!(lossgain_mtx.nzval,0.)
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
        lossgain_mtx_T=sparse(js,is,vs)
    end#lossgain_mtx_T=transpose(lossgain_mtx)#num_eqns*num_reactants
    @timeit to "Dydt part3" begin
        for reactant_ind in 1:num_reactants
            dydt[reactant_ind]=sum(nonzeros(lossgain_mtx_T[:,reactant_ind]))*(-1)#dydt negative for reactants, positive for products 
        end
    end
    nothing
end

function dydt!(dydt,reactants,p::Dict,t::Real)
    rate_values,stoich_mtx,stoich_list,reactants_list,lossgain_mtx,num_eqns,num_reactants=
        [p[ind] for ind in 
            ["rate_values","stoich_mtx","stoich_list","reactants_list","lossgain_mtx",
             "num_eqns","num_reactants"]
        ]
    @timeit to "Dydt eval" loss_gain!(num_reactants,num_eqns,reactants,lossgain_mtx,stoich_mtx,stoich_list,reactants_list,rate_values,dydt)
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

function prepare_gas(mechfile,initfile)
    println("Parsing Reactants")
    rate_values,stoich_mtx,reactants_mtx,num_eqns,num_reactants,reactants2ind=parse_reactants(mechfile)
    reactants_init=parse_init(initfile,reactants2ind,num_reactants)
    reactants_list=mk_reactants_list(num_reactants,num_eqns,reactants_mtx)
    stoich_list=mk_reactants_list(num_reactants,num_eqns,stoich_mtx)
    @printf("num_eqns: %d, num_reactants: %d\n",num_eqns,num_reactants)

    param_dict=Dict("rate_values"=>rate_values,"stoich_mtx"=>stoich_mtx,"lossgain_mtx"=>similar(stoich_mtx,Real),#"dydt"=>dydt,
                    "stoich_list"=>stoich_list,"reactants_list"=>reactants_list,
                    "num_eqns"=>num_eqns,"num_reactants"=>num_reactants)
    return reactants_init,param_dict,reactants2ind
end

function gen_simulation_gas(mechfile,initfile)
    simulation_time= 10.0 # seconds
    tspan=(0,simulation_time)
    reactants_init,param_dict,reactants2ind=prepare_gas(mechfile,initfile)
    sparsity=mk_sparsity_pattern(param_dict["num_reactants"],param_dict["num_eqns"],param_dict["stoich_mtx"],param_dict["stoich_list"],param_dict["reactants_list"])
    println("Coloring sparse Jacobian")
    colorvec=matrix_colors(sparsity)
    println("Producing ODEstruct")
    odefun=ODEFunction(dydt!; jac=gas_jac!)
    odefun_sp=ODEFunction(dydt!; colorvec=colorvec, jac_prototype=sparsity)
    prob=ODEProblem{true}(dydt!,reactants_init,tspan,param_dict)
    prob_jac=ODEProblem{true}(odefun,reactants_init,tspan,param_dict)
    prob_sp=ODEProblem{true}(odefun_sp,reactants_init,tspan,param_dict)
    #sol = solve(prob,solver,dense=false)#,reltol=1e-8,abstol=1.0e-8
                #callback=PositiveDomain(reactants_init,abstol=1.0e-3)
                #isoutofdomain=(u,p,t) -> any(x -> x < 0, u)
    return prob,prob_jac,prob_sp
end
