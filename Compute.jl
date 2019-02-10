#include("JlBoxModule.jl")
using ..Parse_eqn:parse_reactants,gen_evaluate_rates
using ..Optimize:constant_folding!,extract_constants!,generate_loss_gain,mk_stoich_list

using StaticArrays
using SparseArrays
using Printf


function loss_gain!(num_reactants::Int,num_eqns::Int,
                   reactants::Array{Float64,1},#num_reactants
                   stoich_mtx::SparseMatrixCSC{Float64,Int64},#num_reactants*num_eqns
                   stoich_list::Array{Tuple{Int8,SVector{15,Int8},SVector{16,Int64}},1},#num_eqns
                   rate_values::Array{Float64,1},#num_eqns
                   dydt::Array{Float64,1}#num_reactants
                   )
    lossgain_mtx=spzeros(num_reactants,num_eqns)
    for eqn_ind in 1:num_eqns
        prod=rate_values[eqn_ind]
        #reactant_inds=findn(stoich_mtx[:,eqn_ind])
        num_stoichs,stoichvec,indvec=stoich_list[eqn_ind]
        for i in 1:num_stoichs
            reactant_ind=indvec[i]
            stoich=stoichvec[i]
            #stoich=stoich_mtx[reactant_ind,eqn_ind]
            if stoich>0
                prod*=reactants[reactant_ind]^stoich
            end
        end
        #lossgain_mtx[:,eqn_ind]=stoich_mtx[:,eqn_ind]*prod #!!!!! HUGE PERFORMANCE COST
        #for reactant_ind in reactant_inds
        for i in 1:num_stoichs
            reactant_ind=indvec[i]
            lossgain_mtx[reactant_ind,eqn_ind]=stoich_mtx[reactant_ind,eqn_ind]*prod
        end
    end
    lossgain_mtx=transpose(lossgain_mtx)#num_eqns*num_reactants
    for reactant_ind in 1:num_reactants
        dydt[reactant_ind]=sum(nonzeros(lossgain_mtx[:,reactant_ind]))*(-1)#dydt negative for reactants, positive for products 
    end #*reactants[reactant_ind]=>wrong!!!
    return dydt
end

dydt_expr=quote function dydt!(reactants,p,t)#::Array{Float64,1}
    @printf("t=%e\n",t)
    dy,rate_values,J,stoich_mtx,stoich_list,RO2_inds,num_eqns,num_reactants=p
    #dy,rate_values,rate_prods,J,RO2_inds,num_eqns,num_reactants=p
    time_of_day_seconds=start_time+t
    RO2=sum(reactants[RO2_inds])
    evaluate_rates!(time_of_day_seconds,RO2,H2O,temp,rate_values,J)# =>ratevalues
    #Base.invokelatest(evaluate_rates!,time_of_day_seconds,RO2,H2O,temp,rate_values,J)
    loss_gain!(num_reactants,num_eqns,reactants,stoich_mtx,stoich_list,rate_values,dy)
    #loss_gain_static!(num_reactants,num_eqns,reactants,rate_values,rate_prods,dy)
    return dy
end
end

function generate_codes()
    println("Parsing Reactants")
    stoich_mtx,RO2_inds,num_eqns,num_reactants,reactants2ind=parse_reactants(file)
    stoich_list=mk_stoich_list(num_reactants,num_eqns,stoich_mtx)
    #stoich_mtx_gpu=CudaSparseMatrixCSR(transpose(stoich_mtx))

    reactants_initial=zeros(Float64,num_reactants)
    @printf("num_eqns: %d, num_reactants: %d\n",num_eqns,num_reactants)
    for (k,v) in reactants_initial_dict
        reactants_initial[reactants2ind[k]]=v*Cfactor#pbb to molcules/cc
    end
    println("Generating evaluate_rates()")
    evaluate_rates_expr=gen_evaluate_rates(file)
    #println("Generating loss_gain_static()")
    #loss_gain_expr=generate_loss_gain(num_reactants,num_eqns,stoich_mtx)
    rate_values=zeros(Float64,num_eqns)
    rate_prods=zeros(Float64,num_eqns)
    J=zeros(Float64,62)
    dy=zeros(Float64,num_reactants)
    println("Done Generation")
    println("Performing constant folding")
    constant_folding!(evaluate_rates_expr,constantdict,rate_values);
    extract_constants!(evaluate_rates_expr);
    println("Evaluating evaluate_rates&loss_gain codes")
    open("generated_code.jl", "w") do f
        write(f,repr(evaluate_rates_expr.args[2])[3:end-1])# repr(expr)->:(code)--str[3:end-1]->code
        write(f,"\n")
        write(f,repr(dydt_expr.args[2])[3:end-1])
    end
    #include("generated_code.jl") #FATAL ERROR: ... AT ENDS OF VERY LONG LINES
    #@eval($evaluate_rates_expr)#function evaluate_rates!(ttime::Float64,
                             #RO2::Float64,H2O::Float64,temp::Float64,
                             #rate_values::Array{Float64,1},J::Array{Float64,1})
    #eval(loss_gain_expr)#function loss_gain_static!(num_reactants::Int,num_eqns::Int,
                                  #reactants::Array{Float64,1},#num_reactants
                                  #rate_values::Array{Float64,1},#num_eqns
                                  #rate_prods::Array{Float64,1},#num_eqns k*[A]^a*[B]^b
                                  #dydt::Array{Float64,1})#num_reactants
    return reactants2ind,reactants_initial,dy,rate_values,J,stoich_mtx,stoich_list,RO2_inds,num_eqns,num_reactants
end
