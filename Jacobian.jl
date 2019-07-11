using StaticArrays
using SparseArrays
function loss_gain_jac!(num_reactants::Int,num_eqns::Int,
                       reactants::Array{<:Real,1},#num_reactants
                       stoich_mtx::SparseMatrixCSC{<:Real,Int64},#num_reactants*num_eqns
                       stoich_list::Array{Tuple{Int8,SVector{15,Int8},SVector{16,Int64}},1},#num_eqns, both reac and prod
                       reactants_list::Array{Tuple{Int8,SVector{15,Int8},SVector{16,Int64}},1},#num_eqns, only reac
                       rate_values::Array{<:Real,1},#num_eqns
                       lossgain_jac_mtx#::Array{Float64,2}#SparseMatrixCSC{Float64,Int64},#num_output(dydt)*num_input(y)
                       )
    #lossgain_jac_mtx=spzeros(num_reactants,num_reactants)#num_output(dydt)*num_input(y)
    for eqn_ind in 1:num_eqns
        num_reacs,stoichvec,indvec=reactants_list[eqn_ind]
        num_stoichs,_,stoich_indvec=stoich_list[eqn_ind]
        for y_ind in 1:num_reacs
            prod=rate_values[eqn_ind]
            for i in [i for i in 1:num_reacs if i!=y_ind]
                reactant_ind=indvec[i]
                stoich=stoichvec[i]
                prod*=reactants[reactant_ind]^stoich#reactants_list come from reactants_mtx (for catalyse A+B=A+C)
            end
            reactant_y_ind=indvec[y_ind]
            stoich_y::Integer=stoichvec[y_ind]
            prod*=stoich_y*reactants[reactant_y_ind]^(stoich_y-1)
            for i in 1:num_stoichs
                reactant_ind=stoich_indvec[i]
                lossgain_jac_mtx[reactant_ind,reactant_y_ind]+=stoich_mtx[reactant_ind,eqn_ind]*prod*(-1)
            end
        end
    end
    return lossgain_jac_mtx
end

function mk_sparsity_pattern(num_reactants::Int,num_eqns::Int,
    stoich_mtx::SparseMatrixCSC{<:Real,Int64},#num_reactants*num_eqns
    stoich_list::Array{Tuple{Int8,SVector{15,Int8},SVector{16,Int64}},1},#num_eqns, both reac and prod
    reactants_list::Array{Tuple{Int8,SVector{15,Int8},SVector{16,Int64}},1}#num_eqns, only reac
    )
    lossgain_jac_mtx=spzeros(num_reactants,num_reactants)
    for eqn_ind in 1:num_eqns
        num_reacs,stoichvec,indvec=reactants_list[eqn_ind]
        num_stoichs,_,stoich_indvec=stoich_list[eqn_ind]
        for y_ind in 1:num_reacs
            reactant_y_ind=indvec[y_ind]
            for i in 1:num_stoichs
                reactant_ind=stoich_indvec[i]
                lossgain_jac_mtx[reactant_ind,reactant_y_ind]=1
            end
        end
    end
    return lossgain_jac_mtx
end
